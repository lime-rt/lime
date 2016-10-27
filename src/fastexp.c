/*
 *  fastexp.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */

#include "lime.h"

void fillErfTable() {
  int i;
  for (i=0;i<ERF_TABLE_SIZE;i++) {
    ERF_TABLE[i]=(SPI/2.)*erf(i*ERF_TABLE_LIMIT/(ERF_TABLE_SIZE-1.));
  }
}


int factorial(const int n){
  int i, result;

  result = 1.0;
  for(i=1;i<n+1;i++){
    result *= i;
  }
  return result;
}

double taylor(const int maxOrder, const float x){
  double result=1.0;
  int i;

  for (i=maxOrder;i>0;i--){
    result = 1.0 + x*result/(double)i;
  }

  return result;
}


void calcFastExpRange(const int maxTaylorOrder, const int maxNumBitsPerMantField, int *numMantissaFields, int *lowestExponent, int *numExponentsUsed){
  /*
We want to approximate an exponential call exp(-|x|) (where x is a standard, i.e. IEEE 754-format, float) by a look-up table. The scheme proposed makes use of two arrays: a 2D one of shape (J=2^{B-B0},L) and a 3D one of shape (J=2^B,K-1,L), where B is the input argument maxNumBitsPerMantField. The purpose of the present function is to calculate and return B0, K and L, as well as the value of lowestExponent (explained in section 2 below).

1) The fast-exp lookup algorithm.
=================================
The algorithm makes use of the fact that exp(a+b+c) = exp(a)*exp(b)*exp(c) to return a full precision exponential value via the product of several separate lower-precision lookups. This is done by dividing the mantissa of |x| into K contiguous fields, where K is calculated from

	K = ceiling(23/B).

Note that the IEEE 754 standard specifies that the mantissa of a float is 23 bits long. Each of these mantissa fields, considered as an integer, can generate a value for the first index j of the appropriate array. (The 0th field lookup values are in the 2D array, but every other field accesses values in the 3D array.) For field numbers k>0, the second index of the 3D array is found from k-1; and the third index l is taken by treating the bits of the float exponent (bits 1 to 9 in the IEEE 754 standard) as in integer, suitably offset as described below.

The detailed working of the algorithm is illustrated with an example. We construct an example float x equal to -1.3a07b2x * 2^3. This should have a IEEE 754 bit representation as follows:

	    seee eeee  emmm mmmm  mmmm mmmm  mmmm mmmm
	x = 1100 0001  0001 1101  0000 0011  1101 1001

or, in hex,

	0xc11d03d9

If we set B=8, that gives the number of fields K=3. We thus divide the mantissa into three fields a, b and c as follows:

	00111010000001111011001
	aaaaaaabbbbbbbbcccccccc

Note that field a has only 7 bits, not 8: because 8 does not divide evenly into 23. It is this shorter k=0 field which necessitates the use of a separate array for this field. The value B0 is the number of 'missing' bits in this first field - i.e. =1 in the present case.

The example above will yield

	a = 0x1d
	b = 0x03
	c = 0xd9.

Suppose we construct a new IEEE 754 bit representation for each of these. Starting with 'a', this gives

	     seee eeee  emmm mmmm  mmmm mmmm  mmmm mmmm
	A0 = 1100 0001  0001 1101  0000 0000  0000 0000

which is -1.3ax * 2^3 as desired. For 'b' however, if we naively shift the 'b' bit pattern 7 bits to the left, and decrement the exponent by 7, we get

	     seee eeee  emmm mmmm  mmmm mmmm  mmmm mmmm
	A1 = 1011 1101  1000 0001  1000 0000  0000 0000

According to the IEEE 754 rules for encoding floats, this is not -0.03x * 2^-4 as we need, but -1.03x * 2^-4. Similar will hold for 'c'. I.e. we left-shift the 'c' bit pattern by 15 bits, decrementing the exponent also by 15, to give

	     seee eeee  emmm mmmm  mmmm mmmm  mmmm mmmm
	A2 = 1011 1001  1110 1100  1000 0000  0000 0000

This is not -0.d9x * 2^-12 as we need, but -1.d9x * 2^-12. Thus we must add 1 * 2^{e-7} and 1 * 2^{e-15} to the total afterwards, where e is the original exponent (in the present example = 3). Thus we have

	x = A0 + (A1 + 1*2^{e-7}) + (A2 + 1*2^{e-15})

and thus

	exp(x) = exp(A0)*exp(A1 + 1*2^{e-7})*exp(A2 + 1*2^{e-15}).

2) Calculating the range of exponents.
======================================
Values of |x| which are smaller than a given cutoff are calculated via a Pth-order Taylor expansion; values greater than another cutoff return a value of zero. Ideally, we set the low-x cutoff to the point at which the absolute value of the next-higher-order Taylor term equals the floating-point precision, and the high-x cutoff to the point where exp(-|x|) reaches the same precision.

The fractional precision of the lookup table is clearly epsilon=1/2^23, the same precision as an ordinary floating-point number. This will also be approximately the absolute precision sigma for small values of |x| where exp(-|x|) ~ 1. We should choose the low cutoff value x_lo such that the next higher term (P+1th term) of the Taylor series is equal to sigma, i.e.:

	 |x_lo|^{P+1}      1
	-------------- = ------.
	    {P+1}!        2^23

Having calculated |x_lo|, we want to find that exponent N_lo of 2 such that 2^N_lo < |x_lo| < 2^{N_lo+1}.

The second part of the calculation is to find N_hi such that exp(-2^N_hi) < sigma < exp(-2^{N_hi-1}). The number of exponents L which needed is thus 1+N_hi-N_lo.
  */

  int ieee754NumMantBits=23;
  double sigma, xLo;
  int nLo, nHi;

  *numMantissaFields = 1+floor(ieee754NumMantBits/(float)maxNumBitsPerMantField);

  sigma = 1.0/pow(2.,ieee754NumMantBits);
  xLo = pow((double)factorial(maxTaylorOrder+1)*sigma, 1/(double)(maxTaylorOrder+1));
  nLo = floor(log(xLo)/log(2.));

  /*
	exp(-2^N_hi) < sigma
thus
	2^N_hi > -ln(sigma)
thus
	N_hi > ln(-ln(sigma))/ln(2).
  */
  nHi = 1+floor(log(-log(sigma))/log(2.));

  *lowestExponent = nLo;
  *numExponentsUsed = 1+nHi-nLo;
}

void calcTableEntries(const int maxTaylorOrder, const int maxNumBitsPerMantField){
  /*
See description of the lookup algorithm in function calcFastExpRange().
  */

  int ieee754ExpOffset=127,ieee754NumMantBits=23;
  int negativeSignMask=0x80000000;
  int numJs,numJs0,numMantissaFields,lowestExponent,numExponentsUsed,exponentOffset,mantShift,bitOffset0,fieldBitOffset,fieldI,j,k,l,exponentMask;
  float argOffset;
  union
  {
    float f;
    int m;
  } floPo;

  // Should raise an exception here #ifndef FASTEXP?

  calcFastExpRange(maxTaylorOrder, maxNumBitsPerMantField, &numMantissaFields, &lowestExponent, &numExponentsUsed);

  exponentOffset = ieee754ExpOffset + lowestExponent;
  mantShift = ieee754NumMantBits - maxNumBitsPerMantField;
  bitOffset0 = maxNumBitsPerMantField*numMantissaFields - ieee754NumMantBits;
  numJs  = (int)pow(2.,maxNumBitsPerMantField);
  numJs0 = (int)pow(2.,maxNumBitsPerMantField-bitOffset0);

  fieldBitOffset = 0.0;
  for (l=0;l<numExponentsUsed;l++){
    argOffset = 0.0;
    exponentMask = (l+exponentOffset+fieldBitOffset)<<ieee754NumMantBits;

    for (j=0;j<numJs0;j++){
      floPo.m = negativeSignMask | exponentMask | (j<<(mantShift+bitOffset0));
      EXP_TABLE_2D[j][l] = exp(floPo.f + argOffset);
    }
  }

  for (fieldI=1;fieldI<numMantissaFields;fieldI++){
    k = fieldI-1;
    fieldBitOffset = bitOffset0 - fieldI*maxNumBitsPerMantField;

    for (l=0;l<numExponentsUsed;l++){
      argOffset = pow(2.,l+lowestExponent+fieldBitOffset);
      exponentMask = (l+exponentOffset+fieldBitOffset)<<ieee754NumMantBits;

      for (j=0;j<numJs;j++){
        floPo.m = negativeSignMask | exponentMask | (j<<mantShift);
        EXP_TABLE_3D[j][k][l] = exp(floPo.f + argOffset);
      }
    }
  }

  /*We also construct the table of 1/i to be used for faster calculation of the Taylor approximation.*/
  oneOver_i[0]=0.0;
  for (j=1;j<=FAST_EXP_MAX_TAYLOR;j++) oneOver_i[j]=1.0/(1.0*j);
}

double
FastExp(const float negarg){
  /*
See description of the lookup algorithm in function calcFastExpRange(). ****NOTE!**** Most numbers here are hard-wired for the sake of speed. If need be, they can be verified (or recalculated for different conditions) via calcTableEntries().
  */
  int exponentMask=0x7f800000,ieee754NumMantBits=23;
  int exponentOffset=122,numExponentsUsed=10;
  /*
This value should be calculated from 127+lowestExponent, where 127 is the offset for an exponent of zero laid down in the IEEE 754 standard, and both lowestExponent and numExponentsUsed can be calculated via calcFastExpRange().

  exponentOffset = ieee754ExpOffset + lowestExponent;
  */

  int mantMask0=0x007f0000, mantMask1=0x0000ff00, mantMask2=0x000000ff;
  int mantOffset0=16, mantOffset1=8, mantOffset2=0;
  int i,j0,j1,j2,l;
  union
  {
    float f;
    int m;
  } floPo;
  double result;

  // Should raise an exception here #ifndef FASTEXP?

  if (negarg<0.0) return exp(-negarg);
  if (negarg==0.0) return 1.0;

  floPo.f = negarg;
  l = ((floPo.m & exponentMask)>>ieee754NumMantBits)-exponentOffset;

  if (l<0){ // do the Taylor approximation.
    result = 1.0;
    for (i=FAST_EXP_MAX_TAYLOR;i>0;i--){
      result = 1.0 - negarg*result*oneOver_i[i];
    }
    return result;

  }else if(l>=numExponentsUsed){
    return 0.0;
  }

  j0 = (floPo.m & mantMask0)>>mantOffset0;
  j1 = (floPo.m & mantMask1)>>mantOffset1;
  j2 = (floPo.m & mantMask2)>>mantOffset2;

  return (EXP_TABLE_2D[j0]   [l]*
          EXP_TABLE_3D[j1][0][l]*
          EXP_TABLE_3D[j2][1][l]);
}



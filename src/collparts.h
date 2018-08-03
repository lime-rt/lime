/*
 *  collparts.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  See ../COPYRIGHT
 *
 */

#ifndef COLLPARTS_H
#define COLLPARTS_H

/* Collision partner ID numbers from LAMDA */
#define CP_H2			1
#define CP_p_H2			2
#define CP_o_H2			3
#define CP_e			4
#define CP_H			5
#define CP_He			6
#define CP_Hplus		7

/* Bit codes for par->collPartUserSetFlags */
#define CPF_BIT_ids        0
#define CPF_BIT_weights    1
#define CPF_BIT_names      2
#define CPF_BIT_MolWeights 3

#endif /* COLLPARTS_H */

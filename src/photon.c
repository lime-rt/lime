/*
 *  photon.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */


#include "lime.h"


/*....................................................................*/
double
geterf(const double x0, const double x1) {
  /* table lookup erf thingy */

  double val0=0.,val1=0.;
  
  if (fabs(x0)>=ERF_TABLE_LIMIT) val0=(SPI/2.);
  else {
    int index = (int)(fabs(x0*IBIN_WIDTH));
    double inter_coeff = (fabs(x0*IBIN_WIDTH)-index);
    val0=(1-inter_coeff)*ERF_TABLE[index]+inter_coeff*ERF_TABLE[index+1];
  }
  if (x0<0.) val0=-val0;
 
  if (fabs(x1)>=ERF_TABLE_LIMIT) val1=(SPI/2.);
  else {
    int index = (int)(fabs(x1*IBIN_WIDTH));
    double inter_coeff = (fabs(x1*IBIN_WIDTH)-index);
    val1=(1-inter_coeff)*ERF_TABLE[index]+inter_coeff*ERF_TABLE[index+1];
  }
  if (x1<0.) val1=-val1;

  return fabs((val1-val0)/(x1-x0));
}

/*....................................................................*/
double
gaussline(const double v, const double oneOnSigma){
  double val;
  val = v*v*oneOnSigma*oneOnSigma;
#ifdef FASTEXP
  return FastExp(val);
#else
  return exp(-val);
#endif
}

/*....................................................................*/
int
getNextEdge(double *inidir, const int startGi, const int presentGi\
  , struct grid *gp, const gsl_rng *ran){
  /*
The idea here is to select for the next grid point, that one which lies closest (with a little randomizing jitter) to the photon track, while requiring the direction of the edge to be in the 'forward' hemisphere of the photon direction.
  */
  int i,ni,niOfSmallest,niOfNextSmallest;
  double dirCos,distAlongTrack,dirFromStart[3],coord,distToTrackSquared,smallest,nextSmallest;
  const static double scatterReduction = 0.4;
  /*
This affects the ratio of N_2/N_1, where N_2 is the number of times the edge giving the 2nd-smallest distance from the photon track is chosen and N_1 ditto the smallest. Some ratio values obtained from various values of scatterReduction:

	scatterReduction	<N_2/N_1>
		1.0		  0.42
		0.5		  0.75
		0.4		  0.90
		0.2		  1.52

Note that the equivalent ratio value produced by the 1.6 code was 0.91.
  */

  i = 0;
  for(ni=0;ni<gp[presentGi].numNeigh;ni++){
    dirCos = dotProduct3D(inidir, gp[presentGi].dir[ni].xn);

    if(dirCos<=0.0)
  continue; /* because the edge points in the backward direction. */

    dirFromStart[0] = gp[presentGi].neigh[ni]->x[0] - gp[startGi].x[0];
    dirFromStart[1] = gp[presentGi].neigh[ni]->x[1] - gp[startGi].x[1];
    dirFromStart[2] = gp[presentGi].neigh[ni]->x[2] - gp[startGi].x[2];
    distAlongTrack = dotProduct3D(inidir, dirFromStart);

    coord = dirFromStart[0] - distAlongTrack*inidir[0];
    distToTrackSquared  = coord*coord;
    coord = dirFromStart[1] - distAlongTrack*inidir[1];
    distToTrackSquared += coord*coord;
    coord = dirFromStart[2] - distAlongTrack*inidir[2];
    distToTrackSquared += coord*coord;

    if(i==0){
      smallest = distToTrackSquared;
      niOfSmallest = ni;
    }else{
      if(distToTrackSquared<smallest){
        nextSmallest = smallest;
        niOfNextSmallest = niOfSmallest;
        smallest = distToTrackSquared;
        niOfSmallest = ni;
      }else if(i==1 || distToTrackSquared<nextSmallest){
        nextSmallest = distToTrackSquared;
        niOfNextSmallest = ni;
      }
    }

    i++;
  }

  /* Choose the edge to follow.
  */
  if(i>1){ /* then nextSmallest, niOfNextSmallest should exist. */
    if((smallest + scatterReduction*nextSmallest)*gsl_rng_uniform(ran)<smallest){
      return niOfNextSmallest;
    }else{
      return niOfSmallest;
    }
  }else if(i>0){
    return niOfSmallest;
  }else{
    if(!silent)
      bail_out("Photon propagation error - no valid edges.");
    exit(1);
  }
}

/*....................................................................*/
void calcLineAmpPWLin(struct grid *g, const int id, const int k\
  , const int molI, const double deltav, double *inidir, double *vfac_in, double *vfac_out){

  /* convolution of a Gaussian with a box */
  double binv_this, binv_next, v[5];

  binv_this=g[id].mol[molI].binv;
  binv_next=(g[id].neigh[k])->mol[molI].binv;
  v[0]=deltav-dotProduct3D(inidir,g[id].vel);
  v[1]=deltav-dotProduct3D(inidir,&(g[id].v1[3*k]));
  v[2]=deltav-dotProduct3D(inidir,&(g[id].v2[3*k]));
  v[3]=deltav-dotProduct3D(inidir,&(g[id].v3[3*k]));
  v[4]=deltav-dotProduct3D(inidir,g[id].neigh[k]->vel);

  /* multiplying by the appropriate binv changes from velocity to doppler widths(?) */
  /* if the values were be no more than 2 erf table bins apart, we just take a single Gaussian */

  /*
  vfac_out is the lineshape for the part of the edge in the current Voronoi cell,
  vfac_in is for the part in the next cell
  */

  if (fabs(v[1]-v[0])*binv_this>(2.0*BIN_WIDTH)) {
     *vfac_out=0.5*geterf(v[0]*binv_this,v[1]*binv_this);
  } else *vfac_out=0.5*gaussline(0.5*(v[0]+v[1]),binv_this);
  if (fabs(v[2]-v[1])*binv_this>(2.0*BIN_WIDTH)) {
    *vfac_out+=0.5*geterf(v[1]*binv_this,v[2]*binv_this);
  } else *vfac_out+=0.5*gaussline(0.5*(v[1]+v[2]),binv_this);

  if (fabs(v[3]-v[2])*binv_next>(2.0*BIN_WIDTH)) {
     *vfac_in=0.5*geterf(v[2]*binv_next,v[3]*binv_next);
  } else *vfac_in=0.5*gaussline(0.5*(v[2]+v[3]),binv_next);
  if (fabs(v[4]-v[3])*binv_next>(2.0*BIN_WIDTH)) {
    *vfac_in+=0.5*geterf(v[3]*binv_next,v[4]*binv_next);
  } else *vfac_in+=0.5*gaussline(0.5*(v[3]+v[4]),binv_next);
}


void calcLineAmpLin(struct grid *g, const int id, const int k\
  , const int molI, const double deltav, double *inidir, double *vfac_in, double *vfac_out){


  /* convolution of a Gaussian with a box */
  double binv_this, binv_next, v[3];

  binv_this=g[id].mol[molI].binv;
  binv_next=(g[id].neigh[k])->mol[molI].binv;
  v[0]=deltav-dotProduct3D(inidir,g[id].vel);
  v[2]=deltav-dotProduct3D(inidir,g[id].neigh[k]->vel);
  v[1]=0.5*(v[0]+v[2]);

  if (fabs(v[1]-v[0])*binv_this>(2.0*BIN_WIDTH)) {
     *vfac_out=geterf(v[0]*binv_this,v[1]*binv_this);
  } else *vfac_out+=gaussline(0.5*(v[0]+v[1]),binv_this);

  if (fabs(v[2]-v[1])*binv_next>(2.0*BIN_WIDTH)) {
     *vfac_in=geterf(v[1]*binv_next,v[2]*binv_next);
  } else *vfac_in+=gaussline(0.5*(v[1]+v[2]),binv_next);
}

/*....................................................................*/
void calcSourceFn(double dTau, const configInfo *par, double *remnantSnu, double *expDTau){
  /*
  The source function S is defined as j_nu/alpha, which is clearly not
  defined for alpha==0. However S is used in the algorithm only in the
  term (1-exp[-alpha*ds])*S, which is defined for all values of alpha.
  The present function calculates this term and returns it in the
  argument remnantSnu. For values of abs(alpha*ds) less than a pre-
  calculated cutoff supplied in configInfo, a Taylor approximation is
  used.

  Note that the same cutoff condition holds for replacement of
  exp(-dTau) by its Taylor expansion to 3rd order.
  */

#ifdef FASTEXP
  *expDTau = FastExp(dTau);
  if (fabs(dTau)<par->taylorCutoff){
    *remnantSnu = 1. - dTau*(1. - dTau*(1./3.))*(1./2.);
  } else {
    *remnantSnu = (1.-(*expDTau))/dTau;
  }
#else
  if (fabs(dTau)<par->taylorCutoff){
    *remnantSnu = 1. - dTau*(1. - dTau*(1./3.))*(1./2.);
    *expDTau = 1. - dTau*(*remnantSnu);
  } else {
    *expDTau = exp(-dTau);
    *remnantSnu = (1.-(*expDTau))/dTau;
  }
#endif
}

/*....................................................................*/
void
photon(int id, struct grid *gp, molData *md, const gsl_rng *ran\
  , configInfo *par, const int nlinetot, struct blendInfo blends\
  , gridPointData *mp, double *halfFirstDs){

  int iphot,iline,here,there,firststep,neighI;
  int nextMolWithBlend, nextLineWithBlend, molI, lineI, molJ, lineJ, bi;
  double segment,vblend_in,vblend_out,dtau,expDTau,ds_in=0.0,ds_out=0.0,pt_theta,pt_z,semiradius;
  double deltav[par->nSpecies],vfac_in[par->nSpecies],vfac_out[par->nSpecies],vfac_inprev[par->nSpecies];
  double expTau[nlinetot],inidir[3];
  double remnantSnu,velProj;

  for(iphot=0;iphot<gp[id].nphot;iphot++){
    firststep=1;
    iline = 0;
    for(molI=0;molI<par->nSpecies;molI++){
      for(lineI=0;lineI<md[molI].nline;lineI++){
        mp[molI].phot[lineI+iphot*md[molI].nline]=0.;
        expTau[iline]=1.;
        iline++;
      }
    }

    /* Choose random initial photon direction (the distribution used here is even over the surface of a sphere of radius 1).
    */
    pt_theta=gsl_rng_uniform(ran)*2*PI;
    pt_z=2*gsl_rng_uniform(ran)-1;
    semiradius = sqrt(1.-pt_z*pt_z);
    inidir[0]=semiradius*cos(pt_theta);
    inidir[1]=semiradius*sin(pt_theta);
    inidir[2]=pt_z;

    /* Choose the photon frequency/velocity offset.
    */
    segment=gsl_rng_uniform(ran)-0.5;
    /*
    Values of segment should be evenly distributed (considering the
    entire ensemble of photons) between -0.5 and +0.5.
    */

    for (molI=0;molI<par->nSpecies;molI++){
      /* Is factor 4.3=[-2.15,2.15] enough?? */
      deltav[molI]=4.3*segment*gp[id].mol[molI].dopb+dotProduct3D(inidir,gp[id].vel);
      /*
      This is the local (=evaluated at a grid point, not averaged over the local cell) lineshape.
      We store this for later use in ALI loops.
      */
      mp[molI].vfac_loc[iphot]=gaussline(deltav[molI]-dotProduct3D(inidir,gp[id].vel),gp[id].mol[molI].binv);
    }

    here = gp[id].id;

    /* Photon propagation loop */
    while(!gp[here].sink){ /* Testing for sink at loop start is redundant for the first step, since we only start photons from non-sink points, but it makes for simpler code. */
      neighI = getNextEdge(inidir,id,here,gp,ran);

      there=gp[here].neigh[neighI]->id;

      if(firststep){
        firststep=0;
        ds_out=0.5*gp[here].ds[neighI]*dotProduct3D(inidir,gp[here].dir[neighI].xn);
        halfFirstDs[iphot]=ds_out;

        for(molI=0;molI<par->nSpecies;molI++){
          if(!par->doPregrid) {
            calcLineAmpPWLin(gp,here,neighI,molI,deltav[molI],inidir,&vfac_in[molI],&vfac_out[molI]);
         } else
            calcLineAmpLin(gp,here,neighI,molI,deltav[molI],inidir,&vfac_in[molI],&vfac_out[molI]);

          mp[molI].vfac[iphot]=vfac_out[molI];
        }
        /*
        Contribution of the local cell to emission and absorption is done in getjbar.
        We only store the vfac for the local cell for use in ALI loops.
        */
        here=there;
    continue;
      }

      /* If we've got to here, we have progressed beyond the first edge. Length of the new "in" edge is the length of the previous "out".
      */
      ds_in=ds_out;
      ds_out=0.5*gp[here].ds[neighI]*dotProduct3D(inidir,gp[here].dir[neighI].xn);

      for(molI=0;molI<par->nSpecies;molI++){
        vfac_inprev[molI]=vfac_in[molI];
        if(!par->doPregrid)
          calcLineAmpPWLin(gp,here,neighI,molI,deltav[molI],inidir,&vfac_in[molI],&vfac_out[molI]);
        else
          calcLineAmpLin(gp,here,neighI,molI,deltav[molI],inidir,&vfac_in[molI],&vfac_out[molI]);
      }

      nextMolWithBlend = 0;
      iline = 0;
      for(molI=0;molI<par->nSpecies;molI++){
        nextLineWithBlend = 0;
        for(lineI=0;lineI<md[molI].nline;lineI++){
          double jnu_line_in=0., jnu_line_out=0., jnu_cont=0., jnu_blend=0.;
          double alpha_line_in=0., alpha_line_out=0., alpha_cont=0., alpha_blend=0.;

          sourceFunc_line(&md[molI],vfac_inprev[molI],&(gp[here].mol[molI]),lineI,&jnu_line_in,&alpha_line_in);
          sourceFunc_line(&md[molI],vfac_out[molI],&(gp[here].mol[molI]),lineI,&jnu_line_out,&alpha_line_out);
          sourceFunc_cont(gp[here].mol[molI].cont[lineI],&jnu_cont,&alpha_cont);

          /* cont and blend could use the same alpha and jnu counter, but maybe it's clearer this way */

          /* Line blending part.
          */
          if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI\
          && lineI==blends.mols[nextMolWithBlend].lines[nextLineWithBlend].lineI){

            for(bi=0;bi<blends.mols[nextMolWithBlend].lines[nextLineWithBlend].numBlends;bi++){
              molJ  = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].molJ;
              lineJ = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].lineJ;
              velProj = deltav[molI] - blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].deltaV;
	      /*  */
              if(!par->doPregrid)
                calcLineAmpPWLin(gp,here,neighI,molJ,velProj,inidir,&vblend_in,&vblend_out);
              else
                calcLineAmpLin(gp,here,neighI,molJ,velProj,inidir,&vblend_in,&vblend_out);

	      /* we should use also the previous vblend_in, but I don't feel like writing the necessary code now */
              sourceFunc_line(&md[molJ],vblend_out,&(gp[here].mol[molJ]),lineJ,&jnu_blend,&alpha_blend);
              /* note that sourceFunc* increment jnu and alpha, they don't overwrite it  */
            }

            nextLineWithBlend++;
            if(nextLineWithBlend>=blends.mols[nextMolWithBlend].numLinesWithBlends){
              nextLineWithBlend = 0;
              /* The reason for doing this is as follows. Firstly, we only enter the present IF block if molI has at least 1 line which is blended with others; and further, if we have now processed all blended lines for that molecule. Thus no matter what value lineI takes for the present molecule, it won't appear as blends.mols[nextMolWithBlend].lines[i].lineI for any i. Yet we will still test blends.mols[nextMolWithBlend].lines[nextLineWithBlend], thus we want nextLineWithBlend to at least have a sensible value between 0 and blends.mols[nextMolWithBlend].numLinesWithBlends-1. We could set nextLineWithBlend to any number in this range in safety, but zero is simplest. */
            }
          }
          /* End of line blending part */

	  /* as said above, out-in split should be done also for blended lines... */

	  dtau=(alpha_line_out+alpha_cont+alpha_blend)*ds_out;
          if(dtau < -30.) dtau = -30.;
          calcSourceFn(dtau, par, &remnantSnu, &expDTau);
          remnantSnu *= (jnu_line_out+jnu_cont+jnu_blend)*ds_out;
          mp[molI].phot[lineI+iphot*md[molI].nline]+=expTau[iline]*remnantSnu;
	  expTau[iline]*=expDTau;

	  dtau=(alpha_line_in+alpha_cont+alpha_blend)*ds_in;
          if(dtau < -30.) dtau = -30.;
          calcSourceFn(dtau, par, &remnantSnu, &expDTau);
          remnantSnu *= (jnu_line_in+jnu_cont+jnu_blend)*ds_in;
          mp[molI].phot[lineI+iphot*md[molI].nline]+=expTau[iline]*remnantSnu;
	  expTau[iline]*=expDTau;

          if(expTau[iline] > exp(30.)){
            if(!silent) warning("Maser warning: optical depth has dropped below -30");
            expTau[iline]=exp(30.);
          }

          iline++;
        } /* Next line this molecule. */

        if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI)
          nextMolWithBlend++;
      }

      here=there;
    };

    /* Add cmb contribution.
    */
    iline = 0;
    for(molI=0;molI<par->nSpecies;molI++){
      for(lineI=0;lineI<md[molI].nline;lineI++){
        mp[molI].phot[lineI+iphot*md[molI].nline]+=expTau[iline]*md[molI].cmb[lineI];
        iline++;
      }
    }
  }
}

/*....................................................................*/
void
getjbar(int posn, molData *md, struct grid *gp, const int molI\
  , configInfo *par, struct blendInfo blends, int nextMolWithBlend\
  , gridPointData *mp, double *halfFirstDs){

  int lineI,iphot,bi,molJ,lineJ,nextLineWithBlend;
  double dtau,expDTau,remnantSnu,vsum=0.;
  
  for(lineI=0;lineI<md[molI].nline;lineI++) mp[molI].jbar[lineI]=0.;

  for(iphot=0;iphot<gp[posn].nphot;iphot++){
    if(mp[molI].vfac_loc[iphot]>0){
      nextLineWithBlend = 0;
      for(lineI=0;lineI<md[molI].nline;lineI++){
        double jnu=0.0;
        double alpha=0;

        sourceFunc_line(&md[molI],mp[molI].vfac[iphot],&(gp[posn].mol[molI]),lineI,&jnu,&alpha);
        sourceFunc_cont(gp[posn].mol[molI].cont[lineI],&jnu,&alpha);

        /* Line blending part.
        */
        if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI\
        && lineI==blends.mols[nextMolWithBlend].lines[nextLineWithBlend].lineI){
          for(bi=0;bi<blends.mols[nextMolWithBlend].lines[nextLineWithBlend].numBlends;bi++){
            molJ  = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].molJ;
            lineJ = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].lineJ;
            /*
            The next line is not quite correct, because vfac may be different for other molecules due to different values of binv. Unfortunately we don't necessarily have vfac for molJ available yet.
            */
            sourceFunc_line(&md[molJ],mp[molI].vfac[iphot],&(gp[posn].mol[molJ]),lineJ,&jnu,&alpha);
	    /* note that sourceFunc* increment jnu and alpha, they don't overwrite it  */
          }

          nextLineWithBlend++;
          if(nextLineWithBlend>=blends.mols[nextMolWithBlend].numLinesWithBlends){
            nextLineWithBlend = 0;
            /* The reason for doing this is as follows. Firstly, we only enter the present IF block if molI has at least 1 line which is blended with others; and further, if we have now processed all blended lines for that molecule. Thus no matter what value lineI takes for the present molecule, it won't appear as blends.mols[nextMolWithBlend].lines[i].lineI for any i. Yet we will still test blends.mols[nextMolWithBlend].lines[nextLineWithBlend], thus we want nextLineWithBlend to at least have a sensible value between 0 and blends.mols[nextMolWithBlend].numLinesWithBlends-1. We could set nextLineWithBlend to any number in this range in safety, but zero is simplest. */
          }
        }
        /* End of line blending part */

        dtau=alpha*halfFirstDs[iphot];
        calcSourceFn(dtau, par, &remnantSnu, &expDTau);
        remnantSnu *= jnu*halfFirstDs[iphot];
        mp[molI].jbar[lineI]+=mp[molI].vfac_loc[iphot]*(expDTau*mp[molI].phot[lineI+iphot*md[molI].nline]+remnantSnu);

      }
      vsum+=mp[molI].vfac_loc[iphot];
    }
  }
  for(lineI=0;lineI<md[molI].nline;lineI++) mp[molI].jbar[lineI] /= vsum;
}


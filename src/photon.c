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
veloproject(const double dx[3], const double *vel){
  return dx[0]*vel[0]+dx[1]*vel[1]+dx[2]*vel[2];
}

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
getNextEdge_A(double *inidir, const int startGi, const int presentGi, struct grid *gp, const gsl_rng *ran){
  int i,ni,di,niOfLargest,niOfNextLargest;
  double trackCos,dirCos,largest,nextLargest,mytest,dirsFromStart[DIM],lenSquared,norm;
  /*
The idea here is to select for the next grid point, that one which lies closest (with a little randomizing jitter) to the photon track, while requiring the direction of the edge to be in the 'forward' hemisphere of the photon direction.
  */

  /* Calculate dot products between inidir and delta Rs to all the possible next points. Store the largest of these and the next largest.
  */
  i = 0;
  for(ni=0;ni<gp[presentGi].numNeigh;ni++){
    dirCos=( inidir[0]*gp[presentGi].dir[ni].xn[0]
            +inidir[1]*gp[presentGi].dir[ni].xn[1]
            +inidir[2]*gp[presentGi].dir[ni].xn[2]);

    if(dirCos<=0.0)
  continue; /* because the edge points in the backward direction. */

    lenSquared = 0.0;
    for(di=0;di<DIM;di++){
      dirsFromStart[di] = gp[presentGi].neigh[ni]->x[di] - gp[startGi].x[di];
      lenSquared += dirsFromStart[di]*dirsFromStart[di];
    }

    if(lenSquared<=0.0)
  continue;

    norm = 1.0/sqrt(lenSquared);

    trackCos = 0.0;
    for(di=0;di<DIM;di++)
      trackCos += inidir[di]*dirsFromStart[di]*norm;

    if(i==0){
      largest = trackCos;
      niOfLargest = ni;
    }else{
      if(trackCos>largest){
        nextLargest = largest;
        niOfNextLargest = niOfLargest;
        largest = trackCos;
        niOfLargest = ni;
      }else if(i==1 || trackCos>nextLargest){
        nextLargest = trackCos;
        niOfNextLargest = ni;
      }
    }

    i++;
  }

  /* Choose the edge to follow.
  */
  if(i>1){ /* then nextLargest, niOfNextLargest should exist. */
    mytest = (1.0 + nextLargest)/(2.0 + nextLargest + largest);
    /* The addition of the scalars here is I think essentially arbitrary - they just serve to make the choices a bit more even, which tends to scatter the photon a bit more. */
    if(gsl_rng_uniform(ran)<mytest)
      return niOfNextLargest;
    else
      return niOfLargest;
  }else if(i>0){
    return niOfLargest;
  }else{
    if(!silent)
      bail_out("Photon propagation error - no valid edges.");
    exit(1);
  }
}

/*....................................................................*/
int
getNextEdge_B(double *inidir, const int presentGi, struct grid *gp\
  , const gsl_rng *ran, double (*deltaRHats)[DIM], double *deltaRInvMags\
  , const double lastDeltaRInvMag){
  /*
The idea here is to select for the next grid point, that one which lies closest (with a little randomizing jitter) to the photon track, while requiring it to be further from the start point than the 'presentGi'.
  */

  int i,ni,di,niOfLargest,niOfNextLargest,testGi;
  double trackCos,dirCos,largest,nextLargest,mytest;

  /* Calculate dot products between inidir and delta Rs to all the possible next points. Store the largest of these and the next largest.
  */
  i = 0;
  for(ni=0;ni<gp[presentGi].numNeigh;ni++){
    testGi = gp[presentGi].neigh[ni]->id;

    if(lastDeltaRInvMag>0.0){
      if(deltaRInvMags[testGi]>lastDeltaRInvMag)
  continue;

    }else{ /* still at start */
      dirCos=( inidir[0]*gp[presentGi].dir[ni].xn[0]
              +inidir[1]*gp[presentGi].dir[ni].xn[1]
              +inidir[2]*gp[presentGi].dir[ni].xn[2]);

      if(dirCos<=0.0)
  continue; /* because this edge points in the backward direction. */
    }

    trackCos = 0.0;
    for(di=0;di<DIM;di++)
      trackCos += inidir[di]*deltaRHats[testGi][di];

    if(i==0){
      largest = trackCos;
      niOfLargest = ni;
    }else{
      if(trackCos>largest){
        nextLargest = largest;
        niOfNextLargest = niOfLargest;
        largest = trackCos;
        niOfLargest = ni;
      }else if(i==1 || trackCos>nextLargest){
        nextLargest = trackCos;
        niOfNextLargest = ni;
      }
    }

    i++;
  }

  /* Choose the edge to follow.
  */
  if(i>1){ /* then nextLargest, niOfNextLargest should exist. */
    mytest = (1.0 + nextLargest)/(2.0 + nextLargest + largest);
    /* The addition of the scalars here is I think essentially arbitrary - they just serve to make the choices a bit more even, which tends to scatter the photon a bit more. */
    if(gsl_rng_uniform(ran)<mytest)
      return niOfNextLargest;
    else
      return niOfLargest;
  }else if(i>0){
    return niOfLargest;
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
  v[0]=deltav-veloproject(inidir,g[id].vel);
  v[1]=deltav-veloproject(inidir,&(g[id].v1[3*k]));
  v[2]=deltav-veloproject(inidir,&(g[id].v2[3*k]));
  v[3]=deltav-veloproject(inidir,&(g[id].v3[3*k]));
  v[4]=deltav-veloproject(inidir,g[id].neigh[k]->vel);

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
  v[0]=deltav-veloproject(inidir,g[id].vel);
  v[2]=deltav-veloproject(inidir,g[id].neigh[k]->vel);
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
photon(int id, struct grid *gp, molData *md, int iter, const gsl_rng *ran\
  , configInfo *par, const int nlinetot, struct blendInfo blends\
  , gridPointData *mp, double *halfFirstDs){

  int iphot,iline,here,there,firststep,neighI,np_per_line,ip_at_line,j,di;
  int nextMolWithBlend, nextLineWithBlend, molI, lineI, molJ, lineJ, bi;
  double segment,vblend_in,vblend_out,dtau,expDTau,ds_in=0.0,ds_out=0.0,pt_theta,pt_z,semiradius;
  double deltav[par->nSpecies],vfac_in[par->nSpecies],vfac_out[par->nSpecies],vfac_inprev[par->nSpecies];
  double expTau[nlinetot],inidir[3];
  double remnantSnu,velProj;
  double (*deltaRHats)[DIM],*deltaRInvMags,lenSquared,lastDeltaRInvMag;

_Bool doPreCalc=0;

  np_per_line=(int) gp[id].nphot/gp[id].numNeigh; // Works out to be equal to ininphot. :-/

  if(doPreCalc){
    deltaRHats    = malloc(sizeof(*deltaRHats)   *par->ncell);
    deltaRInvMags = malloc(sizeof(*deltaRInvMags)*par->ncell);

    /* Pre-calculate displacement vectors (also 1/length) for all the grid points:
    */
    for(j=0;j<par->ncell;j++){
      lenSquared = 0.0;
      for(di=0;di<DIM;di++){
        deltaRHats[j][di] = gp[j].x[di] - gp[id].x[di];
        lenSquared += deltaRHats[j][di]*deltaRHats[j][di];
      }
      if(lenSquared>0.0)
        deltaRInvMags[j] = 1.0/sqrt(lenSquared);
      else
        deltaRInvMags[j] = 0.0;

      for(di=0;di<DIM;di++)
        deltaRHats[j][di] *= deltaRInvMags[j];
    }
  }

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
    iter=(int) (gsl_rng_uniform(ran)*(double)N_RAN_PER_SEGMENT); /* can have values in [0,1,..,N_RAN_PER_SEGMENT-1]*/
    ip_at_line=(int) iphot/gp[id].numNeigh;
    segment=(N_RAN_PER_SEGMENT*(ip_at_line-np_per_line*0.5)+iter)/(double)(np_per_line*N_RAN_PER_SEGMENT);
    /*
    Values of segment should be evenly distributed (considering the
    entire ensemble of photons) between -0.5 and +0.5, and are chosen
    from a sequence of possible values separated by
    1/(N_RAN_PER_SEGMENT*ininphot).

    */

    for (molI=0;molI<par->nSpecies;molI++){
      /* Is factor 4.3=[-2.15,2.15] enough?? */
      deltav[molI]=4.3*segment*gp[id].mol[molI].dopb+veloproject(inidir,gp[id].vel);
      /*
      This is the local (=evaluated at a grid point, not averaged over the local cell) lineshape.
      We store this for later use in ALI loops.
      */
      mp[molI].vfac_loc[iphot]=gaussline(deltav[molI]-veloproject(inidir,gp[id].vel),gp[id].mol[molI].binv);
    }

    here = gp[id].id;
    if(doPreCalc)
      lastDeltaRInvMag = 0.0;

    /* Photon propagation loop */
    while(!gp[here].sink){ /* Testing for sink at loop start is redundant for the first step, since we only start photons from non-sink points, but it makes for simpler code. */
      if(doPreCalc){
        neighI=getNextEdge_B(inidir,here,gp,ran,deltaRHats,deltaRInvMags,lastDeltaRInvMag);
      }else{
        neighI=getNextEdge_A(inidir,id,here,gp,ran);
      }

      there=gp[here].neigh[neighI]->id;
      if(doPreCalc)
        lastDeltaRInvMag = deltaRInvMags[there];

      if(firststep){
        firststep=0;
        ds_out=0.5*gp[here].ds[neighI]*veloproject(inidir,gp[here].dir[neighI].xn);
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
      ds_out=0.5*gp[here].ds[neighI]*veloproject(inidir,gp[here].dir[neighI].xn);

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

  if(doPreCalc){
    free(deltaRInvMags);
    free(deltaRHats);
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


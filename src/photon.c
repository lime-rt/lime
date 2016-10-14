/*
 *  photon.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
 *
 */


#include "lime.h"

double veloproject(const double dx[3], const double *vel){
  return dx[0]*vel[0]+dx[1]*vel[1]+dx[2]*vel[2];
}

double geterf(const double x0, const double x1) {
  /* table lookup erf thingy */

  double val0=0.,val1=0.;
  
  if (fabs(x0)>=6.) val0=(SPI/2.);
  else {
    int index = (int)(fabs(x0*1024.));
    double inter_coeff = (fabs(1024.*x0)-index);
    val0=(1-inter_coeff)*ERF_TABLE[index]+inter_coeff*ERF_TABLE[index+1];
  }
  if (x0<0.) val0=-val0;
 
  if (fabs(x1)>=6.) val1=(SPI/2.);
  else {
    int index = (int)(fabs(x1*1024.));
    double inter_coeff = (fabs(1024.*x1)-index);
    val1=(1-inter_coeff)*ERF_TABLE[index]+inter_coeff*ERF_TABLE[index+1];
  }
  if (x1<0.) val1=-val1;

  return fabs((val1-val0)/(x1-x0));
}

double gaussline(const double v, const double oneOnSigma){
  double val;
  val = v*v*oneOnSigma*oneOnSigma;
#ifdef FASTEXP
  return FastExp(val);
#else
  return exp(-val);
#endif
}


int getNextEdge(double *inidir, int id, struct grid *gp, const gsl_rng *ran){
  int i,iOfLargest,iOfNextLargest,numPositive;
  double cosAngle,largest,nextLargest,mytest;

  /* Calculate dot products between inidir and all the edges. Store the largest of these and the next largest.
  */
  numPositive = 0;
  for(i=0;i<gp[id].numNeigh;i++){
    cosAngle=( inidir[0]*gp[id].dir[i].xn[0]
              +inidir[1]*gp[id].dir[i].xn[1]
              +inidir[2]*gp[id].dir[i].xn[2]);

    if(cosAngle>0.0)
      numPositive++;

    if(i==0){
      largest = cosAngle;
      iOfLargest = i;
    }else if(i==1){
      if(cosAngle>largest){
        nextLargest = largest;
        iOfNextLargest = iOfLargest;
        largest = cosAngle;
        iOfLargest = i;
      }else{
        nextLargest = cosAngle;
        iOfNextLargest = i;
      }
    }else{
      if(cosAngle>largest){
        nextLargest = largest;
        iOfNextLargest = iOfLargest;
        largest = cosAngle;
        iOfLargest = i;
      }else if(cosAngle>nextLargest){
        nextLargest = cosAngle;
        iOfNextLargest = i;
      }
    }
  }

  if(!silent && numPositive<=0)
    warning("Photon propagation error - there are no forward-going edges.");

  /* Choose the edge to follow.
  */
  mytest = (1.0 + nextLargest)/(2.0 + nextLargest + largest);
  /* The addition of the scalars here is I think essentially arbitrary - they just serve to make the choices a bit more even, which tends to scatter the photon a bit more. */
  if(gsl_rng_uniform(ran)<mytest)
    return iOfNextLargest;
  else
    return iOfLargest;
}


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
  vfac_out is the the lineshape for the part of the edge in the current Voronoi cell,
  vfac_in is for the part in the next cell
  */

  if (fabs(v[1]-v[0])*binv_this>1./512.) {
     *vfac_out=0.5*geterf(v[0]*binv_this,v[1]*binv_this);
  } else *vfac_out=0.5*gaussline(0.5*(v[0]+v[1]),binv_this);
  if (fabs(v[2]-v[1])*binv_this>1./512.) {
    *vfac_out+=0.5*geterf(v[1]*binv_this,v[2]*binv_this);
  } else *vfac_out+=0.5*gaussline(0.5*(v[1]+v[2]),binv_this);

  if (fabs(v[3]-v[2])*binv_next>1./512.) {
     *vfac_in=0.5*geterf(v[2]*binv_next,v[3]*binv_next);
  } else *vfac_in=0.5*gaussline(0.5*(v[2]+v[3]),binv_next);
  if (fabs(v[4]-v[3])*binv_next>1./512.) {
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

  if (fabs(v[1]-v[0])*binv_this>1./512.) {
     *vfac_out=geterf(v[0]*binv_this,v[1]*binv_this);
  } else *vfac_out+=gaussline(0.5*(v[0]+v[1]),binv_this);

  if (fabs(v[2]-v[1])*binv_next>1./512.) {
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
void photon(int id, struct grid *gp, molData *md, int iter, const gsl_rng *ran\
  , configInfo *par, const int nlinetot, struct blendInfo blends\
  , gridPointData *mp, double *halfFirstDs){

  int iphot,iline,here,there,firststep,neighI,np_per_line,ip_at_line;
  int nextMolWithBlend, nextLineWithBlend, molI, lineI, molJ, lineJ, bi;
  double segment,vblend_in,vblend_out,dtau,expDTau,jnu,alpha,ds_in=0.0,ds_out=0.0,pt_theta,pt_z,semiradius;
  double deltav[par->nSpecies],vfac_in[par->nSpecies],vfac_out[par->nSpecies],vfac_inprev[par->nSpecies];
  double *tau,*expTau,inidir[3];
  double remnantSnu, velProj;

  tau    = malloc(sizeof(*tau)   *nlinetot);
  expTau = malloc(sizeof(*expTau)*nlinetot);  
  np_per_line=(int) gp[id].nphot/gp[id].numNeigh; // Works out to be equal to ininphot. :-/

  for(iphot=0;iphot<gp[id].nphot;iphot++){
    firststep=1;
    iline = 0;
    for(molI=0;molI<par->nSpecies;molI++){
      for(lineI=0;lineI<md[molI].nline;lineI++){
        mp[molI].phot[lineI+iphot*md[molI].nline]=0.;
        tau[iline]=0.;
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

    here=gp[id].id;

    /* Photon propagation loop */
    do{
      neighI=getNextEdge(inidir,here,gp,ran);
      there=gp[here].neigh[neighI]->id;

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
        continue;
      } else {
        /* length of the new "in" edge is the length of the previous "out"*/
        ds_in=ds_out;
        ds_out=0.5*gp[here].ds[neighI]*veloproject(inidir,gp[here].dir[neighI].xn);

        for(molI=0;molI<par->nSpecies;molI++){
          vfac_inprev[molI]=vfac_in[molI];
          if(!par->doPregrid)
            calcLineAmpPWLin(gp,here,neighI,molI,deltav[molI],inidir,&vfac_in[molI],&vfac_out[molI]);
          else
            calcLineAmpLin(gp,here,neighI,molI,deltav[molI],inidir,&vfac_in[molI],&vfac_out[molI]);
        }
      }

      nextMolWithBlend = 0;
      iline = 0;
      for(molI=0;molI<par->nSpecies;molI++){
        nextLineWithBlend = 0;
        for(lineI=0;lineI<md[molI].nline;lineI++){
          double jnu_in=0., jnu_out=0.;
          double alpha_in=0., alpha_out=0.;
          double dtau_all_lines, expDTau_all_lines;
          alpha=0.;
          jnu=0.;

          sourceFunc_line(md[molI],vfac_inprev[molI],gp[here].mol[molI],lineI,&jnu_in,&alpha_in);
          sourceFunc_line(md[molI],vfac_out[molI],gp[here].mol[molI],lineI,&jnu_out,&alpha_out);
          sourceFunc_cont(gp[here].mol[molI].cont[lineI],&jnu,&alpha);

          dtau=alpha_in*ds_in+alpha_out*ds_out+alpha*(ds_in+ds_out);

          if(dtau < -30) dtau = -30;
          calcSourceFn(dtau, par, &remnantSnu, &expDTau);
          remnantSnu *= jnu*(ds_in+ds_out)+jnu_in*ds_in+jnu_out*ds_out;
          mp[molI].phot[lineI+iphot*md[molI].nline]+=expTau[iline]*remnantSnu;

          dtau_all_lines=dtau;
          expDTau_all_lines=expDTau;

          if(tau[iline] < -30.){
            if(!silent) warning("Maser warning: optical depth has dropped below -30");
            tau[iline]= -30.;
            expTau[iline]=exp(-tau[iline]);
          }

          /* Line blending part.
          */
          if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI\
          && lineI==blends.mols[nextMolWithBlend].lines[nextLineWithBlend].lineI){
            jnu=0.;
            alpha=0.;
            for(bi=0;bi<blends.mols[nextMolWithBlend].lines[nextLineWithBlend].numBlends;bi++){
              molJ  = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].molJ;
              lineJ = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].lineJ;
              velProj = deltav[molI] - blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].deltaV;
	      /*  */
              if(!par->doPregrid)
                calcLineAmpPWLin(gp,here,neighI,molJ,velProj,inidir,&vblend_in,&vblend_out);
              else
                calcLineAmpLin(gp,here,neighI,molJ,velProj,inidir,&vblend_in,&vblend_out);

	      /* we should really use the previous vblend_in, but I don't feel like writing the necessary code now... */
              sourceFunc_line(md[molJ],vblend_out,gp[here].mol[molJ],lineJ,&jnu,&alpha);
              dtau=alpha*(ds_in+ds_out);
              if(dtau < -30) dtau = -30;
              calcSourceFn(dtau, par, &remnantSnu, &expDTau);
              remnantSnu *= jnu*(ds_in+ds_out);

              mp[molI].phot[lineI+iphot*md[molI].nline]+=expTau[iline]*remnantSnu;

              dtau_all_lines+=dtau;
              expDTau_all_lines*=expDTau;

              if(tau[iline] < -30.){
                if(!silent) warning("Optical depth has dropped below -30");
                tau[iline]= -30.;
                expTau[iline]=exp(-tau[iline]);
              }
            }

            nextLineWithBlend++;
            if(nextLineWithBlend>=blends.mols[nextMolWithBlend].numLinesWithBlends){
              nextLineWithBlend = 0;
              /* The reason for doing this is as follows. Firstly, we only enter the present IF block if molI has at least 1 line which is blended with others; and further, if we have now processed all blended lines for that molecule. Thus no matter what value lineI takes for the present molecule, it won't appear as blends.mols[nextMolWithBlend].lines[i].lineI for any i. Yet we will still test blends.mols[nextMolWithBlend].lines[nextLineWithBlend], thus we want nextLineWithBlend to at least have a sensible value between 0 and blends.mols[nextMolWithBlend].numLinesWithBlends-1. We could set nextLineWithBlend to any number in this range in safety, but zero is simplest. */
            }
          }
          /* End of line blending part */
          tau[iline]+=dtau_all_lines;
          expTau[iline]*=expDTau_all_lines;

          iline++;
        } /* Next line this molecule. */

        if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI)
          nextMolWithBlend++;
      }

      here=there;
    } while(!gp[here].sink);
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
  free(expTau);
  free(tau);
}

/*....................................................................*/
void getjbar(int posn, molData *md, struct grid *gp, const int molI\
  , configInfo *par, struct blendInfo blends, int nextMolWithBlend\
  , gridPointData *mp, double *halfFirstDs){

  int lineI,iphot,bi,molJ,lineJ,nextLineWithBlend;
  double tau,expTau,remnantSnu,vsum=0.,jnu,alpha;
  
  for(lineI=0;lineI<md[molI].nline;lineI++) mp[molI].jbar[lineI]=0.;

  for(iphot=0;iphot<gp[posn].nphot;iphot++){
    if(mp[molI].vfac[iphot]>0){
      nextLineWithBlend = 0;
      for(lineI=0;lineI<md[molI].nline;lineI++){
        jnu=0.;
        alpha=0.;

        sourceFunc_line(md[molI],mp[molI].vfac[iphot],gp[posn].mol[molI],lineI,&jnu,&alpha);
        sourceFunc_cont(gp[posn].mol[molI].cont[lineI],&jnu,&alpha);

        tau=alpha*halfFirstDs[iphot];
        calcSourceFn(tau, par, &remnantSnu, &expTau);
        remnantSnu *= jnu*halfFirstDs[iphot];

        mp[molI].jbar[lineI]+=mp[molI].vfac_loc[iphot]*(expTau*mp[molI].phot[lineI+iphot*md[molI].nline]+remnantSnu);

        /* Line blending part.
        */
        if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI\
        && lineI==blends.mols[nextMolWithBlend].lines[nextLineWithBlend].lineI){
          jnu=0.;
          alpha=0.;
          for(bi=0;bi<blends.mols[nextMolWithBlend].lines[nextLineWithBlend].numBlends;bi++){
            molJ  = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].molJ;
            lineJ = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].lineJ;
            /*
            The next line is not quite correct, because vfac may be different for other molecules due to different values of binv. Unfortunately we don't necessarily have vfac for molJ available yet.
            */
            sourceFunc_line(md[molJ],mp[molI].vfac[iphot],gp[posn].mol[molJ],lineJ,&jnu,&alpha);
            tau=alpha*halfFirstDs[iphot];
            calcSourceFn(tau, par, &remnantSnu, &expTau);
            remnantSnu *= jnu*halfFirstDs[iphot];

            mp[molI].jbar[lineI]+=mp[molI].vfac_loc[iphot]*(expTau*mp[molI].phot[lineI+iphot*md[molI].nline]+remnantSnu);
          }

          nextLineWithBlend++;
          if(nextLineWithBlend>=blends.mols[nextMolWithBlend].numLinesWithBlends){
            nextLineWithBlend = 0;
            /* The reason for doing this is as follows. Firstly, we only enter the present IF block if molI has at least 1 line which is blended with others; and further, if we have now processed all blended lines for that molecule. Thus no matter what value lineI takes for the present molecule, it won't appear as blends.mols[nextMolWithBlend].lines[i].lineI for any i. Yet we will still test blends.mols[nextMolWithBlend].lines[nextLineWithBlend], thus we want nextLineWithBlend to at least have a sensible value between 0 and blends.mols[nextMolWithBlend].numLinesWithBlends-1. We could set nextLineWithBlend to any number in this range in safety, but zero is simplest. */
          }
        }
        /* End of line blending part */
      }
      vsum+=mp[molI].vfac_loc[iphot];
    }
  }
  for(lineI=0;lineI<md[molI].nline;lineI++) mp[molI].jbar[lineI] /= vsum;
}


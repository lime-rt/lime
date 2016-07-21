/*
 *  photon.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2016 The LIME development team
 *
 */

#include "lime.h"

int getNextEdge(double *inidir, int id, struct grid *g, const gsl_rng *ran){
  int i,iOfLargest,iOfNextLargest,numPositive;
  double cosAngle,largest,nextLargest,mytest;

  /* Calculate dot products between inidir and all the edges. Store the largest of these and the next largest.
  */
  numPositive = 0;
  for(i=0;i<g[id].numNeigh;i++){
    cosAngle=( inidir[0]*g[id].dir[i].xn[0]
              +inidir[1]*g[id].dir[i].xn[1]
              +inidir[2]*g[id].dir[i].xn[2]);

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

/*....................................................................*/
void calcLineAmpSpline(struct grid *g, const int id, const int k\
  , const double binv, const double deltav, double *vfac){

  int nspline,ispline,naver,iaver;
  double v1,v2,s1,s2,sd,v,vfacsub,d;
  
  v1=deltav-veloproject(g[id].dir[k].xn,g[id].vel);
  v2=deltav-veloproject(g[id].dir[k].xn,g[id].neigh[k]->vel);

  nspline=(fabs(v1-v2)*binv < 1) ? 1 : (int)(fabs(v1-v2)*binv);
  *vfac=0.;
  s2=0;
  v2=v1;
  
  for(ispline=0;ispline<nspline;ispline++){
    s1=s2;
    s2=((double)(ispline+1))/(double)nspline;
    v1=v2;
    d=s2*g[id].ds[k];
    v2=deltav-((((g[id].a4[k]*d+g[id].a3[k])*d+g[id].a2[k])*d+g[id].a1[k])*d+g[id].a0[k]);
    naver=(1 > fabs(v1-v2)*binv) ? 1 : (int)(fabs(v1-v2)*binv);
    for(iaver=0;iaver<naver;iaver++){
      sd=s1+(s2-s1)*((double)iaver+0.5)/(double)naver;
      d=sd*g[id].ds[k];
      v=deltav-((((g[id].a4[k]*d+g[id].a3[k])*d+g[id].a2[k])*d+g[id].a1[k])*d+g[id].a0[k]);
      vfacsub=gaussline(v,binv);
      *vfac+=vfacsub/(double)naver;
    }
  }
  *vfac= *vfac/(double)nspline;
  return;
}

/*....................................................................*/
void calcLineAmpLinear(struct grid *g, const int id, const int k\
  , const double binv, const double deltav, double *vfac){

  int nspline,ispline,naver,iaver;
  double v1,v2,s1,s2,sd,v,vfacsub,d;
  
  v1=deltav-veloproject(g[id].dir[k].xn,g[id].vel);
  v2=deltav-veloproject(g[id].dir[k].xn,g[id].neigh[k]->vel);

  nspline=(fabs(v1-v2)*binv < 1) ? 1 : (int)(fabs(v1-v2)*binv);
  *vfac=0.;
  s2=0;
  v2=v1;
  
  for(ispline=0;ispline<nspline;ispline++){
    s1=s2;
    s2=((double)(ispline+1))/(double)nspline;
    v1=v2;
    d=s2*g[id].ds[k];
    v2=deltav-(g[id].a1[k]*d+g[id].a0[k]);
    naver=(1 > fabs(v1-v2)*binv) ? 1 : (int)(fabs(v1-v2)*binv);
    for(iaver=0;iaver<naver;iaver++){
      sd=s1+(s2-s1)*((double)iaver+0.5)/(double)naver;
      d=sd*g[id].ds[k];
      v=deltav-(g[id].a1[k]*d+g[id].a0[k]);
      vfacsub=gaussline(v,binv);
      *vfac+=vfacsub/(double)naver;
    }
  }
  *vfac= *vfac/(double)nspline;
  return;
}

double veloproject(double dx[3], double *vel){
  return dx[0]*vel[0]+dx[1]*vel[1]+dx[2]*vel[2];
}


double gaussline(double v, double oneOnSigma){
  double val;

  val = v*v*oneOnSigma*oneOnSigma;
#ifdef FASTEXP
  return FastExp(val);
#else
  return exp(-val);
#endif
}


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
    *remnantSnu = 1. - dTau*(1. - dTau/3.)/2.;
  } else {
    *remnantSnu = (1.-(*expDTau))/dTau;
  }
#else
  if (fabs(dTau)<par->taylorCutoff){
    *remnantSnu = 1. - dTau*(1. - dTau/3.)/2.;
    *expDTau = 1. - dTau*(*remnantSnu);
  } else {
    *expDTau = exp(-dTau);
    *remnantSnu = (1.-(*expDTau))/dTau;
  }
#endif
}

void
photon(int id, struct grid *g, molData *m, int iter, const gsl_rng *ran\
  , configInfo *par, const int nlinetot, struct blendInfo blends\
  , gridPointData *mp, double *halfFirstDs){

  int iphot,iline,here,there,firststep,neighI,np_per_line,ip_at_line;
  int nextMolWithBlend, nextLineWithBlend, molI, lineI, molJ, lineJ, bi;
  double deltav,segment,vblend,dtau,expDTau,jnu,alpha,ds,vfac[par->nSpecies],pt_theta,pt_z,semiradius;
  double *tau,*expTau,inidir[3];
  double remnantSnu, velProj;

  tau    = malloc(sizeof(*tau)   *nlinetot);
  expTau = malloc(sizeof(*expTau)*nlinetot);
  
  np_per_line=(int) g[id].nphot/g[id].numNeigh; // Works out to be equal to ininphot. :-/

  for(iphot=0;iphot<g[id].nphot;iphot++){
    firststep=1;
    iline = 0;
    for(molI=0;molI<par->nSpecies;molI++){
      for(lineI=0;lineI<m[molI].nline;lineI++){
        mp[molI].phot[lineI+iphot*m[molI].nline]=0.;
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
    ip_at_line=(int) iphot/g[id].numNeigh;
    segment=(N_RAN_PER_SEGMENT*(ip_at_line-np_per_line*0.5)+iter)/(double)(np_per_line*N_RAN_PER_SEGMENT);
    /*
    Values of segment should be evenly distributed (considering the
    entire ensemble of photons) between -0.5 and +0.5, and are chosen
    from a sequence of possible values separated by
    1/(N_RAN_PER_SEGMENT*ininphot).
    */

    here=g[id].id;
    deltav=segment*4.3*g[id].dopb+veloproject(inidir,g[id].vel);

    /* Photon propagation loop */
    do{
      neighI=getNextEdge(inidir,here,g,ran);
      there=g[here].neigh[neighI]->id;

      if(firststep){
        firststep=0;				
        ds=g[here].ds[neighI]*0.5;
        halfFirstDs[iphot]=ds;

        for(molI=0;molI<par->nSpecies;molI++){
          if(!par->doPregrid)
            calcLineAmpSpline(g,here,neighI,g[id].mol[molI].binv,deltav,&vfac[molI]);
          else
            calcLineAmpLinear(g,here,neighI,g[id].mol[molI].binv,deltav,&vfac[molI]);
          mp[molI].vfac[iphot]=vfac[molI];
        }
      } else {
        ds=g[here].ds[neighI];
      
        for(molI=0;molI<par->nSpecies;molI++){
          if(!par->doPregrid)
            calcLineAmpSpline(g,here,neighI,g[id].mol[molI].binv,deltav,&vfac[molI]);
          else
            calcLineAmpLinear(g,here,neighI,g[id].mol[molI].binv,deltav,&vfac[molI]);
        }
      }

      nextMolWithBlend = 0;
      iline = 0;
      for(molI=0;molI<par->nSpecies;molI++){
        nextLineWithBlend = 0;
        for(lineI=0;lineI<m[molI].nline;lineI++){
          jnu=0.;
          alpha=0.;

          sourceFunc_line(m[molI],vfac[molI],g[here].mol[molI],lineI,&jnu,&alpha);
          sourceFunc_cont(g[here].mol[molI],lineI,&jnu,&alpha);

          dtau=alpha*ds;
          if(dtau < -30) dtau = -30;
          calcSourceFn(dtau, par, &remnantSnu, &expDTau);
          remnantSnu *= jnu*m[molI].norminv*ds;

          mp[molI].phot[lineI+iphot*m[molI].nline]+=expTau[iline]*remnantSnu;
          tau[iline]+=dtau;
          expTau[iline]*=expDTau;
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
              velProj = deltav - blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].deltaV;

              if(!par->doPregrid)
                calcLineAmpSpline(g,here,neighI,g[id].mol[molJ].binv,velProj,&vblend);
              else
                calcLineAmpLinear(g,here,neighI,g[id].mol[molJ].binv,velProj,&vblend);

              sourceFunc_line(m[molJ],vblend,g[here].mol[molJ],lineJ,&jnu,&alpha);
              dtau=alpha*ds;
              if(dtau < -30) dtau = -30;
              calcSourceFn(dtau, par, &remnantSnu, &expDTau);
              remnantSnu *= jnu*m[molJ].norminv*ds;

              mp[molI].phot[lineI+iphot*m[molI].nline]+=expTau[iline]*remnantSnu;
              tau[iline]+=dtau;
              expTau[iline]*=expDTau;
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

          iline++;
        }

        if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI)
          nextMolWithBlend++;
      }
      
      here=there;
    } while(!g[here].sink);
    
    /* Add cmb contribution.
    */
    iline = 0;
    for(molI=0;molI<par->nSpecies;molI++){
      for(lineI=0;lineI<m[molI].nline;lineI++){
        mp[molI].phot[lineI+iphot*m[molI].nline]+=expTau[iline]*m[molI].cmb[lineI];
        iline++;
      }
    }
  }
  free(expTau);
  free(tau);
}

void
getjbar(int posn, molData *m, struct grid *g, const int molI\
  , configInfo *par, struct blendInfo blends, int nextMolWithBlend\
  , gridPointData *mp, double *halfFirstDs){

  int lineI,iphot,bi,molJ,lineJ,nextLineWithBlend;
  double tau, expTau, remnantSnu, vsum=0., jnu, alpha;
  
  for(lineI=0;lineI<m[molI].nline;lineI++) mp[molI].jbar[lineI]=0.;

  for(iphot=0;iphot<g[posn].nphot;iphot++){
    if(mp[molI].vfac[iphot]>0){
      nextLineWithBlend = 0;
      for(lineI=0;lineI<m[molI].nline;lineI++){
        jnu=0.;
        alpha=0.;

        sourceFunc_line(m[molI],mp[molI].vfac[iphot],g[posn].mol[molI],lineI,&jnu,&alpha);
        sourceFunc_cont(g[posn].mol[molI],lineI,&jnu,&alpha);
        tau=alpha*halfFirstDs[iphot];
        calcSourceFn(tau, par, &remnantSnu, &expTau);
        remnantSnu *= jnu*m[molI].norminv*halfFirstDs[iphot];

        mp[molI].jbar[lineI]+=mp[molI].vfac[iphot]*(expTau*mp[molI].phot[lineI+iphot*m[molI].nline]+remnantSnu);

        /* Line blending part.
        */
        if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI\
        && lineI==blends.mols[nextMolWithBlend].lines[nextLineWithBlend].lineI){
          jnu=0.;
          alpha=0.;
          for(bi=0;bi<blends.mols[nextMolWithBlend].lines[nextLineWithBlend].numBlends;bi++){
            molJ  = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].molJ;
            lineJ = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].lineJ;

            sourceFunc_line(m[molJ],mp[molI].vfac[iphot],g[posn].mol[molJ],lineJ,&jnu,&alpha);
            tau=alpha*halfFirstDs[iphot];
            calcSourceFn(tau, par, &remnantSnu, &expTau);
            remnantSnu *= jnu*m[molJ].norminv*halfFirstDs[iphot];

            mp[molI].jbar[lineI]+=mp[molI].vfac[iphot]*(expTau*mp[molI].phot[lineI+iphot*m[molI].nline]+remnantSnu);
          }

          nextLineWithBlend++;
          if(nextLineWithBlend>=blends.mols[nextMolWithBlend].numLinesWithBlends){
            nextLineWithBlend = 0;
            /* The reason for doing this is as follows. Firstly, we only enter the present IF block if molI has at least 1 line which is blended with others; and further, if we have now processed all blended lines for that molecule. Thus no matter what value lineI takes for the present molecule, it won't appear as blends.mols[nextMolWithBlend].lines[i].lineI for any i. Yet we will still test blends.mols[nextMolWithBlend].lines[nextLineWithBlend], thus we want nextLineWithBlend to at least have a sensible value between 0 and blends.mols[nextMolWithBlend].numLinesWithBlends-1. We could set nextLineWithBlend to any number in this range in safety, but zero is simplest. */
          }
        }
        /* End of line blending part */
      }
      vsum+=mp[molI].vfac[iphot];
    }
  }
  for(lineI=0;lineI<m[molI].nline;lineI++) mp[molI].jbar[lineI] *= m[molI].norm/vsum;
}


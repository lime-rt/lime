/*
 *  photon.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 15/11/06.
 *  Copyright 2006-2014, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *  All rights reserved.
 *
 */

#include "lime.h"

int
sortangles(double *inidir, int id, struct grid *g, const gsl_rng *ran) {
  int i,n[2];
  double angle,exitdir[2];

  exitdir[0]=1e30;
  exitdir[1]=1e31;
  n[0]=-1;
  n[1]=-1;
  for(i=0;i<g[id].numNeigh;i++){
    angle=( inidir[0]*g[id].dir[i].xn[0]
           +inidir[1]*g[id].dir[i].xn[1]
           +inidir[2]*g[id].dir[i].xn[2]);
    if(angle<exitdir[0]){
      exitdir[1]=exitdir[0];
      n[1]=n[0];
      exitdir[0]=angle;
      n[0]=i;
    } else if(angle<exitdir[1]) {
      exitdir[1]=angle;
      n[1]=i;
    }
  }
  if(gsl_rng_uniform(ran)<1./((1-exitdir[0])/(1-exitdir[1])+1) ) {
    if(n[0]==-1){
      if(!silent) bail_out("Photon propagation error");
      exit(1);
    }
    return n[0];
  } else {
    if(n[1]==-1){
      if(!silent) bail_out("Photon propagation error");
      exit(1);
    }
    return n[1];
  }
}



void
velocityspline(struct grid *g, int id, int k, double binv, double deltav, double *vfac){
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
      sd=s1+(s2-s1)*((double)iaver-0.5)/(double)naver;
      d=sd*g[id].ds[k];
      v=deltav-((((g[id].a4[k]*d+g[id].a3[k])*d+g[id].a2[k])*d+g[id].a1[k])*d+g[id].a0[k]);
      vfacsub=gaussline(v,binv);
      *vfac+=vfacsub/(double)naver;
    }
  }
  *vfac= *vfac/(double)nspline;
  return;
}


void
velocityspline_lin(struct grid *g, int id, int k, double binv, double deltav, double *vfac){
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
      sd=s1+(s2-s1)*((double)iaver-0.5)/(double)naver;
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
  return exp(-v*v*oneOnSigma*oneOnSigma);
}


void calcSourceFn(double dTau, const inputPars *par, double *remnantSnu, double *expDTau){
  /*
  The source function S is defined as j_nu/alpha, which is clearly not
  defined for alpha==0. However S is used in the algorithm only in the
  term (1-exp[-alpha*ds])*S, which is defined for all values of alpha.
  The present function calculates this term and returns it in the
  argument remnantSnu. For values of abs(alpha*ds) less than a pre-
  calculated cutoff supplied in inputPars, a Taylor approximation is
  used.

  Note that the same cutoff condition holds for replacement of
  exp(-dTau) by its Taylor expansion to 3rd order.
  */

  if (fabs(dTau)<par->taylorCutoff){
    *remnantSnu = 1. - dTau*(1. - dTau/3.)/2.;
    *expDTau = 1. - dTau*(*remnantSnu);
  } else {
    *expDTau = exp(-dTau);
    *remnantSnu = (1.-(*expDTau))/dTau;
  }
}


void
photon(int id, struct grid *g, molData *m, int iter, const gsl_rng *ran,inputPars *par,blend *matrix){
  int iphot,iline,jline,here,there,firststep,dir,np_per_line,ip_at_line,l;
  int *counta, *countb,nlinetot;
  double deltav,segment,vblend,dtau,expDTau,jnu,alpha,ds,vfac[par->nSpecies],pt_theta,pt_z,semiradius;
  double *tau,*expTau,vel[3],x[3], inidir[3];
  double remnantSnu;

  lineCount(par->nSpecies, m, &counta, &countb, &nlinetot);
  tau=malloc(sizeof(double)*nlinetot);
  expTau=malloc(sizeof(double)*nlinetot);
  velocity(g[id].x[0],g[id].x[1],g[id].x[2],vel);
  
  np_per_line=(int) g[id].nphot/g[id].numNeigh; // Works out to be equal to ininphot. :-/

  for(iphot=0;iphot<g[id].nphot;iphot++){
    firststep=1;
    for(iline=0;iline<nlinetot;iline++){
      m[0].phot[iline+iphot*m[0].nline]=0.;
      tau[iline]=0.;
      expTau[iline]=1.;
    }
    
    /* Initial velocity, direction and frequency offset  */		
    pt_theta=gsl_rng_uniform(ran)*2*PI;
    pt_z=2*gsl_rng_uniform(ran)-1;
    semiradius = sqrt(1.-pt_z*pt_z);
    inidir[0]=semiradius*cos(pt_theta);
    inidir[1]=semiradius*sin(pt_theta);
    inidir[2]=pt_z;
    
    iter=(int) (gsl_rng_uniform(ran)*(double)N_RAN_PER_SEGMENT); // can have values in [0,1,..,N_RAN_PER_SEGMENT-1]
    ip_at_line=(int) iphot/g[id].numNeigh;
    segment=(N_RAN_PER_SEGMENT*(ip_at_line-np_per_line/2.)+iter)/(double)(np_per_line*N_RAN_PER_SEGMENT);
    /*
    Values of segment should be evenly distributed (considering the
    entire ensemble of photons) between -0.5 and +0.5, and are chosen
    from a sequence of possible values separated by
    1/(N_RAN_PER_SEGMENT*ininphot).
    */
    
    dir=sortangles(inidir,id,g,ran);
    here=g[id].id;
    there=g[here].neigh[dir]->id;
    deltav=segment*4.3*g[id].dopb+veloproject(g[id].dir[dir].xn,vel);
    
    /* Photon propagation loop */
    do{
      if(firststep){
        firststep=0;				
        ds=g[here].ds[dir]/2.;
        for(l=0;l<par->nSpecies;l++){
          if(!par->doPregrid) velocityspline(g,here,dir,g[id].mol[l].binv,deltav,&vfac[l]);
          else velocityspline_lin(g,here,dir,g[id].mol[l].binv,deltav,&vfac[l]);
          m[l].vfac[iphot]=vfac[0];
          m[l].ds[iphot]=ds;
        }
        for(l=0;l<3;l++) x[l]=g[here].x[l]+(g[here].dir[dir].xn[l] * g[id].ds[dir]/2.);
      } else {
        ds=g[here].ds[dir];
        for(l=0;l<3;l++) x[l]=g[here].x[l];
      }
      
      for(l=0;l<par->nSpecies;l++){
        if(!par->doPregrid) velocityspline(g,here,dir,g[id].mol[l].binv,deltav,&vfac[l]);
        else velocityspline_lin(g,here,dir,g[id].mol[l].binv,deltav,&vfac[l]);
      }
      
      for(iline=0;iline<nlinetot;iline++){
        jnu=0.;
        alpha=0.;
        
        sourceFunc_line(&jnu,&alpha,m,vfac[counta[iline]],g,here,counta[iline],countb[iline]);
        sourceFunc_cont(&jnu,&alpha,g,here,counta[iline],countb[iline]);

        dtau=alpha*ds;
        if(dtau < -30) dtau = -30;
        calcSourceFn(dtau, par, &remnantSnu, &expDTau);
        remnantSnu *= jnu*m[0].norminv*ds;

        m[0].phot[iline+iphot*m[0].nline]+=expTau[iline]*remnantSnu;
        tau[iline]+=dtau;
        expTau[iline]*=expDTau;
        if(tau[iline] < -30.){
          if(!silent) warning("Maser warning: optical depth has dropped below -30");
          tau[iline]= -30.; 
          expTau[iline]=exp(-tau[iline]);
        }
        
        /* Line blending part */
        if(par->blend){
          jnu=0.;
          alpha=0.;
          for(jline=0;jline<sizeof(matrix)/sizeof(blend);jline++){
            if(matrix[jline].line1 == jline || matrix[jline].line2 == jline){	
              if(!par->doPregrid) velocityspline(g,here,dir,g[id].mol[counta[jline]].binv,deltav-matrix[jline].deltav,&vblend);
              else velocityspline_lin(g,here,dir,g[id].mol[counta[jline]].binv,deltav-matrix[jline].deltav,&vblend);	
              sourceFunc_line(&jnu,&alpha,m,vblend,g,here,counta[jline],countb[jline]);
              dtau=alpha*ds;
              if(dtau < -30) dtau = -30;
              calcSourceFn(dtau, par, &remnantSnu, &expDTau);
              remnantSnu *= jnu*m[0].norminv*ds;

              m[0].phot[jline+iphot*m[0].nline]+=expTau[jline]*remnantSnu;
              tau[jline]+=dtau;
              expTau[jline]*=expDTau;
              if(tau[jline] < -30.){
                if(!silent) warning("Optical depth has dropped below -30");
                tau[jline]= -30.; 
                expTau[jline]=exp(-tau[jline]);
              }
            }
          }
        }
        /* End of line blending part */
      }
      
      dir=sortangles(inidir,there,g,ran);
      here=there;
      there=g[here].neigh[dir]->id;
    } while(!g[there].sink);
    
    /* Add cmb contribution */
    if(m[0].cmb[0]>0.){
      for(iline=0;iline<nlinetot;iline++){
        m[0].phot[iline+iphot*m[0].nline]+=expTau[iline]*m[counta[iline]].cmb[countb[iline]];
      }
    }
  }
  free(tau);
  free(counta);
  free(countb);
}

void
getjbar(int posn, molData *m, struct grid *g, inputPars *par){
  int iline,iphot;
  double tau, expTau, remnantSnu, vsum=0., jnu, alpha;
  int *counta, *countb,nlinetot;

  lineCount(par->nSpecies, m, &counta, &countb, &nlinetot);
  
  for(iline=0;iline<m[0].nline;iline++) m[0].jbar[iline]=0.;
  for(iphot=0;iphot<g[posn].nphot;iphot++){
    if(m[0].vfac[iphot]>0){
      for(iline=0;iline<m[0].nline;iline++){
        jnu=0.;
        alpha=0.;
        
        sourceFunc_line(&jnu,&alpha,m,m[0].vfac[iphot],g,posn,counta[iline],countb[iline]);
        sourceFunc_cont(&jnu,&alpha,g,posn,counta[iline],countb[iline]);
        tau=alpha*m[0].ds[iphot];
        calcSourceFn(tau, par, &remnantSnu, &expTau);
        remnantSnu *= jnu*m[0].norminv*m[0].ds[iphot];

        m[0].jbar[iline]+=m[0].vfac[iphot]*(expTau*m[0].phot[iline+iphot*m[0].nline]+remnantSnu);
      }
      vsum+=m[0].vfac[iphot];
    }
  }
  for(iline=0;iline<m[0].nline;iline++) m[0].jbar[iline]=m[0].norm*m[0].jbar[iline]/vsum;
  free(counta);
  free(countb);
}


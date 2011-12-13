/*
 *  photon.c
 *  LIME, The versatile 3D line modeling environment 
 *
 *  Created by Christian Brinch on 15/11/06.
 *  Copyright 2006-2011, Christian Brinch, 
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *  All rights reserved.
 *
 */

#include "lime.h"

int
sortangles(int id, int dir, int newid, struct grid *g, const gsl_rng *ran) {
	int i,j,itmp,val;
	/* int b,l,k; */
	double *angle,tmp,best;
	int *n;

	angle=malloc(sizeof(double)*g[newid].numNeigh);	
	n=malloc(sizeof(int)*g[newid].numNeigh);	

	best=1e30;
	for(i=0;i<g[newid].numNeigh;i++){
	  angle[i]=( g[id].dir[dir].xn[0]*g[newid].dir[i].xn[0]
       		    +g[id].dir[dir].xn[1]*g[newid].dir[i].xn[1]
	    		+g[id].dir[dir].xn[2]*g[newid].dir[i].xn[2]);
	  n[i]=i;
/*	  if(angle[i]<best){
	 	best=angle[i];
        b=i;
      }*/
	}
	
	for(i=g[newid].numNeigh-1; i>=0; i--) {
	  for(j=1; j<=i; j++) {
	    if(angle[j-1]<angle[j]) {
		  tmp=angle[j-1];
		  angle[j-1]=angle[j];
		  angle[j]=tmp;

		  itmp=n[j-1];
		  n[j-1]=n[j];
		  n[j]=itmp;
	    }
	  }	
	}	  

  	if(gsl_rng_uniform(ran)<1./((1-angle[0])/(1-angle[1])+1) ) {
  	  val=n[0];
	} else {
	  val=n[1];	
	}
	free(angle);
	free(n);
	return val;

/*    return angle[b]; */
}



void
velocityspline(struct grid *g, int id, int k, double binv, double deltav, double *vfac){
  int nspline,ispline,naver,iaver;
  double v1,v2,s1,s2,sd,v,vfacsub,d;
	
  v1=deltav-veloproject(g[id].dir[k].xn,g[id].vel);
  v2=deltav-veloproject(g[id].dir[k].xn,g[g[id].neigh[k]].vel);

  nspline=(fabs(v1-v2)*binv < 1) ? 1 : (int)(fabs(v1-v2)*binv);
  *vfac=0.;
  s2=0;
  v2=v1;
	
  for(ispline=0;ispline<nspline;ispline++){
	  s1=s2;
    s2=((double)(ispline+1))/(double)nspline;					
    v1=v2;
    d=s2*g[id].ds[k];
    v2=deltav-(g[id].a4[k]*pow(d,4)+g[id].a3[k]*pow(d,3)+g[id].a2[k]*pow(d,2)+g[id].a1[k]*d+g[id].a0[k]);
		naver=(1 > fabs(v1-v2)*binv) ? 1 : (int)(fabs(v1-v2)*binv);
		for(iaver=0;iaver<naver;iaver++){
	    sd=s1+(s2-s1)*((double)iaver-0.5)/(double)naver;
      d=sd*g[id].ds[k];
	    v=deltav-(g[id].a4[k]*pow(d,4)+g[id].a3[k]*pow(d,3)+g[id].a2[k]*pow(d,2)+g[id].a1[k]*d+g[id].a0[k]);
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
  v2=deltav-veloproject(g[id].dir[k].xn,g[g[id].neigh[k]].vel);
	
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



void
velocityspline2(double x[3], double dx[3], double ds, double binv, double deltav, double *vfac){
	int nspline,ispline,naver,iaver;
	double v1,v2,s1,s2,sd,v,vfacsub,vel[3];
	
	velocity(x[0],x[1],x[2],vel);
	v1=deltav-veloproject(dx,vel);
	velocity(x[0]+(dx[0]*ds),x[1]+(dx[1]*ds),x[2]+(dx[2]*ds),vel);
	v2=deltav-veloproject(dx,vel);

	nspline=(fabs(v1-v2)*binv < 1) ? 1 : (int)(fabs(v1-v2)*binv);
	*vfac=0.;
	s2=0;
	v2=v1;
	
	for(ispline=0;ispline<nspline;ispline++){
		s1=s2;
		s2=((double)(ispline+1))/(double)nspline;					
		v1=v2;
		velocity(x[0]+(s2*dx[0]*ds),x[1]+(s2*dx[1]*ds),x[2]+(s2*dx[2]*ds),vel);
		v2=deltav-veloproject(dx,vel);
		naver=(1 > fabs(v1-v2)*binv) ? 1 : (int)(fabs(v1-v2)*binv);
		for(iaver=0;iaver<naver;iaver++){
			sd=s1+(s2-s1)*((double)iaver-0.5)/(double)naver;
			velocity(x[0]+(sd*dx[0]*ds),x[1]+(sd*dx[1]*ds),x[2]+(sd*dx[2]*ds),vel);
			v=deltav-veloproject(dx,vel);
			vfacsub=gaussline(v,binv);
			*vfac+=vfacsub/(double)naver;
		}
	}
	*vfac= *vfac/(double)nspline;

	return;
}


double veloproject(double dx[3], double *vel)
{
	return dx[0]*vel[0]+dx[1]*vel[1]+dx[2]*vel[2];
}


double gaussline(double v, double sigma)
{
	int maxgau=101,maxsig=4,ival;
	double fac,val;
	
	fac=(maxgau-1)/maxsig;
	ival=(int)(fac*(fabs(v)*sigma))+1;
	if((ival-1)>=maxgau) return 0.;
	val=(ival*ival)/(fac*fac);
	return exp(-val);
}


void
photon(int id, struct grid *g, molData *m, int iter, const gsl_rng *ran,inputPars *par,blend *matrix){
	int iphot,iline,jline,here,there,firststep,inidir,dir,np_per_line,ip_at_line,l;
	int *counta, *countb,nlinetot;
	double deltav,segment,vblend,snu,dtau,jnu,alpha,ds,vfac[par->nSpecies];
	double *tau,vel[3],x[3];
				
	lineCount(par->nSpecies, m, &counta, &countb, &nlinetot);	

	tau=malloc(sizeof(double)*nlinetot);
	
	velocity(g[id].x[0],g[id].x[1],g[id].x[2],vel);
	
	for(iphot=0;iphot<g[id].nphot;iphot++){
		firststep=1;
		for(iline=0;iline<nlinetot;iline++){
			m[0].phot[iline+iphot*m[0].nline]=0.;
			tau[iline]=0.;
		}

/* Initial velocity, direction and frequency offset  */		
		inidir=iphot%g[id].numNeigh;

		iter=(int) (gsl_rng_uniform(ran)*3.);
		np_per_line=(int) g[id].nphot/g[id].numNeigh;
		ip_at_line=(int) iphot/g[id].numNeigh;
		segment=1/(2.*np_per_line)*(2*ip_at_line-np_per_line+iter);
		deltav=segment*4.3*g[id].dopb+veloproject(g[id].dir[inidir].xn,vel);

		dir=inidir;
		here=g[id].id;
		there=g[g[here].neigh[dir]].id;

/* Photon propagation loop */
		do{
		  if(firststep){
			firststep=0;				
			ds=g[here].ds[dir]/2.;
			for(l=0;l<par->nSpecies;l++){
              if(!par->pregrid) velocityspline(g,here,dir,g[id].mol[l].binv,deltav,&vfac[l]);
              else velocityspline_lin(g,here,dir,g[id].mol[l].binv,deltav,&vfac[l]);
			  m[l].vfac[iphot]=vfac[0];
              m[l].weight[iphot]=g[here].w[inidir];	
		    }
			for(l=0;l<3;l++) x[l]=g[here].x[l]+(g[here].dir[dir].xn[l] * g[id].ds[inidir]/2.);
		  } else {
			ds=g[here].ds[dir];
			for(l=0;l<3;l++) x[l]=g[here].x[l];
		  }

		  for(l=0;l<par->nSpecies;l++){
		    if(!par->pregrid) velocityspline(g,here,dir,g[id].mol[l].binv,deltav,&vfac[l]);
		    else velocityspline_lin(g,here,dir,g[id].mol[l].binv,deltav,&vfac[l]);
		  }
		
		  for(iline=0;iline<nlinetot;iline++){
			jnu=0.;
			alpha=0.;

			sourceFunc_line(&jnu,&alpha,m,vfac[counta[iline]],g,here,counta[iline],countb[iline]);
			sourceFunc_cont(&jnu,&alpha,g,here,counta[iline],countb[iline]);
			if(fabs(alpha)>0.){
			  snu=(jnu/alpha)*m[0].norminv;
		  	  dtau=alpha*ds;
			}
		    m[0].phot[iline+iphot*m[0].nline]+=exp(-tau[iline])*(1.-exp(-dtau))*snu;	
			tau[iline]+=dtau;
				
			if(par->blend){
			  jnu=0.;
			  alpha=0.;
  			  for(jline=0;jline<sizeof(matrix)/sizeof(blend);jline++){
	  			if(matrix[jline].line1 == jline || matrix[jline].line2 == jline){	
				  if(!par->pregrid) velocityspline(g,here,dir,g[id].mol[counta[jline]].binv,deltav-matrix[jline].deltav,&vblend);
		  		  else velocityspline_lin(g,here,dir,g[id].mol[counta[jline]].binv,deltav-matrix[jline].deltav,&vblend);	
  				  sourceFunc_line(&jnu,&alpha,m,vblend,g,here,counta[jline],countb[jline]);
				  if(fabs(alpha)>0.){
				    snu=(jnu/alpha)*m[0].norminv;
				    dtau=alpha*ds;
				  }
				  m[0].phot[jline+iphot*m[0].nline]+=exp(-tau[jline])*(1.-exp(-dtau))*snu;
				  tau[jline]+=dtau;
				  if(tau[jline]<-30) tau[iline]=30;
			    }
			  }
			}
		  }
			
	  	  dir=sortangles(id,inidir,there,g,ran);
		  here=there;
		  there=g[g[here].neigh[dir]].id;
		} while(!g[there].sink);
		
		/* Add cmb contribution */
		if(m[0].cmb[0]>0.){
			for(iline=0;iline<nlinetot;iline++){
			  m[0].phot[iline+iphot*m[0].nline]+=exp(-tau[iline])*m[counta[iline]].cmb[countb[iline]];
			}
		}
				

	}
	free(tau);
	free(counta);
	free(countb);
}



void getjbar(int posn, molData *m, struct grid *g, inputPars *par){
  int iline,iphot;	      
  double tau, snu, vsum=0., jnu, alpha;
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
		if(fabs(alpha)>0.){
		  snu=(jnu/alpha)*m[0].norminv;
	  	  tau=alpha*m[0].ds[iphot];
		}
		m[0].jbar[iline]+=m[0].vfac[iphot]*(exp(-tau)*m[0].phot[iline+iphot*m[0].nline]+(1.-exp(-tau))*snu)*m[0].weight[iphot];
	  }
	  vsum+=m[0].vfac[iphot]*m[0].weight[iphot];	
	}
  } 
  for(iline=0;iline<m[0].nline;iline++) m[0].jbar[iline]=m[0].norm*m[0].jbar[iline]/vsum;
}




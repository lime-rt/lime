/*
 *  blowpops.c
 *  LIME, The versatile 3D line modeling environment 
 *
 *  Created by Christian Brinch on 14/11/07.
 *  Copyright 2006-2013, Christian Brinch, 
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */
#include "lime.h"

void
popsout(inputPars *par, struct grid *g, molData *m){
	FILE *fp;
	int j,k,l;
	double dens;
	/* int i,mi,c,q=0,best; */
	/* double vel[3],ra[100],rb[100],za[100],zb[100],min; */
	
    if((fp=fopen(par->outputfile, "w"))==NULL){
      if(!silent) bail_out("Error writing output populations file!");
      exit(1);
    }
	fprintf(fp,"# Column definition: x, y, z, H2 density, kinetic gas temperature, molecular abundance, convergence flag, pops_0...pops_n\n"); 
	for(j=0;j<par->ncell-par->sinkPoints;j++){
	    dens=0.;
	 	for(l=0;l<par->collPart;l++) dens+=g[j].dens[l];
	    fprintf(fp,"%e %e %e %e %e %e %d ", g[j].x[0], g[j].x[1], g[j].x[2], dens, g[j].t[0], g[j].nmol[0]/dens, g[j].conv);
		for(k=0;k<m[0].nlev;k++) fprintf(fp,"%e ",g[j].mol[0].pops[k]);
		fprintf(fp,"\n");
	}	
	fclose(fp);


/*
	fp=fopen(par->outputfile,"w");
	fprintf(fp,"rmax=%e\n", par->radius);
	fprintf(fp,"zmax=%e\n", par->radius);
	fprintf(fp,"ncell=100\n");
    fprintf(fp,"tcmb=2.728\n");
	fprintf(fp,"gas:dust=100\n");
	fprintf(fp,"columns=id,ra,rb,za,zb,nh,tk,nm,vr,vz,va,db,td,lp\n");
	fprintf(fp,"molfile=/disks/chem19/brinch/ratran/molec/p-h2o.dat\n");
	fprintf(fp,"velo=/disks/chem19/ratran/velocity/vgrid_2d_sky.f\n");
	fprintf(fp,"kappa=jena,bare,e7\n");
	fprintf(fp,"@\n");
	j=0;

	for(k=0;k<1;k++){
	  for(i=0;i<100;i++){
	    if(i%((int) sqrt(100))>((int) sqrt(900))/2-1 || i>(900/2-1) || k==0){
	      ra[j]=(i%(int) sqrt(100))  * (pow(2,k)*par->radius/(5*pow(2,1)));
	      rb[j]=(i%(int) sqrt(100)+1)* (pow(2,k)*par->radius/(5*pow(2,1)));
	      za[j]=(i/(int) sqrt(100))  * (pow(2,k)*par->radius/(5*pow(2,1)));
	      zb[j]=(i/(int) sqrt(100)+1)* (pow(2,k)*par->radius/(5*pow(2,1)));
		  j++;
		}
	  }
	}
	c=0;

	
	for(k=0;k<100;k++){	
	  min=1e30;
	  mi=-1;
	  for(j=0;j<par->pIntensity;j++){
	    if(pow((ra[k]+rb[k])/2.-sqrt(g[j].x[0]*g[j].x[0]+g[j].x[1]*g[j].x[1]),2) +
		   pow((za[k]+zb[k])/2.-fabs(g[j].x[2]),2)	< min){
		  min=pow((ra[k]+rb[k])/2.-sqrt(g[j].x[0]*g[j].x[0]+g[j].x[1]*g[j].x[1]),2) +
			  pow((za[k]+zb[k])/2.-fabs(g[j].x[2]),2);
		  mi=j;
		}
      }
  	
	  
	  fprintf(fp,"%d %e %e %e %e ", ++c, ra[k],rb[k], za[k], zb[k]);
	  if(mi==-1) fprintf(fp,"%e %e %e %e %e %e %e %e ",0,10,0,0.0,0.0,0.0,0.2,10); 
	  else fprintf(fp,"%e %e %e %e %e %e %e %e ",g[mi].dens[0]/1e6,g[mi].t[0],g[mi].nmol[0]/1e6, 0.0,0.0,0.0, 0.2,g[mi].t[1]);
	  for(i=0;i<m[0].nlev;i++) fprintf(fp,"%e ",g[mi].mol[0].pops[i]);
	  fprintf(fp,"\n");
	}
		 
	
	fclose(fp);
*/
}

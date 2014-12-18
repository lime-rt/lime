/*
 *  statistics.c
 *  LIME, The versatile 3D line modeling environment
 *
 *  Created by Christian Brinch on 18/12/08.
 *  Copyright 2006-2014, Christian Brinch,
 *  <brinch@nbi.dk>
 *  Niels Bohr institutet
 *  University of Copenhagen
 *	All rights reserved.
 *
 */


#include "lime.h"

void statistics(int id, molData *m, struct grid *g, int *exceed, double *opops, double *oopops, int *conv){
  int ilev;
  double avepops,var=0.,snr1=fixset,snr2=fixset;

  for(ilev=0;ilev<m[0].nlev;ilev++){
    avepops=(g[id].mol[0].pops[ilev]+opops[ilev+id*m[0].nlev]+oopops[ilev+id*m[0].nlev])/3.;
    if(avepops>=0.01){
      var=gsl_max(fabs(g[id].mol[0].pops[ilev]-avepops)/avepops,
                  gsl_max(fabs( opops[ilev+id*m[0].nlev]-avepops)/avepops,
                          fabs(oopops[ilev+id*m[0].nlev]-avepops)/avepops));
      snr1 = (snr1>var) ? snr1 : var;
    } else if(avepops>=minpop && avepops<0.01){
      var=gsl_max(fabs(g[id].mol[0].pops[ilev]-avepops)/avepops,
                  gsl_max(fabs( opops[ilev+id*m[0].nlev]-avepops)/avepops,
                          fabs(oopops[ilev+id*m[0].nlev]-avepops)/avepops));
      snr2 = (snr2>var) ? snr2 : var;
    }

  }

  snr1=1./snr1;
  snr2=1./snr2;

  if(snr1>=goal && snr2>=5) {
    *conv+=1;
  }else {
    g[id].nphot=g[id].nphot+10;
    if(g[id].nphot>(max_phot*200)){
      g[id].nphot=max_phot*200;
      *exceed+=1;
      if(!silent) warning("Warning: limiting nphot reached in a grid point");
    }
  }
}

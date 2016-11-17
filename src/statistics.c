/*
 *  statistics.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2016 The LIME development team
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
    if(g[id].nphot>(MAX_RAYS_PER_POINT)){
      g[id].nphot=MAX_RAYS_PER_POINT;
      *exceed+=1;
      if(!silent) warning("Warning: limiting nphot reached in a grid point");
    }
  }
}

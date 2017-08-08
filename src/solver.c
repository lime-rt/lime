/*
 *  solver.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
TODO:
  - The test to run calculateJBar() etc in levelPops just tests dens[0]. This is a bit sloppy.
 */

#include "lime.h"
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_errno.h>

/* Data concerning a single grid vertex which is passed from calculateJBar() to solveStatEq(). This data needs to be thread-safe. */
typedef struct {
  double *jbar,*phot,*vfac,*vfac_loc;
} gridPointData;

struct blend{
  int molJ, lineJ;
  double deltaV;
};

struct lineWithBlends{
  int lineI, numBlends;
  struct blend *blends;
};

struct molWithBlends{
  int molI, numLinesWithBlends;
  struct lineWithBlends *lines;
};

struct blendInfo{
  int numMolsWithBlends;
  struct molWithBlends *mols;
};

/*....................................................................*/
void
freeMolsWithBlends(struct molWithBlends *mols, const int numMolsWithBlends){
  int mi, li;

  if(mols != NULL){
    for(mi=0;mi<numMolsWithBlends;mi++){
      if(mols[mi].lines != NULL){
        for(li=0;li<mols[mi].numLinesWithBlends;li++)
          free(mols[mi].lines[li].blends);
        free(mols[mi].lines);
      }
    }
    free(mols);
  }
}

/*....................................................................*/
void
freeGridPointData(const int nSpecies, gridPointData *mol){
  int i;
  if(mol!= NULL){
    for(i=0;i<nSpecies;i++){
      free(mol[i].jbar);
      free(mol[i].phot);
      free(mol[i].vfac);
      free(mol[i].vfac_loc);
    }
  }
}

/*....................................................................*/
void lineBlend(molData *m, configInfo *par, struct blendInfo *blends){
  /*
This obtains information on all the lines of all the radiating species which have other lines within some cutoff velocity separation.

A variable of type 'struct blendInfo' has a nested structure which can be illustrated diagrammaticaly as follows.

  Structs:	blendInfo		molWithBlends		lineWithBlends		blend

  Variables:	blends
		  .numMolsWithBlends     ____________________
		  .*mols--------------->|.molI               |
		                        |.numLinesWithBlends |   ___________
		                        |.*lines--------------->|.lineI     |
		                        |____________________|  |.numBlends |           ________
		                        |        etc         |  |.*blends------------->|.molJ   |
		                                                |___________|          |.lineJ  |
		                                                |    etc    |          |.deltaV |
		                                                                       |________|
		                                                                       |   etc  |

Pointers are indicated by a * before the attribute name and an arrow to the memory location pointed to.
  */
  int molI, lineI, molJ, lineJ;
  int nmwb, nlwb, numBlendsFound, li, bi;
  double deltaV;
  struct blend *tempBlends=NULL;
  struct lineWithBlends *tempLines=NULL;

  /* Dimension blends.mols first to the total number of species, then realloc later if need be.
  */
  (*blends).mols = malloc(sizeof(struct molWithBlends)*par->nSpecies);
  (*blends).numMolsWithBlends = 0;

  nmwb = 0;
  for(molI=0;molI<par->nSpecies;molI++){
    tempBlends = malloc(sizeof(struct blend)*m[molI].nline);
    tempLines  = malloc(sizeof(struct lineWithBlends)*m[molI].nline);

    nlwb = 0;
    for(lineI=0;lineI<m[molI].nline;lineI++){
      numBlendsFound = 0;
      for(molJ=0;molJ<par->nSpecies;molJ++){
        for(lineJ=0;lineJ<m[molJ].nline;lineJ++){
          if(!(molI==molJ && lineI==lineJ)){
            deltaV = (m[molJ].freq[lineJ] - m[molI].freq[lineI])*CLIGHT/m[molI].freq[lineI];
            if(fabs(deltaV)<maxBlendDeltaV){
              tempBlends[numBlendsFound].molJ   = molJ;
              tempBlends[numBlendsFound].lineJ  = lineJ;
              tempBlends[numBlendsFound].deltaV = deltaV;
              numBlendsFound++;
            }
          }
        }
      }

      if(numBlendsFound>0){
        tempLines[nlwb].lineI = lineI;
        tempLines[nlwb].numBlends = numBlendsFound;
        tempLines[nlwb].blends = malloc(sizeof(struct blend)*numBlendsFound);
        for(bi=0;bi<numBlendsFound;bi++)
          tempLines[nlwb].blends[bi] = tempBlends[bi];

        nlwb++;
      }
    }

    if(nlwb>0){
      (*blends).mols[nmwb].molI = molI;
      (*blends).mols[nmwb].numLinesWithBlends = nlwb;
      (*blends).mols[nmwb].lines = malloc(sizeof(struct lineWithBlends)*nlwb);
      for(li=0;li<nlwb;li++){
        (*blends).mols[nmwb].lines[li].lineI     = tempLines[li].lineI;
        (*blends).mols[nmwb].lines[li].numBlends = tempLines[li].numBlends;
        (*blends).mols[nmwb].lines[li].blends = malloc(sizeof(struct blend)*tempLines[li].numBlends);
        for(bi=0;bi<tempLines[li].numBlends;bi++)
          (*blends).mols[nmwb].lines[li].blends[bi] = tempLines[li].blends[bi];
      }

      nmwb++;
    }

    free(tempLines);
    free(tempBlends);
  }

  (*blends).numMolsWithBlends = nmwb;
  if(nmwb>0){
    if(!par->blend)
      if(!silent) warning("There are blended lines, but line blending is switched off.");

    (*blends).mols = realloc((*blends).mols, sizeof(struct molWithBlends)*nmwb);
  }else{
    if(par->blend)
      if(!silent) warning("Line blending is switched on, but no blended lines were found.");

    free((*blends).mols);
    (*blends).mols = NULL;
  }
}

/*....................................................................*/
void calcGridCollRates(configInfo *par, molData *md, struct grid *gp){
  int i,id,ipart,itrans,itemp,tnint=-1;
  struct cpData part;
  double fac;

  for(i=0;i<par->nSpecies;i++){
    for(id=0;id<par->ncell;id++){
      gp[id].mol[i].partner = malloc(sizeof(struct rates)*md[i].npart);
    }

    for(ipart=0;ipart<md[i].npart;ipart++){
      part = md[i].part[ipart];
      for(id=0;id<par->ncell;id++){
        for(itrans=0;itrans<part.ntrans;itrans++){
          if((gp[id].t[0]>part.temp[0])&&(gp[id].t[0]<part.temp[part.ntemp-1])){
            for(itemp=0;itemp<part.ntemp-1;itemp++){
              if((gp[id].t[0]>part.temp[itemp])&&(gp[id].t[0]<=part.temp[itemp+1])){
                tnint=itemp;
              }
            }
            fac=(gp[id].t[0]-part.temp[tnint])/(part.temp[tnint+1]-part.temp[tnint]);
            gp[id].mol[i].partner[ipart].t_binlow = tnint;
            gp[id].mol[i].partner[ipart].interp_coeff = fac;

	  } else if(gp[id].t[0]<=part.temp[0]) {
	    gp[id].mol[i].partner[ipart].t_binlow = 0;
	    gp[id].mol[i].partner[ipart].interp_coeff = 0.0;
	  } else {
	    gp[id].mol[i].partner[ipart].t_binlow = part.ntemp-2;
	    gp[id].mol[i].partner[ipart].interp_coeff = 1.0;
	  }
        } /* End loop over transitions. */
      } /* End loop over grid points. */
    } /* End loop over collision partners. */
  } /* End loop over radiating molecules. */
}

/*....................................................................*/
void calcGridLinesDustOpacity(configInfo *par, molData *md, double *lamtab\
  , double *kaptab, const int nEntries, struct grid *gp){

  int iline,id,si;
  double *kappatab,gtd;
  gsl_spline *spline = NULL;
  gsl_interp_accel *acc = NULL;
  double *knus=NULL, *dusts=NULL;

  if(par->dust != NULL){
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline,nEntries);
    gsl_spline_init(spline,lamtab,kaptab,nEntries);
  }

  for(id=0;id<par->ncell;id++){
    for(si=0;si<par->nSpecies;si++){
      free(gp[id].mol[si].cont);
      gp[id].mol[si].cont = malloc(sizeof(*(gp[id].mol[si].cont))*md[si].nline);
    }
  }

  for(si=0;si<par->nSpecies;si++){
    kappatab = malloc(sizeof(*kappatab)*md[si].nline);
    knus     = malloc(sizeof(*knus)    *md[si].nline);
    dusts    = malloc(sizeof(*dusts)   *md[si].nline);

    if(par->dust == NULL){
      for(iline=0;iline<md[si].nline;iline++)
        kappatab[iline] = 0.;
    }else{
      for(iline=0;iline<md[si].nline;iline++)
        kappatab[iline] = interpolateKappa(md[si].freq[iline]\
                        , lamtab, kaptab, nEntries, spline, acc);
    }

    for(id=0;id<par->ncell;id++){
      gasIIdust(gp[id].x[0],gp[id].x[1],gp[id].x[2],&gtd);
      calcDustData(par, gp[id].dens, md[si].freq, gtd, kappatab, md[si].nline, gp[id].t, knus, dusts);
      for(iline=0;iline<md[si].nline;iline++){
        gp[id].mol[si].cont[iline].knu  = knus[iline];
        gp[id].mol[si].cont[iline].dust = dusts[iline];
      }
    }

    free(kappatab);
    free(knus);
    free(dusts);
  }

  if(par->dust != NULL){
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
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

/*....................................................................*/
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
void
calculateJBar(int id, struct grid *gp, molData *md, const gsl_rng *ran\
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
    pt_theta=gsl_rng_uniform(ran)*2*M_PI;
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
          if(par->edgeVelsAvailable) {
            calcLineAmpPWLin(gp,here,neighI,molI,deltav[molI],inidir,&vfac_in[molI],&vfac_out[molI]);
         } else
            calcLineAmpLin(gp,here,neighI,molI,deltav[molI],inidir,&vfac_in[molI],&vfac_out[molI]);

          mp[molI].vfac[iphot]=vfac_out[molI];
        }
        /*
        Contribution of the local cell to emission and absorption is done in updateJBar.
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
        if(par->edgeVelsAvailable)
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
              if(par->edgeVelsAvailable)
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
updateJBar(int posn, molData *md, struct grid *gp, const int molI\
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

/*....................................................................*/
void
getFixedMatrix(molData *md, int ispec, struct grid *gp, int id, gsl_matrix *colli, configInfo *par){
  int ipart,k,l,ti;

  /* Initialize matrix with zeros */
  if(md[ispec].nlev<=0){
    if(!silent) bail_out("Matrix initialization error in solveStatEq");
    exit(1);
  }
  gsl_matrix_set_zero(colli);

  /* Populate matrix with collisional transitions */
  for(ipart=0;ipart<md[ispec].npart;ipart++){
    struct cpData part = md[ispec].part[ipart];
    double *downrates = part.down;
    int di = md[ispec].part[ipart].densityIndex;
    if (di<0) continue;

    for(ti=0;ti<part.ntrans;ti++){
      int coeff_index = ti*part.ntemp + gp[id].mol[ispec].partner[ipart].t_binlow;
      double down = downrates[coeff_index]\
                  + gp[id].mol[ispec].partner[ipart].interp_coeff*(downrates[coeff_index+1]\
                  - downrates[coeff_index]);
      double up = down*md[ispec].gstat[part.lcu[ti]]/md[ispec].gstat[part.lcl[ti]]\
                *exp(-HCKB*(md[ispec].eterm[part.lcu[ti]]-md[ispec].eterm[part.lcl[ti]])/gp[id].t[0]);

      gsl_matrix_set(colli, part.lcu[ti], part.lcl[ti], gsl_matrix_get(colli, part.lcu[ti], part.lcl[ti]) - down*gp[id].dens[di]);
      gsl_matrix_set(colli, part.lcl[ti], part.lcu[ti], gsl_matrix_get(colli, part.lcl[ti], part.lcu[ti]) - up*gp[id].dens[di]);
    }

  }

  /* Does this work with >1 coll. part? */
  double *ctot  = malloc(sizeof(double)*md[ispec].nlev);
  for(k=0;k<md[ispec].nlev;k++){     
    ctot[k]=0.0;
    for(l=0;l<md[ispec].nlev;l++)
      ctot[k] += gsl_matrix_get(colli,k,l);
    gsl_matrix_set(colli,k,k,gsl_matrix_get(colli,k,k) - ctot[k]);
  }
  free(ctot);

  /* For some reason collision matrix as constructed above is transposed.
  Someone who is not this lazy could fix the loops instead of using gsl_matrix_transpose. */
  gsl_matrix_transpose(colli);

  double *girtot = malloc(sizeof(double)*md[ispec].nlev);
  if(par->girdatfile!=NULL){
    for(k=0;k<md[ispec].nlev;k++){
      girtot[k] = 0;
      for(l=0;l<md[ispec].nlev;l++)
        girtot[k] += md[ispec].gir[k*md[ispec].nlev+l];
    }
    for(k=0;k<md[ispec].nlev;k++){
      gsl_matrix_set(colli,k,k,gsl_matrix_get(colli,k,k)+girtot[k]);
      for(l=0;l<md[ispec].nlev;l++){
        if(k!=l){
          if(par->girdatfile!=NULL)
            gsl_matrix_set(colli,k,l,gsl_matrix_get(colli,k,l)-md[ispec].gir[l*md[ispec].nlev+k]);
        }
      }
    }
  }
  free(girtot);
}

/*....................................................................*/
void
getMatrix(gsl_matrix *matrix, molData *md, int ispec, gridPointData *mp, gsl_matrix *colli){
  int k,l,li;

  /* Initialize matrix by copying the fixed part */
  gsl_matrix_memcpy(matrix, colli);

  /* Populate matrix with radiative transitions */
  for(li=0;li<md[ispec].nline;li++){
    k=md[ispec].lau[li];
    l=md[ispec].lal[li];
    gsl_matrix_set(matrix, k, k, gsl_matrix_get(matrix, k, k)+md[ispec].beinstu[li]*mp[ispec].jbar[li]+md[ispec].aeinst[li]);
    gsl_matrix_set(matrix, l, l, gsl_matrix_get(matrix, l, l)+md[ispec].beinstl[li]*mp[ispec].jbar[li]);
    gsl_matrix_set(matrix, k, l, gsl_matrix_get(matrix, k, l)-md[ispec].beinstl[li]*mp[ispec].jbar[li]);
    gsl_matrix_set(matrix, l, k, gsl_matrix_get(matrix, l, k)-md[ispec].beinstu[li]*mp[ispec].jbar[li]-md[ispec].aeinst[li]);
  }
}

/*....................................................................*/
void
lteOnePoint(molData *md, const int ispec, const double temp, double *pops){
  int ilev;
  double sum;

  sum = 0.0;
  for(ilev=0;ilev<md[ispec].nlev;ilev++){
    pops[ilev] = md[ispec].gstat[ilev]*exp(-HCKB*md[ispec].eterm[ilev]/temp);
    sum += pops[ilev];
  }
  for(ilev=0;ilev<md[ispec].nlev;ilev++)
    pops[ilev] /= sum;
}

/*....................................................................*/
void
LTE(configInfo *par, struct grid *gp, molData *md){
  int id,ispec;

  for(id=0;id<par->pIntensity;id++){
    for(ispec=0;ispec<par->nSpecies;ispec++){
      lteOnePoint(md, ispec, gp[id].t[0], gp[id].mol[ispec].pops);
    }
  }
  if(par->outputfile) popsout(par,gp,md);
}

/*....................................................................*/
void
solveStatEq(int id, struct grid *gp, molData *md, const int ispec, configInfo *par\
  , struct blendInfo blends, int nextMolWithBlend, gridPointData *mp\
  , double *halfFirstDs, _Bool *luWarningGiven){

  int t,s,iter,status;
  double *opop,*oopop,*tempNewPop=NULL;
  double diff;
  const double minpop_for_convergence_check = 1.e-6;
  char errStr[80];

  gsl_matrix *colli  = gsl_matrix_alloc(md[ispec].nlev, md[ispec].nlev);
  gsl_matrix *matrix = gsl_matrix_alloc(md[ispec].nlev, md[ispec].nlev);
  gsl_vector *newpop = gsl_vector_alloc(md[ispec].nlev);
  gsl_vector *rhVec  = gsl_vector_alloc(md[ispec].nlev);
  gsl_permutation *p = gsl_permutation_alloc (md[ispec].nlev);

  opop       = malloc(sizeof(*opop)      *md[ispec].nlev);
  oopop      = malloc(sizeof(*oopop)     *md[ispec].nlev);
  tempNewPop = malloc(sizeof(*tempNewPop)*md[ispec].nlev);

  for(t=0;t<md[ispec].nlev;t++){
    opop[t]=0.;
    oopop[t]=0.;
    gsl_vector_set(rhVec,t,0.);
  }
  gsl_vector_set(rhVec,md[ispec].nlev-1,1.);
  diff=1;
  iter=0;

  getFixedMatrix(md,ispec,gp,id,colli,par);

  while((diff>TOL && iter<MAXITER) || iter<5){
    updateJBar(id,md,gp,ispec,par,blends,nextMolWithBlend,mp,halfFirstDs);

    getMatrix(matrix,md,ispec,mp,colli);

    /* this could also be done in getFixedMatrix */ 
    for(s=0;s<md[ispec].nlev;s++){
      gsl_matrix_set(matrix,md[ispec].nlev-1,s,1.);
    }

    status = gsl_linalg_LU_decomp(matrix,p,&s);
    if(status){
      if(!silent){
        sprintf(errStr, "LU decomposition failed for point %d, iteration %d (GSL error %d).", id, iter, status);
        bail_out(errStr);
      }
      exit(1);
    }

    status = gsl_linalg_LU_solve(matrix,p,rhVec,newpop);
    if(status){
      if(!silent && !(*luWarningGiven)){
        *luWarningGiven = 1;
        sprintf(errStr, "LU solver failed for point %d, iteration %d (GSL error %d).", id, iter, status);
        warning(errStr);
        warning("Doing LSE for this point. NOTE that no further warnings will be issued.");
      }
      lteOnePoint(md, ispec, gp[id].t[0], tempNewPop);
      for(s=0;s<md[ispec].nlev;s++)
        gsl_vector_set(newpop,s,tempNewPop[s]);
    }

    diff=0.;
    for(t=0;t<md[ispec].nlev;t++){
      gsl_vector_set(newpop,t,gsl_max(gsl_vector_get(newpop,t),EPS));
      oopop[t]=opop[t];
      opop[t]=gp[id].mol[ispec].pops[t];

      /* Removed pragma omp critical; no two threads have the same cell id. */
      {
        gp[id].mol[ispec].pops[t]=gsl_vector_get(newpop,t);
      }

      if(gsl_min(gp[id].mol[ispec].pops[t],gsl_min(opop[t],oopop[t]))>minpop_for_convergence_check){
        diff=gsl_max(fabs(gp[id].mol[ispec].pops[t]-opop[t])/gp[id].mol[ispec].pops[t]\
            ,gsl_max(fabs(gp[id].mol[ispec].pops[t]-oopop[t])/gp[id].mol[ispec].pops[t],diff));
      }
    }
    iter++;
  }

  gsl_matrix_free(colli);
  gsl_matrix_free(matrix);
  gsl_vector_free(rhVec);
  gsl_vector_free(newpop);
  gsl_permutation_free(p);
  free(tempNewPop);
  free(opop);
  free(oopop);
}

/*....................................................................*/
void gridPopsInit(configInfo *par, molData *md, struct grid *gp){
  int i,id,ilev;

  for(i=0;i<par->nSpecies;i++){
    /* Allocate space for populations etc */
    for(id=0;id<par->ncell; id++){
      gp[id].mol[i].pops = malloc(sizeof(double)*md[i].nlev);
      for(ilev=0;ilev<md[i].nlev;ilev++)
        gp[id].mol[i].pops[ilev] = 0.0;
    }
  }
}

/*....................................................................*/
void
levelPops(molData *md, configInfo *par, struct grid *gp, int *popsdone, double *lamtab, double *kaptab, const int nEntries){
  int id,iter,ilev,ispec,c=0,n,i,threadI,nVerticesDone,nItersDone,nlinetot;
  double percent=0.,*median,result1=0,result2=0,snr,delta_pop;
  int nextMolWithBlend;
  struct statistics { double *pop, *ave, *sigma; } *stat;
  const gsl_rng_type *ranNumGenType = gsl_rng_ranlxs2;
  struct blendInfo blends;
  _Bool luWarningGiven=0;
  gsl_error_handler_t *defaultErrorHandler=NULL;
  int RNG_seeds[par->nThreads];

  nlinetot = 0;
  for(ispec=0;ispec<par->nSpecies;ispec++)
    nlinetot += md[ispec].nline;

  gridPopsInit(par,md,gp);

  if(par->lte_only){
    LTE(par,gp,md);
    if(par->outputfile) popsout(par,gp,md);

  }else{ /* Non-LTE */
    stat=malloc(sizeof(struct statistics)*par->pIntensity);

    /* Random number generator */
    gsl_rng *ran = gsl_rng_alloc(ranNumGenType);
#ifdef TEST
    gsl_rng_set(ran, 1237106) ;
#else 
    gsl_rng_set(ran,time(0));
#endif

    gsl_rng **threadRans;
    threadRans = malloc(sizeof(gsl_rng *)*par->nThreads);

    for (i=0;i<par->nThreads;i++){
      threadRans[i] = gsl_rng_alloc(ranNumGenType);
      if (par->resetRNG==1) RNG_seeds[i] = (int)(gsl_rng_uniform(ran)*1e6);
      else gsl_rng_set(threadRans[i],(int)(gsl_rng_uniform(ran)*1e6));
    }

    calcGridCollRates(par,md,gp);
    calcGridLinesDustOpacity(par, md, lamtab, kaptab, nEntries, gp);

    /* Check for blended lines */
    lineBlend(md, par, &blends);

    if(par->init_lte) LTE(par,gp,md);

    for(id=0;id<par->pIntensity;id++){
      stat[id].pop=malloc(sizeof(double)*md[0].nlev*5);
      stat[id].ave=malloc(sizeof(double)*md[0].nlev);
      stat[id].sigma=malloc(sizeof(double)*md[0].nlev);
      for(ilev=0;ilev<md[0].nlev;ilev++) {
        for(iter=0;iter<5;iter++) stat[id].pop[ilev+md[0].nlev*iter]=gp[id].mol[0].pops[ilev];
      }
    }

    if(par->outputfile) popsout(par,gp,md);

    /* Initialize convergence flag */
    for(id=0;id<par->ncell;id++){
      gp[id].conv=0;
    }

    defaultErrorHandler = gsl_set_error_handler_off();
    /*
This is done to allow proper handling of errors which may arise in the LU solver within solveStatEq(). It is done here because the GSL documentation does not recommend leaving the error handler at the default within multi-threaded code.

While this is off however, other gsl_* etc calls will not exit if they encounter a problem. We may need to pay some attention to trapping their errors.
    */

    nItersDone = par->nSolveItersDone;
    while(nItersDone < par->nSolveIters){ /* Not a 'for' loop because we will probably later want to add a convergence criterion. */
      if(!silent) progressbar2(par, 0, nItersDone, 0, result1, result2);

      for(id=0;id<par->pIntensity;id++){
        for(ilev=0;ilev<md[0].nlev;ilev++) {
          for(iter=0;iter<4;iter++) stat[id].pop[ilev+md[0].nlev*iter]=stat[id].pop[ilev+md[0].nlev*(iter+1)];
          stat[id].pop[ilev+md[0].nlev*4]=gp[id].mol[0].pops[ilev];
        }
      }
      calcGridMolSpecNumDens(par,md,gp);

      nVerticesDone=0;
      omp_set_dynamic(0);
#pragma omp parallel private(id,ispec,threadI,nextMolWithBlend) num_threads(par->nThreads)
      {
        threadI = omp_get_thread_num();

        if (par->resetRNG==1) gsl_rng_set(threadRans[threadI],RNG_seeds[threadI]);
        /* Declare and allocate thread-private variables */
        gridPointData *mp;	// Could have declared them earlier
        double *halfFirstDs;	// and included them in private() I guess.
        mp=malloc(sizeof(gridPointData)*par->nSpecies);

#pragma omp for
        for(id=0;id<par->pIntensity;id++){
#pragma omp atomic
          ++nVerticesDone;

          for (ispec=0;ispec<par->nSpecies;ispec++){
            mp[ispec].jbar = malloc(sizeof(double)*md[ispec].nline);
            mp[ispec].phot = malloc(sizeof(double)*md[ispec].nline*gp[id].nphot);
            mp[ispec].vfac = malloc(sizeof(double)*                gp[id].nphot);
            mp[ispec].vfac_loc = malloc(sizeof(double)*            gp[id].nphot);
          }
          halfFirstDs = malloc(sizeof(*halfFirstDs)*gp[id].nphot);

          if (threadI == 0){ /* i.e., is master thread. */
            if(!silent) progressbar(nVerticesDone/(double)par->pIntensity,10);
          }
          if(gp[id].dens[0] > 0 && gp[id].t[0] > 0){
            calculateJBar(id,gp,md,threadRans[threadI],par,nlinetot,blends,mp,halfFirstDs);
            nextMolWithBlend = 0;
            for(ispec=0;ispec<par->nSpecies;ispec++){
              solveStatEq(id,gp,md,ispec,par,blends,nextMolWithBlend,mp,halfFirstDs,&luWarningGiven);
              if(par->blend && blends.mols!=NULL && ispec==blends.mols[nextMolWithBlend].molI)
                nextMolWithBlend++;
            }
          }
          if (threadI == 0){ /* i.e., is master thread */
            if(!silent) warning("");
          }
          freeGridPointData(par->nSpecies, mp);
          free(halfFirstDs);
        }
        free(mp);
      } /* end parallel block. */

      for(id=0;id<par->pIntensity;id++){
        snr=0;
        n=0;
        for(ilev=0;ilev<md[0].nlev;ilev++) {
          stat[id].ave[ilev]=0;
          for(iter=0;iter<5;iter++) stat[id].ave[ilev]+=stat[id].pop[ilev+md[0].nlev*iter];
          stat[id].ave[ilev]=stat[id].ave[ilev]/5.;
          stat[id].sigma[ilev]=0;
          for(iter=0;iter<5;iter++) {
            delta_pop = stat[id].pop[ilev+md[0].nlev*iter]-stat[id].ave[ilev];
            stat[id].sigma[ilev]+=delta_pop*delta_pop;
          }
          stat[id].sigma[ilev]=sqrt(stat[id].sigma[ilev])/5.;
          if(gp[id].mol[0].pops[ilev] > 1e-12) c++;

          if(gp[id].mol[0].pops[ilev] > 1e-12 && stat[id].sigma[ilev] > 0.){
            snr+=gp[id].mol[0].pops[ilev]/stat[id].sigma[ilev];
            n++;
          }
        }
        if(n>0) snr=snr/n;
        else if(n==0) snr=1e6;
        if(snr > 3.) gp[id].conv=2;
        if(snr <= 3 && gp[id].conv==2) gp[id].conv=1;
      }

      median=malloc(sizeof(*median)*gsl_max(c,1));
      c=0;
      for(id=0;id<par->pIntensity;id++){
        for(ilev=0;ilev<md[0].nlev;ilev++){
          if(gp[id].mol[0].pops[ilev] > 1e-12) median[c++]=gp[id].mol[0].pops[ilev]/stat[id].sigma[ilev];
        }
      }

      gsl_sort(median, 1, c);
      if(nItersDone>par->nSolveItersDone+1){
        result1=median[0];
        result2 =gsl_stats_median_from_sorted_data(median, 1, c);
      }
      free(median);

      if(!silent) progressbar2(par, 1, nItersDone, percent, result1, result2);
      if(par->outputfile != NULL) popsout(par,gp,md);
      nItersDone++;
    }
    gsl_set_error_handler(defaultErrorHandler);

    freeMolsWithBlends(blends.mols, blends.numMolsWithBlends);

    for (i=0;i<par->nThreads;i++){
      gsl_rng_free(threadRans[i]);
    }
    free(threadRans);
    gsl_rng_free(ran);

    for(id=0;id<par->pIntensity;id++){
      free(stat[id].pop);
      free(stat[id].ave);
      free(stat[id].sigma);
    }
    free(stat);
  }

  par->dataFlags |= (1 << DS_bit_populations);

  if(par->binoutputfile != NULL) binpopsout(par,gp,md);

  *popsdone=1;
}


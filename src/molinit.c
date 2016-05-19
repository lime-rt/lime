/*
 *  molinit.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015 The LIME development team
 *
 */

#include "lime.h"

void
kappa(molData *md, struct grid *gp, inputPars *par, int s){
  FILE *fp;
  char string[80];
  int i=0,k,j,iline,id;
  double loglam, *lamtab, *kaptab, *kappatab, gtd;
  gsl_spline *spline;

  kappatab   	 = malloc(sizeof(*kappatab)*md[s].nline);
  md[s].cmb	 = malloc(sizeof(double)*md[s].nline);
  md[s].local_cmb = malloc(sizeof(double)*md[s].nline);

  if(par->dust == NULL){
    for(i=0;i<md[s].nline;i++) kappatab[i]=0.;
  } else {
    gsl_interp_accel *acc=gsl_interp_accel_alloc();
    if((fp=fopen(par->dust, "r"))==NULL){
      if(!silent) bail_out("Error opening dust opacity data file!");
      exit(1);
    }
    while(fgetc(fp) != EOF){
      fgets(string,80,fp);
      i++;
    }
    rewind(fp);
    if(i>0){
      lamtab=malloc(sizeof(*lamtab)*i);
      kaptab=malloc(sizeof(*kaptab)*i);
    } else {
      if(!silent) bail_out("No opacities read");
      exit(1);
    }
    for(k=0;k<i;k++){
      fscanf(fp,"%lf %lf\n", &lamtab[k], &kaptab[k]);
      lamtab[k]=log10(lamtab[k]/1e6);
      kaptab[k]=log10(kaptab[k]);
    }
    fclose(fp);
    spline=gsl_spline_alloc(gsl_interp_cspline,i);
    gsl_spline_init(spline,lamtab,kaptab,i);
    for(j=0;j<md[s].nline;j++) {
      loglam=log10(CLIGHT/md[s].freq[j]);
      if(loglam < lamtab[0]){
        kappatab[j]=0.1*pow(10.,kaptab[0] + (loglam-lamtab[0]) * (kaptab[1]-kaptab[0])/(lamtab[1]-lamtab[0]));
      } else if(loglam > lamtab[i-1]){
        kappatab[j]=0.1*pow(10.,kaptab[i-2] + (loglam-lamtab[i-2]) * (kaptab[i-1]-kaptab[i-2])/(lamtab[i-1]-lamtab[i-2]));
      } else kappatab[j]=0.1*pow(10.,gsl_spline_eval(spline,loglam,acc));
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    free(kaptab);
    free(lamtab);
  }

  for(iline=0;iline<md[s].nline;iline++){
    for(id=0;id<par->ncell;id++){
      gasIIdust(gp[id].x[0],gp[id].x[1],gp[id].x[2],&gtd);
      gp[id].mol[s].knu[iline]=kappatab[iline]*2.4*AMU/gtd*gp[id].dens[0];
      //Check if input model supplies a dust temperature. Otherwise use the kinetic temperature
      if(gp[id].t[1]==-1) {
        gp[id].mol[s].dust[iline]=planckfunc(iline,gp[id].t[0],md,s);
      } else {
        gp[id].mol[s].dust[iline]=planckfunc(iline,gp[id].t[1],md,s);
      }
    }
    // fix the normalization at 230GHz
    md[s].norm=planckfunc(0,par->tcmb,md,0);
    md[s].norminv=1./md[s].norm;
    if(par->tcmb>0.) md[s].cmb[iline]=planckfunc(iline,par->tcmb,md,s)/md[s].norm;
    else md[s].cmb[iline]=0.;
    md[s].local_cmb[iline]=planckfunc(iline,2.728,md,s)/md[s].norm;
  }

  free(kappatab);
  return;
}


double
planckfunc(int iline, double temp, molData *md, int s){
  double bb=10.,wn;
  if(temp<eps) bb = 0.0;
  else {
    wn=md[s].freq[iline]/CLIGHT;
    if (HPLANCK*md[s].freq[iline]>100.*KBOLTZ*temp) 
      bb=2.*HPLANCK*wn*wn*md[s].freq[iline]*exp(-HPLANCK*md[s].freq[iline]/KBOLTZ/temp);
    else 
      bb=2.*HPLANCK*wn*wn*md[s].freq[iline]/(exp(HPLANCK*md[s].freq[iline]/KBOLTZ/temp)-1);
  }
  return bb;
}

void
molinit(molData *md, inputPars *par, struct grid *gp, int i){
  int id, ilev, iline, itrans, ispec, itemp, *ntemp, tnint=-1, idummy, ipart, *count,flag=0;
  char *collpartname[] = {"H2","p-H2","o-H2","electrons","H","He","H+"}; /* definition from LAMDA */
  double fac, uprate, downrate=0, dummy, amass;
  struct data { double *colld, *temp; } *part;

  char string[200], specref[90], partstr[90];
  FILE *fp;

  if((fp=fopen(par->moldatfile[i], "r"))==NULL) {
    if(!silent) bail_out("Error opening molecular data file");
    exit(1);
  }

  /* Read the header of the data file */
  fgets(string, 80, fp);
  fgets(specref, 90, fp);
  fgets(string, 80, fp);
  fscanf(fp, "%lf\n", &amass);
  fgets(string, 80, fp);
  fscanf(fp, "%d\n", &md[i].nlev);
  fgets(string, 80, fp);

  md[i].eterm=malloc(sizeof(double)*md[i].nlev);
  md[i].gstat=malloc(sizeof(double)*md[i].nlev);

  /* Read the level energies and statistical weights */
  for(ilev=0;ilev<md[i].nlev;ilev++){
    fscanf(fp, "%d %lf %lf", &idummy, &md[i].eterm[ilev], &md[i].gstat[ilev]);
    fgets(string, 80, fp);
  }

  /* Read the number of transitions and allocate array space */
  fgets(string, 80, fp);
  fscanf(fp, "%d\n", &md[i].nline);
  fgets(string, 80, fp);

  md[i].lal     = malloc(sizeof(int)*md[i].nline);
  md[i].lau     = malloc(sizeof(int)*md[i].nline);
  md[i].aeinst  = malloc(sizeof(double)*md[i].nline);
  md[i].freq    = malloc(sizeof(double)*md[i].nline);
  md[i].beinstu = malloc(sizeof(double)*md[i].nline);
  md[i].beinstl = malloc(sizeof(double)*md[i].nline);

  /* Read transitions, Einstein A, and frequencies */
  for(iline=0;iline<md[i].nline;iline++){
    fscanf(fp, "%d %d %d %lf %lf %lf\n", &idummy, &md[i].lau[iline], &md[i].lal[iline], &md[i].aeinst[iline], &md[i].freq[iline], &dummy);
    md[i].freq[iline]*=1e9;
    md[i].lau[iline]-=1;
    md[i].lal[iline]-=1;
  }

  /* Calculate Einsten B's */
  for(iline=0;iline<md[i].nline;iline++){
    /*		md[i].freq[iline]=(md[i].eterm[md[i].lau[iline]]-md[i].eterm[md[i].lal[iline]])*100*CLIGHT; */
    md[i].beinstu[iline]=md[i].aeinst[iline]*(CLIGHT/md[i].freq[iline])*(CLIGHT/md[i].freq[iline])/(HPLANCK*md[i].freq[iline])/2.;
    md[i].beinstl[iline]=md[i].gstat[md[i].lau[iline]]/md[i].gstat[md[i].lal[iline]]*md[i].beinstu[iline];
  }

  /* Calculate Doppler and thermal line broadening */
  amass*=AMU;
  for(id=0;id<par->ncell;id++) {
    gp[id].mol[i].dopb=sqrt(gp[id].dopb*gp[id].dopb+2.*KBOLTZ/amass*gp[id].t[0]);
    gp[id].mol[i].binv=1./gp[id].mol[i].dopb;
  }

  /* Collision rates below here */
  if(par->lte_only==0){
    fgets(string, 80, fp);
    fscanf(fp,"%d\n", &md[i].npart);
    count=malloc(sizeof(*count)*md[i].npart);
    /* collision partner sanity check */

    if(md[i].npart > par->collPart) flag=1;
    if(md[i].npart < par->collPart){
      if(!silent) bail_out("Too many density profiles defined");
      exit(1);
    }


    md[i].ntrans = malloc(sizeof(int)*md[i].npart);
    ntemp = malloc(sizeof(*ntemp)*md[i].npart);
    part = malloc(sizeof(struct data) * md[i].npart);

    for(ipart=0;ipart<md[i].npart;ipart++){
      fgets(string, 80, fp);
      fscanf(fp,"%d\n", &count[ipart]);
      fgets(string, 80, fp);
      fgets(string, 80, fp);
      fscanf(fp,"%d\n", &md[i].ntrans[ipart]);
      fgets(string, 80, fp);
      fscanf(fp,"%d\n", &ntemp[ipart]);
      fgets(string, 80, fp);

      part[ipart].temp=malloc(sizeof(double)*ntemp[ipart]);

      if(ipart==0){
        md[i].lcl = malloc(sizeof(int)*md[i].ntrans[ipart]);
        md[i].lcu = malloc(sizeof(int)*md[i].ntrans[ipart]);
      }

      for(itemp=0;itemp<ntemp[ipart];itemp++){
        fscanf(fp, "%lf", &part[ipart].temp[itemp]);
      }

      fscanf(fp,"\n");
      fgets(string, 80, fp);

      part[ipart].colld=malloc(sizeof(double)*md[i].ntrans[ipart]*ntemp[ipart]);

      for(itrans=0;itrans<md[i].ntrans[ipart];itrans++){
        fscanf(fp, "%d %d %d", &idummy, &md[i].lcu[itrans], &md[i].lcl[itrans]);
        md[i].lcu[itrans]-=1;
        md[i].lcl[itrans]-=1;
        for(itemp=0;itemp<ntemp[ipart];itemp++){
          fscanf(fp, "%lf", &part[ipart].colld[itrans*ntemp[ipart]+itemp]);
          part[ipart].colld[itrans*ntemp[ipart]+itemp]/=1.e6;
        }
        fscanf(fp,"\n");
      }
    }

    /* Print out collision partner information */
    strcpy(partstr, collpartname[count[0]-1]);
    for(ipart=1;ipart<md[i].npart;ipart++){
      strcat( partstr, ", ");
      strcat( partstr, collpartname[count[ipart]-1]);
    }
    if(!silent) {
      collpartmesg(specref, md[i].npart);
      collpartmesg2(partstr, ipart);
      collpartmesg3(par->collPart, flag);
    }

    /* Calculate molecular density */
    for(id=0;id<par->ncell; id++){
      for(ispec=0;ispec<par->nSpecies;ispec++){
        if(md[i].npart == 1 && (count[0] == 1 || count[0] == 2 || count[0] == 3)){
          gp[id].nmol[ispec]=gp[id].abun[ispec]*gp[id].dens[0];
        } else if(md[i].npart == 2 && (count[0] == 2 || count[0] == 3) && (count[1] == 2 || count[1] == 3)){
          if(!flag){
            gp[id].nmol[ispec]=gp[id].abun[ispec]*(gp[id].dens[0]+gp[id].dens[1]);
          } else {
            gp[id].nmol[ispec]=gp[id].abun[ispec]*gp[id].dens[0];
            if(!silent) warning("Calculating molecular density with respect to first collision partner only");
          }
        } else if(md[i].npart > 2 && !flag){
          gp[id].nmol[ispec]=gp[id].abun[ispec]*(gp[id].dens[0]+gp[id].dens[1]);
          if(!silent) warning("Calculating molecular density with respect first and second collision partner");
        }
      }
    }

    for(id=0;id<par->ncell;id++){
      gp[id].mol[i].partner=malloc(sizeof(struct rates)*md[i].npart);
      for(ipart=0;ipart<md[i].npart;ipart++){
        gp[id].mol[i].partner[ipart].up = malloc(sizeof(double)*md[i].ntrans[ipart]);
        gp[id].mol[i].partner[ipart].down = malloc(sizeof(double)*md[i].ntrans[ipart]);
      }
    }

    for(id=0;id<par->ncell;id++){
      for(ipart=0;ipart<md[i].npart;ipart++){
        for(itrans=0;itrans<md[i].ntrans[ipart];itrans++){
          if((gp[id].t[0]>part[ipart].temp[0])&&(gp[id].t[0]<part[ipart].temp[ntemp[ipart]-1])){
            for(itemp=0;itemp<ntemp[ipart]-1;itemp++){
              if((gp[id].t[0]>part[ipart].temp[itemp])&&(gp[id].t[0]<=part[ipart].temp[itemp+1])){
                tnint=itemp;
              }
            }
            fac=(gp[id].t[0]-part[ipart].temp[tnint])/(part[ipart].temp[tnint+1]-part[ipart].temp[tnint]);
            downrate=part[ipart].colld[itrans*ntemp[ipart]+tnint]+fac*(part[ipart].colld[itrans*ntemp[ipart]+tnint+1]-part[ipart].colld[itrans*ntemp[ipart]+tnint]);
          } else {
            if(gp[id].t[0]<=part[ipart].temp[0]) downrate=part[ipart].colld[itrans*ntemp[ipart]];
            if(gp[id].t[0]>=part[ipart].temp[ntemp[ipart]-1]) downrate=part[ipart].colld[itrans*ntemp[ipart]+ntemp[ipart]-1];
          }
          uprate=md[i].gstat[md[i].lcu[itrans]]/md[i].gstat[md[i].lcl[itrans]]*downrate*exp(-HCKB*(md[i].eterm[md[i].lcu[itrans]]-md[i].eterm[md[i].lcl[itrans]])/gp[id].t[0]);
          gp[id].mol[i].partner[ipart].up[itrans]=uprate;
          gp[id].mol[i].partner[ipart].down[itrans]=downrate;
        }
      }
    }
    for(ipart=0;ipart<md[i].npart;ipart++){
      free(part[ipart].colld);
      free(part[ipart].temp);
    }
    free(ntemp);
    free(part);
    free(count);
  }
  /* End of collision rates */

  fclose(fp);

  /* Allocate space for populations and opacities */
  for(id=0;id<par->ncell; id++){
    gp[id].mol[i].pops = malloc(sizeof(double)*md[i].nlev);
    gp[id].mol[i].dust = malloc(sizeof(double)*md[i].nline);
    gp[id].mol[i].knu  = malloc(sizeof(double)*md[i].nline);
    for(ilev=0;ilev<md[i].nlev;ilev++) gp[id].mol[i].pops[ilev]=0.0;
  }

  /* Get dust opacities */
  kappa(md,gp,par,i);
}



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
kappa(molData *m, struct grid *g, inputPars *par, int s){
  FILE *fp;
  char string[80];
  int i=0,k,j,iline,id;
  double loglam, *lamtab, *kaptab, *kappatab, gtd;
  gsl_spline *spline;

  kappatab   	 = malloc(sizeof(*kappatab)*m[s].nline);
  m[s].cmb	 = malloc(sizeof(double)*m[s].nline);
  m[s].local_cmb = malloc(sizeof(double)*m[s].nline);

  if(par->dust == NULL){
    for(i=0;i<m[s].nline;i++) kappatab[i]=0.;
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
    for(j=0;j<m[s].nline;j++) {
      loglam=log10(CLIGHT/m[s].freq[j]);
      if(loglam < lamtab[0]){
        kappatab[j]=0.1*pow(10.,kaptab[0] + (loglam-lamtab[0]) * (kaptab[1]-kaptab[0])/(lamtab[1]-lamtab[0]));
      } else if(loglam > lamtab[i-1]){
        kappatab[j]=0.1*pow(10.,kaptab[i-2] + (loglam-lamtab[i-2]) * (kaptab[i-1]-kaptab[i-2])/(lamtab[i-1]-lamtab[i-2]));
      } else
        kappatab[j]=0.1*pow(10.,gsl_spline_eval(spline,loglam,acc));
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    free(kaptab);
    free(lamtab);
  }

  for(iline=0;iline<m[s].nline;iline++){
    for(id=0;id<par->ncell;id++){
      gasIIdust(g[id].x[0],g[id].x[1],g[id].x[2],&gtd);
      g[id].mol[s].knu[iline]=kappatab[iline]*2.4*AMU/gtd*g[id].dens[0];

      /* Check if input model supplies a dust temperature. Otherwise use the kinetic temperature. */
      if(g[id].t[1]==-1) {
        g[id].mol[s].dust[iline]=planckfunc(iline,g[id].t[0],m,s);
      } else {
        g[id].mol[s].dust[iline]=planckfunc(iline,g[id].t[1],m,s);
      }
    }
    /* Fix the normalization at 230GHz. */
    m[s].norm=planckfunc(0,par->tcmb,m,0);
    m[s].norminv=1./m[s].norm;
    if(par->tcmb>0.) m[s].cmb[iline]=planckfunc(iline,par->tcmb,m,s)/m[s].norm;
    else m[s].cmb[iline]=0.;
    m[s].local_cmb[iline]=planckfunc(iline,2.728,m,s)/m[s].norm;
  }

  free(kappatab);
  return;
}


double
planckfunc(int iline, double temp, molData *m,int s){
  double bb=10.,wn;
  if(temp<eps) bb = 0.0;
  else {
    wn=m[s].freq[iline]/CLIGHT;
    if (HPLANCK*m[s].freq[iline]>100.*KBOLTZ*temp) 
      bb=2.*HPLANCK*wn*wn*m[s].freq[iline]*exp(-HPLANCK*m[s].freq[iline]/KBOLTZ/temp);
    else 
      bb=2.*HPLANCK*wn*wn*m[s].freq[iline]/(exp(HPLANCK*m[s].freq[iline]/KBOLTZ/temp)-1);
  }
  return bb;
}

void
molinit(molData *m, inputPars *par, struct grid *g,int i){
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
  fscanf(fp, "%d\n", &m[i].nlev);
  fgets(string, 80, fp);

  m[i].eterm=malloc(sizeof(double)*m[i].nlev);
  m[i].gstat=malloc(sizeof(double)*m[i].nlev);

  /* Read the level energies and statistical weights */
  for(ilev=0;ilev<m[i].nlev;ilev++){
    fscanf(fp, "%d %lf %lf", &idummy, &m[i].eterm[ilev], &m[i].gstat[ilev]);
    fgets(string, 80, fp);
  }

  /* Read the number of transitions and allocate array space */
  fgets(string, 80, fp);
  fscanf(fp, "%d\n", &m[i].nline);
  fgets(string, 80, fp);

  m[i].lal     = malloc(sizeof(int)*m[i].nline);
  m[i].lau     = malloc(sizeof(int)*m[i].nline);
  m[i].aeinst  = malloc(sizeof(double)*m[i].nline);
  m[i].freq    = malloc(sizeof(double)*m[i].nline);
  m[i].beinstu = malloc(sizeof(double)*m[i].nline);
  m[i].beinstl = malloc(sizeof(double)*m[i].nline);

  /* Read transitions, Einstein A, and frequencies */
  for(iline=0;iline<m[i].nline;iline++){
    fscanf(fp, "%d %d %d %lf %lf %lf\n", &idummy, &m[i].lau[iline], &m[i].lal[iline], &m[i].aeinst[iline], &m[i].freq[iline], &dummy);
    m[i].freq[iline]*=1e9;
    m[i].lau[iline]-=1;
    m[i].lal[iline]-=1;
  }

  /* Calculate Einsten B's */
  for(iline=0;iline<m[i].nline;iline++){
    /*		m[i].freq[iline]=(m[i].eterm[m[i].lau[iline]]-m[i].eterm[m[i].lal[iline]])*100*CLIGHT; */
    m[i].beinstu[iline]=m[i].aeinst[iline]*(CLIGHT/m[i].freq[iline])*(CLIGHT/m[i].freq[iline])/(HPLANCK*m[i].freq[iline])/2.;
    m[i].beinstl[iline]=m[i].gstat[m[i].lau[iline]]/m[i].gstat[m[i].lal[iline]]*m[i].beinstu[iline];
  }

  /* Calculate Doppler and thermal line broadening */
  amass*=AMU;
  for(id=0;id<par->ncell;id++) {
    g[id].mol[i].dopb=sqrt(g[id].dopb*g[id].dopb+2.*KBOLTZ/amass*g[id].t[0]);
    g[id].mol[i].binv=1./g[id].mol[i].dopb;
  }

  /* Collision rates below here */
  if(par->lte_only==0){
    fgets(string, 80, fp);
    fscanf(fp,"%d\n", &m[i].npart);
    count=malloc(sizeof(*count)*m[i].npart);
    /* collision partner sanity check */

    if(m[i].npart > par->collPart) flag=1;
    if(m[i].npart < par->collPart){
      if(!silent) bail_out("Error: Too many density profiles defined");
      exit(1);
    }


    m[i].ntrans = malloc(sizeof(int)*m[i].npart);
    ntemp = malloc(sizeof(*ntemp)*m[i].npart);
    part = malloc(sizeof(struct data) * m[i].npart);

    for(ipart=0;ipart<m[i].npart;ipart++){
      fgets(string, 80, fp);
      fscanf(fp,"%d\n", &count[ipart]);
      fgets(string, 80, fp);
      fgets(string, 80, fp);
      fscanf(fp,"%d\n", &m[i].ntrans[ipart]);
      fgets(string, 80, fp);
      fscanf(fp,"%d\n", &ntemp[ipart]);
      fgets(string, 80, fp);

      part[ipart].temp=malloc(sizeof(double)*ntemp[ipart]);

      if(ipart==0){
        m[i].lcl = malloc(sizeof(int)*m[i].ntrans[ipart]);
        m[i].lcu = malloc(sizeof(int)*m[i].ntrans[ipart]);
      }

      for(itemp=0;itemp<ntemp[ipart];itemp++){
        fscanf(fp, "%lf", &part[ipart].temp[itemp]);
      }

      fscanf(fp,"\n");
      fgets(string, 80, fp);

      part[ipart].colld=malloc(sizeof(double)*m[i].ntrans[ipart]*ntemp[ipart]);

      for(itrans=0;itrans<m[i].ntrans[ipart];itrans++){
        fscanf(fp, "%d %d %d", &idummy, &m[i].lcu[itrans], &m[i].lcl[itrans]);
        m[i].lcu[itrans]-=1;
        m[i].lcl[itrans]-=1;
        for(itemp=0;itemp<ntemp[ipart];itemp++){
          fscanf(fp, "%lf", &part[ipart].colld[itrans*ntemp[ipart]+itemp]);
          part[ipart].colld[itrans*ntemp[ipart]+itemp]/=1.e6;
        }
        fscanf(fp,"\n");
      }
    }
    fclose(fp);

    /* Print out collision partner information */
    strcpy(partstr, collpartname[count[0]-1]);
    for(ipart=1;ipart<m[i].npart;ipart++){
      strcat( partstr, ", ");
      strcat( partstr, collpartname[count[ipart]-1]);
    }
    if(!silent) {
      collpartmesg(specref, m[i].npart);
      collpartmesg2(partstr, ipart);
      collpartmesg3(par->collPart, flag);
    }

    /* Calculate molecular density */
    for(id=0;id<par->ncell; id++){
      for(ispec=0;ispec<par->nSpecies;ispec++){
        if(m[i].npart == 1 && (count[0] == 1 || count[0] == 2 || count[0] == 3)){
          g[id].nmol[ispec]=g[id].abun[ispec]*g[id].dens[0];
        } else if(m[i].npart == 2 && (count[0] == 2 || count[0] == 3) && (count[1] == 2 || count[1] == 3)){
          if(!flag){
            g[id].nmol[ispec]=g[id].abun[ispec]*(g[id].dens[0]+g[id].dens[1]);
          } else {
            g[id].nmol[ispec]=g[id].abun[ispec]*g[id].dens[0];
            if(!silent) warning("Calculating molecular density with respect to first collision partner only");
          }
        } else if(m[i].npart > 2 && !flag){
          g[id].nmol[ispec]=g[id].abun[ispec]*(g[id].dens[0]+g[id].dens[1]);
          if(!silent) warning("Calculating molecular density with respect first and second collision partner");
        }
      }
    }

    for(id=0;id<par->ncell;id++){
      g[id].mol[i].partner=malloc(sizeof(struct rates)*m[i].npart);
      for(ipart=0;ipart<m[i].npart;ipart++){
        g[id].mol[i].partner[ipart].up = malloc(sizeof(double)*m[i].ntrans[ipart]);
        g[id].mol[i].partner[ipart].down = malloc(sizeof(double)*m[i].ntrans[ipart]);
      }
    }

    for(id=0;id<par->ncell;id++){
      for(ipart=0;ipart<m[i].npart;ipart++){
        for(itrans=0;itrans<m[i].ntrans[ipart];itrans++){
          if((g[id].t[0]>part[ipart].temp[0])&&(g[id].t[0]<part[ipart].temp[ntemp[ipart]-1])){
            for(itemp=0;itemp<ntemp[ipart]-1;itemp++){
              if((g[id].t[0]>part[ipart].temp[itemp])&&(g[id].t[0]<=part[ipart].temp[itemp+1])){
                tnint=itemp;
              }
            }
            fac=(g[id].t[0]-part[ipart].temp[tnint])/(part[ipart].temp[tnint+1]-part[ipart].temp[tnint]);
            downrate=part[ipart].colld[itrans*ntemp[ipart]+tnint]+fac*(part[ipart].colld[itrans*ntemp[ipart]+tnint+1]-part[ipart].colld[itrans*ntemp[ipart]+tnint]);
          } else {
            if(g[id].t[0]<=part[ipart].temp[0]) downrate=part[ipart].colld[itrans*ntemp[ipart]];
            if(g[id].t[0]>=part[ipart].temp[ntemp[ipart]-1]) downrate=part[ipart].colld[itrans*ntemp[ipart]+ntemp[ipart]-1];
          }
          uprate=m[i].gstat[m[i].lcu[itrans]]/m[i].gstat[m[i].lcl[itrans]]*downrate*exp(-HCKB*(m[i].eterm[m[i].lcu[itrans]]-m[i].eterm[m[i].lcl[itrans]])/g[id].t[0]);
          g[id].mol[i].partner[ipart].up[itrans]=uprate;
          g[id].mol[i].partner[ipart].down[itrans]=downrate;
        }
      }
    }
    for(ipart=0;ipart<m[i].npart;ipart++){
      free(part[ipart].colld);
      free(part[ipart].temp);
    }
    free(ntemp);
    free(part);
    free(count);
  }
  /* End of collision rates */

  /* Allocate space for populations and opacities */
  for(id=0;id<par->ncell; id++){
    g[id].mol[i].pops = malloc(sizeof(double)*m[i].nlev);
    g[id].mol[i].dust = malloc(sizeof(double)*m[i].nline);
    g[id].mol[i].knu  = malloc(sizeof(double)*m[i].nline);
    for(ilev=0;ilev<m[i].nlev;ilev++) g[id].mol[i].pops[ilev]=0.0;
  }

  /* Get dust opacities */
  kappa(m,g,par,i);
}



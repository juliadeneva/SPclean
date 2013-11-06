#include "spclean.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


const int IDMCUTOFF_AND = 50;
const float TBINW_AND = 1.0;
const int MAXPENALTY = 2;
extern int NFILES;

void ANDclean(char** infiles, char** outfiles, char* stampfile, char* dfile, float dt, int ndm, float tobs, int mode)
{
  FILE* delayfile, *pulsefile, *outfile, *stfile;
  //char outfilename[128];
  int delaysamp[ndm];
  float ttol = TBINW_AND/2.0;
  int nbins = floor(tobs/TBINW_AND)+1;
  float snrsum[NFILES][nbins], t_corr;
  int penalties[NFILES][nbins], penalty;
  long int npulses[NFILES], npulses_good[NFILES];
  double stamp_corrections[NFILES], stampref;
  Pulsus p;
  int i,j,k,idummy;

  fprintf(stderr,"ANDclean starting\n");

  //Zero out arrays for penalties and SNR sum bins
  for(j=0;j<NFILES;j++)
    {
      for(i=0;i<nbins;i++)
	{
	  snrsum[j][i] = 0.0;
	  penalties[j][i] = 0;
	}
      stamp_corrections[j] = 0;
    }
  
  
  fprintf(stderr,"Reading delay file %s\n",dfile);
  //Read the delay file
  if((delayfile = fopen(dfile,"r")) == NULL)
    {
      fprintf(stderr,"Couldn't open delay file %s, exiting.\n",dfile);
      exit(1);
    }

  for(i=0;i<ndm;i++)
    fscanf(delayfile,"%d %d %d",&idummy,&idummy,&delaysamp[i]);
  fclose(delayfile);


  //Read the file with beam timestamps
  if((stfile = fopen(stampfile,"r")) == NULL)
    fprintf(stderr,"Couldn't open timestamp file %s, assuming equal beam time stamps.\n",stampfile);
  else
    {
      fprintf(stderr,"stampfile: %s\n",stampfile);
      stampref = -1;
      i=0;
      for(i=0;i<NFILES;i++)
	{
	  fscanf(stfile, "%lf",&stamp_corrections[i]);
	  //fprintf(stderr,"%d %lf\n",i,stamp_corrections[i]);
	  
	  if(stamp_corrections[i] > 0 && stampref < 0)
	    stampref = stamp_corrections[i];
	  if(stamp_corrections[i] != 0)
	    stamp_corrections[i] = (stamp_corrections[i] - stampref)*24*3600;
	  fprintf(stderr,"beam: %d stamp correction: %f\n",i,stamp_corrections[i]);
	}
      fclose(stfile);
    }

  //fprintf(stderr,"Before SNR hist loop\n");
  //Read pulse files and construct SNR sums binned in time
  for(i=0;i<NFILES;i++)
    {
      //printf("Opening infile: %s \n",infiles[i]);
      if((pulsefile = fopen(infiles[i],"r")) == NULL)
	{
	  fprintf(stderr,"Couldn't open pulse file %s, skipping.\n",infiles[i]);
	  continue;
	}

      npulses[i] = 0;
      while(!feof(pulsefile))
	{
	  if(mode == 1) //For sift.dat files
	    fscanf(pulsefile,"%d %d %d %d %f %f %f %f %d %f %f",&p.idm,&p.ns,&p.nsm,&p.isamp,&p.snr,&p.mean,&p.rms,&p.suma,&p.nsum,&p.wgp,&p.diff);
	  else if(mode == 0) //For match.best files
	    fscanf(pulsefile,"%d %d %d %d %f %f %f",&p.idm,&p.ns,&p.nsm,&p.isamp,&p.snr,&p.mean,&p.rms);
	  else
	    {
	      fprintf(stderr,"Unrecognized mode: %d\n",mode);
	      fprintf(stderr,"Use 0 for match.best files or 1 for sift.dat files\n");
	      exit(1);
	    }
	      
	  //printf("idm: %d isamp: %d\n",p.idm,p.isamp);
	  npulses[i]++;

	  if(p.idm < IDMCUTOFF_AND)
	    {
	      p.isamp_corr = p.isamp - (int)(delaysamp[p.idm]/2.0+0.5) + (int)(stamp_corrections[i]/dt);
	      t_corr = p.isamp_corr*dt;
	      j = (int)floor(t_corr/TBINW_AND);
	      if(j > nbins-1)
		{
		  fprintf(stderr,"Trying to access %d-th of %d nbins, make sure Tobs is correct in caller.\n",j+1,nbins);
		  //exit(1);
		}
	      else
		snrsum[i][j] += p.snr;
	    }
	}
      fclose(pulsefile);
    }

  //fprintf(stderr,"Before penalties loop\n");
  //Assign penalties to time bins depending on how many beams an event was
  //detected in and how adjacent they are
  for(i=0;i<NFILES;i++)
    for(j=i+1;j<NFILES;j++)
      {
	if(i == 0)
	  penalty = 1;
	else if(j-i > 3)
	  penalty = 6-(j-i);
	else
	  penalty = j-i;

	//fprintf(stderr,"i: %d j: %d penalty: %d\n",i,j,penalty);

	for(k=0;k<nbins;k++)
	  if(snrsum[i][k]>0 && snrsum[j][k]>0)
	    {
	      penalties[i][k] = penalties[i][k] + penalty;
	      penalties[j][k] = penalties[j][k] + penalty;
	    }
	 
      }

  //fprintf(stderr,"Before output loop\n");
  //Read pulse files again, and zap events within 'time zones' with high penalties
  for(i=0;i<NFILES;i++)
    {
      printf("Opening infile: %s \n",infiles[i]);
      if((pulsefile = fopen(infiles[i],"r")) == NULL)
	{
	  fprintf(stderr,"Couldn't open pulse file %s, skipping.\n",infiles[i]);
	  continue;
	}

      //sprintf(outfilename,"%s%s",infiles[i],".andclean");
      if((outfile = fopen(outfiles[i],"w")) == NULL)
	{
	  fprintf(stderr,"Couldn't open outfile %s, skipping.\n",outfiles[i]);
	  continue;
	}
      
      npulses_good[i] = 0;
      while(!feof(pulsefile))
	{
	  if(mode == 1) //For sift.dat files
	    fscanf(pulsefile,"%d %d %d %d %f %f %f %f %d %f %f",&p.idm,&p.ns,&p.nsm,&p.isamp,&p.snr,&p.mean,&p.rms,&p.suma,&p.nsum,&p.wgp,&p.diff);
	  else if(mode == 0) //For match.best files
	    fscanf(pulsefile,"%d %d %d %d %f %f %f",&p.idm,&p.ns,&p.nsm,&p.isamp,&p.snr,&p.mean,&p.rms);

	  p.isamp_corr = p.isamp - (int)(delaysamp[p.idm]/2.0+0.5) + (int)(stamp_corrections[i]/dt);
	  t_corr = p.isamp_corr*dt;

	  if(penalties[i][(int)floor(t_corr/TBINW_AND)] <= MAXPENALTY)
	    {
	      if(mode == 1) //For sift.dat files
		fprintf(outfile,"%d %d %d %d %f %f %f %f %d %f %f\n",p.idm,p.ns,p.nsm,p.isamp,p.snr,p.mean,p.rms,p.suma,p.nsum,p.wgp,p.diff);
	      else if(mode == 0) //For match.best files
		fprintf(outfile,"%d %d %d %d %f %f %f\n",p.idm,p.ns,p.nsm,p.isamp,p.snr,p.mean,p.rms);
	      npulses_good[i]++;
	    }
	}
      fprintf(stderr,"Pulses total: %ld good: %ld\n",npulses[i],npulses_good[i]);
      fclose(pulsefile);
      fclose(outfile);
    }

  fprintf(stderr,"ANDclean exiting\n");
}

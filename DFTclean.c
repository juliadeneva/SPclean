#include "spclean.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <sys/types.h>
#include <sys/stat.h>


const int IDMCUTOFF_DFT = 10;
const float TBINW_DFT = 0.01;
const float TTOL = 1.0;
const int NSIGMA = 7;
const int NSIGMADFT = 3;
const int NIT = 2;
const float DF = 0.0001; //Hz
const float MAXF = 0.1; //Hz
const float MINF = 0.01;//Hz
const float FTOL = 0.002; //Hz
const int NMAXBAD = 3; //if at least this many beams have radar, clean all
extern int NFILES;

float rms(float* v, int n);
float mean(float* v, int n);
int get_file_size(char *path,off_t *size);

void DFTclean(char** infiles, char** outfiles, char* stampfile, char* dfile, float dt, int ndm, float tobs, float fradar, int mode)
{
  FILE* delayfile, *pulsefile, *outfile, *specfile, *histfile, *stfile;
  //char outfilename[128],specfilename[128],histfilename[128];
  char specfilename[128],histfilename[128];
  int delaysamp[ndm];
  int nbins = floor(tobs/TBINW_DFT)+1;
  float snrsum[nbins], t_corr, sigma, f_curr, tmpsum, tmpsumsq, p_frac;
  float dft2peak,fpeak,amp,tradar,dft2[nbins],dft2rms;
  Pulsus p;
  int ii,i,j,k,idummy,nbadbeams;
  long int npulses[NFILES],npulses_good[NFILES];
  complex dft[nbins],dftpeak;
  float fpeaks[NFILES],phases[NFILES];
  double stamp_corrections[NFILES], stampref;

  fprintf(stderr,"DFTclean starting\n");

  for(i=0;i<NFILES;i++)
    {
      fpeaks[i] = 0;
      phases[i] = 0;
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

  nbadbeams = 0;
  for(ii=0;ii<NFILES;ii++)
    {
      printf("\nOpening infile: %s \n",infiles[ii]);
      if((pulsefile = fopen(infiles[ii],"r")) == NULL)
	{
	  fprintf(stderr,"Couldn't open pulse file %s, skipping.\n",infiles[ii]);
	  continue;
	}

      for(i=0;i<nbins;i++)
	snrsum[i] = 0.0;
      
      //Bin the pulses in time, and add the SNRs of pulses in each bin
      npulses[ii] = 0;
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
	  p.isamp_corr = p.isamp - (int)(delaysamp[p.idm]/2.0+0.5);
	  
	  if(p.idm < IDMCUTOFF_DFT)
	    {
	      t_corr = p.isamp_corr*dt;
	      //snrsum[(int)floor(t_corr/TBINW_DFT)] += p.snr;
	      snrsum[(int)floor(t_corr/TBINW_DFT)] += 1.0;
	    }
	  npulses[ii]++;
	}
      fclose(pulsefile);

      /*
	sprintf(histfilename,"%s%s",infile,".hist");
	if((histfile = fopen(histfilename,"w")) == NULL)
	{
	fprintf(stderr,"Couldn't open file %s, exiting.\n",histfilename);
	exit(1);
	}
	for(j=0;j<nbins;j++)
	fprintf(histfile,"%f %f\n",TBINW_DFT*(j+0.5),snrsum[j]);
	fclose(histfile);
      */
      
      //Find sigma of the resulting "binned time series", removing peaks
      sigma = rms(snrsum,nbins);
      //printf("Sigma0: %f\n",sigma);
      for(i=0; i<NIT; i++)
	{
	  tmpsum = 0;
	  tmpsumsq = 0;
	  k = 0;
	  for(j=0;j<nbins;j++)
	    {
	      if(snrsum[j]/sigma < NSIGMA)
		{
		  tmpsum = tmpsum + snrsum[j];
		  tmpsumsq = tmpsumsq + snrsum[j]*snrsum[j];
		  k++;
		}
	    }
	  sigma = sqrt(tmpsumsq/k - tmpsum*tmpsum/(k*k));
	  //printf("Sigma%d: %f\n",i+1,sigma);
	}
      
      /*
	sprintf(specfilename,"%s%s",infile,".spec");
	if((specfile = fopen(specfilename,"w")) == NULL)
	{
	fprintf(stderr,"Couldn't open file %s, exiting.\n",specfilename);
	exit(1);
	}
      */
      
      //DFT the peaks
      f_curr = MINF;
      i = 0;
      while(f_curr <= MAXF)
	{
	  for(j=0;j<nbins;j++)
	    {
	      if(snrsum[j]/sigma > NSIGMA)
		dft[i] = dft[i] + cexpf(-2.0*M_PI*I*(j+0.5)*TBINW_DFT*f_curr);
	    }
	  dft2[i] = crealf(dft[i]*conjf(dft[i]));
	  dft2[i] = dft2[i]/dft2[0];
	  
	  //fprintf(specfile,"%f %f\n",f_curr,dft2[i]);
	  f_curr = f_curr + DF;
	  i++;
	}
      //fclose(specfile);
      //printf("After DFTing\n");
      

      //Find the power spectrum rms, and where the peak is
      dft2rms = rms(dft2,i);
      f_curr = MINF;
      dft2peak = 0;
      i = 0;
      while(f_curr <= MAXF)
	{
	  if(dft2[i]/dft2rms > dft2peak)
	    {
	      dftpeak = dft[i];
	      dft2peak = dft2[i]/dft2rms;
	      fpeak = f_curr;
	    }
	  f_curr = f_curr + DF;
	  i++;
	}
      fprintf(stderr,"Peak f: %f peak SNR: %f\n",fpeak,dft2peak);
      //fprintf(stderr,"fabs(fpeak-fradar): %f\n",fabs(fpeak-fradar));
  
      if(fabs(fpeak-fradar) < FTOL && dft2peak > NSIGMADFT)
	{
	  fpeaks[ii] = fpeak;
	  amp = sqrt(dft2peak);
	  phases[ii] = crealf(-clogf(dftpeak/amp)/I);
	  nbadbeams++;
	  fprintf(stderr,"RADAR peak found in spectrum.\n");
	}
      else
	fprintf(stderr,"NO radar peak found in spectrum.\n");
    } //END read loop over infiles

  //If too many beams have radar,clean it in all because it's probably 
  //there anyway. Use the average of radar frequencies and phases from 
  //the bad beams to clean the ones where no radar was detected in the DFT
  if(nbadbeams > NMAXBAD)
    {
      tmpsum = 0; //use this for the frequency
      tmpsumsq = 0; //use this for the phase
      j = 0;
      for(i=0;i<NFILES;i++)
	if(fpeaks[i] > 0)
	  {
	    tmpsum += fpeaks[i];
	    tmpsumsq += phases[i];
	    j++;
	  }
      tmpsum = tmpsum/(float)j;
      tmpsumsq = tmpsumsq/(float)j;
      for(i=0;i<NFILES;i++)
	if(fpeaks[i] == 0)
	  {
	    fpeaks[i] = tmpsum;
	    phases[i] = tmpsumsq;
	  }
    }
  
  //Now loop over input files again and do the cleaning
  for(ii=0;ii<NFILES;ii++)
    {
      printf("\nOpening infile: %s \n",infiles[ii]);
      if((pulsefile = fopen(infiles[ii],"r")) == NULL)
	{
	  fprintf(stderr,"Couldn't open pulse file %s, skipping.\n",infiles[ii]);
	  continue;
	}

      //sprintf(outfilename,"%s%s",infiles[ii],".dftclean");
      if((outfile = fopen(outfiles[ii],"w")) == NULL)
	{
	  fprintf(stderr,"Couldn't open outfile %s, exiting.\n",outfiles[ii]);
	  exit(1);
	}
      
      if(fpeaks[ii] != 0)
	{
	  tradar = phases[ii]/(fpeaks[ii]*2.0*M_PI); //t of 1st radar blast
	  printf("tradar0: %f fpeak: %f\n",tradar,fpeaks[ii]);
	  if(tradar < 0)
	    {
	      tradar = tradar + 1.0/fpeaks[ii];
	      printf("tradar: %f\n",tradar);
	    }
	}

      npulses_good[ii] = 0;
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

	  p.isamp_corr = p.isamp - (int)(delaysamp[p.idm]/2.0+0.5);

	  if(fpeaks[ii] != 0)
	    {
	      t_corr = p.isamp_corr*dt;
	      p_frac = fabs(t_corr-tradar)*fpeaks[ii]-floor(fabs(t_corr-tradar)*fpeaks[ii]);	      
	      if(p_frac > TTOL*fpeaks[ii] && 1.0-p_frac > TTOL*fpeaks[ii])
		{
		  //Print good event to outfile
		  if(mode == 1) //For sift.dat files
		    fprintf(outfile,"%d %d %d %d %f %f %f %f %d %f %f\n",p.idm,p.ns,p.nsm,p.isamp,p.snr,p.mean,p.rms,p.suma,p.nsum,p.wgp,p.diff);
		  else if(mode == 0) //For sift.dat files
		    fprintf(outfile,"%d %d %d %d %f %f %f\n",p.idm,p.ns,p.nsm,p.isamp,p.snr,p.mean,p.rms);
		  else
		    {
		      fprintf(stderr,"Unrecognized mode: %d\n",mode);
		      fprintf(stderr,"Use 0 for match.best files or 1 for sift.dat files\n");
		      exit(1);
		    }

		  npulses_good[ii]++;
		}
	    }
	  else
	    {
	      if(mode == 1) //For sift.dat files
		fprintf(outfile,"%d %d %d %d %f %f %f %f %d %f %f\n",p.idm,p.ns,p.nsm,p.isamp,p.snr,p.mean,p.rms,p.suma,p.nsum,p.wgp,p.diff);
	      else if(mode == 0)//For match.best files
		fprintf(outfile,"%d %d %d %d %f %f %f\n",p.idm,p.ns,p.nsm,p.isamp,p.snr,p.mean,p.rms);
	      else
		{
		  fprintf(stderr,"Unrecognized mode: %d\n",mode);
		  fprintf(stderr,"Use 0 for match.best files or 1 for sift.dat files\n");
		  exit(1);
		}
	      npulses_good[ii]++;
	    }
	}
      fprintf(stderr,"Pulses total: %ld good: %ld\n",npulses[ii],npulses_good[ii]);
      fclose(outfile);
      fclose(pulsefile);
    } //END write loop over beam pulse files

  fprintf(stderr,"DFTclean exiting\n\n");
}



float mean(float* v, int n)
{
  int i;
  float sum = 0;

  for(i=0;i<n;i++)
    sum = sum + v[i];
  return sum/n;
}

float rms(float* v, int n)
{
  int i;
  float sum = 0;
  float av;

  av = mean(v,n);
  for(i=0;i<n;i++)
    sum = sum + v[i]*v[i];
  
  return sqrt(sum/n - av*av);
}

int get_file_size(char *path,off_t *size)
{
  struct stat file_stats;

  if(stat(path,&file_stats))
    return -1;
  *size = file_stats.st_size;

  //printf("File size: %ld\n",(long int)(*size));
  return 0;
}

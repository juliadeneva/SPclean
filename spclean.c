#include "spclean.h"
#include <stdio.h>
#include <stdlib.h>

const int NFILES = 7;
const float dt = 0.000064;
const float tobs = 268;
const int ndm = 1272;

void ANDclean(char** infiles, char** outfiles, char* stampfile, char* dfile, float dt, int ndm, float tobs, int mode);
void DFTclean(char** infiles, char** outfiles, char* stampfile, char* dfile, float dt, int ndm, float tobs, float fradar, int mode);

int main(int argc, char** argv)
{
  char* infiles[7], *outfiles[7], *outfiles2[7];
  int i,mode;

  if(argc != 10)
    {
      fprintf(stderr, "Usage: spclean <7 single pulse files> <timestamp file> <mode>\n");
      fprintf(stderr, "(Note1: Pulse files must be in beam ascending order, enter dummy name for missing file.)\n");
      fprintf(stderr, "(Note2: If timestamp file is dummy, beam timestamps are assumed equal.)\n");
      fprintf(stderr, "(Note3: Use mode=0 for match.best files or mode=1 for sift.dat files.)\n");
      exit(1);
    }
  
  mode = atoi(argv[argc-1]);

  for(i=0; i<7; i++)
    {
      infiles[i] = malloc(128);
      sprintf(infiles[i],"%s",argv[i+1]);
      //printf("%s\n",infiles[i]);

      outfiles[i] = malloc(128);
      sprintf(outfiles[i],"%s%d",argv[i+1],1);

      outfiles2[i] = malloc(128);
      sprintf(outfiles2[i],"%s%d",argv[i+1],2);
    }

  //Clean radar
  DFTclean(infiles, outfiles, argv[argc-2], "delays", dt, ndm, tobs, 1.0/12.0, mode);
  
  //Clean multiple-beam RFI
  ANDclean(outfiles, outfiles2, argv[argc-2], "delays", dt, ndm, tobs, mode);
    
        
  return 0;
}

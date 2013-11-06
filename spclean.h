#ifndef _SPCLEAN
#define _SPCLEAN

#define _FILE_OFFSET_BITS 64  //for opening files larger than 2GB

typedef struct plinfo {
  int idm;
  int ns;
  int nsm;
  int isamp;
  int isamp_corr;
  float snr;
  float mean;
  float rms;
  float suma;
  int nsum;
  float wgp;
  float diff;
} Pulsus;

/*
#define IDMCUTOFF_DFT  10;
#define TBINW_DFT  0.01;
#define TTOL  0.2;
#define NSIGMA  7;
#define NIT  2;
#define DF  0.0001; //Hz
#define MAXF  0.1; //Hz
#define MINF  0.01;//Hz
#define FTOL  0.001; //Hz

#define IDMCUTOFF_AND  50;
#define TBINW_AND  1.0;
#define MAXPENALTY  2;
#define NFILES  7;
*/

/*
const int IDMCUTOFF_DFT = 10;
const float TBINW_DFT = 0.01;
const float TTOL = 0.2;
const int NSIGMA = 7;
const int NIT = 2;
const float DF = 0.0001; //Hz
const float MAXF = 0.1; //Hz
const float MINF = 0.01;//Hz
const float FTOL = 0.001; //Hz


const int IDMCUTOFF_AND = 50;
const float TBINW_AND = 1.0;
const int MAXPENALTY = 2;
const int NFILES = 7;
*/

#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "config.h"
#include "segmenter.h"
#include "se.h"
#include "seeds.h"
#include "map.h"
// #include "graphCuts.h"

//static const char *data_file_path = "/home/zappd/SIP_Challenge_Data/StanfordMemorial_rat1_rot90_crop_bandCrop.img";
//static const char *header_file_path = "/home/zappd/SIP_Challenge_Data/StanfordMemorial_rat1_rot90_crop_bandCrop.img.hdr";

static const char *data_file_path = "StanfordMemorial_rat1_rot90_crop_bandCrop.img";
static const char *header_file_path = "StanfordMemorial_rat1_rot90_crop_bandCrop.img.hdr";
static const char *png_extension = ".png";

/* Global configuration, initialized to defaults */
struct gconf gconf = {
  .inputData="data.data",
  .nx=64,
  .ny=64,
  .nw=64,
  .sigma=0.1,
  .evcrit=0.333,
  .t=600,
  .maxevfact=8 ,
  .kcent=16,
  .krep=8,
  .kiter=200,
  .verbosity=1,
  .taucardinality=0.89,
  .kelbw=0.002,
    .cardmin=10
};

void readE(float * data)
{
  FILE * fp;
    
  if((fp=fopen(gconf.inputData,"r"))==NULL)
    fprintf(stderr,"Could not open file %s\n",gconf.inputData);

  if(fread(data,sizeof(float),gconf.nx*gconf.ny*gconf.nw,fp)!=
      (gconf.nx*gconf.ny*gconf.nw))
  {
    fprintf(stderr,"Could not find enough data in the file: %s\n",gconf.inputData);
  }

  fclose(fp);
}

int main(int argc, char* argv[])
{
  readConfig("seconfig");

  gconf.evcrit=strtof(argv[1],NULL);
  gconf.maxevfact=strtof(argv[2],NULL);

  float * A;    /* Nonzero elements       */
  int   * JA;   /* Column indices         */
  int   * IA;   /* Row pointer            */
  int     nnz;  /* Number nonzero elements */

  float * D;    /* Normalization vector to change A into a transition 
                   probability matrix */

  float * E;    /* Eigenvalues */
  float * U;    /* Eigenvectors */
  int     M0;   /* max number of eigenvalues */
  int     M;    /* Number of eigenvalues found */
  float * data=malloc(sizeof(float)*gconf.nx*gconf.ny);
  int * O   =malloc(sizeof(int)*gconf.nx*gconf.ny);
  seed_region_t * id  =malloc(sizeof(seed_region_t)*gconf.kcent);
  readE(data);

  printConfig();
  M0 =(int)(gconf.nx*gconf.ny)/(gconf.maxevfact);

  A  =calloc(gconf.nx*gconf.ny*9,sizeof(float));
  D  =calloc(gconf.nx*gconf.ny  ,sizeof(float));
  JA =calloc(gconf.nx*gconf.ny*9,sizeof(int));
  IA =calloc(gconf.nx*gconf.ny+1,sizeof(int));
  U  =calloc(2*2*M0*gconf.ny*gconf.nx,sizeof(float)); /* complex */
  E  =calloc(2*M0,sizeof(float));                     /* complex */

  nnz=setAffinityGaussian(data,A,JA,IA);
  setNormalizingVector(A,IA,D);
  setMarkovMatrix(A,JA,D,nnz);
  M=setEigenPairs(A,JA,IA,E,U);

  {
    fprintf(stderr,"Found %d eigenvalues\n",M);
  }

  /* Throw out imaginary components (should be 0) and left eigenvectors */
  condenseEigenvectors(U,M);
  condenseEigenvalues(E,M);
  normalizeEigenvectors(U,M);
  
  float * tmp;
  tmp=realloc(U,M*gconf.ny*gconf.nx*sizeof(float));
  U=tmp;
  tmp=realloc(E,M*sizeof(float));
  E=tmp;

  FILE * fp1=fopen("A.data","w");
  FILE * fp2=fopen("JA.data","w");
  FILE * fp3=fopen("IA.data","w");
  FILE * fp4=fopen("E.data","w");
  FILE * fp5=fopen("U.data","w");

  fwrite(A,sizeof(float),nnz,fp1);
  fwrite(JA,sizeof(int),nnz,fp2);
  fwrite(IA,sizeof(int),gconf.nx*gconf.ny+1,fp3);
  fwrite(E,sizeof(float),M,fp4);
  fwrite(U,sizeof(float),M*gconf.nx*gconf.ny,fp5);
    
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
  fclose(fp5);

    return 0;
}

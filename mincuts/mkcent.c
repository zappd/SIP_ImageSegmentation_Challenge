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

void readE(float * data, char * name, size_t count)
{
  FILE * fp;
    
  if((fp=fopen(name,"r"))==NULL)
    fprintf(stderr,"Could not open file %s\n",name);

  if(fread(data,sizeof(float),count,fp)!=
      (gconf.nx*gconf.ny*gconf.nw))
  {
    fprintf(stderr,"Could not find enough data in the file: %s\n",gconf.inputData);
  }

  fclose(fp);
}

int main(int argc, char* argv[])
{
  readConfig("seconfig");


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

  printConfig();

  A  =calloc(gconf.nx*gconf.ny*9,sizeof(float));
  D  =calloc(gconf.nx*gconf.ny  ,sizeof(float));
  JA =calloc(gconf.nx*gconf.ny*9,sizeof(int));
  IA =calloc(gconf.nx*gconf.ny+1,sizeof(int));

  M=(int)strtol(argv[1],NULL,10);
  gconf.kcent=(int)strtol(argv[2],NULL,10);

  readE(data,gconf.inputData,gconf.nx*gconf.ny);
  
  int nk;
  float tau;

  nnz=setAffinityGaussian(data,A,JA,IA);
  setNormalizingVector(A,IA,D);
  setMarkovMatrix(A,JA,D,nnz);
  
  float * W   =calloc(M*gconf.ny*gconf.nx,sizeof(float)); 
  float * Z   =calloc(M*gconf.ny*gconf.nx,sizeof(float)); 
  float * R   =calloc(gconf.nx*gconf.ny  ,sizeof(float));
  double * Cent =calloc(M*gconf.kcent,sizeof(double));
  int * Card =calloc(gconf.kcent  ,sizeof(int));
  float * F =calloc(gconf.kcent*gconf.kcent  ,sizeof(int));
  U  =calloc(M*gconf.ny*gconf.nx,sizeof(float)); /* complex */
  E  =calloc(M,sizeof(float));                     /* complex */

  readE(U,"U.data",M*gconf.ny*gconf.nx);
  readE(E,"E.data",M);

  setSpectralWeights(W,E,U,D,M);
  setCoordinateTransform(Z,W,E,U,D,M);
  setTransformCenters(Z,R,O,M,Cent);

  tau=findThresholdOverlay(R,O);
  setThresholdOverlay(R,O,tau);
  setMinimumCardinality(R,O);
  setCenterCardinality(O,Card);
  setF(F,Card,Cent,M,tau);


  nk=setUniqueOrderedIndices(O);

  float * Cx=calloc(nk,sizeof(float));
  float * Cy=calloc(nk,sizeof(float));

  setSeedCentroids(O,Cx,Cy,nk);

  FILE * fp6=fopen("O.data","w");
  FILE * fp7=fopen("Z.data","w");
  FILE * fp8=fopen("D.data","w");
  FILE * fp9=fopen("R.data","w");
  FILE * fp10=fopen("W.data","w");
  FILE * fp11=fopen("Cx.data","w");
  FILE * fp12=fopen("Cy.data","w");
  fwrite(O,sizeof(int),gconf.nx*gconf.ny,fp6);
  fwrite(R,sizeof(float),gconf.nx*gconf.ny,fp9);
  fwrite(Z,sizeof(float),M*gconf.nx*gconf.ny,fp7);
  fwrite(D,sizeof(float),gconf.nx*gconf.ny,fp8);
  fwrite(W,sizeof(float),M*gconf.nx*gconf.ny,fp10);
  fwrite(Cx,sizeof(float),nk,fp11);
  fwrite(Cy,sizeof(float),nk,fp12);
  fclose(fp6);
  fclose(fp7);
  fclose(fp8);
  fclose(fp9);
  fclose(fp10);
  fclose(fp11);
  fclose(fp12);

  free(A);
  free(Cx);
  free(Cy);
  free(D);
  free(JA);
  free(IA);
  free(W);
  free(Z);
  free(R);

    return 0;
}

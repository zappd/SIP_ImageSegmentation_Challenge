#include <math.h>
#include <stdio.h>


#include "feast.h"
#include "feast_sparse.h"
#include "mkl_cblas.h"
#include "mkl_trans.h"
#include "mkl_lapacke.h"
#include "se.h"
#include "config.h"

#include "vl/kmeans.h"
#include "vl/mathop.h"
static void setCentroids();
static void triangualateCentroids();
static void setNeighbors();


static void writeIntermediate( void * data, size_t size, const char * name)
{

  FILE * fp;
    
  if((fp=fopen(name,"w"))==NULL)
    fprintf(stderr,"Could not open file %s\n",gconf.inputData);

  if(fwrite(data,1,size,fp)!=(size))
  {
    fprintf(stderr,"Error while writing to file %s\n",name);
  }

  fclose(fp);
}

static void readIntermediate( void * data, size_t size,const char * name)
{

  FILE * fp;
    
  if((fp=fopen(name,"r"))==NULL)
    fprintf(stderr,"Could not open file %s\n",gconf.inputData);

  if(fread(data,1,size,fp)!=(size))
  {
    fprintf(stderr,"Error while writing to file %s\n",name);
  }

  fclose(fp);
}

void setSpectralSeeds(float * data, int * O)
{ 
  /* data is the gconf.nx X gconf.ny intensity map
   * O is the cluster overlay map which maps to integers 1 to C, where 0 
   * indicates there is no seed cluster */


                /* Sparse affinity matrix */
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

  float * W;    /* Weights of pixels per each blur kernel */
  float * Z;    /* Projected values onto blur kernels */

  float * R;    /* Overlay with distances from cluster centers*/

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

  if(gconf.verbosity==1)
  {
    fprintf(stderr,"Found %d eigenvalues\n",M);
  }

  /* Throw out imaginary components (should be 0) and left eigenvectors */
  condenseEigenvectors(U,M);
  condenseEigenvalues(E,M);

  /* Purge excess space*/
  float * tmp;
  tmp=realloc(U,M*gconf.ny*gconf.nx*sizeof(float));
  U=tmp;
  tmp=realloc(E,M*sizeof(float));
  E=tmp;


  /* Does not raise high water mark due to purging */
  W   =calloc(M*gconf.ny*gconf.nx,sizeof(float)); 
  Z   =calloc(M*gconf.ny*gconf.nx,sizeof(float)); 
  R   =calloc(gconf.nx*gconf.ny  ,sizeof(float));

  setSpectralWeights(W,E,U,D,M);
  setCoordinateTransform(Z,W,E,U,D,M);
  
  /* Z is in row major to play nice with vlfeat library */
  setTransformCenters(Z,R,O,M);

  /* now threshold */
  thresholdOverlay(R,O);

  free(A);
  free(D);
  free(JA);
  free(IA);
  free(U);
  free(E);
  free(W);
  free(Z);
  free(R);
}

static void thresholdOverlay(float * R, float * O)
{
int i,j,k;
int card;

const int div=100;
float thresh[div];


thresh[0]=0;
for(i=1;i<div;i++);
{
  thresh[i]=thresh[i-1]+1/div;
}

for (i=0;i<div;i++)
{
  for (j=0;i<gconf.nx;j++)
  {
    for (k=0;k<gconf.ny;k++)
    {


    }
  }
}

}


static void squareRootMatrix(float * A, int N)
{
  /* gets the squareroot of an NxN positive semidefinite matrix */
#define Idx(i,j) (N*i+j)
  int M=N;
  float * Wtmp = malloc(sizeof(float)*N);
  float * Ztmp = malloc(sizeof(float)*N*N);
  float * Btmp = malloc(sizeof(float)*N*N);
  lapack_int * Itmp = malloc(sizeof(lapack_int)*2*N);
  int info;
  int i,j;
  float prec=LAPACKE_slamch('S');

  /* http://www.seehuhn.de/pages/matrixfn */
  info=LAPACKE_ssyevr(LAPACK_COL_MAJOR,'V','A','L',N,A,N,0,0,0,0,prec,&M,Wtmp,Ztmp,N,Itmp);

  if((N!=M)||info!=0)
    fprintf(stderr,"Failures when trying to find square root of matrix");

  /* B=Z*D */
  for(i=0;i<N;i++)
  {
    for(j=0;j<N;j++)
    {
      Btmp[Idx(i,j)]=Ztmp[Idx(i,j)]*sqrt(Wtmp[i]);
    }
  }

  /* get square root */
  cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0, Btmp, N, Ztmp, N, 0, A, N);


  free(Wtmp);
  free(Ztmp);
  free(Btmp);
  free(Itmp);


#undef Idx
}
  static void normalizeEigenvectors(float * U, int M)
{
  /* Feast doesn not normalize EV, do it here after condense */
#define Idx(i,j) ((i)*gconf.nx*gconf.ny+(j))
    int i,j;
    double sumsq;
    for (i=0;i<M;i++)
    {
      sumsq=0;
      for(j=0;j<gconf.nx*gconf.ny;j++)
      {
        sumsq+=U[Idx(i,j)]*U[Idx(i,j)];
      }
      for(j=0;j<gconf.nx*gconf.ny;j++)
      {
        U[Idx(i,j)]/=sqrt(sumsq);
      }

    }
#undef Idx

}

static void setTransformCenters(float * Z, float * R, int * O, int M)
{
  /* Z is coordinate transofrmed,  R is overlay distance from center, O is overlap map*/
  /* Z is in row major for the kmeans library */
  double energy;
  double * centers;
  double gmin;
  double sumsq;
  int i,j,imin;

  /* krep, kiter, kcent */
  vl_uint32 * assignments = vl_malloc(sizeof(vl_uint32)*gconf.nx*gconf.ny);

  /* normalize */
#define Jdx(i,j) (i*M+j)
  for(i=0;i<gconf.nx*gconf.ny;i++)
  {
    sumsq=0;
    for(j=0;j<M;j++)
    {
      sumsq+=Z[Jdx(i,j)]*Z[Jdx(i,j)];
    }

    for(j=0;j<M;j++)
    {
      Z[Jdx(i,j)]/=sqrt(sumsq);
    }
  }
#undef Jdx


  VlKMeans * kmeans;
  kmeans = vl_kmeans_new(VlDistanceL2, VL_TYPE_FLOAT);
  if(gconf.verbosity==1)
  {
  vl_kmeans_set_verbosity(kmeans,1);
  }

  vl_kmeans_set_algorithm(kmeans,VlKMeansLloyd);
  vl_kmeans_init_centers_plus_plus(kmeans,Z,M,gconf.nx*gconf.ny,gconf.kcent);
  vl_kmeans_set_max_num_iterations(kmeans,gconf.kiter);
  vl_kmeans_set_min_energy_variation(kmeans,10e-5);
  vl_kmeans_set_num_repetitions(kmeans,gconf.krep);
  vl_kmeans_cluster(kmeans,Z,M,gconf.nx*gconf.ny,gconf.kcent);
  vl_kmeans_refine_centers(kmeans,Z,gconf.nx*gconf.ny);
  energy=vl_kmeans_get_energy(kmeans);
  centers = (double *)vl_kmeans_get_centers(kmeans);

  vl_kmeans_quantize(kmeans,assignments,R,Z,gconf.nx*gconf.ny);


  for(i=0;i<gconf.nx*gconf.ny;i++)
  {
    O[i]=(int)assignments[i];
  }

  vl_kmeans_delete(kmeans);
  vl_free(assignments);
}


static void setCoordinateTransform(float * Z, float * W, float * E, float * U, float * D, int M)
{
  /* computes Q=U^T D U, z=Q^(1/2)W*/
  /* Q is MxM */
#define Idx(i,j) ((i)*gconf.nx*gconf.ny+(j))
  int i,j;
  int K=gconf.nx*gconf.ny;

  float * Q = malloc(sizeof(float)*M*K);

  /*Z=D U*/
  for (i=0;i<M;i++)
  {
    for(j=0;j<K;j++)
    {
      Z[Idx(i,j)]=D[j]*U[Idx(i,j)];
    }
  }

  /*Z=U^T Z*/
  cblas_sgemm (CblasColMajor, CblasTrans, CblasNoTrans, M, M, K, 1.0, U, K, Z, K, 0.0, Q, M);

  /* Q is formed */
  memcpy(Z,Q,sizeof(float)*M*M);

  squareRootMatrix(Z,M);

  cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, M, K, M, 1.0, Z, M, W, K, 0.0, Q, M);
  memcpy(Z,Q,sizeof(float)*M*K);

  /*Z is formed, size (MxK)*/


  /* Kmeans algo need row-major data, no not transpose */
  //mkl_simatcopy('C','T',M,K,1.0,Z,M,K);
  
  free(Q);
#undef Idx
}

static void setSpectralWeights(float * W, float * E, float * U, float * D, int M)
{
  /* Computes transpose of  W=E^t U^T D^(-1/2), eqn (6),
   * as W^T=D^(-1/2)^T U E^t */
#define Idx(i,j) ((i)*gconf.nx*gconf.ny+(j))
  int i,j;
  for (i=0;i<M;i++)
  {
    for(j=0;j<gconf.nx*gconf.ny;j++)
    {
      W[Idx(i,j)]=1/sqrt(D[j])*U[Idx(i,j)]*pow(E[i],gconf.t);
    }
  }
#undef Idx
}

static void condenseEigenvalues(float * E,int m)
{
  /* get rid of imaginary components (should be zero) */
  int i;
  for(i=1;i<m;i++)
  {
    E[i]=E[2*i];
  }

}

static void condenseEigenvectors(float * U, int m)
{
  /* Eigenvectors are returned from feast to include imaginary components
   * followed by the left eigenvector, this function aligns real components
   * sequentially.
   *
   * Sequentializes real components then memcpy each EV
   *
   * On exit, U only has relevant information on gconf.xy*gconf.ny*m elements
   */

#define DST(i,j) ((i)*gconf.nx*gconf.ny*2+(j))
#define SRC(i,j) ((i)*gconf.nx*gconf.ny*2+(j)*2)

  int i,j;
  int src;
  int des;

  for(i=0;i<m;i++)
  {
    for(j=1;j<gconf.nx*gconf.ny;j++)
    {
      U[DST(i,j)]=U[SRC(i,j)];
    } 

    if(i!=0) 
    {
      memmove(&U[DST(i,0)/2],&U[DST(i,0)],sizeof(float)*gconf.nx*gconf.ny);
    }
  }

#undef DST
#undef SRC
}

static int setEigenPairs(float * A, int * JA, int * IA, float * E, float * U)
{
  /* Feast in */
  int fpm[64];    /* Feast configuration */
  int loop;       /* Number of subspace iterations */
  int M0;         /* Maximum EV count */
  float Emid[2];  /* Contour ellipse (complex) center */
  float r;        /* ellipse radius */
  int N,nnz;      /* Problem size*/

  /* Feast out */
  float epsout;  /* Relative error (out) */
  int M;          /* # of EV found in interval */
  int info;       /* Error handling */
  float * res;    /*residual*/

  /* Set problem parameters */
  N=gconf.nx*gconf.ny;

  r=(1-gconf.evcrit)/2;
  Emid[0]=gconf.evcrit+r;
  r*=1.0;

  M0=(int)(gconf.nx*gconf.ny)/gconf.maxevfact;

  /* Configure feast */
  fpm[0]=1;       /* verbose output */
  fpm[3]=600;       /* refinement */
  fpm[5]=1;       /* error type */
  // fpm[6]=5;       /* 10^-x convergence criteria */
  fpm[17]=5;     /* ellipse contour ratio (im/real) */
  //fpm[17]=10-;     /* ellipse contour ratio (im/real) */

  feastinit(fpm);

  res=calloc(2*M0,sizeof(float));

  if(gconf.verbosity==1)
  {
  fprintf(stderr,"Feast Initialized\n");
  }

  sfeast_gcsrev(&N,A,IA,JA,fpm,&epsout,&loop,Emid,&r,&M0,E,U,&M,res,&info);

  if(gconf.verbosity==1)
  {
  fprintf(stderr,"Feast Finished\n");
  }

  free(res);

  return M;

}

static inline float affinityFunctionGaussian(float I1, float I2)
{
  return exp(-(I1-I2)*(I1-I2)/(2*gconf.sigma*gconf.sigma));
}

static void setMarkovMatrix(float * A, int * JA, float * D, int nnz)
{
  /* Divide each column by its normalizing factor */
  int i,j;
#define FortIndex 1
  for(i=0;i<nnz;i++)
  {
    j=JA[i]-FortIndex;
    A[i]/=D[j];
  }
#undef FortIndex
}

static void setNormalizingVector(float * A, int * IA, float * D)
{
  /* sum along each column (or row, its symmetric), put in D */
  long i,j;
  double sum;

#define FortIndex 1
  for(i=0;i<gconf.nx*gconf.ny;i++)
  {
    D[i]=0;
    for(j=IA[i]-FortIndex;j<IA[i+1]-FortIndex;j++)
    {  
      D[i]+=A[j];
    }
  }
#undef FortIndex
}

/* Put affinity matrix into sparse CSR format */
static int setAffinityGaussian(float * data, float * A, int * JA, int * IA)
{
#define Data(x,y,z) data[gconf.nx*gconf.ny*(z)+gconf.nx*(y)+(x)]
#define Idx(x,y,z)  (gconf.nx*gconf.ny*(z)+gconf.nx*(y)+(x))
#define SubX(n)     (((n)%(gconf.nx*gconf.ny))%gconf.nx)
#define SubY(n)     (((n)-SubX(n)%(gconf.nx*gconf.ny))/gconf.nx)
#define SubZ(n)     ((((n)-SubX(n)-gconf.nx*SubY(n))/(gconf.nx*gconf.ny)))
#define FortIndex 1 /* Use fortran indices for eigenvalue solver */
  //#define FortIndex

  int i,j,k;
  int ni,nj;
  int ak=0,jk=0,ik=0;

  /* Doing one row at a time*/
  for(ni=0;ni<gconf.nx*gconf.ny;ni++)
  {
    /* Row pointer */
    IA[ik++]=ak+FortIndex;

    /* ni is the current pixel */
    k=SubZ(ni);
    j=SubY(ni);
    i=SubX(ni);

    /* SOUTH WEST */
    if( ((j-1) >= 0) && ((i-1)>=0) )
    {
      A[ak++]=affinityFunctionGaussian(Data(i,j,0),Data(i-1,j-1,0));
      JA[jk++]=Idx(i-1,j-1,0)+FortIndex;
    }

    /* SOUTH */
    if( ((j-1)>=0) )
    {
      A[ak++]=affinityFunctionGaussian(Data(i,j,0),Data(i,j-1,0));
      JA[jk++]=Idx(i,j-1,0)+FortIndex;
    }

    /* SOUTH EAST */
    if( ((j-1) >= 0) && ((i+1)<gconf.nx) )
    {
      A[ak++]=affinityFunctionGaussian(Data(i,j,0),Data(i+1,j-1,0));
      JA[jk++]=Idx(i+1,j-1,0)+FortIndex;
    }

    /* WEST */
    if( ((i-1)>=0) )
    {
      A[ak++]=affinityFunctionGaussian(Data(i,j,0),Data(i-1,j,0));
      JA[jk++]=Idx(i-1,j,0)+FortIndex;
    }

    /* CENTRAL */
    {
      A[ak++]=affinityFunctionGaussian(Data(i,j,0),Data(i,j,0));
      JA[jk++]=Idx(i,j,0)+FortIndex;
    }

    /* EAST */
    if( ((i+1)<gconf.nx) )
    {
      A[ak++]=affinityFunctionGaussian(Data(i,j,0),Data(i+1,j,0));
      JA[jk++]=Idx(i+1,j,0)+FortIndex;
    }

    /* NORTH WEST */
    if( ((j+1) < gconf.ny) && ((i-1) >=0) )
    {
      A[ak++]=affinityFunctionGaussian(Data(i,j,0),Data(i-1,j+1,0));
      JA[jk++]=Idx(i-1,j+1,0)+FortIndex;
    }

    /* NORTH */
    if( ((j+1) < gconf.ny) )
    {
      A[ak++]=affinityFunctionGaussian(Data(i,j,0),Data(i,j+1,0));
      JA[jk++]=Idx(i,j+1,0)+FortIndex;
    }

    /* NORTH EAST */
    if( ((j+1) < gconf.ny) && ((i+1)<gconf.nx) )
    {
      A[ak++]=affinityFunctionGaussian(Data(i,j,0),Data(i+1,j+1,0));
      JA[jk++]=Idx(i+1,j+1,0)+FortIndex;
    }

  }
  IA[ik++]=ak+FortIndex;
  return ak;

#undef Data
#undef Idx
#undef SubZ
#undef SubY
#undef SubX
#undef FortIndex
}

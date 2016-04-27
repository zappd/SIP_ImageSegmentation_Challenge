#if !defined(_SE_H)
#define _SE_H

void setSpectralSeeds(float * data, int * O);

static void writeIntermediate(float *, size_t, const char *);
static void readIntermediate(float *, size_t, const char *);

static int setAffinityGaussian(float * , float *, int *, int *);
static void setMarkovMatrix(float * , int * , float * , int);
static void setNormalizingVector(float * , int * , float * );
static int setEigenPairs(float * , int * , int * , float * , float * );

static void condenseEigenvectors(float *, int );
static void condenseEigenvalues(float *, int );
static void normalizeEigenvectors(float *,int);

static void setSpectralWeights(float * , float * , float * , float * , int);
static void setCoordinateTransform(float * , float * , float * , float * , float * , int );
static void setTransformCenters(float * , float * , int * , int );


#endif

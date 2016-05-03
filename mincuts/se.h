#if !defined(_SE_H)
#define _SE_H
#include "seeds.h"

void setSpectralSeeds(float * data, int * O, seed_region_t * );

 void setMinimumCardinality(float *, int *);

 float findThresholdOverlay(float * , int * );
 void setThresholdOverlay(float * , int * , float);
 void setSeedCentroids(int * , float * , float * , int );
 int setUniqueOrderedIndices(int *);

 void writeIntermediate(float *, size_t, const char *);
 void readIntermediate(float *, size_t, const char *);

 int setAffinityGaussian(float * , float *, int *, int *);
 void setMarkovMatrix(float * , int * , float * , int);
 void setNormalizingVector(float * , int * , float * );
 int setEigenPairs(float * , int * , int * , float * , float * );

 void condenseEigenvectors(float *, int );
 void condenseEigenvalues(float *, int );
 void normalizeEigenvectors(float *,int);

 void setSpectralWeights(float * , float * , float * , float * , int);
 void setCoordinateTransform(float * , float * , float * , float * , float * , int );
 void setTransformCenters(float * , float * , int * , int );


#endif

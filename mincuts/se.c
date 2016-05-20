#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#include "feast.h"
#include "feast_sparse.h"
#include "cblas.h"
#include "trans.h"
#include "lapacke.h"
#include "se.h"

#include "vl/kmeans.h"
#include "vl/mathop.h"


/********************************
 * Implementation
 ********************************/

void setSpectralSeeds(float *data, int *O, seed_region_t *seeds)
{
    /* data is the global_config.nx X global_config.ny intensity map
     * O is the cluster overlay map which maps to integers 1 to C, where 0
     * indicates there is no seed cluster */

    /* Sparse affinity matrix */
    float *A;  /* Nonzero elements */
    int   *JA; /* Column indices */
    int   *IA; /* Row pointer */
    int   nnz; /* Number nonzero elements */

    float *D; /* Normalization vector to change A into a transition probability matrix */

    float *E; /* Eigenvalues */
    float *U; /* Eigenvectors */
    int   M0; /* max number of eigenvalues */
    int   M;  /* Number of eigenvalues found */

    float *W; /* Weights of pixels per each blur kernel */
    float *Z; /* Projected values onto blur kernels */

    float *R; /* Overlay with distances from cluster centers*/

    float tau; /* Procedurally determined threshold for seed regions */

    M0 = (int) ((global_config.nx * global_config.ny) / (global_config.maxevfact));

    A  = calloc((size_t) (global_config.nx * global_config.ny * 9), sizeof(float));
    D  = calloc((size_t) (global_config.nx * global_config.ny), sizeof(float));
    JA = calloc((size_t) (global_config.nx * global_config.ny * 9), sizeof(int));
    IA = calloc((size_t) (global_config.nx * global_config.ny + 1), sizeof(int));
    U  = calloc((size_t) (2 * 2 * M0 * global_config.ny * global_config.nx), sizeof(float)); /* complex */
    E  = calloc((size_t) (2 * M0), sizeof(float));                     /* complex */

    nnz = setAffinityGaussian(data, A, JA, IA);
    setNormalizingVector(A, IA, D);
    setMarkovMatrix(A, JA, D, nnz);
    M = setEigenPairs(A, JA, IA, E, U);

    //if(global_config.verbosity==1)
    {
        fprintf(stderr, "Found %d eigenvalues\n", M);
    }

    /* Throw out imaginary components (should be 0) and left eigenvectors */
    condenseEigenvectors(U, M);
    condenseEigenvalues(E, M);
    normalizeEigenvectors(U, M);

    /* Purge excess space*/
    float *tmp;
    tmp = realloc(U, M * global_config.ny * global_config.nx * sizeof(float));
    U   = tmp;
    tmp = realloc(E, M * sizeof(float));
    E   = tmp;


    /* Does not raise high water mark due to purging */
    W = calloc((size_t) (M * global_config.ny * global_config.nx), sizeof(float));
    Z = calloc((size_t) (M * global_config.ny * global_config.nx), sizeof(float));
    R = calloc((size_t) (global_config.nx * global_config.ny), sizeof(float));

    setSpectralWeights(W, E, U, D, M);
    setCoordinateTransform(Z, W, E, U, D, M);

    /* Z is in row major to play nice with vlfeat library */
    setTransformCenters(Z, R, O, M);

    int nk;
    /* now threshold */
    tau = findThresholdOverlay(R, O);
    setThresholdOverlay(R, O, tau);
    setMinimumCardinality(R, O);
    nk = setUniqueOrderedIndices(O);

    float *Cx = calloc((size_t) nk, sizeof(float));
    float *Cy = calloc((size_t) nk, sizeof(float));

    setSeedCentroids(O, Cx, Cy, nk);

    FILE *fp1  = fopen("A.data", "w");
    FILE *fp2  = fopen("JA.data", "w");
    FILE *fp3  = fopen("IA.data", "w");
    FILE *fp4  = fopen("E.data", "w");
    FILE *fp5  = fopen("U.data", "w");
    FILE *fp6  = fopen("O.data", "w");
    FILE *fp7  = fopen("Z.data", "w");
    FILE *fp8  = fopen("D.data", "w");
    FILE *fp9  = fopen("R.data", "w");
    FILE *fp10 = fopen("W.data", "w");
    FILE *fp11 = fopen("Cx.data", "w");
    FILE *fp12 = fopen("Cy.data", "w");

    fwrite(A, sizeof(float), (size_t) nnz, fp1);
    fwrite(JA, sizeof(int), (size_t) nnz, fp2);
    fwrite(IA, sizeof(int), (size_t) (global_config.nx * global_config.ny + 1), fp3);
    fwrite(E, sizeof(float), (size_t) M, fp4);
    fwrite(U, sizeof(float), (size_t) (M * global_config.nx * global_config.ny), fp5);
    fwrite(O, sizeof(int), (size_t) (global_config.nx * global_config.ny), fp6);
    fwrite(R, sizeof(float), (size_t) (global_config.nx * global_config.ny), fp9);
    fwrite(Z, sizeof(float), (size_t) (M * global_config.nx * global_config.ny), fp7);
    fwrite(D, sizeof(float), (size_t) (global_config.nx * global_config.ny), fp8);
    fwrite(W, sizeof(float), (size_t) (M * global_config.nx * global_config.ny), fp10);
    fwrite(Cx, sizeof(float), (size_t) nk, fp11);
    fwrite(Cy, sizeof(float), (size_t) nk, fp12);

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
    fclose(fp5);
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
    free(U);
    free(E);
    free(W);
    free(Z);
    free(R);
}

void setSeedCentroids(int *O, float *Cx, float *Cy, int countkey)
{
#define SubX(n)     (((n)%(global_config.nx*global_config.ny))%global_config.nx)
#define SubY(n)     (((n)-SubX(n)%(global_config.nx*global_config.ny))/global_config.nx)
    int i;
    int j;
    int x;
    int y;
    int count;

    for (i = 1; i <= countkey; i++)
    {
        count  = 0;
        for (j = 0; j < global_config.nx * global_config.ny; j++)
        {
            if (O[j] == i)
            {
                x = SubX(j);
                y = SubY(j);


                Cx[i - 1] += (float) x;
                Cy[i - 1] += (float) y;

                count++;
            }
        }
        Cx[i - 1] /= (float) count;
        Cy[i - 1] /= (float) count;
    }
#undef SubX
#undef Suby
}

int findKeyInArray(int key, int *array, int asize)
{
    /* returns index >=0 with index in array where found, <0 if not found */
    int i;
    int ret = -99;

    for (i = 0; i < asize; i++)
    {
        if (array[i] == key)
        {
            ret = i;
            break;
        }
    }

    return ret;
}


int setUniqueOrderedIndices(int *O)
{
    int *newkeys = malloc(sizeof(int) * global_config.kcent + 1);
    int countkey = 1;
    int oldkey;
    int i;
    newkeys[0] = 0;

    for (i = 1; i <= global_config.kcent; i++)
    {
        if (findKeyInArray(i, O, global_config.nx * global_config.ny) >= 0)
        {
            newkeys[i] = countkey;
            countkey++;
        }
        else
        {
            newkeys[i] = 0;
        }
    }

    for (i = 0; i < global_config.nx * global_config.ny; i++)
    {
        oldkey = O[i];
        O[i] = newkeys[oldkey];
    }

    return countkey;
}

void setF(float *F, int *Card, float *Cent, int M, float tau)
{
    int   i, j, k;
    float r;
    double *Centact = calloc((size_t) (M * global_config.kcent), sizeof(double));

#define Idx(i, j) (i*M+j) /* i is center, j is dimension */
    for (i = 1; i <= global_config.kcent; i++)
    {
        for (j = 0; j < 50; j++)
        {
            r = (float)j * 0.02f;

            for (k = 0; k < M; k++)
            {
                Centact[Idx(i - 1, k)];
            }
        }
    }

#undef Idx
}

void setCenterCardinality(int *O, int *Card)
{
    /* gets rid of seed regions with cardinality < global_config.cardmin */
    int i, k;
    int card;

    for (k = 1; k <= global_config.kcent; k++)
    {
        Card[k - 1] = 0;

        for (i = 0; i < global_config.nx * global_config.ny; i++)
        {
            if (O[i] == k)
            {
                Card[k - 1]++;
            }
        }
    }
}

void setMinimumCardinality(float *R, int *O)
{
    /* gets rid of seed regions with cardinality < global_config.cardmin */
    int i, k;
    int card;

    for (k = 1; k <= global_config.kcent; k++)
    {
        card = 0;

        for (i = 0; i < global_config.nx * global_config.ny; i++)
        {
            if (O[i] == k)
            {
                card++;
            }
        }

        if (card < global_config.cardmin)
        {
            for (i = 0; i < global_config.nx * global_config.ny; i++)
            {
                if (O[i] == k)
                {
                    O[i] = 0;
                }
            }
        }

    }
}

void setThresholdOverlay(float *R, int *O, float tau)
{
    /* Set ID to zero if above threshold, add one to each other ID */
    int i;

    for (i = 0; i < global_config.nx * global_config.ny; i++)
    {
        O[i] += 1;

        if (R[i] > tau)
        {
            O[i] = 0;
        }
    }
}

float findThresholdOverlay(float *R, int *O)
{
    int   i, j;
    int   card;
    float tau_ratio;

    const int div = 100;
    float     thresh[div];

    if (global_config.verbosity == 1)
    {
        fprintf(stderr, "Finding correct Tau");
    }
    thresh[0] = 0;
    for (i = 1; i < 100; i++)
    {
        thresh[i] = thresh[i - 1] + 1 / (float) div;
    }

    for (i = 0; i < div; i++)
    {
        card   = 0;
        for (j = 0; j < global_config.nx * global_config.ny; j++)
        {
            if (R[j] < thresh[i])
            {
                card++;
            }
        }

        tau_ratio = (float) card / ((float) global_config.nx * global_config.ny);
        if (global_config.verbosity == 1)
        {
            fprintf(stderr, "tau: %d %e %e\n", card, tau_ratio, thresh[i]);
        }

        if (tau_ratio > global_config.taucardinality)
        {
            break;
        }

    }

    return thresh[i];

}


void squareRootMatrix(float *A, int N)
{
    /* gets the squareroot of an NxN positive semidefinite matrix */
#define Idx(i, j) (N*i+j)
    int        M     = N;
    float      *Wtmp = malloc(sizeof(float) * N);
    float      *Ztmp = malloc(sizeof(float) * N * N);
    float      *Btmp = malloc(sizeof(float) * N * N);
    lapack_int *Itmp = malloc(sizeof(lapack_int) * 2 * N);
    int        info;
    int        i, j;
    float      prec  = LAPACKE_slamch('S');

    /* http://www.seehuhn.de/pages/matrixfn */
    info = LAPACKE_ssyevr(LAPACK_COL_MAJOR, 'V', 'A', 'L', N, A, N, 0, 0, 0, 0, prec, &M, Wtmp, Ztmp, N, Itmp);

    if ((N != M) || info != 0)
    {
        fprintf(stderr, "Failures when trying to find square root of matrix");
    }

    /* B=Z*D */
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            Btmp[Idx(i, j)] = (float) (Ztmp[Idx(i, j)] * sqrt(Wtmp[i]));
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

void normalizeEigenvectors(float *U, int M)
{
    /* Feast does not normalize EV, do it here after condense */
#define Idx(i, j) ((i)*global_config.nx*global_config.ny+(j))

    int    i, j;
    double sumsq;

    for (i = 0; i < M; i++)
    {
        sumsq = 0;

        for (j = 0; j < global_config.nx * global_config.ny; j++)
        {
            sumsq += U[Idx(i, j)] * U[Idx(i, j)];
        }

        for (j = 0; j < global_config.nx * global_config.ny; j++)
        {
            U[Idx(i, j)] /= sqrt(sumsq);
        }

    }

#undef Idx
}

double setKmeans(float *Z, float *R, int *O, int M, vl_size K, double *Cent)
{
    double energy;
    double *centers;

    int i;

    VlKMeans *kmeans;
    kmeans = vl_kmeans_new(VlDistanceL2, VL_TYPE_FLOAT);

    if (global_config.verbosity == 1)
    {
        vl_kmeans_set_verbosity(kmeans, 1);
    }

    vl_kmeans_set_algorithm(kmeans, VlKMeansLloyd);
    vl_kmeans_init_centers_plus_plus(kmeans, Z, (vl_size) M, (vl_size) (global_config.nx * global_config.ny), K);
    vl_kmeans_set_max_num_iterations(kmeans, (vl_size) global_config.kiter);
    vl_kmeans_set_min_energy_variation(kmeans, 10e-5);
    vl_kmeans_set_num_repetitions(kmeans, (vl_size) global_config.krep);
    vl_kmeans_cluster(kmeans, Z, (vl_size) M, (vl_size) (global_config.nx * global_config.ny), K);

    energy = vl_kmeans_refine_centers(kmeans, Z, (vl_size) (global_config.nx * global_config.ny));

    centers = (double *) vl_kmeans_get_centers(kmeans);

    for (i = 1; i < M * global_config.kcent; i++)
    {
        Cent[i] = centers[i];
    }

    /* krep, kiter, kcent */
    vl_uint32 *assignments = vl_malloc(sizeof(vl_uint32) * global_config.nx * global_config.ny);

    vl_kmeans_quantize(kmeans, assignments, R, Z, (vl_size) (global_config.nx * global_config.ny));

    for (i = 0; i < global_config.nx * global_config.ny; i++)
    {
        O[i] = (int) assignments[i];
    }

    vl_kmeans_delete(kmeans);
    vl_free(assignments);

    return energy;
}

void setTransformCenters(float *Z, float *R, int *O, int M)
{
    /* Z is coordinate transofrmed,  R is overlay distance from center, O is overlap map*/
    /* Z is in row major for the kmeans library */

    double sumsq;
    int    i, j;

    /* normalize */
#define Jdx(i, j) (i*M+j)

    for (i = 0; i < global_config.nx * global_config.ny; i++)
    {
        sumsq  = 0;
        for (j = 0; j < M; j++)
        {
            sumsq += Z[Jdx(i, j)] * Z[Jdx(i, j)];
        }

        for (j = 0; j < M; j++)
        {
            Z[Jdx(i, j)] /= sqrt(sumsq);
        }
    }

#undef Jdx
}


void setCoordinateTransform(float *Z, float *W, float *E, float *U, float *D, int M)
{
    /* computes Q=U^T D U, z=Q^(1/2)W*/
    /* Q is MxM */
#define Idx(i, j) ((i)*global_config.nx*global_config.ny+(j))

    int i, j;
    int K = global_config.nx * global_config.ny;

    float *Q = malloc(sizeof(float) * M * K);

    /*Z=D U*/
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < K; j++)
        {
            Z[Idx(i, j)] = D[j] * U[Idx(i, j)];
        }
    }

    /*Z=U^T Z*/
    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, M, M, K, 1.0, U, K, Z, K, 0.0, Q, M);

    /* Q is formed */
    memcpy(Z, Q, sizeof(float) * M * M);

    squareRootMatrix(Z, M);

    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, M, K, M, 1.0, Z, M, W, K, 0.0, Q, M);
    memcpy(Z, Q, sizeof(float) * M * K);

    /*Z is formed, size (MxK)*/

    free(Q);
#undef Idx
}

void setSpectralWeights(float *W, float *E, float *U, float *D, int M)
{
    /* Computes transpose of  W=E^t U^T D^(-1/2), eqn (6),
     * as W^T=D^(-1/2)^T U E^t */
#define Idx(i, j) ((i)*global_config.nx*global_config.ny+(j))

    int i, j;
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < global_config.nx * global_config.ny; j++)
        {
            W[Idx(i, j)] = (float) (1 / sqrt(D[j]) * U[Idx(i, j)] * pow(E[i], global_config.t));
        }
    }

#undef Idx
}

void condenseEigenvalues(float *E, int m)
{
    /* get rid of imaginary components (should be zero) */
    int i;

    for (i = 1; i < m; i++)
    {
        E[i] = E[2 * i];
    }
}

void condenseEigenvectors(float *U, int m)
{
    /* Eigenvectors are returned from feast to include imaginary components
     * followed by the left eigenvector, this function aligns real components
     * sequentially.
     *
     * Sequentializes real components then memcpy each EV
     *
     * On exit, U only has relevant information on global_config.xy*global_config.ny*m elements
     */

#define DST(i, j) ((i)*global_config.nx*global_config.ny*2+(j))
#define SRC(i, j) ((i)*global_config.nx*global_config.ny*2+(j)*2)

    int i, j;

    for (i = 0; i < m; i++)
    {
        for (j = 1; j < global_config.nx * global_config.ny; j++)
        {
            U[DST(i, j)] = U[SRC(i, j)];
        }

        if (i != 0)
        {
            memmove(&U[DST(i, 0) / 2], &U[DST(i, 0)], sizeof(float) * global_config.nx * global_config.ny);
        }
    }

#undef DST
#undef SRC
}

int setEigenPairs(float *A, int *JA, int *IA, float *E, float *U)
{
    /* Feast in */
    int   fpm[64];    /* Feast configuration */
    int   loop;       /* Number of subspace iterations */
    int   M0;         /* Maximum EV count */
    float Emid[2];  /* Contour ellipse (complex) center */
    float r;        /* ellipse radius */
    int   N;      /* Problem size*/

    /* Feast out */
    float epsout;  /* Relative error (out) */
    int   M;          /* # of EV found in interval */
    int   info;       /* Error handling */
    float *res;    /*residual*/

    /* Set problem parameters */
    N = global_config.nx * global_config.ny;

    r = (1 - global_config.evcrit) / 2;
    Emid[0] = global_config.evcrit + r;
    r *= 1.0;

    M0 = (int) ((global_config.nx * global_config.ny) / global_config.maxevfact);

    /* Configure feast */
    fpm[0]  = 0;   /* verbose output */
    fpm[1]  = 16;
    fpm[3]  = 600; /* refinement */
    fpm[5]  = 1;   /* error type */
    fpm[17] = 5;   /* ellipse contour ratio (im/real) */

    feastinit(fpm);

    res = calloc((size_t) (2 * M0), sizeof(float));

    if (global_config.verbosity == 1)
    {
        fprintf(stderr, "Feast Initialized, assuming %d ev\n", M0);
    }

    sfeast_gcsrev(&N, A, IA, JA, fpm, &epsout, &loop, Emid, &r, &M0, E, U, &M, res, &info);

    if (global_config.verbosity == 1)
    {
        fprintf(stderr, "Feast Finished\n");
    }

    free(res);

    return M;

}

inline float affinityFunctionGaussian(float I1, float I2)
{
    return (float) exp(-(I1 - I2) * (I1 - I2) / (2 * global_config.sigma * global_config.sigma));
}

void setMarkovMatrix(float *A, int *JA, float *D, int nnz)
{
    /* Divide each column by its normalizing factor */
    int i, j;
#define FortIndex 1

    for (i = 0; i < nnz; i++)
    {
        j = JA[i] - FortIndex;
        A[i] /= D[j];
    }

#undef FortIndex
}

void setNormalizingVector(float *A, int *IA, float *D)
{
    /* sum along each column (or row, its symmetric), put in D */
    long i, j;
#define FortIndex 1

    for (i = 0; i < global_config.nx * global_config.ny; i++)
    {
        D[i] = 0;
        for (j = IA[i] - FortIndex; j < IA[i + 1] - FortIndex; j++)
        {
            D[i] += A[j];
        }
    }

#undef FortIndex
}

/* Put affinity matrix into sparse CSR format */
int setAffinityGaussian(float *data, float *A, int *JA, int *IA)
{
#define Data(x, y, z) data[global_config.nx*global_config.ny*(z)+global_config.nx*(y)+(x)]
#define Idx3(x, y, z)  (global_config.nx*global_config.ny*(z)+global_config.nx*(y)+(x))
#define SubX(n)     (((n)%(global_config.nx*global_config.ny))%global_config.nx)
#define SubY(n)     (((n)-SubX(n)%(global_config.nx*global_config.ny))/global_config.nx)
#define SubZ(n)     ((((n)-SubX(n)-global_config.nx*SubY(n))/(global_config.nx*global_config.ny)))
#define FortIndex 1 /* Use fortran indices for eigenvalue solver */
    //#define FortIndex

    int i, j, k;
    int ni, nj;
    int ak = 0, jk = 0, ik = 0;

    /* Doing one row at a time*/
    for (ni = 0; ni < global_config.nx * global_config.ny; ni++)
    {
        /* Row pointer */
        IA[ik++] = ak + FortIndex;

        /* ni is the current pixel */
        k = SubZ(ni);
        j = SubY(ni);
        i = SubX(ni);

        /* SOUTH WEST */
        if (((j - 1) >= 0) && ((i - 1) >= 0))
        {
            A[ak++]  = affinityFunctionGaussian(Data(i, j, 0), Data(i - 1, j - 1, 0));
            JA[jk++] = Idx3(i - 1, j - 1, 0) + FortIndex;
        }

        /* SOUTH */
        if (((j - 1) >= 0))
        {
            A[ak++]  = affinityFunctionGaussian(Data(i, j, 0), Data(i, j - 1, 0));
            JA[jk++] = Idx3(i, j - 1, 0) + FortIndex;
        }

        /* SOUTH EAST */
        if (((j - 1) >= 0) && ((i + 1) < global_config.nx))
        {
            A[ak++]  = affinityFunctionGaussian(Data(i, j, 0), Data(i + 1, j - 1, 0));
            JA[jk++] = Idx3(i + 1, j - 1, 0) + FortIndex;
        }

        /* WEST */
        if (((i - 1) >= 0))
        {
            A[ak++]  = affinityFunctionGaussian(Data(i, j, 0), Data(i - 1, j, 0));
            JA[jk++] = Idx3(i - 1, j, 0) + FortIndex;
        }

        /* CENTRAL */
        {
            A[ak++]  = affinityFunctionGaussian(Data(i, j, 0), Data(i, j, 0));
            JA[jk++] = Idx3(i, j, 0) + FortIndex;
        }

        /* EAST */
        if (((i + 1) < global_config.nx))
        {
            A[ak++]  = affinityFunctionGaussian(Data(i, j, 0), Data(i + 1, j, 0));
            JA[jk++] = Idx3(i + 1, j, 0) + FortIndex;
        }

        /* NORTH WEST */
        if (((j + 1) < global_config.ny) && ((i - 1) >= 0))
        {
            A[ak++]  = affinityFunctionGaussian(Data(i, j, 0), Data(i - 1, j + 1, 0));
            JA[jk++] = Idx3(i - 1, j + 1, 0) + FortIndex;
        }

        /* NORTH */
        if (((j + 1) < global_config.ny))
        {
            A[ak++]  = affinityFunctionGaussian(Data(i, j, 0), Data(i, j + 1, 0));
            JA[jk++] = Idx3(i, j + 1, 0) + FortIndex;
        }

        /* NORTH EAST */
        if (((j + 1) < global_config.ny) && ((i + 1) < global_config.nx))
        {
            A[ak++]  = affinityFunctionGaussian(Data(i, j, 0), Data(i + 1, j + 1, 0));
            JA[jk++] = Idx3(i + 1, j + 1, 0) + FortIndex;
        }

    }
    IA[ik++] = ak + FortIndex;
    return ak;

#undef Data
#undef Idx3
#undef SubZ
#undef SubY
#undef SubX
#undef FortIndex
}

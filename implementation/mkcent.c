#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "se.h"

/* Global configuration, initialized to defaults */
global_config_t global_config = {
        .inputData="data.data",
        .nx=64,
        .ny=64,
        .nw=64,
        .sigma=0.1,
        .evcrit=0.333,
        .t=600,
        .maxevfact=8,
        .kcent=16,
        .krep=8,
        .kiter=200,
        .verbosity=1,
        .taucardinality=0.89,
        .kelbw=0.002,
        .cardmin=10
};

void readE(float *data, char *name, size_t count)
{
    FILE *fp;

    if ((fp = fopen(name, "r")) == NULL)
        fprintf(stderr, "Could not open file %s\n", name);

    if (fread(data, sizeof(float), count, fp) !=
        (global_config.nx * global_config.ny * global_config.nw))
    {
        fprintf(stderr, "Could not find enough data in the file: %s\n", global_config.inputData);
    }

    fclose(fp);
}

int main(int argc, char *argv[])
{
    readConfig("seconfig");

    float *A;  /* Nonzero elements */
    int   *JA; /* Column indices */
    int   *IA; /* Row pointer */
    int   nnz; /* Number nonzero elements */

    float *D;    /* Normalization vector to change A into a transition probability matrix */

    float *E; /* Eigenvalues */
    float *U; /* Eigenvectors */
    int   M0; /* max number of eigenvalues */
    int   M;  /* Number of eigenvalues found */
    float *data = malloc(sizeof(float) * global_config.nx * global_config.ny);
    int   *O    = malloc(sizeof(int) * global_config.nx * global_config.ny);

    seed_region_t *id = malloc(sizeof(seed_region_t) * global_config.kcent);

    printConfig();

    A  = calloc((size_t) (global_config.nx * global_config.ny * 9), sizeof(float));
    D  = calloc((size_t) (global_config.nx * global_config.ny), sizeof(float));
    JA = calloc((size_t) (global_config.nx * global_config.ny * 9), sizeof(int));
    IA = calloc((size_t) (global_config.nx * global_config.ny + 1), sizeof(int));

    M = (int) strtol(argv[1], NULL, 10);
    global_config.kcent = (int) strtol(argv[2], NULL, 10);

    readE(data, global_config.inputData, (size_t) (global_config.nx * global_config.ny));

    int   nk;
    float tau;

    nnz = setAffinityGaussian(data, A, JA, IA);
    setNormalizingVector(A, IA, D);
    setMarkovMatrix(A, JA, D, nnz);

    float  *W    = calloc((size_t) (M * global_config.ny * global_config.nx), sizeof(float));
    float  *Z    = calloc((size_t) (M * global_config.ny * global_config.nx), sizeof(float));
    float  *R    = calloc((size_t) (global_config.nx * global_config.ny), sizeof(float));
    double *Cent = calloc((size_t) (M * global_config.kcent), sizeof(double));
    int    *Card = calloc((size_t) global_config.kcent, sizeof(int));
    float  *F    = calloc((size_t) (global_config.kcent * global_config.kcent), sizeof(int));
    U = calloc((size_t) (M * global_config.ny * global_config.nx), sizeof(float)); /* complex */
    E = calloc((size_t) M, sizeof(float));                     /* complex */

    readE(U, "U.data", (size_t) (M * global_config.ny * global_config.nx));
    readE(E, "E.data", (size_t) M);

    setSpectralWeights(W, E, U, D, M);
    setCoordinateTransform(Z, W, E, U, D, M);
    setTransformCenters(Z, R, O, M);

    tau = findThresholdOverlay(R, O);
    setThresholdOverlay(R, O, tau);
    setMinimumCardinality(R, O);
    setCenterCardinality(O, Card);
    setF(F, Card, Cent, M, tau);

    nk = setUniqueOrderedIndices(O);

    float *Cx = calloc(nk, sizeof(float));
    float *Cy = calloc(nk, sizeof(float));

    setSeedCentroids(O, Cx, Cy, nk);

    FILE *fp6  = fopen("O.data", "w");
    FILE *fp7  = fopen("Z.data", "w");
    FILE *fp8  = fopen("D.data", "w");
    FILE *fp9  = fopen("R.data", "w");
    FILE *fp10 = fopen("W.data", "w");
    FILE *fp11 = fopen("Cx.data", "w");
    FILE *fp12 = fopen("Cy.data", "w");
    fwrite(O, sizeof(int), (size_t) (global_config.nx * global_config.ny), fp6);
    fwrite(R, sizeof(float), (size_t) (global_config.nx * global_config.ny), fp9);
    fwrite(Z, sizeof(float), (size_t) (M * global_config.nx * global_config.ny), fp7);
    fwrite(D, sizeof(float), (size_t) (global_config.nx * global_config.ny), fp8);
    fwrite(W, sizeof(float), (size_t) (M * global_config.nx * global_config.ny), fp10);
    fwrite(Cx, sizeof(float), (size_t) nk, fp11);
    fwrite(Cy, sizeof(float), (size_t) nk, fp12);
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
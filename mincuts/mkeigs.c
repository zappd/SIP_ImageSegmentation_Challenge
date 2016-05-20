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

void readE(float *data)
{
    FILE *fp;

    if ((fp = fopen(global_config.inputData, "r")) == NULL)
        fprintf(stderr, "Could not open file %s\n", global_config.inputData);

    if (fread(data, sizeof(float), (size_t) (global_config.nx * global_config.ny * global_config.nw), fp) !=
        (global_config.nx * global_config.ny * global_config.nw))
    {
        fprintf(stderr, "Could not find enough data in the file: %s\n", global_config.inputData);
    }

    fclose(fp);
}

int main(int argc, char *argv[])
{
    readConfig("seconfig");

    global_config.evcrit    = strtof(argv[1], NULL);
    global_config.maxevfact = strtof(argv[2], NULL);

    float *A;  /* Nonzero elements */
    int   *JA; /* Column indices */
    int   *IA; /* Row pointer */
    int   nnz; /* Number nonzero elements */

    float *D; /* Normalization vector to change A into a transition probability matrix */

    float *E; /* Eigenvalues */
    float *U; /* Eigenvectors */
    int   M0; /* max number of eigenvalues */
    int   M;  /* Number of eigenvalues found */

    float *data = malloc(sizeof(float) * global_config.nx * global_config.ny);

    readE(data);

    printConfig();
    M0 = (int) ((global_config.nx * global_config.ny) / (global_config.maxevfact));

    A  = calloc((size_t) (global_config.nx * global_config.ny * 9), sizeof(float));
    D  = calloc((size_t) (global_config.nx * global_config.ny), sizeof(float));
    JA = calloc((size_t) (global_config.nx * global_config.ny * 9), sizeof(int));
    IA = calloc((size_t) (global_config.nx * global_config.ny + 1), sizeof(int));
    U  = calloc((size_t) (2 * 2 * M0 * global_config.ny * global_config.nx), sizeof(float)); /* complex */
    E  = calloc((size_t) (2 * M0), sizeof(float)); /* complex */

    nnz = setAffinityGaussian(data, A, JA, IA);
    setNormalizingVector(A, IA, D);
    setMarkovMatrix(A, JA, D, nnz);
    M = setEigenPairs(A, JA, IA, E, U);

    {
        fprintf(stderr, "Found %d eigenvalues\n", M);
    }

    /* Throw out imaginary components (should be 0) and left eigenvectors */
    condenseEigenvectors(U, M);
    condenseEigenvalues(E, M);
    normalizeEigenvectors(U, M);

    float *tmp;
    tmp = realloc(U, M * global_config.ny * global_config.nx * sizeof(float));
    U   = tmp;
    tmp = realloc(E, M * sizeof(float));
    E   = tmp;

    FILE *fp1 = fopen("A.data", "w");
    FILE *fp2 = fopen("JA.data", "w");
    FILE *fp3 = fopen("IA.data", "w");
    FILE *fp4 = fopen("E.data", "w");
    FILE *fp5 = fopen("U.data", "w");

    fwrite(A, sizeof(float), (size_t) nnz, fp1);
    fwrite(JA, sizeof(int), (size_t) nnz, fp2);
    fwrite(IA, sizeof(int), (size_t) (global_config.nx * global_config.ny + 1), fp3);
    fwrite(E, sizeof(float), (size_t) (M), fp4);
    fwrite(U, sizeof(float), (size_t) (M * global_config.nx * global_config.ny), fp5);

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
    fclose(fp5);

    return 0;
}

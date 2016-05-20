#include <math.h>

#include "map.h"
#include "config.h"


/********************************
 * Implementation
 ********************************/

void collapseHyperspectralMeans(float *img, float ***specimg)
{
    /* Naively collapses spectral field to scalar field by taking the mean across
     * each band */
    int    i, j, k;
    int    ir;
    double mean;

    for (i = 0, ir = 0; i < global_config.nx; i++)
    {
        for (j = 0; j < global_config.ny; j++, ir++)
        {
            mean   = 0;
            for (k = 0; k < global_config.nx; k++)
            {
                mean += specimg[i][j][k];
            }
            img[ir] = (float) (mean / (float) global_config.nw);
        }
    }
}

void normalizeMap(float *map)
{
    int    i;
    double sum = 0, sumsq = 0;

    for (i = 0; i < global_config.nx * global_config.ny; i++)
    {
        sum += map[i] / (global_config.nx * global_config.ny);
        sumsq += map[i] * map[i] / (global_config.nx * global_config.ny);
    }

    for (i = 0; i < global_config.nx * global_config.ny; i++)
    {
        map[i] -= sum;
        map[i] /= sqrt(sumsq - sum * sum);
    }

}

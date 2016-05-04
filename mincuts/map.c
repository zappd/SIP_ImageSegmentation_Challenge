#include <math.h>
#include <stdint.h>

#include "map.h"


void collapseHyperspectralMeans(float *img, float ***specimg, image_info_t *info)
{
    /* Naively collapses spectral field to scalar field by taking the mean across
     * each band */
    int    i, j, k;
    int    ir;
    double mean;

    for (i = 0, ir = 0; i < info->lines; i++)
    {
        for (j = 0; j < info->samples; j++, ir++)
        {
            mean   = 0;
            for (k = 0; k < info->bands; k++)
            {
                mean += specimg[i][j][k];
            }
            img[ir] = mean / (float) info->bands;
        }
    }
}

void normalizeMap(float *map, image_info_t *image_info)
{
    int    i;
    double sum = 0, sumsq = 0;

    for (i = 0; i < image_info->lines * image_info->samples; i++)
    {
        sum += map[i] / (image_info->lines * image_info->samples);
        sumsq += map[i] * map[i] / (image_info->lines * image_info->samples);
    }

    for (i = 0; i < image_info->lines * image_info->samples; i++)
    {
        map[i] -= sum;
        map[i] /= sqrt(sumsq - sum * sum);
    }

}

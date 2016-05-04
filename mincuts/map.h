#if !defined(_MAP_H)
#define _MAP_H

#include "readENVI.h"

void collapseHyperspectralMeans(float *img, float ***specimg, image_info_t *info);

void normalizeMap(float *map, image_info_t *image_info);

#endif

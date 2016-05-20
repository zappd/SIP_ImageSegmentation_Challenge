#ifndef MAP_H
#define MAP_H

#include "readENVI.h"

void collapseHyperspectralMeans(float *img, float ***specimg);

void normalizeMap(float *map);

#endif // MAP_H

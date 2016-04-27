#include "map.h"

void collapseHyperspectralMeans( float * img, float *** specimg, image_info_t * info)
{
  /* Naively collapses spectral field to scalar field by taking the mean across
   * each band */
  int i,j,k;
  int ir;
  double mean;

  for(i=0,ir=0;i<info->lines;i++)
  {
    for(j=0;j<info->samples;j++,ir++)
    {
      mean=0;
      for(k=0;k<info->bands;k++)
      {
        mean+=specimg[i][j][k];
      }
      img[ir]=mean/(float)info->bands;
    }
  }
}


#include <math.h>
#include "map.h"
#include "config.h"


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

void normalizeMap(float * map)
{
  int i;
  double sum=0,sumsq=0;

  for(i=0;i<gconf.nx*gconf.ny;i++)
  {
    sum+=map[i]/(gconf.nx*gconf.ny);
    sumsq+=map[i]*map[i]/(gconf.nx*gconf.ny);
  }

  for(i=0;i<gconf.nx*gconf.ny;i++)
  {
    map[i]-=sum;
    map[i]/=sqrt(sumsq-sum*sum);
  }

}

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "config.h"
#include "segmenter.h"
#include "se.h"
#include "seeds.h"
#include "map.h"

static const char *data_file_path = "StanfordMemorial_rat1_rot90_crop_bandCrop.img";
static const char *header_file_path = "StanfordMemorial_rat1_rot90_crop_bandCrop.img.hdr";
static const char *png_extension = ".png";

/* Global configuration, initialized to defaults */
struct gconf gconf = {
  .inputData="data.data",
  .nx=64,
  .ny=64,
  .nw=64,
  .sigma=0.1,
  .evcrit=0.333,
  .t=600,
  .maxevfact=8 ,
  .kcent=16,
  .krep=8,
  .kiter=200,
  .verbosity=1,
  .taucardinality=0.89,
  .kelbw=0.002,
  .cardmin=10
};

void readE(float * data)
{
  FILE * fp;
    
  if((fp=fopen(gconf.inputData,"r"))==NULL)
    fprintf(stderr,"Could not open file %s\n",gconf.inputData);

  if(fread(data,sizeof(float),gconf.nx*gconf.ny*gconf.nw,fp)!=
      (gconf.nx*gconf.ny*gconf.nw))
  {
    fprintf(stderr,"Could not find enough data in the file: %s\n",gconf.inputData);
  }

  fclose(fp);
}

int main(int argc, char* argv[])
{
  readConfig("seconfig");

	image_info_t image_info;
	char *segment_image_path = "./";

	float ***image_cube = readImageCube(data_file_path, header_file_path, &image_info);

  gconf.nx=image_info.lines;
  gconf.ny=image_info.samples;
  gconf.nw=image_info.bands;


  float * data=malloc(sizeof(float)*gconf.nx*gconf.ny);
  int * O   =malloc(sizeof(int)*gconf.nx*gconf.ny);
  seed_region_t * id  =malloc(sizeof(seed_region_t)*gconf.kcent);


	segment_image_path = (char *)malloc((strlen(data_file_path) + strlen(png_extension) + 1) * sizeof(char));
	strcpy(segment_image_path, data_file_path);
	strcat(segment_image_path, png_extension);

  collapseHyperspectralSingle(data,image_cube,&image_info,95);
  gconf.nw=1;
  normalizeMap(data);

  FILE * fp4=fopen("Map2.data","w");

  fwrite(data,sizeof(int),gconf.nx*gconf.ny,fp4);

  fclose(fp4);

  freeImageCube(image_cube, &image_info);

  return 0;
}

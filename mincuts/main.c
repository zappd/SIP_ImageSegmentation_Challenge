#include <stdio.h>

#include "config.h"
#include "segmenter.h"
// #include "graphCuts.h"

static const char *data_file_path = "/home/zappd/SIP_Challenge_Data/StanfordMemorial_rat1_rot90_crop_bandCrop.img";
static const char *header_file_path = "/home/zappd/SIP_Challenge_Data/StanfordMemorial_rat1_rot90_crop_bandCrop.img.hdr";

/* Global configuration, initialized to defaults */
struct gconf gconf = {
  .inputData="data.data",
  .nx=64,
  .ny=64,
  .nw=64
};

int main(int argc, char* argv[])
{
	image_info_t image_info;

	float ***image_cube = readImageCube(data_file_path, header_file_path, &image_info);

	printf("Read Finished\n\n");

	recursivelySegment(image_cube, &image_info);

	freeImageCube(image_cube, &image_info);

    // readConfig("seconfig");
    /*printConfig();*/
    
    return 0;
}
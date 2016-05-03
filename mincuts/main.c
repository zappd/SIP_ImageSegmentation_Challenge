#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "seeds.h"
#include "config.h"
#include "segmenter.h"
// #include "graphCuts.h"

//static const char *data_file_path   = "/home/zappd/SIP_Challenge_Data/StanfordMemorial_rat1_rot90_crop_bandCrop.img";
//static const char *header_file_path = "/home/zappd/SIP_Challenge_Data/StanfordMemorial_rat1_rot90_crop_bandCrop.img.hdr";
static const char *data_file_path   = "StanfordMemorial_rat1_rot90_crop_bandCrop.img";
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
        .maxevfact=8,
        .kcent=16,
        .krep=8,
        .kiter=200
};

int main(int argc, char *argv[])
{
    image_info_t image_info;
    char         *segment_image_path;
    uint32_t ** segment_id_map;
    seed_region_t * seed_region_arr;

    readConfig("seconfig");

    float ***image_cube = readImageCube(data_file_path, header_file_path, &image_info);

    printf("Read Finished\n\n");

    segment_image_path = (char *) malloc((strlen(data_file_path) + strlen(png_extension) + 1) * sizeof(char));
    segment_id_map=malloc(image_info.lines*sizeof(uint32_t *));


    for (uint32_t line = 0; line < image_info.lines; line++)
    {
        segment_id_map[line] = (uint32_t *) malloc(image_info.samples * sizeof(uint32_t));
    }

    strcpy(segment_image_path, data_file_path);
    strcat(segment_image_path, png_extension);

    getSeedIDFromFile("binid.data",segment_id_map);
    seed_region_arr=getSeedDataFromFile("Fnumadj_int.data","Fadjacent_int.data","Ftype_int.data");


    //cutSeedRegions(image_cube, image_info, segment_id_map, seed_region_arr, gconf.kcent);

    //recursivelySegment(image_cube, &image_info, segment_image_path);

    freeImageCube(image_cube, &image_info);
    freeSeeds(seed_region_arr);

    // readConfig("seconfig");
    /*printConfig();*/

    return 0;
}

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "seeds.h"

static const char *data_file_path   = "StanfordMemorial_rat1_rot90_crop_bandCrop.img";
static const char *header_file_path = "StanfordMemorial_rat1_rot90_crop_bandCrop.img.hdr";

static const char *segment_id_file_path      = "binid.data";
static const char *neighbor_length_file_path = "Fnumadj_int.data";
static const char *neighbor_file_path        = "Fadjacent_int.data";
static const char *neighbor_type_file_path   = "Ftype_int.data";

static const char *global_config_file = "seconfig";

static const char *png_extension = ".png";

/* Global configuration, initialized to defaults */
global_config_t global_config = {
        .inputData="data.data",
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
    char  *segment_image_path;
    float ***image_cube;

    image_info_t  image_info;
    seed_region_t *seed_region_arr;
    uint32_t      **segment_id_map;

    // load config values into global struct
    readConfig(global_config_file);

    // read in image cube from the specific path;
    image_cube = readImageCubeFromFile(data_file_path, header_file_path, &image_info);

    // make some space for the output PNG file path
    segment_image_path = (char *) malloc((strlen(data_file_path) + strlen(png_extension) + 1) * sizeof(char));

    // create the output PNG file path as a function of the image input
    strcpy(segment_image_path, data_file_path);
    strcat(segment_image_path, png_extension);

    // read in segment ids from binary file
    readSegmentIdsFromFile(segment_id_file_path, &segment_id_map, &image_info);

    // read in all seed region data from file
    seed_region_arr = getSeedDataFromFile(neighbor_length_file_path,
                                          neighbor_file_path,
                                          neighbor_type_file_path);

    // cut those seed regions!!!
    cutSeedRegions(image_cube, &image_info, segment_id_map, seed_region_arr, (uint32_t) global_config.kcent);

    freeImageCube(image_cube, &image_info);

    freeSeeds(seed_region_arr);

    return 0;
}

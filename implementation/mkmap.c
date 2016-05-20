#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "config.h"
#include "se.h"
#include "map.h"

static const char *data_file_path   = "StanfordMemorial_rat1_rot90_crop_bandCrop.img";
static const char *header_file_path = "StanfordMemorial_rat1_rot90_crop_bandCrop.img.hdr";

static const char *png_extension    = ".png";

/* Global configuration, initialized to defaults */
global_config_t global_config = {
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
        .kiter=200,
        .verbosity=1,
        .taucardinality=0.89,
        .kelbw=0.002,
        .cardmin=10
};

void readE(float *data)
{
    FILE *fp;

    if ((fp = fopen(global_config.inputData, "r")) == NULL)
        fprintf(stderr, "Could not open file %s\n", global_config.inputData);

    if (fread(data, sizeof(float), (size_t) (global_config.nx * global_config.ny * global_config.nw), fp) !=
        (global_config.nx * global_config.ny * global_config.nw))
    {
        fprintf(stderr, "Could not find enough data in the file: %s\n", global_config.inputData);
    }

    fclose(fp);
}

int main(int argc, char *argv[])
{
    readConfig("seconfig");

    image_info_t image_info;
    char         *segment_image_path = "./";

    float ***image_cube = readImageCubeFromFile(data_file_path, header_file_path, &image_info);

    global_config.nx = image_info.lines;
    global_config.ny = image_info.samples;
    global_config.nw = image_info.bands;


    float         *data = malloc(sizeof(float) * global_config.nx * global_config.ny);
    int           *O    = malloc(sizeof(int) * global_config.nx * global_config.ny);
    seed_region_t *id   = malloc(sizeof(seed_region_t) * global_config.kcent);


    segment_image_path = (char *) malloc((strlen(data_file_path) + strlen(png_extension) + 1) * sizeof(char));
    strcpy(segment_image_path, data_file_path);
    strcat(segment_image_path, png_extension);

    collapseHyperspectralSingle(data, image_cube, &image_info, 95);
    global_config.nw = 1;
    normalizeMap(data);

    FILE *fp4 = fopen("Map2.data", "w");

    fwrite(data, sizeof(int), global_config.nx * global_config.ny, fp4);

    fclose(fp4);

    freeImageCube(image_cube, &image_info);

    return 0;
}

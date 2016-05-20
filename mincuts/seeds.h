#include <stdint.h>

#include "config.h"
#include "readENVI.h"

#ifndef SEEDS_H
#define SEEDS_H

/********************************
 * Typedef Enums
 ********************************/

typedef enum
{
    AMBIGUOUS = 0,
    SINK      = 1,
    SOURCE    = 2
} neighbor_type_t;


/********************************
 * Typedef Structs
 ********************************/

typedef struct segment_bounds
{
    int32_t x_start;
    int32_t x_end;
    int32_t y_start;
    int32_t y_end;
} segment_bounds_t;

typedef struct neighbor
{
    uint32_t        segment_id;
    neighbor_type_t neighbor_type;
} neighbor_t;

typedef struct seed_region
{
    uint32_t   segment_id;
    uint32_t   neighbor_count;
    neighbor_t *neighbors;
} seed_region_t;


/********************************
 * Functions
 ********************************/

// TODO - make sure that seed region array comes ordered from lowest to highest segment id
void cutSeedRegions(float ***image, image_info_t *image_info, uint32_t **segment_id_map,
                    seed_region_t *seed_region_arr, uint32_t number_of_seed_regions);

void readSegmentIdsFromFile(const char *segment_id_file_path, uint32_t ***segment_id_map, image_info_t *image_info);

seed_region_t *getSeedDataFromFile(const char *neighbor_length_file_path, const char *neighbor_file_path,
                                   const char *neighbor_type_file_path);

void freeSeeds(seed_region_t *seeds_regions_to_free);

#endif // SEEDS_H

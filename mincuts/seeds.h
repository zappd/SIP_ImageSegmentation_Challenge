#ifndef MINCUTS_SEEDS_H
#define MINCUTS_SEEDS_H
#include <stdint.h>

#include "config.h"
#include "readENVI.h"




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
    neighbor_t *neighbors;
} seed_region_t;

/********************************
 * Functions
 ********************************/

// TODO - make sure that seed region array comes ordered from lowest to highest segment id
void cutSeedRegions(float ***image, image_info_t *image_info, uint32_t **segment_id_map, seed_region_t *seed_region_arr, uint32_t segments);

void getSeedIDFromFile(const char *, uint32_t **);


#endif // MINCUTS_SEEDS_H





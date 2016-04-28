
#include <stdint.h>
//
// Created by root on 4/27/16.
//

#ifndef MINCUTS_SEEDS_H
#define MINCUTS_SEEDS_H

typedef enum
{
    SOURCE,
    SINK,
    AMBIGUOUS
} neighbor_type_t;

typedef struct neighbor
{
    uint32_t segment_id;
    neighbor_type_t neighbor_type;
} neighbor_t;

typedef struct seed_region
{
    uint32_t segment_id;
    neighbor_t *neighbors;
} seed_region_t;

#endif //MINCUTS_SEEDS_H



#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

#include "graphCuts.h"
#include "math.h"
#include "config.h"
#include "seeds.h"
#include "segmentPngIO.h"


/********************************
 * Definitions
 ********************************/

#define NUM_OFFSETS 4
#define NUM_NEIGHBORS 8
#define INTERSTITIAL_SEGMENT_ID 0
#define UNUSED_NODE_ID -1
#define DEFAULT_SEGMENT_BOUNDS (segment_bounds_t){-1, -1, (samples + 1), -1}
#define MAXIMUM_NUM_SEGMENTS 10000

/********************************
 * Typedef Structs
 ********************************/

typedef struct segment_info
{
    uint32_t         segment_id;
    int32_t          background_segment_id;
    uint32_t         pixel_count;
    neighbor_type_t  neighbor_type;
    segment_bounds_t segment_bounds;
} segment_info_t;

typedef struct {
    uint32_t begin_line;
    uint32_t end_line;
    float *source_centroid;
    float *sink_centroid;
} terminal_weight_info_t;


/********************************
 * Functions
 ********************************/

static void normalizeImageCube();

static void calculateEdgeWeights();

static void calculateBandVariances();

static void initialize();

static void deallocate();

static void segmentImage();

static void addSegmentToNodeMap(uint32_t segment_id);

static void removeSegmentFromNodeMap(uint32_t segment_id);

static void setAllEdgeWeights();

static void calculateSegmentBounds();

static void resetNodeIdMap();

static uint32_t probeForSegmentId();

static void storeCut();

static void addEdgesFromNodeAt(uint32_t u_line, uint32_t u_sample);

static void removeEdgesFromNodeAt(uint32_t u_line, uint32_t u_sample);

/* Vector Functions */
static inline float hyperspectralDistanceFromOffset(int32_t line, int32_t sample, int32_t line_offset, int32_t sample_offset);

static inline float hyperspectralDistanceFromVectors(float *vector_1, float *vector_2);

static inline float addVectors(float *sum_vector, float *addend_vector);

/* Terminal Calculation */
static void setAllTerminalWeights();

static void *computeSubsetOfTerminalWeights(void *wi);

static void parallelComputeAmbiguousTerminalWeights(float *source_centroid, float *sink_centroid);

/********************************
 * Global Variables
 ********************************/

static const int X_OFFSETS[NUM_OFFSETS] = {1, 1, 1, 0};
static const int Y_OFFSETS[NUM_OFFSETS] = {-1, 0, 1, 1};

static const int X_NEIGHBORS[NUM_NEIGHBORS] = {0, 1, 1, 1, 0, -1, -1, -1};
static const int Y_NEIGHBORS[NUM_NEIGHBORS] = {-1, -1, 0, 1, 1, 1, 0, -1};

static uint32_t      **segment_ids;
static seed_region_t *seed_regions;

static float   ***image_cube;
static float   ***edge_weights;
static int32_t **node_ids;

static float *band_variances;

static uint32_t minimum_node_id = 0;
static uint32_t number_of_edges = 0;

static uint32_t number_of_segments;
static uint32_t lines;
static uint32_t samples;
static uint32_t bands;

static segment_info_t *segments_info;

static uint32_t **new_segment_ids;
static uint32_t number_of_new_segments;

static float ***ambiguous_terminal_storage;


/********************************
 * Implementation
 ********************************/


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ~ Seed Reading
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void readSegmentIdsFromFile(const char *segment_id_file_path, uint32_t ***segment_id_map, image_info_t *image_info)
{
    FILE *file_pointer = fopen(segment_id_file_path, "r");

    uint32_t number_of_pixels = image_info->lines * image_info->samples;
    *segment_id_map = (uint32_t **) malloc(image_info->lines * sizeof(uint32_t *));

    for (uint32_t line = 0; line < image_info->lines; line++)
    {
        fseek(file_pointer, (image_info->samples * sizeof(uint32_t) * line), SEEK_SET);

        (*segment_id_map)[line] = (uint32_t *) malloc(image_info->samples * sizeof(uint32_t));

        if ((fread((*segment_id_map)[line], sizeof(uint32_t), (size_t) image_info->samples, file_pointer)) !=
            image_info->samples)
        {
            fprintf(stderr, "Could not read all %d segment ids from file '%s'\n", number_of_pixels,
                    segment_id_file_path);
        }
    }

    fclose(file_pointer);
}

void freeSeeds(seed_region_t *seeds_regions_to_free)
{
    for (uint32_t i = 0; i < global_config.kcent; i++)
    {
        free(seeds_regions_to_free[i].neighbors);
    }

    free(seeds_regions_to_free);
}

seed_region_t *getSeedDataFromFile(const char *neighbor_length_file_path, const char *neighbor_file_path,
                                   const char *neighbor_type_file_path)
{
    FILE          *file_pointer;
    seed_region_t *seed_region_arr;
    uint32_t      *neighbors;
    uint32_t      *neighbor_types;
    uint32_t      sum_of_lengths = 0;

    // store the number of neighbors per node
    uint32_t *neigbors_per_segment = (uint32_t *) malloc(sizeof(uint32_t) * global_config.kcent);

    // open the neighbor lengths file
    file_pointer = fopen(neighbor_length_file_path, "r");

    if ((fread(neigbors_per_segment, sizeof(int), (size_t) global_config.kcent, file_pointer)) != global_config.kcent)
    {
        fprintf(stderr, "Could not read %d seed regions from file '%s'\n",
                global_config.kcent, neighbor_length_file_path);
    }

    fclose(file_pointer);

    for (uint32_t i = 0; i < global_config.kcent; i++)
    {
        sum_of_lengths += neigbors_per_segment[i];
    }

    // allocate memory for neighbor types
    neighbor_types = malloc(sizeof(uint32_t) * sum_of_lengths);

    file_pointer = fopen(neighbor_type_file_path, "r");

    if ((fread(neighbor_types, sizeof(int), sum_of_lengths, file_pointer)) != sum_of_lengths)
    {
        fprintf(stderr, "Could not read %d neighbor types from file '%s'\n", sum_of_lengths, neighbor_type_file_path);
    }

    fclose(file_pointer);


    // allocate memory for neighbors
    neighbors = malloc(sizeof(uint32_t) * sum_of_lengths);

    file_pointer = fopen(neighbor_file_path, "r");

    if ((fread(neighbors, sizeof(int), sum_of_lengths, file_pointer)) != sum_of_lengths)
    {
        fprintf(stderr, "Could not read %d neighbor types from file '%s'\n", sum_of_lengths, neighbor_file_path);
    }

    fclose(file_pointer);

    // allocate space for seed regions
    seed_region_arr = malloc(sizeof(seed_region_t) * global_config.kcent);

    for (uint32_t i = 0; i < global_config.kcent; i++)
    {
        seed_region_arr[i].neighbors      = (neighbor_t *) malloc(sizeof(neighbor_t) * neigbors_per_segment[i]);
        seed_region_arr[i].neighbor_count = neigbors_per_segment[i];
    }

    int           neighbor_index = 0;
    for (uint32_t i              = 0; i < global_config.kcent; i++)
    {
        seed_region_arr[i].segment_id = (uint32_t) i + 1;

        for (uint32_t j = 0; j < neigbors_per_segment[i]; j++, neighbor_index++)
        {
            seed_region_arr[i].neighbors[j].segment_id = (uint32_t) neighbors[neighbor_index];
            if (neighbor_types[neighbor_index] == 0)
            {
                seed_region_arr[i].neighbors[j].neighbor_type = AMBIGUOUS;
            }
            else if (neighbor_types[neighbor_index] == 1)
            {
                seed_region_arr[i].neighbors[j].neighbor_type = SINK;
            }
            else
            {
                seed_region_arr[i].neighbors[j].neighbor_type = SOURCE;
            }

        }
    }

    free(neigbors_per_segment);
    free(neighbor_types);
    free(neighbors);

    return seed_region_arr;
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ~ Cuts
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void cutSeedRegions(float ***image, image_info_t *image_info, uint32_t **segment_id_map,
                    seed_region_t *seed_region_arr,
                    uint32_t number_of_seed_regions)
{
    image_cube = image;

    lines   = image_info->lines;
    samples = image_info->samples;
    bands   = image_info->bands;

    segment_ids  = segment_id_map;
    seed_regions = seed_region_arr;

    number_of_segments     = number_of_seed_regions + 1;
    number_of_new_segments = number_of_segments;

    initialize();

    segmentImage();

    deallocate();
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ~ Private
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

static void initialize()
{
    // init band variances
    band_variances = (float *) calloc(bands, sizeof(float));

    // init node ids
    node_ids = (int32_t **) malloc(lines * sizeof(uint32_t *));
    for (uint32_t line = 0; line < lines; line++)
    {
        // initializes all pixels to node 0
        node_ids[line] = (int32_t *) calloc(samples, sizeof(uint32_t));
    }

    // init edge weights
    edge_weights = (float ***) malloc(lines * sizeof(float **));

    for (uint32_t line = 0; line < lines; line++)
    {
        edge_weights[line] = (float **) malloc(samples * sizeof(float *));

        for (uint32_t sample = 0; sample < samples; sample++)
        {
            edge_weights[line][sample] = (float *) malloc(NUM_OFFSETS * sizeof(float));
        }
    }

    // init segment info, add one info struct at the beginning to
    // hold info about the interstitial segment, segment 0
    segments_info = (segment_info_t *) calloc(MAXIMUM_NUM_SEGMENTS, sizeof(segment_info_t));
    for (uint32_t i = 0; i < MAXIMUM_NUM_SEGMENTS; i++)
    {
        segments_info[i].segment_id = i;
    }

    // init new segment ids and copy over old values
    new_segment_ids = (uint32_t **) malloc(lines * sizeof(uint32_t *));
    for (uint32_t line = 0; line < lines; line++)
    {
        // initializes all pixels to node 0
        new_segment_ids[line] = (uint32_t *) malloc(samples * sizeof(uint32_t));
        memcpy(new_segment_ids[line], segment_ids[line], samples * sizeof(uint32_t));
    }

    // init space for holding terminal values
    ambiguous_terminal_storage = (float ***) malloc(lines * sizeof(float **));

    for (uint32_t line = 0; line < lines; line++)
    {
        ambiguous_terminal_storage[line] = (float **) malloc(samples * sizeof(float *));

        for (uint32_t sample = 0; sample < samples; sample++)
        {
            ambiguous_terminal_storage[line][sample] = (float *) malloc(2 * sizeof(float));
        }
    }


    normalizeImageCube();

    calculateEdgeWeights();

    resetNodeIdMap();

    calculateSegmentBounds();
}

static void deallocate()
{
    // free band variances
    free(band_variances);

    // free node ids
    for (uint32_t line = 0; line < lines; line++)
    {
        free(node_ids[line]);
    }
    free(node_ids);

    // free edge weights
    for (uint32_t line = 0; line < lines; line++)
    {

        for (uint32_t sample = 0; sample < samples; sample++)
        {
            free(edge_weights[line][sample]);
        }
        free(edge_weights[line]);

    }
    free(edge_weights);

    // free segment info
    free(segments_info);

    // free new segment ids
    for (uint32_t line = 0; line < lines; line++)
    {
        free(new_segment_ids[line]);
    }
    free(new_segment_ids);

    // free ambiguous terminal storage
    for (uint32_t line = 0; line < lines; line++)
    {

        for (uint32_t sample = 0; sample < samples; sample++)
        {
            free(ambiguous_terminal_storage[line][sample]);
        }
        free(ambiguous_terminal_storage[line]);
    }
    free(ambiguous_terminal_storage);

}

static void normalizeImageCube()
{
    float maximum = 0.0f;

    // find maximum value
    for (uint32_t line = 0; line < lines; line++)
    {
        for (uint32_t sample = 0; sample < samples; sample++)
        {
            for (uint32_t band = 0; band < bands; band++)
            {
                if (image_cube[line][sample][band] > maximum)
                {
                    maximum = image_cube[line][sample][band];
                }
            }
        }
    }

    // normalize
    for (uint32_t line = 0; line < lines; line++)
    {
        for (uint32_t sample = 0; sample < samples; sample++)
        {
            for (uint32_t band = 0; band < bands; band++)
            {
                image_cube[line][sample][band] = image_cube[line][sample][band] / (maximum);
            }
        }
    }
}

static void calculateEdgeWeights()
{
    calculateBandVariances();

    for (int32_t line = 0; line < lines; line++)
    {
        for (int32_t sample = 0; sample < samples; sample++)
        {
            for (uint32_t offset_index = 0; offset_index < NUM_OFFSETS; offset_index++)
            {
                if ((line + X_OFFSETS[offset_index]) >= lines ||
                    (sample + Y_OFFSETS[offset_index]) >= samples ||
                    (line + X_OFFSETS[offset_index]) < 0 ||
                    (sample + Y_OFFSETS[offset_index]) < 0)
                {
                    continue;
                }

                edge_weights[line][sample][offset_index] =
                        hyperspectralDistanceFromOffset(line, sample, X_OFFSETS[offset_index], Y_OFFSETS[offset_index]);
            }
        }
    }
}

static void calculateBandVariances()
{
    float minus_mean;
    float minimum = FLT_MAX;

    float nodes          = (float) lines * (float) samples;
    float *band_averages = (float *) calloc(bands, sizeof(float));

    // calculate band averages
    for (uint32_t line = 0; line < lines; line++)
    {
        for (uint32_t sample = 0; sample < samples; sample++)
        {
            for (uint32_t band = 0; band < bands; band++)
            {
                band_averages[band] += image_cube[line][sample][band];
            }
        }
    }

    for (uint32_t band = 0; band < bands; band++)
    {
        band_averages[band] = (band_averages[band] / nodes);
    }

    // calculate band variances
    for (uint32_t line = 0; line < lines; line++)
    {
        for (uint32_t sample = 0; sample < samples; sample++)
        {
            for (uint32_t k = 0; k < bands; k++)
            {
                minus_mean = (image_cube[line][sample][k] - band_averages[k]);
                band_variances[k] += (minus_mean * minus_mean);
            }
        }
    }

    for (uint32_t band = 0; band < bands; band++)
    {
        band_variances[band] = band_variances[band] / nodes;
        minimum = fminf(minimum, band_variances[band]);
    }

    // create an offset for the logarithm to not go negative
    // minimum -= 10.1;

    for (uint32_t band = 0; band < bands; band++)
    {
        // take natural log of function, then prep it for division
        band_variances[band] = band_variances[band] / minimum;
    }
}

static void resetNodeIdMap()
{
    minimum_node_id = 0;
    number_of_edges = 0;

    segments_info[0].neighbor_type = AMBIGUOUS;

    for (uint32_t line = 0; line < lines; line++)
    {
        for (uint32_t sample = 0; sample < samples; sample++)
        {
            if (segment_ids[line][sample] == INTERSTITIAL_SEGMENT_ID)
            {
                node_ids[line][sample] = minimum_node_id++;
                addEdgesFromNodeAt(line, sample);
            }
            else
            {
                node_ids[line][sample] = UNUSED_NODE_ID;
            }
        }
    }
}

static void calculateSegmentBounds()
{
    uint32_t segment_id;

    // we are going to be accessing segments a lot
    // let's not make a map or anything crazy like that,
    // but pre-process to decrease the range of the search

    // init segment bounds to defaults
    for (uint32_t segment = 0; segment < MAXIMUM_NUM_SEGMENTS; segment++)
    {
        segments_info[segment].segment_bounds = DEFAULT_SEGMENT_BOUNDS;
    }

    for (int32_t line = 0; line < lines; line++)
    {
        for (int32_t sample = 0; sample < samples; sample++)
        {
            segment_id = segment_ids[line][sample];

            segment_bounds_t *segment_bounds = &segments_info[segment_id].segment_bounds;

            // test x start
            if (segment_bounds->x_start == -1)
            {
                segment_bounds->x_start = line;
            }

            // set x end
            segment_bounds->x_end = line;

            // test y start
            if (sample < segment_bounds->y_start)
            {
                segment_bounds->y_start = sample;
            }

            // test y end
            if (sample > segment_bounds->y_end)
            {
                segment_bounds->y_end = sample;
            }
        }
    }
}

static void segmentImage()
{
    neighbor_t *neighbors;

    clock_t start;

    // ignore segment zero because it is not a seed region
    for (uint32_t seed_index = 1; seed_index < number_of_segments; seed_index++)
    {
        neighbors = seed_regions[seed_index - 1].neighbors;
        uint32_t number_of_neighbors = sizeof(neighbors) / sizeof(neighbor_t);

        printf("Processing Segment %d\n", seed_index);

        for (uint32_t seed = 0; seed < number_of_neighbors; seed++)
        {
            // clean slate
            resetNodeIdMap();

            // add this seed region to the node map
            segments_info[seed_index].neighbor_type = SOURCE;
            addSegmentToNodeMap(seed_index);

            // we always want all sink nodes included.
            // we want to try all combinations of including surrounding
            // source nodes.
            // all ambiguous nodes should always be left included

            // find number of source neighbors
            uint32_t source_neighbors = 0;
            uint32_t *source_indices  = (uint32_t *) calloc(number_of_neighbors, sizeof(uint32_t));

            for (uint32_t i = 0; i < number_of_neighbors; i++)
            {
                if (neighbors[i].neighbor_type == SOURCE)
                {
                    // find all of the indices at which source nodes are
                    // stored in the neighbors array
                    source_indices[source_neighbors++] = i;
                }
                else
                {
                    // find all other nodes in the neighbor list and assign them node
                    // ids in the node id map
                    addSegmentToNodeMap(neighbors[i].segment_id);
                }

                // create a global reference for the neighbor type of this segment relative to the
                // segment currently under investigation
                segments_info[neighbors[i].segment_id].neighbor_type = neighbors[i].neighbor_type;
            }

            // assumes there are fewer than 33 neighbors
            assert(source_neighbors < 33);

            uint32_t mask_max = (uint32_t) (1 << source_neighbors);

            for (uint32_t neighbor_mask = 0x00000000; neighbor_mask < mask_max; neighbor_mask++)
            {
                // add all of the included source regions
                for (uint32_t i = 0; i < source_neighbors; i++)
                {
                    if (neighbor_mask & (1 << source_neighbors))
                    {
                        addSegmentToNodeMap(neighbors[source_indices[i]].segment_id);
                    }
                }

                /* START: Perform Algorithm */

                initializeGraph(minimum_node_id, number_of_edges);

                start = clock();
                setAllTerminalWeights();
                printf("SATW: %ld ms\n", (clock() - start));

                setAllEdgeWeights();

                computeMaximumFlow(false, NULL, 0);

                storeCut();

                writeSegmentImage(new_segment_ids,
                                  "/media/zappd/Data/Dropbox/School/ECE697SP/Challenge/SIP_ImageSegmentation_Challenge/mincuts/node_ids_dbg.png",
                                  lines, samples);

                /* END: Perform Algorithm */

                // remove all included source regions
                for (uint32_t i = 0; i < source_neighbors; i++)
                {
                    if (neighbor_mask & (1 << source_neighbors))
                    {
                        removeSegmentFromNodeMap(neighbors[source_indices[i]].segment_id);
                    }
                }

            }

            free(source_indices);
        }
    }

    writeSegmentImage(new_segment_ids,
                      "/media/zappd/Data/Dropbox/School/ECE697SP/Challenge/SIP_ImageSegmentation_Challenge/mincuts/StanfordMemorial_rat1_rot90_crop_bandCrop.png",
                      lines, samples);

    // we are done!
    printf("Done\n");
}


static void addSegmentToNodeMap(uint32_t segment_id)
{
    segment_bounds_t bounds = segments_info[segment_id].segment_bounds;

    assert(bounds.x_start >= 0);
    assert(bounds.x_end >= 0);
    assert(bounds.y_start <= samples);
    assert(bounds.y_end >= 0);

    for (uint32_t line = (uint32_t) bounds.x_start; line < bounds.x_end; line++)
    {
        for (uint32_t sample = (uint32_t) bounds.y_start; sample < bounds.y_end; sample++)
        {
            if (segment_ids[line][sample] == segment_id)
            {
                node_ids[line][sample] = minimum_node_id++;
                addEdgesFromNodeAt(line, sample);
            }
        }
    }
}

static void removeSegmentFromNodeMap(uint32_t segment_id)
{
    segment_bounds_t bounds = segments_info[segment_id].segment_bounds;

    assert(bounds.x_start >= 0);
    assert(bounds.x_end >= 0);
    assert(bounds.y_start <= samples);
    assert(bounds.y_end >= 0);

    for (uint32_t line = (uint32_t) bounds.x_start; line < bounds.x_end; line++)
    {
        for (uint32_t sample = (uint32_t) bounds.y_start; sample < bounds.y_end; sample++)
        {
            if (segment_ids[line][sample] == segment_id)
            {
                node_ids[line][sample] = UNUSED_NODE_ID;

                // this is a somewhat sketchy thing to do...
                // it's efficacy stems from the constraint that whenever we are removing nodes from the map
                // we do it in bulk, so that after all nodes that will be removed have been, there has
                // been no fragmentation in node id assignment
                minimum_node_id--;

                removeEdgesFromNodeAt(line, sample);
            }
        }
    }
}

static void addEdgesFromNodeAt(uint32_t u_line, uint32_t u_sample)
{
    int32_t line   = (int32_t) u_line;
    int32_t sample = (int32_t) u_sample;

    for (uint32_t neighbor_index = 0; neighbor_index < NUM_NEIGHBORS; neighbor_index++)
    {
        if ((line + X_NEIGHBORS[neighbor_index]) >= lines ||
            (sample + Y_NEIGHBORS[neighbor_index]) >= samples ||
            (line + X_NEIGHBORS[neighbor_index]) < 0 ||
            (sample + Y_NEIGHBORS[neighbor_index]) < 0)
        {
            continue;
        }

        if (node_ids[line + X_NEIGHBORS[neighbor_index]][sample + Y_NEIGHBORS[neighbor_index]] >= 0)
        {
            number_of_edges++;
        }
    }
}

static void removeEdgesFromNodeAt(uint32_t u_line, uint32_t u_sample)
{
    int32_t line   = (int32_t) u_line;
    int32_t sample = (int32_t) u_sample;

    for (uint32_t neighbor_index = 0; neighbor_index < NUM_NEIGHBORS; neighbor_index++)
    {
        if ((line + X_NEIGHBORS[neighbor_index]) >= lines ||
            (sample + Y_NEIGHBORS[neighbor_index]) >= samples ||
            (line + X_NEIGHBORS[neighbor_index]) < 0 ||
            (sample + Y_NEIGHBORS[neighbor_index]) < 0)
        {
            continue;
        }

        if (node_ids[line + X_NEIGHBORS[neighbor_index]][sample + Y_NEIGHBORS[neighbor_index]] >= 0)
        {
            number_of_edges--;
        }
    }
}


static void setAllEdgeWeights()
{
    for (int32_t line = 0; line < lines; line++)
    {
        for (int32_t sample = 0; sample < samples; sample++)
        {
            if (node_ids[line][sample] >= 0)
            {
                for (uint32_t offset_index = 0; offset_index < NUM_OFFSETS; offset_index++)
                {
                    if ((line + X_OFFSETS[offset_index]) >= lines ||
                        (sample + Y_OFFSETS[offset_index]) >= samples ||
                        (line + X_OFFSETS[offset_index]) < 0 ||
                        (sample + Y_OFFSETS[offset_index]) < 0)
                    {
                        continue;
                    }

                    if (node_ids[line + X_OFFSETS[offset_index]][sample + Y_OFFSETS[offset_index]] >= 0)
                    {
                        setEdgeWeight((uint32_t) node_ids[line][sample],
                                      (uint32_t) node_ids[line + X_OFFSETS[offset_index]][sample +
                                                                                          Y_OFFSETS[offset_index]],
                                      edge_weights[line][sample][offset_index],
                                      edge_weights[line][sample][offset_index]);
                    }
                }
            }
        }
    }
}

static uint32_t probeForSegmentId()
{
    for (uint32_t i = number_of_segments; i < MAXIMUM_NUM_SEGMENTS; i++)
    {
        if (segments_info[i].pixel_count == 0)
        {
            return i;
        }
    }

    // if this value gets returned, we have bigger problem anyway
    return 0;
}

static void storeCut()
{
    uint32_t fore = 0, back = 0;

    // reinitialize background segment ids for all segments
    for (uint32_t i = 0; i < MAXIMUM_NUM_SEGMENTS; i++)
    {
        segments_info[i].background_segment_id = -1;
        segments_info[i].pixel_count           = 0;
    }

    for (uint32_t line = 0; line < lines; line++)
    {
        for (uint32_t sample = 0; sample < samples; sample++)
        {
            if (node_ids[line][sample] >= 0 &&
                segments_info[segment_ids[line][sample]].neighbor_type == AMBIGUOUS)
            {
                if (getTerminal((uint32_t) node_ids[line][sample]) == BACKGROUND)
                {
                    if (segments_info[new_segment_ids[line][sample]].background_segment_id == -1)
                    {
                        segments_info[new_segment_ids[line][sample]].background_segment_id = probeForSegmentId();

                        assert(number_of_new_segments < MAXIMUM_NUM_SEGMENTS);
                    }

                    new_segment_ids[line][sample] = (uint32_t) segments_info[new_segment_ids[line][sample]].background_segment_id;
                    back++;
                }
                else
                {
                    fore++;
                }
            }

            segments_info[new_segment_ids[line][sample]].pixel_count++;
        }
    }
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ~ Terminal Calculation
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#define SOURCE_INDEX 0
#define SINK_INDEX 1
#define NUM_TERMINAL_WEIGHT_THREADS 8

static void setAllTerminalWeights()
{
    float *source_centroid = (float *) calloc(bands, sizeof(float));
    float *sink_centroid   = (float *) calloc(bands, sizeof(float));

    uint32_t source_count = 0;
    uint32_t sink_count   = 0;

    clock_t start;

    start = clock();
    // set all source and sink weights, and tally up vectors for
    for (uint32_t line = 0; line < lines; line++)
    {
        for (uint32_t sample = 0; sample < samples; sample++)
        {
            if (node_ids[line][sample] >= 0)
            {
                switch (segments_info[segment_ids[line][sample]].neighbor_type)
                {
                    case SOURCE:
                        addVectors(source_centroid, image_cube[line][sample]);
                        source_count++;

                        setTerminalWeights((uint32_t) node_ids[line][sample], FLT_MAX, 0);
                        break;

                    case SINK:
                        addVectors(sink_centroid, image_cube[line][sample]);
                        sink_count++;

                        setTerminalWeights((uint32_t) node_ids[line][sample], 0, FLT_MAX);
                        break;

                    case AMBIGUOUS:
                        break;

                    default:
                        assert(0);
                }
            }
        }
    }

    printf("1: %ld ms\t", (clock() - start));

    start = clock();
    for (uint32_t band = 0; band < bands; band++)
    {
        sink_centroid[band]   = sink_centroid[band] / sink_count;
        source_centroid[band] = source_centroid[band] / source_count;
    }
    printf("2: %ld ms\t", (clock() - start));

    start = clock();

    parallelComputeAmbiguousTerminalWeights(source_centroid, sink_centroid);

    for (uint32_t line = 0; line < lines; line++)
    {
        for (uint32_t sample = 0; sample < samples; sample++)
        {
            if (node_ids[line][sample] >= 0)
            {
                if (segments_info[segment_ids[line][sample]].neighbor_type == AMBIGUOUS)
                {
                    setTerminalWeights((uint32_t) node_ids[line][sample], ambiguous_terminal_storage[line][sample][SOURCE_INDEX], ambiguous_terminal_storage[line][sample][SINK_INDEX]);
                }
            }
        }
    }

    printf("3: %ld ms\t", (clock() - start));
}

static void parallelComputeAmbiguousTerminalWeights(float *source_centroid, float *sink_centroid)
{
    pthread_t threads[NUM_TERMINAL_WEIGHT_THREADS];
    terminal_weight_info_t weights_info[NUM_TERMINAL_WEIGHT_THREADS];

    uint32_t lines_per_thread = lines / NUM_TERMINAL_WEIGHT_THREADS;

    for (uint8_t thread_index = 0; thread_index < NUM_TERMINAL_WEIGHT_THREADS; thread_index++)
    {
        weights_info[thread_index].begin_line = thread_index * lines_per_thread;
        weights_info[thread_index].end_line = (thread_index == (NUM_TERMINAL_WEIGHT_THREADS - 1)) ? lines : ((thread_index + 1) * lines_per_thread);
        weights_info[thread_index].source_centroid = source_centroid;
        weights_info[thread_index].sink_centroid = sink_centroid;

        if(pthread_create(&threads[thread_index], NULL, computeSubsetOfTerminalWeights, (void *)(&weights_info[thread_index])))
        {
            fprintf(stderr, "Error creating thread %d\n", thread_index);
        }
    }

    // wait for all threads to finish
    for (uint8_t thread_index = 0; thread_index < NUM_TERMINAL_WEIGHT_THREADS; thread_index++)
    {
        if(pthread_join(threads[thread_index], NULL))
        {
            fprintf(stderr, "Error joining thread %d\n", thread_index);
        }
    }
}

static void *computeSubsetOfTerminalWeights(void *wi)
{
    terminal_weight_info_t *weight_info = (terminal_weight_info_t *)wi;

    for (uint32_t line = weight_info->begin_line; line < weight_info->end_line; line++)
    {
        for (uint32_t sample = 0; sample < samples; sample++)
        {
            if (node_ids[line][sample] >= 0)
            {
                if (segments_info[segment_ids[line][sample]].neighbor_type == AMBIGUOUS)
                {
                    ambiguous_terminal_storage[line][sample][SOURCE_INDEX] =
                            hyperspectralDistanceFromVectors(weight_info->source_centroid, image_cube[line][sample]);

                    ambiguous_terminal_storage[line][sample][SINK_INDEX] =
                            hyperspectralDistanceFromVectors(weight_info->sink_centroid, image_cube[line][sample]);
                }
            }
        }
    }
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ~ Vector Manipulation
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

static inline float hyperspectralDistanceFromOffset(int32_t line, int32_t sample, int32_t line_offset,
                                                    int32_t sample_offset)
{
    float temp_1, temp_2;
    float distance = 0;

    for (uint32_t band = 0; band < bands; band++)
    {
        temp_1 = image_cube[line + line_offset][sample + sample_offset][band] - image_cube[line][sample][band];
        temp_2 = sqrtf(fabsf(line_offset) + fabsf(sample_offset));

        distance += powf(temp_1 * temp_1 * temp_2, band_variances[band]);
    }

    return expf(((-4) * sqrtf(distance)) / (2 * global_config.sigma * global_config.sigma));
}

static inline float hyperspectralDistanceFromVectors(float *vector_1, float *vector_2)
{
    float temp_1;
    float distance = 0;

    for (uint32_t band = 0; band < bands; band++)
    {
        temp_1 = vector_1[band] - vector_2[band];
        distance += powf(temp_1 * temp_1, band_variances[band]);
    }

    return expf(((-1) * powf(distance, 0.2)) / (2 * global_config.sigma * global_config.sigma));
}

static inline float addVectors(float *sum_vector, float *addend_vector)
{
    for (uint32_t band = 0; band < bands; band++)
    {
        sum_vector[band] += addend_vector[band];
    }
}
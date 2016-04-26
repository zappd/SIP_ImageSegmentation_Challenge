#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "segmenter.h"
#include "graphCuts.h"
#include "segment_png_io.h"

/********************************
 * Definitions
 ********************************/

#define NUM_NEIGHBORS 8

/********************************
 * Functions
 ********************************/

static void initializeSegmenter();
static void deallocateSegmenter();
static void segment(float ***image_cube, uint32_t segment_id);
static void calculateEdgeWeights(float ***image_cube);
static void calculateBandAverages(float ***image_cube, uint32_t segment_id, uint32_t number_of_nodes);


/********************************
 * Global Vars
 ********************************/

static int X_OFFSETS[NUM_NEIGHBORS] = {0, 1, 1, 1, 0, -1, -1, -1};
static int Y_OFFSETS[NUM_NEIGHBORS] = {-1, -1, 0, 1, 1, 1, 0, -1};

static int segment_count = 0;

static uint32_t lines = 0;
static uint32_t samples = 0;
static uint32_t bands = 0;

static float ***edge_weights;
static float *band_averages;
static uint32_t **segment_ids;
static uint32_t **node_ids;

void recursivelySegment(float ***image_cube, image_info_t *image_info, char *segment_image_file_path)
{
    lines = image_info->lines;
    samples = image_info->samples;
    bands = image_info->bands;

    initializeSegmenter();

    calculateEdgeWeights(image_cube);

    segment(image_cube, 0);

    if (segment_image_file_path != NULL)
    {
        writeSegmentImage(segment_ids, segment_image_file_path, lines, samples);
    }

    deallocateSegmenter();
}

static void initializeSegmenter()
{
    // init band averages
    band_averages = (float *)malloc(bands * sizeof(float));

    // init segment_ids
    segment_ids = (uint32_t **)malloc(lines * sizeof(uint32_t *));
    for (uint32_t line = 0; line < lines; line++)
    {
        // initializes all pixels to node 0
        segment_ids[line] = (uint32_t *)calloc(samples, sizeof(uint32_t));
    }

    // init node ids
    node_ids = (uint32_t **)malloc(lines * sizeof(uint32_t*));
    for (uint32_t line = 0; line < lines; line++)
    {
        // initializes all pixels to node 0
        node_ids[line] = (uint32_t *)calloc(samples, sizeof(uint32_t));
    }

    // init edge weights
    edge_weights = (float ***)malloc(lines * sizeof(float **));

    for (uint32_t line = 0; line < lines; line++)
    {
        edge_weights[line] = (float **)malloc(samples * sizeof(float *));

        for (uint32_t sample = 0; sample < samples; sample++)
        {
            edge_weights[line][sample] = (float *)malloc(NUM_NEIGHBORS * sizeof(float));
        }
    }
}

static void deallocateSegmenter()
{
    // free band averages
    free(band_averages);


    // free segment_ids
    for (uint32_t line = 0; line < lines; line++)
    {
        free(segment_ids[line]);
    }
    free(segment_ids);


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
}


static void segment(float ***image_cube, uint32_t segment_id)
{
    uint32_t number_of_nodes = 0;
    uint32_t number_of_edges = 0;

    // START: Get Nodes in segment
    for (uint32_t line = 0; line < lines; line++)
    {
        for (uint32_t sample = 0; sample < samples; sample++)
        {
            if (segment_ids[line][sample] == segment_id)
            {
                node_ids[line][sample] = number_of_nodes++;

                for (uint8_t offset_index = 0; offset_index < NUM_NEIGHBORS; offset_index++)
                {
                    if ((line + X_OFFSETS[offset_index]) >= lines ||
                        (sample + Y_OFFSETS[offset_index]) >= samples ||
                        (line + X_OFFSETS[offset_index]) < 0 ||
                        (sample + Y_OFFSETS[offset_index]) < 0)
                    {
                        continue;
                    }

                    if (segment_ids[line + X_OFFSETS[offset_index]][sample + Y_OFFSETS[offset_index]] != segment_id)
                    {
                        continue;
                    }

                    number_of_edges++;
                }
            }
        }
    }
    // END: Get Nodes in segment

    calculateBandAverages(image_cube, segment_id, number_of_nodes);

    initializeGraph(number_of_nodes, number_of_edges);

    uint32_t node_id = 0;
    double weight = 0, source = 0, sink = 0;

    for (uint32_t line = 0; line < lines; line++)
    {
        for (uint32_t sample = 0; sample < samples; sample++)
        {
            if (segment_ids[line][sample] != segment_id)
            {
                continue;
            }

            weight = 0;
            node_id = node_ids[line][sample];

            for (uint32_t band = 0; band < bands; band++)
            {
                weight += (band_averages[band] - image_cube[line][sample][band]) / band_averages[band];
            }

            source = powf(1.005, weight);
            sink = powf(1.005, (-1) * weight);

            setTerminalWeights(node_id, source, sink);

            // only iterate over the subset of offsets which concern new edges
            for (uint32_t offset_index = 1; offset_index < (NUM_NEIGHBORS - 3); offset_index++)
            {
                if ((line + X_OFFSETS[offset_index]) >= lines ||
                    (sample + Y_OFFSETS[offset_index]) >= samples ||
                    (line + X_OFFSETS[offset_index]) < 0 ||
                    (sample + Y_OFFSETS[offset_index]) < 0)
                {
                    continue;
                }

                if (segment_ids[line + X_OFFSETS[offset_index]][sample + Y_OFFSETS[offset_index]] != segment_id)
                {
                    continue;
                }

                setEdgeWeight(node_id, node_ids[line + X_OFFSETS[offset_index]][sample + Y_OFFSETS[offset_index]],
                              edge_weights[line][sample][offset_index], edge_weights[line][sample][offset_index]);
            }
        }
    }

    computeMaximumFlow(false, NULL, 0);

    if (number_of_nodes != 0) {
        uint32_t fore = 0, back = 0;

        // START: Partition by Fore/Back-ground
        for (uint32_t line = 0; line < lines; line++)
        {
            for (uint32_t sample = 0; sample < samples; sample++)
            {
                if (segment_ids[line][sample] == segment_id)
                {
                    segment_ids[line][sample] = (getTerminal(node_ids[line][sample]) == BACKGROUND) ? (2 * segment_id + 1) : (2 * segment_id + 2);

                    if (segment_ids[line][sample] == (2 * segment_id + 1))
                    {
                        back++;
                    }
                    else
                    {
                        fore++;
                    }
                }
            }

        }
        // END: Partition by Fore/Back-ground

        printf("%d\n", ++segment_count);

        if (fore > 0 && back > 0)
        {
            segment(image_cube, 2 * segment_id + 1);
            segment(image_cube, 2 * segment_id + 2);
        }
    }
}

static void calculateEdgeWeights(float ***image_cube)
{
    float distance = 0;

    for (uint32_t line = 0; line < lines; line++)
    {
        for (uint32_t sample = 0; sample < samples; sample++)
        {

            for (uint32_t offset_index = 0; offset_index < NUM_NEIGHBORS; offset_index++)
            {
                if ((line + X_OFFSETS[offset_index]) >= lines ||
                    (sample + Y_OFFSETS[offset_index]) >= samples ||
                    (line + X_OFFSETS[offset_index]) < 0 ||
                    (sample + Y_OFFSETS[offset_index]) < 0)
                {
                    continue;
                }

                distance = 0;

                for (int band = 0; band < bands; band++)
                {
                    distance += sqrtf( fabsf(X_OFFSETS[offset_index]) + fabsf(Y_OFFSETS[offset_index]) ) * powf(image_cube[line + X_OFFSETS[offset_index]][sample + Y_OFFSETS[offset_index]][band] - image_cube[line][sample][band], 2);
                }

                edge_weights[line][sample][offset_index] = powf( 1.005 , (-1) * powf(distance, 0.2) );
            }
        }
    }
}

static void calculateBandAverages(float ***image_cube, uint32_t segment_id, uint32_t number_of_nodes)
{
    for (uint32_t band = 0; band < bands; band++)
    {
        band_averages[band] = 0;
    }

    for (uint32_t line = 0; line < lines; line++)
    {
        for (uint32_t sample = 0; sample < samples; sample++)
        {
            if (segment_ids[line][sample] == segment_id)
            {
                for (uint32_t band = 0; band < bands; band++)
                {
                    band_averages[band] += image_cube[line][sample][band];
                }
            }
        }
    }

    for (uint32_t band = 0; band < bands; band++)
    {
        band_averages[band] = (band_averages[band] / number_of_nodes);
    }
}
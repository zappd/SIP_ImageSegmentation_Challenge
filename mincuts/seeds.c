#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "graphCuts.h"
#include "math.h"
#include "config.h"
#include "seeds.h"


/********************************
 * Definitions
 ********************************/

#define NUM_OFFSETS 4
#define NUM_NEIGHBORS 8
#define INTERSTITIAL_SEGMENT_ID 0
#define UNUSED_NODE_ID -1
#define DEFAULT_SEGMENT_BOUNDS (segment_bounds_t){-1, -1, (samples + 1), -1}


/********************************
 * Typedef Structs
 ********************************/

typedef struct segment_info
{
    uint32_t         segment_id;
    neighbor_type_t  neighbor_type;
    segment_bounds_t segment_bounds;
} segment_info_t;


/********************************
 * Functions
 ********************************/

static void calculateEdgeWeights();

static void calculateBandVariances();

static void initialize();

static void deallocate();

static void segment();

static void addSegmentToNodeMap(uint32_t segment_id);

static void removeSegmentFromNodeMap(uint32_t segment_id);

static void setAllEdgeWeights();

static void setAllTerminalWeights();

static void calculateSegmentBounds();

static void resetNodeIdMap();

static void storeCut();

static void splitSegment(uint32_t parent_segment_id, uint32_t new_segment_id);

static void flagAffectedSegments();

static void addEdgesFromNodeAt(uint32_t u_line, uint32_t u_sample);

static void removeEdgesFromNodeAt(uint32_t u_line, uint32_t u_sample);

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


/********************************
 * Implementation
 ********************************/

void cutSeedRegions(float ***image, image_info_t *image_info, uint32_t **segment_id_map,
                    seed_region_t *seed_region_arr,
                    uint32_t segments)
{
    image_cube = image;

    lines   = image_info->lines;
    samples = image_info->samples;
    bands   = image_info->bands;

    segment_ids  = segment_id_map;
    seed_regions = seed_region_arr;

    number_of_segments     = segments;
    number_of_new_segments = segments;

    initialize();

    segment();

    deallocate();
}

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

    // init segment info
    segments_info = (segment_info_t *) calloc(number_of_segments, sizeof(segment_info_t));
    for (uint32_t i = 0; i < number_of_segments; i++)
    {
        segments_info[i].segment_id = i;
    }

    // init new segment ids and copy over old values
    new_segment_ids = (uint32_t **) malloc(lines * sizeof(uint32_t *));
    for (uint32_t line = 0; line < lines; line++)
    {
        // initializes all pixels to node 0
        new_segment_ids[line] = (uint32_t *) calloc(samples, sizeof(uint32_t));
        memcpy(new_segment_ids[line], segment_ids[line], samples * sizeof(uint32_t));
    }

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
}

static void calculateEdgeWeights()
{
    float distance = 0, temp_1, temp_2;

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

                distance = 0;

                for (uint32_t band = 0; band < bands; band++)
                {
                    temp_1 = image_cube[line + X_OFFSETS[offset_index]][sample + Y_OFFSETS[offset_index]][band] -
                             image_cube[line][sample][band];
                    temp_2 = sqrtf(fabsf(X_OFFSETS[offset_index]) + fabsf(Y_OFFSETS[offset_index]));
                    distance += temp_1 * temp_1 * temp_2 * band_variances[band];
                }

                edge_weights[line][sample][offset_index] = expf(
                        ((-1) * distance * distance) /
                        (2 * gconf.sigma * gconf.sigma));
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
    minimum -= expf(1.05);

    for (uint32_t band = 0; band < bands; band++)
    {
        // take natural log of function, then prep it for division
        band_variances[band] = 1 / log10f(band_variances[band] - minimum);
    }
}

static void resetNodeIdMap()
{
    minimum_node_id = 0;
    number_of_edges = 0;

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
    // we are going to be accessing segments a lot
    // let's not make a map or anything crazy like that,
    // but pre-process to decrease the range of the search

    // init segment bounds to defaults
    for (uint32_t segment = 0; segment < number_of_segments; segment++)
    {
        segments_info[segment].segment_bounds = DEFAULT_SEGMENT_BOUNDS;
    }

    for (uint32_t line = 0; line < lines; line++)
    {
        for (uint32_t sample = 0; sample < samples; sample++)
        {
            uint32_t segment_id = segment_ids[line][sample];

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

static void segment()
{
    neighbor_t *neighbors;

    for (uint32_t segment = 0; segment < number_of_segments; segment++)
    {
        neighbors = seed_regions[segment].neighbors;
        uint32_t number_of_neighbors = sizeof(neighbors) / sizeof(neighbor_t);

        for (uint32_t seed = 0; seed < number_of_neighbors; seed++)
        {
            // clean slate
            resetNodeIdMap();

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
                    // ids in the node is map
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

                setAllTerminalWeights();

                setAllEdgeWeights();

                computeMaximumFlow(false, NULL, 0);

                storeCut();

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
        }
    }

    // we are done!
    printf("Done");
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
    int32_t line = (int32_t)u_line;
    int32_t sample = (int32_t)u_sample;

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
    int32_t line = (int32_t)u_line;
    int32_t sample = (int32_t)u_sample;

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

static void setAllTerminalWeights()
{
    for (uint32_t line = 0; line < lines; line++)
    {
        for (uint32_t sample = 0; sample < samples; sample++)
        {
            if (node_ids[line][sample] >= 0)
            {
                switch (segments_info[segment_ids[line][sample]].neighbor_type)
                {
                    case SOURCE:
                        setTerminalWeights((uint32_t) node_ids[line][sample], FLT_MAX, 0);
                        break;
                    case SINK:
                        setTerminalWeights((uint32_t) node_ids[line][sample], 0, FLT_MAX);
                        break;
                    case AMBIGUOUS:
                        setTerminalWeights((uint32_t) node_ids[line][sample], 0, 0);
                        break;
                    default:
                        assert(0);
                }
            }
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


#define UNKNOWN_CLASS -1
#define FOREGROUND_CLASS FOREGROUND
#define BACKGROUND_CLASS BACKGROUND
#define SPLIT_CLASS 2

// holds one of the classes defined above
static int8_t *segment_class;

static void storeCut()
{
    segment_class = (int8_t *) malloc(number_of_new_segments * sizeof(int8_t));
    memset(segment_class, UNKNOWN_CLASS, number_of_new_segments);

    uint32_t new_segment_counter = number_of_new_segments;

    // find all segments affected by this cut
    flagAffectedSegments();

    // we now know how many new segments we will be adding, so allocate more memory
    segments_info = realloc(segments_info, number_of_new_segments * sizeof(segment_info_t));

    // for affected segments, add them and update information
    for (uint32_t segment_id = 0; segment_id < number_of_new_segments; segment_id++)
    {
        if (segment_class[segment_id] == 2)
        {
            new_segment_counter++;
            segments_info[new_segment_counter].segment_id = new_segment_counter;

            splitSegment(segment_id, new_segment_counter);
        }
    }

    // update number of segments to reflect the number of new segments created this round
    number_of_new_segments = new_segment_counter;

    free(segment_class);
}


static void flagAffectedSegments()
{
    for (uint32_t line = 0; line < lines; line++)
    {
        for (uint32_t sample = 0; sample < samples; sample++)
        {
            if (node_ids[line][sample] >= 0)
            { // this node was included
                switch (segment_class[new_segment_ids[line][sample]])
                {
                    case UNKNOWN_CLASS:
                        segment_class[new_segment_ids[line][sample]] = getTerminal(
                                (uint32_t) node_ids[line][sample]);
                        break;

                    case BACKGROUND_CLASS:
                    case FOREGROUND_CLASS:
                        if (getTerminal((uint32_t) node_ids[line][sample]) !=
                            segment_class[new_segment_ids[line][sample]])
                        {
                            segment_class[new_segment_ids[line][sample]] = 2;

                            // we have discovered a new segment which will be split:
                            // increment the segment counter
                            number_of_new_segments++;
                        }
                        break;

                    case SPLIT_CLASS:
                        break;

                    default:
                        assert(0);
                }
            }
        }
    }
}


static void splitSegment(uint32_t parent_segment_id, uint32_t new_segment_id)
{
    // copy old segment bounds
    segment_bounds_t old_parent_bounds = segments_info[parent_segment_id].segment_bounds;

    assert(old_parent_bounds.x_start >= 0);
    assert(old_parent_bounds.x_end >= 0);
    assert(old_parent_bounds.y_start <= samples);
    assert(old_parent_bounds.y_end >= 0);

    for (uint32_t line = (uint32_t) old_parent_bounds.x_start;
         line < old_parent_bounds.x_end; line++)
    {
        for (uint32_t sample = (uint32_t) old_parent_bounds.y_start;
             sample < old_parent_bounds.y_end; sample++)
        {
            if (segment_ids[line][sample] == parent_segment_id)
            {
                if (getTerminal((uint32_t) node_ids[line][sample]) == FOREGROUND)
                {
                    // this segment bisected by the most recent cut.
                    // if nodes of this segment are in the background, they maintain their current segment id
                    // if they are in the foreground, then they get the lowest available segment id
                    new_segment_ids[line][sample] = new_segment_id;
                }
            }
        }
    }

    // reset segment bounds
    segments_info[parent_segment_id].segment_bounds = DEFAULT_SEGMENT_BOUNDS;
    segments_info[new_segment_id].segment_bounds    = DEFAULT_SEGMENT_BOUNDS;

    // look within the old segment bounds and recalculate only bounds for the two
    // relevant segments
    for (uint32_t line = (uint32_t) old_parent_bounds.x_start;
         line < old_parent_bounds.x_end; line++)
    {
        for (uint32_t sample = (uint32_t) old_parent_bounds.y_start;
             sample < old_parent_bounds.y_end; sample++)
        {
            if (new_segment_ids[line][sample] != parent_segment_id &&
                new_segment_ids[line][sample] != new_segment_id)
            {
                continue;
            }

            segment_bounds_t *segment_bounds = &segments_info[new_segment_ids[line][sample]].segment_bounds;

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

void getSeedIDFromFile(const char *nm, uint32_t **dest)
{
    FILE *fp;
    int i, j, ir;

    int *src = (int *) malloc(sizeof(int) * gconf.nx * gconf.ny);

    fp = fopen(nm, "r");

    if ((fread(src, sizeof(int), (size_t)(gconf.nx * gconf.ny), fp)) != gconf.nx * gconf.ny)
    {
        fprintf(stderr, "Could not slurp %d ints\n", gconf.nx * gconf.ny);
    }

    fclose(fp);

    for (i = 0, ir=0; i < gconf.nx; i++)
    {
        for (j = 0; j < gconf.ny; j++, ir++)
        {
            dest[j][i] = (uint32_t) src[ir];
        }
    }

    free(src);
}

void freeSeeds( seed_region_t * sr)
{

  for(int i=0;i<gconf.kcent;i++)
  {
    free(sr[i].neighbors);
  }
  free(sr);

}

seed_region_t * getSeedDataFromFile(const char * lengthfile, const char * neighborfile, const char * typefile)
{

  int sum=0,i,j;
  FILE * fp;
  seed_region_t * seed_arr;

  int *lengths=(int *)malloc(sizeof(int)*gconf.kcent); /* number of neighbors per node */

  fp=fopen(lengthfile,"r");

  if((fread(lengths, sizeof(int),gconf.kcent,fp))!=gconf.kcent)
  {
    fprintf(stderr, "Could not slurp %d ints\n", gconf.nx*gconf.ny);
  }

  fclose(fp);

  for(i=0;i<gconf.kcent;i++)
  {
    sum+=lengths[i];
  }

  int *ntype=malloc(sizeof(int)*sum); /* neighbor type */
  int *neigh=malloc(sizeof(int)*sum); /* neighbor index */

  fp=fopen(typefile,"r");

  if((fread(ntype, sizeof(int),sum,fp))!=sum)
  {
    fprintf(stderr, "Could not slurp %d ints\n", gconf.nx*gconf.ny);
  }

  fclose(fp);


  fp=fopen(neighborfile,"r");

  if((fread(neigh, sizeof(int),sum,fp))!=sum)
  {
    fprintf(stderr, "Could not slurp %d ints\n", gconf.nx*gconf.ny);
  }

  fclose(fp);


  seed_arr =malloc(sizeof(seed_region_t)*gconf.kcent);

  for(i=0;i<gconf.kcent;i++)
  {
    seed_arr[i].neighbors=malloc(sizeof(neighbor_t)*lengths[i]);
    seed_arr[i].neighbor_count=(uint32_t)lengths[i];
  }

  int gindex=0;
  for(i=0;i<gconf.kcent;i++)
  {
    seed_arr[i].segment_id =(uint32_t) i+1;
    for(j=0;j<lengths[i];j++,gindex++)
    {
      seed_arr[i].neighbors[j].segment_id=(uint32_t)neigh[gindex];
      if(ntype[gindex]==0)
      {
        seed_arr[i].neighbors[j].neighbor_type=AMBIGUOUS;
      }
      else if(ntype[gindex]==1)
      {
        seed_arr[i].neighbors[j].neighbor_type=SINK;
      }
      else
      {
        seed_arr[i].neighbors[j].neighbor_type=SOURCE;
        }

    }
  }
  free(lengths);
  free(ntype);
  free(neigh);


  return seed_arr;

}



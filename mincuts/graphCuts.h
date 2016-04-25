#include <stdint.h>
#include <stdbool.h>

/********************************
 * Definitions
 ********************************/

#define FOREGROUND 1
#define BACKGROUND 0

#define DELTA 1E-10


/********************************
 * Typedef Structs
 ********************************/

typedef struct node
{
    void *first_outgoing_edge; // edge_t *
    void *parent_edge; // edge_t *

    void *next; // node_t *

    uint32_t timestamp;
    uint32_t distance;

    bool is_in_sink;
    bool is_marked;
    bool is_in_changed_list;

    float residual_capacity;
} node_t;

typedef struct orphan
{
    node_t *this;
    void *previous; // orphan_t *
    void *next; // orphan_t *
} orphan_t;

typedef struct edge
{
    node_t *head;
    void *next; // edge_t *
    void *sister_edge; // edge_t *

    float residual_capacity;
} edge_t;



/********************************
 * Functions
 ********************************/

void initializeGraph(uint32_t num_nodes, uint32_t num_edges);

void setTerminalWeights(uint32_t node_id, float source, float sink);

void setEdgeWeight(uint32_t node_id_1, uint32_t node_id_2, float weight_to, float weight_from);

bool getTerminal(uint32_t node_id);

float computeMaximumFlow(bool reuse_trees, uint32_t *changed_nodes, uint32_t *number_changed_nodes);
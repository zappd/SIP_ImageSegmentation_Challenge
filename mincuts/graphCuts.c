#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "graphCuts.h"


/********************************
 * Functions
 ********************************/

static void setNodeActive(node_t *node);

static node_t *getNextActiveNode();

static void addOrphanAtFront(node_t *node);

static void addOrphanAtBack(node_t *node);

static void clearOrphanList(void);

static void addToChangedList(node_t *node);

static void addOrphanToFrontOfQueue(node_t *node);

static void addOrphanToBackOfQueue(node_t *node);

static node_t *pollOrphanQueue(void);

static void initializeMaxFlow();

static void initializeMaxFlowAndReuseTrees();

static void augment(edge_t *middle_edge);

static void processSourceOrphan(node_t *orphan_node);

static void processSinkOrphan(node_t *orphan_node);

static bool edgeHasResidualCapacity(edge_t *edge);

static bool nodeHasResidualCapacity(node_t *node);

/********************************
 * Global Vars
 ********************************/

static edge_t ORPHAN_EDGE = {0};
static edge_t TERMINAL_EDGE = {0};

// number of nodes in the graph
static uint32_t number_of_nodes = 0;

// number of edges in the graph
static uint32_t number_of_edges = 0;

// counter to keep track of edge creation
static uint32_t current_edge_number = 0;

// the total flow through the whole graph
static float total_flow = 0.0f;

// counter for the numbers of iterations to maxflow
static uint32_t maxflow_iteration = 0;

// counter for iterations of loop
static uint32_t time = 0;

// array of all nodes in the graph
static node_t *nodes = NULL;

// array of all edges in the graph
static edge_t *edges = NULL;

static uint32_t nodes_size = 0;
static uint32_t edges_size = 0;

// head of queue_1 of all active nodes
static node_t *node_queue_1_head = NULL;

// tail of queue_1 of all active nodes
static node_t *node_queue_1_tail = NULL;

// head of queue_2 of all active nodes
static node_t *node_queue_2_head = NULL;

// tail of queue_2 of all active nodes
static node_t *node_queue_2_tail = NULL;



void initializeGraph(uint32_t num_nodes, uint32_t num_edges)
{
    number_of_nodes = num_nodes;
    number_of_edges = num_edges * 2;

    current_edge_number = 0;

    total_flow = 0.0f;

    maxflow_iteration = 0;

    time = 0;

    node_queue_1_head = NULL;
    node_queue_1_tail = NULL;
    node_queue_2_head = NULL;
    node_queue_2_tail = NULL;

    if (nodes == NULL)
    {
        // if this is the first time around, allocate memory
        nodes = (node_t *)malloc(number_of_nodes * sizeof(node_t));

        nodes_size = number_of_nodes;
    }
    else if (number_of_nodes > nodes_size)
    {
        // if this is not the first time around, and we need more memory,
        // reallocate memory
        realloc(nodes, (number_of_edges * 2) * sizeof(edge_t));
        nodes_size = number_of_nodes;
    }


    if (edges == NULL)
    {
        // if this is the first time around, allocate memory
        edges = (edge_t *)malloc((number_of_edges * 2) * sizeof(edge_t));

        edges_size = number_of_edges;
    }
    else if (number_of_edges > edges_size)
    {
        // if this is not the first time around, and we need more memory,
        // reallocate memory
        realloc(edges, (number_of_edges * 2) * sizeof(edge_t));

        edges_size = number_of_edges;
    }

    // clear out the needed memory for both nodes and edges
    memset(nodes, 0, number_of_nodes * sizeof(node_t));
    memset(edges, 0, (number_of_edges * 2) * sizeof(edge_t));
}


/**
 * Set the affinity for one node to belong to the foreground (i.e., source)
 * or background (i.e., sink).
 *
 * @param node_id The number of the node.
 * @param source The affinity of this node to the foreground (i.e., source)
 * @param sink   The affinity of this node to the background (i.e., sink)
 */
void setTerminalWeights(uint32_t node_id, float source, float sink)
{
    assert(node_id >= 0 && node_id < number_of_nodes);

    float delta = nodes[node_id].residual_capacity;

    if (delta > 0)
    {
        source += delta;
    }
    else
    {
        sink -= delta;
    }

    total_flow += (source < sink) ? source : sink;

    nodes[node_id].residual_capacity = (source - sink);
}


/**
 * Set the edge weight of a pair of directed edges between two nodes.
 *
 * Please note that you cannot call any <tt>setEdgeWeight</tt> more often
 * than the number of edges you specified at the time of construction!
 *
 * @param node_id_1    The first node.
 * @param node_id_2    The second node.
 * @param weight_from The weight (i.e., the cost) of the directed edge from node1 to node2.
 * @param weight_to The weight (i.e., the cost) of the directed edge from node2 to node1.
 */
void setEdgeWeight(uint32_t node_id_1, uint32_t node_id_2, float weight_to, float weight_from) {
    assert(node_id_1 >= 0 && node_id_1 < number_of_nodes);
    assert(node_id_2 >= 0 && node_id_2 < number_of_nodes);
    assert(node_id_1 != node_id_2);
    assert(weight_from >= 0);
    assert(weight_to >= 0);
    assert(current_edge_number < (number_of_edges - 2));

    // create new edges
    edge_t *edge_to = &edges[current_edge_number++];
    edge_t *edge_from = &edges[current_edge_number++];

    // get corresponding nodes
    node_t *node_1 = &nodes[node_id_1];
    node_t *node_2 = &nodes[node_id_2];


    // link edges
    edge_to->sister_edge = edge_from;
    edge_from->sister_edge = edge_to;

    // add node1 to edge
    edge_to->next = (edge_t *)node_1->first_outgoing_edge;
    node_1->first_outgoing_edge = edge_to;

    // add node2 to reverseEdge
    edge_from->next = (edge_t *)node_2->first_outgoing_edge;
    node_2->first_outgoing_edge = edge_from;

    // set targets of edges
    edge_to->head = node_2;
    edge_from->head = node_1;

    // set residual capacities
    edge_to->residual_capacity = weight_to;
    edge_from->residual_capacity = weight_from;
}

/**
 * Get the segmentation, i.e., the terminal node that is connected to the
 * specified node. If there are several min-cut solutions, free nodes are
 * assigned to the background.
 *
 * @param node_id the node to check
 * @return Either FOREGROUND or BACKGROUND
 */
bool getTerminal(uint32_t node_id)
{
    assert(node_id >= 0 && node_id < number_of_nodes);

    node_t *node = &nodes[node_id];

    return (node->parent_edge != NULL) ? node->is_in_sink : (bool)BACKGROUND;
}


/**
 * Marks a node as being active and adds it to second queue of active nodes.
 */
static void setNodeActive(node_t *node)
{
    if (node->next == NULL)
    {
        if (node_queue_2_tail != NULL)
        {
            node_queue_2_tail->next = node;
        }
        else
        {
            node_queue_2_head = node;
        }

        node_queue_2_tail = node;
        node->next = node;
    }
}


/**
 * Gets the next active node, that is, the first node of the first queue of
 * active nodes. If this queue is empty, the second queue is used. Returns
 * NULL, if no active node is left.
 */
static node_t *getNextActiveNode()
{
    node_t *node;

    while (true)
    {
        if ((node = node_queue_1_head) == NULL) {

            // queue 0 was empty, try other one
            node = node_queue_2_head;

            // swap queues
            node_queue_1_head = node_queue_2_head;
            node_queue_1_tail  = node_queue_2_tail;
            node_queue_2_head = NULL;
            node_queue_2_tail  = NULL;

            // if other queue was emtpy as well, return null
            if (node == NULL)
            {
                return NULL;
            }
        }

        // remove current node from active list
        if (node->next == node) {
            // this was the last one
            node_queue_1_head = NULL;
            node_queue_1_tail  = NULL;
        }
        else
        {
            node_queue_1_head = node->next;
        }

        // not in any list anymore
        node->next = NULL;

        // return only if it has a parent and is therefore active
        if (node->parent_edge != NULL)
        {
            return node;
        }
    }
}



/**********************************************
 * Orphan Stuff
 **********************************************/

static orphan_t *orphan_head_node = NULL;
static orphan_t *orphan_tail_node = NULL;

/**
 * Mark a node as orphan and add it to the front of the queue.
 */
static void addOrphanAtFront(node_t *node)
{
    node->parent_edge = &ORPHAN_EDGE;
    addOrphanToFrontOfQueue(node);
}

/**
 * Mark a node as orphan and add it to the back of the queue.
 */
static void addOrphanAtBack(node_t *node)
{
    node->parent_edge = &ORPHAN_EDGE;
    addOrphanToBackOfQueue(node);
}

static void clearOrphanList(void)
{
    if (orphan_head_node == NULL)
    {
    	orphan_head_node = (orphan_t *)calloc(1, sizeof(orphan_t));
    	orphan_tail_node = orphan_head_node;
    }

    orphan_t *current_orphan = orphan_head_node->next;
    orphan_t *next_orphan = NULL;

    while (current_orphan != NULL)
    {
        next_orphan = current_orphan->next;
        free(current_orphan);
        current_orphan = next_orphan;
    }
}

/**
 * Add a node to the list of potentially changed nodes.
 */
static void addToChangedList(node_t *node)
{
    node->is_in_changed_list = true;
}

static void addOrphanToFrontOfQueue(node_t *node)
{
    if (orphan_head_node == NULL)
    {
    	orphan_head_node = (orphan_t *)calloc(1, sizeof(orphan_t));
    	orphan_tail_node = orphan_head_node;
    }

    if (orphan_head_node->this == NULL)
    {
        orphan_head_node->this = node;
        return;
    }
    else
    {
        orphan_t *old_head_node = orphan_head_node;
        orphan_head_node = (orphan_t *)calloc(1, sizeof(orphan_t));
        orphan_head_node->this = node;
        orphan_head_node->next = old_head_node;
        old_head_node->previous = orphan_head_node;
    }
}

static void addOrphanToBackOfQueue(node_t *node)
{
    if (orphan_head_node == NULL)
    {
    	orphan_head_node = (orphan_t *)calloc(1, sizeof(orphan_t));
    	orphan_tail_node = orphan_head_node;
    }

    if (orphan_tail_node->this == NULL)
    {
        // list is empty, set tail orphan
        orphan_tail_node->this = node;
        return;
    }
    else
    {
        orphan_t *old_tail_node = orphan_tail_node;
        orphan_tail_node = (orphan_t *)calloc(1, sizeof(orphan_t));
        orphan_tail_node->this = node;
        orphan_head_node->previous = old_tail_node;
        old_tail_node->next = orphan_tail_node;
    }
}

static node_t *pollOrphanQueue(void)
{
    if (orphan_head_node == NULL)
    {
    	orphan_head_node = (orphan_t *)calloc(1, sizeof(orphan_t));
    	orphan_tail_node = orphan_head_node;
    }

    node_t *polled_node = orphan_head_node->this;

    if (orphan_head_node->next != NULL)
    {
        orphan_t *old_head_node = orphan_head_node;
        orphan_head_node = orphan_head_node->next;
        free(old_head_node);
    }
    else
    {
        // preserve the head / tail node, jus clear it out
        *orphan_head_node = (orphan_t){.this = NULL, .previous = NULL, .next = NULL};
    }

    return polled_node;
}


/********************************
 * The Algorithm
 ********************************/

/**
 * Initialise the algorithm.
 *
 * Only called if reuse_trees is false.
 */
static void initializeMaxFlow() 
{
    node_queue_1_head = NULL;
    node_queue_1_tail = NULL;
    node_queue_2_head = NULL;
    node_queue_2_tail = NULL;

    time = 0;

    clearOrphanList();

    uint32_t i;
    node_t *node;
    for (i = 0; i < number_of_nodes; i++)
    {
        node = &nodes[i];

        node->next = NULL;
        node->is_marked = false;
        node->is_in_changed_list = false;
        node->timestamp = time;

        if (node->residual_capacity > DELTA)
        {
            // node is connected to source
            node->is_in_sink = false;
            node->parent_edge = &TERMINAL_EDGE;
            node->distance = 1;

            setNodeActive(node);
        }
        else if (node->residual_capacity < -DELTA)
        {
            node->is_in_sink = true;
            node->parent_edge = &TERMINAL_EDGE;
            node->distance = 1;
            setNodeActive(node);
        }
        else
        {
            node->parent_edge = NULL;
        }
    }
}


/**
 * Initialise the algorithm for reuse of trees
 *
 * Only called if reuse_trees is true.
 */
static void initializeMaxFlowAndReuseTrees() 
{

    node_t *node_1 = NULL;
    node_t *node_2 = NULL;

    node_t *queue_start = node_queue_2_head;

    edge_t *edge = NULL;

    node_queue_1_head = NULL;
    node_queue_1_tail = NULL;
    node_queue_2_head = NULL;
    node_queue_2_tail = NULL;

    time++;

    clearOrphanList();

    while ((node_1 = queue_start) != NULL)
    {
        queue_start = node_1->next;

        if (queue_start == node_1)
        {
            queue_start = NULL;
        }

        node_1->next = NULL;
        node_1->is_marked = false;

        setNodeActive(node_1);

        if (!nodeHasResidualCapacity(node_1))
        {
            if (node_1->parent_edge != NULL)
            {
                addOrphanAtBack(node_1);
            }

            continue;
        }

        if (node_1->residual_capacity > DELTA)
        {
            if (node_1->parent_edge == NULL || node_1->is_in_sink)
            {
                node_1->is_in_sink = false;

                for (edge = node_1->first_outgoing_edge; edge != NULL; edge = edge->next)
                {
                    node_2 = edge->head;

                    if (!node_2->is_marked)
                    {
                        if (node_2->parent_edge == edge->sister_edge)
                        {
                            addOrphanAtBack(node_2);
                        }

                        if (node_2->parent_edge != NULL &&
                                node_2->is_in_sink &&
                                edge->residual_capacity > DELTA)
                        {
                            setNodeActive(node_2);
                        }
                    }
                }

                addToChangedList(node_1);
            }
        }
        else if (node_1->parent_edge == NULL || !node_1->is_in_sink)
        {
            node_1->is_in_sink = true;

            for (edge = node_1->first_outgoing_edge; edge != NULL; edge = edge->next)
            {
                node_2 = edge->head;

                if (!node_2->is_marked)
                {
                    if (node_2->parent_edge == edge->sister_edge)
                    {
                        addOrphanAtBack(node_2);
                    }

                    if (node_2->parent_edge != NULL &&
                            !node_2->is_in_sink &&
                            edge->sister_edge > 0)
                    {
                        setNodeActive(node_2);
                    }
                }
            }

            addToChangedList(node_1);
        }

        node_1->parent_edge = &TERMINAL_EDGE;
        node_1->timestamp = time;
        node_1->distance = 1;
    }

    // adoption

    node_t *current_node = NULL;
    while ((current_node = pollOrphanQueue()) != NULL)
    {

        if (current_node->is_in_sink)
        {
            processSinkOrphan(current_node);
        }
        else
        {
            processSourceOrphan(current_node);
        }
    }
}


/**
 * Perform the augmentation step of the graph cut algorithm.
 *
 * This is done whenever a path between the source and the sink was found.
 */
static void augment(edge_t *middle_edge) 
{
    node_t *node;
    edge_t *edge;

    float bottleneck;

    // 1. find bottleneck capacity

    // 1a - the source tree
    bottleneck = middle_edge->residual_capacity;

    for (node = ((edge_t *)middle_edge->sister_edge)->head; ; node = edge->head)
    {
        edge = node->parent_edge;

        if (edge == &TERMINAL_EDGE)
		{
            break;
    	}
        if (bottleneck > ((edge_t *)edge->sister_edge)->residual_capacity)
        {
            bottleneck = ((edge_t *)edge->sister_edge)->residual_capacity;
        }
    }

    if (bottleneck > node->residual_capacity)
    {
        bottleneck = node->residual_capacity;
    }

    // 1b - the sink tree
    for (node = middle_edge->head; ; node = edge->head) {

        edge = node->parent_edge;

        if (edge == &TERMINAL_EDGE)
        {
            break;
        }

        if (bottleneck > edge->residual_capacity)
        {
            bottleneck = edge->residual_capacity;
        }
    }

    if (bottleneck > (-1.0f * node->residual_capacity))
    {
        bottleneck = (-1.0f * node->residual_capacity);
    }

    // 2. augmenting

    // 2a - the source tree
    ((edge_t *)middle_edge->sister_edge)->residual_capacity += bottleneck;
    middle_edge->residual_capacity -= bottleneck;

    for (node = ((edge_t *)middle_edge->sister_edge)->head; ; node = edge->head) {

        edge = node->parent_edge;

        if (edge == &TERMINAL_EDGE) {
            // end of path
            break;
        }

        edge->residual_capacity += bottleneck;
        ((edge_t *)edge->sister_edge)->residual_capacity -= bottleneck;

        if (!edgeHasResidualCapacity(edge->sister_edge))
        {
            addOrphanAtFront(node);
        }
    }

	node->residual_capacity -= bottleneck;
    
    if (!nodeHasResidualCapacity(node))
    {
        addOrphanAtFront(node);
    }

    // 2b - the sink tree
    for (node = middle_edge->head; ; node = edge->head) 
    {
        edge = node->parent_edge;

        if (edge == &TERMINAL_EDGE)
        {
            // end of path
            break;
        }

	    ((edge_t *)edge->sister_edge)->residual_capacity += bottleneck;
	    edge->residual_capacity -= bottleneck;

        if (!edgeHasResidualCapacity(edge))
        {
            addOrphanAtFront(node);
        }
    }

    node->residual_capacity += bottleneck;

    if (!nodeHasResidualCapacity(node))
    {
        addOrphanAtFront(node);
    }

    total_flow += bottleneck;
}

/**
 * Adopt an orphan.
 */
static void processSourceOrphan(node_t *orphan) 
{
    edge_t *best_edge = NULL;
    edge_t *orphan_edge = NULL;

    uint32_t min_distance = UINT32_MAX;

    for (orphan_edge = orphan->first_outgoing_edge; orphan_edge != NULL; orphan_edge = orphan_edge->next)
    {
        if (edgeHasResidualCapacity(orphan_edge->sister_edge))
        {
            node_t *node = orphan_edge->head;
            edge_t *parent_edge = node->parent_edge;

            if (!node->is_in_sink && parent_edge != NULL) 
            {
                // check the origin of node
                uint32_t distance = 0;
                while (true)
                {
                    if (node->timestamp == time) 
                    {
                        distance += node->distance;
                        break;
                    }
                    
                    parent_edge = node->parent_edge;
                    distance++;

                    if (parent_edge == &TERMINAL_EDGE) {
                        node->timestamp = time;
                        node->distance = 1;
                        break;
                    }

                    if (parent_edge == &ORPHAN_EDGE) {
                        distance = UINT32_MAX;
                        break;
                    }

                    // otherwise, proceed to the next node
                    node = parent_edge->head;
                }

                if (distance < UINT32_MAX)
                { 
                // node originates from the source

                    if (distance < UINT32_MAX) 
                    {
                        best_edge = orphan_edge;
                        min_distance = distance;
                    }

                    // set marks along the path
                    for (node = orphan_edge->head; node->timestamp != time; node = ((edge_t *)node->parent_edge)->head) 
                    {
                        node->timestamp = time;
                        node->distance = distance;
                        distance--;
                    }
                }
            }
        }
    }

    orphan->parent_edge = best_edge;

    if (best_edge != NULL) 
    {
        orphan->timestamp = time;
        orphan->distance = min_distance + 1;
    } 
    else
    {
        // no parent found
        addToChangedList(orphan);

        // process neighbors
        for (orphan_edge = orphan->first_outgoing_edge; orphan_edge != NULL; orphan_edge = orphan_edge->next) {

            node_t *node = orphan_edge->head;
            edge_t *parent_edge = node->parent_edge;

            if (!node->is_in_sink && parent_edge != NULL) 
            {
                if (edgeHasResidualCapacity(orphan_edge->sister_edge))
                {
                    setNodeActive(node);
                }

                if (parent_edge != &TERMINAL_EDGE && parent_edge != &ORPHAN_EDGE && parent_edge->head == orphan)
                {
                    addOrphanAtBack(node);
                }
            }
        }
    }
}

/**
 * Adopt an orphan.
 */
static void processSinkOrphan(node_t *orphan) 
{
    edge_t *best_edge = NULL;
    edge_t *orphan_edge = NULL;

    uint32_t min_distance = UINT32_MAX;

    for (orphan_edge = orphan->first_outgoing_edge; orphan_edge != NULL; orphan_edge = orphan_edge->next)
    {
        if (edgeHasResidualCapacity(orphan_edge)) {

            node_t *node = orphan_edge->head;
            edge_t *parent_edge = node->parent_edge;

            if (node->is_in_sink && parent_edge != NULL) 
            {
                // check the origin of node
                uint32_t distance = 0;

                while (true) 
                {
                    if (node->timestamp == time) 
                    {
                        distance += node->distance;
                        break;
                    }

                    parent_edge = node->parent_edge;
                    distance++;

                    if (parent_edge == &TERMINAL_EDGE)
                    {
                        node->timestamp = time;
                        node->distance = 1;
                        break;
                    }

                    if (parent_edge == &ORPHAN_EDGE)
                    {
                        distance = UINT32_MAX;
                        break;
                    }

                    // otherwise, proceed to the next node
                    node = parent_edge->head;
                }

                if (distance < UINT32_MAX) 
                {
                    // node originates from the sink
                    if (distance < min_distance) 
                    {
                        best_edge = orphan_edge;
                        min_distance = distance;
                    }

                    // set marks along the path
                    for (node = orphan_edge->head; node->timestamp != time; node = ((edge_t *)node->parent_edge)->head) 
                    {
                        node->timestamp = time;
                        node->distance = distance;
                        distance--;
                    }
                }
            }
        }
    }

    orphan->parent_edge = best_edge;

    if (best_edge != NULL) 
    {
        orphan->timestamp = time;
        orphan->distance = min_distance + 1;
    } 
    else 
    {
        // no parent found
        addToChangedList(orphan);

        // process neighbors
    	for (orphan_edge = orphan->first_outgoing_edge; orphan_edge != NULL; orphan_edge = orphan_edge->next)
    	{
            node_t *node = orphan_edge->head;
            edge_t *parent_edge = node->parent_edge;

            if (node->is_in_sink && parent_edge != NULL) 
            {
                if (edgeHasResidualCapacity(orphan_edge))
                {
                    setNodeActive(node);
                }

                if (parent_edge != &TERMINAL_EDGE && parent_edge != &ORPHAN_EDGE && parent_edge->head == orphan)
                {
                    addOrphanAtBack(node);
                }
            }
        }
    }
}


/**
 * Performs the actual max-flow/min-cut computation.
 *
 * @param reuse_trees reuse trees of a previos call
 * @param changed_nodes list of nodes that potentially changed their
 *                      segmentation compared to a previous call, can be set
 *                      to NULL
 * @param number_changed_nodes number of nodes in the changed list
 */
float computeMaximumFlow(bool reuse_trees, uint32_t *changed_nodes, uint32_t *number_changed_nodes)
{
	if (maxflow_iteration == 0)
	{
		reuse_trees = false;
	}

	if (reuse_trees)
	{
		initializeMaxFlowAndReuseTrees();
	}
	else
	{
		initializeMaxFlow();
	}

	node_t *current_node = NULL;
	edge_t *edge = NULL;

	// main loop
	while (true) 
	{
		node_t *active_node = current_node;

		if (active_node != NULL) 
		{
			// remove active flag
			active_node->next = NULL;

			if (active_node->parent_edge == NULL)
			{
				active_node = NULL;
			}
		}
		if (active_node == NULL) 
		{
			if ((active_node = getNextActiveNode()) == NULL)
			{
				// no more active nodes - we're done here
				break;
			}
		}

		// groth
		if (!active_node->is_in_sink) 
		{
			// grow source tree
			for (edge = active_node->first_outgoing_edge; edge != NULL; edge = edge->next) 
			{
				if (edgeHasResidualCapacity(edge))
				{
					node_t *head_node = edge->head;

					if (head_node->parent_edge == NULL) 
					{
						// free node found, add to source tree
						head_node->is_in_sink = false;
						head_node->parent_edge = edge->sister_edge;
						head_node->timestamp = active_node->timestamp;

						head_node->distance = active_node->distance + 1;

						setNodeActive(head_node);
						addToChangedList(head_node);

					}
					else if (head_node->is_in_sink) 
					{
						// node is not free and belongs to other tree - path
						// via edge found
						break;
					} 
					else if (head_node->timestamp <= active_node->timestamp &&
					         head_node->distance > active_node->distance) 
					{
						// node is not free and belongs to our tree - try to
						// shorten its distance to the source
						head_node->parent_edge = edge->sister_edge;
						head_node->timestamp = active_node->timestamp;
						head_node->distance = active_node->distance + 1;
					}
				}
			}
		} 
		else
		{
			// active_node is in sink, grow sink tree
			for (edge = active_node->first_outgoing_edge; edge != NULL; edge = edge->next) 
			{
				if (edgeHasResidualCapacity(edge->sister_edge))
				{
					node_t *head_node = edge->head;

					if (head_node->parent_edge == NULL) 
					{
						// free node found, add to sink tree
						head_node->is_in_sink = true;
						head_node->parent_edge = edge->sister_edge; // n.setParent(edge.getsister_edge());
						head_node->timestamp = active_node->timestamp; // .setTimestamp(active_node.getTimestamp());
						head_node->distance = active_node->distance + 1; // .setDistance(active_node.getDistance() + 1);
						setNodeActive(head_node);
						addToChangedList(head_node);

					} 
					else if (!head_node->is_in_sink) 
					{
						// node is not free and belongs to other tree - path
						// via edge's sister_edge found
						edge = edge->sister_edge;
						break;

					} 
					else if (head_node->timestamp <= active_node->timestamp &&
					         head_node->distance > active_node->distance) 
					{
						// node is not free and belongs to our tree - try to
						// shorten its distance to the sink
						head_node->parent_edge = edge->sister_edge;
						head_node->timestamp = active_node->timestamp;
						head_node->distance = active_node->distance + 1;
					}
				}
			}
		}

		time++;

		if (edge != NULL) 
		{
			// we found a path via edge

			// set active flag
			active_node->next = active_node;
			current_node = active_node;

			// augmentation
			augment(edge);

			// adoption
			node_t *orphan_node;
			while ((orphan_node = pollOrphanQueue()) != NULL) 
			{
				if (orphan_node->is_in_sink)
				{
					processSinkOrphan(orphan_node);
				}
				else
				{
					processSourceOrphan(orphan_node);
				}
			}
		} 
		else 
		{
			// no path found
			current_node = NULL;
		}
	}

	maxflow_iteration++;

	// create list of changed nodes
	if (changed_nodes != NULL)
	{
		free(changed_nodes);

		*number_changed_nodes = 0;

		uint32_t i, changed_index;

		for (i = 0; i < number_of_nodes; i++)
		{
			(*number_changed_nodes)++;
		}

		changed_nodes = (uint32_t *)calloc((*number_changed_nodes), sizeof(uint32_t));

		for (i = 0, changed_index = 0; i < number_of_nodes; i++)
		{
			if (nodes[i].is_in_changed_list)
			{
				changed_nodes[changed_index++] = i;
			}
		}
	}

	return total_flow;
}

static bool edgeHasResidualCapacity(edge_t *edge)
{
    return (edge->residual_capacity > DELTA || edge->residual_capacity < -DELTA);
}

static bool nodeHasResidualCapacity(node_t *node)
{
    return (node->residual_capacity > DELTA || node->residual_capacity < -DELTA);
}
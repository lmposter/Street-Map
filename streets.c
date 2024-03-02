#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "streets.h"

// Structure definition for a map node.
struct node {
    int id;           // Unique identifier for the node.
    double lat;       // Latitude of the node.
    double lon;       // Longitude of the node.
    int num_ways;     // Number of ways this node is part of.
    int *way_ids;     // Dynamic array of way IDs this node is part of.
};

// Structure definition for a way (road).
struct way {
    int id;           // Unique identifier for the way.
    char *name;       // Name of the way, dynamically allocated.
    float maxspeed;   // Maximum speed limit on the way.
    bool oneway;      // Flag indicating if the way is one-way.
    int num_nodes;    // Number of nodes in this way.
    int *node_ids;    // Dynamic array of node IDs that form this way.
};

// Structure representing the entire map.
struct ssmap {
    struct node *nodes;  // Dynamic array of nodes.
    int nr_nodes;        // Total number of nodes.
    struct way *ways;    // Dynamic array of ways.
    int nr_ways;         // Total number of ways.
};

// Creates and initializes a map with specified numbers of nodes and ways.
// Parameters: nr_nodes - number of nodes, nr_ways - number of ways.
// Returns: Pointer to the newly created map, or NULL on failure.
struct ssmap * 
ssmap_create(int nr_nodes, int nr_ways)
{
    if (nr_nodes <= 0 || nr_ways <= 0) return NULL;

    struct ssmap *m = malloc(sizeof(struct ssmap));

    if (m == NULL) return NULL; // Out of memory check

    m->nodes = malloc(nr_nodes * sizeof(struct node));
    m->ways = malloc(nr_ways * sizeof(struct way));

    if (m->nodes == NULL || m->ways == NULL) {
        // Cleanup if any allocation fails
        free(m->nodes); // Safe even if NULL
        free(m->ways); // Safe even if NULL
        free(m);
        return NULL;
    }

    m->nr_nodes = nr_nodes;
    m->nr_ways = nr_ways;
    return m;
}

// Optional function that was not used but included in the starter file.
bool
ssmap_initialize(struct ssmap * m) 
{
    return true;
}

// Destroys a map and frees all allocated resources.
// Parameters: m - pointer to the map to destroy.
// No return value.
void
ssmap_destroy(struct ssmap * m)
{
    for (int i = 0; i < m->nr_nodes; i++) {
        free(m->nodes[i].way_ids);
    }
    free(m->nodes);

    for (int i = 0; i < m->nr_ways; i++) {
        free(m->ways[i].name);
        free(m->ways[i].node_ids);
    }
    free(m->ways);

    free(m);
}

// Adds a way to the map with specified attributes.
// Parameters: m - map to add to, id - way ID, name - name of the way, maxspeed - max speed limit,
//             oneway - one-way flag, num_nodes - number of nodes, node_ids - array of node IDs.
// Returns: Pointer to the added way, or NULL on failure.
struct way * 
ssmap_add_way(struct ssmap * m, int id, const char * name, float maxspeed, bool oneway, 
              int num_nodes, const int node_ids[num_nodes])
{
    if (id < 0 || id >= m->nr_ways) return NULL; // ID check

    struct way *w = &m->ways[id];
    w->id = id;
    w->name = strdup(name);
    w->maxspeed = maxspeed;
    w->oneway = oneway;
    w->num_nodes = num_nodes;
    w->node_ids = malloc(num_nodes * sizeof(int));

    if (w->name == NULL || w->node_ids == NULL) {
        // Cleanup on failure
        free(w->name);
        free(w->node_ids);
        return NULL;
    }

    memcpy(w->node_ids, node_ids, num_nodes * sizeof(int));
    return w;
}

// Adds a node to the map with specified attributes.
// Parameters: m - map to add to, id - node ID, lat - latitude, lon - longitude,
//             num_ways - number of ways, way_ids - array of way IDs.
// Returns: Pointer to the added node, or NULL on failure.
struct node * 
ssmap_add_node(struct ssmap * m, int id, double lat, double lon, 
               int num_ways, const int way_ids[num_ways])
{
    struct node *n = &m->nodes[id];
    n->id = id;
    n->lat = lat;
    n->lon = lon;
    n->num_ways = num_ways;
    n->way_ids = malloc(num_ways * sizeof(int));

    if (n->way_ids == NULL) {
        // Cleanup on failure
        free(n->way_ids);
        return NULL;
    }

    memcpy(n->way_ids, way_ids, num_ways * sizeof(int));
    return n;
}

// Prints information about a way with a given ID.
// Parameters: m - map to search in, id - ID of the way to print.
// No return value.
void
ssmap_print_way(const struct ssmap * m, int id)
{
    if (id < 0 || id >= m->nr_ways) {
        printf("error: way %d does not exist.\n", id);
        return;
    }

    struct way w = m->ways[id];
    printf("Way %d: %s\n", w.id, w.name);
}

// Prints information about a node with a given ID.
// Parameters: m - map to search in, id - ID of the node to print.
// No return value.
void
ssmap_print_node(const struct ssmap * m, int id)
{
    if (id < 0 || id >= m->nr_nodes) {
        printf("error: node %d does not exist.\n", id);
        return;
    }

    struct node n = m->nodes[id];
    printf("Node %d: (%.7lf, %.7lf)\n", n.id, n.lat, n.lon);
}

// Finds ways by name and prints their IDs.
// Parameters: m - map to search in, name - name to search for.
// No return value.
void 
ssmap_find_way_by_name(const struct ssmap * m, const char * name)
{
    for (int i = 0; i < m->nr_ways; i++) {
        if (strstr(m->ways[i].name, name) != NULL) {
            printf("%d ", m->ways[i].id);
        }
    }
        printf("\n");
}

// Finds nodes connected by ways with specified names.
// Parameters: m - map to search in, name1 - first name to search for, name2 - second name to search for.
// No return value.
void 
ssmap_find_node_by_names(const struct ssmap * m, const char * name1, const char * name2)
{
    for (int i = 0; i < m->nr_nodes; i++) {
        bool foundFirstMatch = false, foundSecondMatch = false;

        for (int j = 0; j < m->nodes[i].num_ways; j++) {
            struct way w = m->ways[m->nodes[i].way_ids[j]];

            // If name2 is NULL or same as name1, we're looking for two distinct ways matching name1
            if (name2 == NULL || strcmp(name1, name2) == 0){
                if (strstr(w.name, name1) != NULL) {
                    if (!foundFirstMatch) {
                        foundFirstMatch = true;
                    } else {
                        foundSecondMatch = true; // Found a second distinct way matching name1
                    }
                }
            } else {
                // For two different names, check each match separately
                if (!foundFirstMatch && strstr(w.name, name1) != NULL) {
                    foundFirstMatch = true; // Found a way that matches name1
                } else {
                    if (!foundSecondMatch && strstr(w.name, name2) != NULL) {
                        foundSecondMatch = true; // Found a way that matches name2
                    }
                }
            }       
        }
        if (foundFirstMatch && (name2 == NULL || foundSecondMatch)) {
            printf("%d ", m->nodes[i].id);
        }
    }
    printf("\n");
}

/**
 * Converts from degree to radian
 *
 * @param deg The angle in degrees.
 * @return the equivalent value in radian
 */
#define d2r(deg) ((deg) * M_PI/180.)

/**
 * Calculates the distance between two nodes using the Haversine formula.
 *
 * @param x The first node.
 * @param y the second node.
 * @return the distance between two nodes, in kilometre.
 */
static double
distance_between_nodes(const struct node * x, const struct node * y) {
    double R = 6371.;       
    double lat1 = x->lat;
    double lon1 = x->lon;
    double lat2 = y->lat;
    double lon2 = y->lon;
    double dlat = d2r(lat2-lat1); 
    double dlon = d2r(lon2-lon1); 
    double a = pow(sin(dlat/2), 2) + cos(d2r(lat1)) * cos(d2r(lat2)) * pow(sin(dlon/2), 2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a)); 
    return R * c; 
}

// Finds a way connecting two nodes.
// Parameters: m - map to search in, node_a_id - ID of the first node, node_b_id - ID of the second node.
// Returns: Pointer to the way connecting the nodes, or NULL if not found.
static struct way * find_way_between_nodes(const struct ssmap * m, int node_a_id, int node_b_id) {
    for (int i = 0; i < m->nr_ways; ++i) {
        struct way *w = &m->ways[i];
        // Check if both nodes exist in this way
        bool node_a_found = false, node_b_found = false;
        for (int j = 0; j < w->num_nodes; ++j) {
            if (w->node_ids[j] == node_a_id) node_a_found = true;
            if (w->node_ids[j] == node_b_id) node_b_found = true;
            if(node_a_found && node_b_found) return w;
        }
    }
    return NULL; // No way found
}

// Calculates the total travel time for a path through specified nodes.
// Parameters: m - map containing the nodes, size - number of nodes in the path, node_ids - array of node IDs.
// Returns: Total travel time in minutes, or -1.0 on error.
double 
ssmap_path_travel_time(const struct ssmap * m, int size, int node_ids[size])
{
    double total_time = 0.0; // Total travel time in minutes
    struct node *nodes = m->nodes; // Shortcut to the nodes array

    // Used to check if a node appears more than once
    bool node_visited[m->nr_nodes];
    memset(node_visited, 0, sizeof(node_visited));

    for (int i = 0; i < size - 1; i++) {
        int node_id = node_ids[i];
        int next_node_id = node_ids[i + 1];

        // Check if nodes exist
        if (node_id < 0 || node_id >= m->nr_nodes || next_node_id < 0 || next_node_id >= m->nr_nodes) {
            printf("error: node %d does not exist.\n", node_id < 0 || node_id >= m->nr_nodes ? node_id : next_node_id);
            return -1.0;
        }

        // Check for repeated nodes
        if (node_visited[node_id]) {
            printf("error: node %d appeared more than once.\n", node_id);
            return -1.0;
        }
        node_visited[node_id] = true;

        struct node current_node = nodes[node_id];
        struct node next_node = nodes[next_node_id];
        struct way *connecting_way = find_way_between_nodes(m, node_id, next_node_id);

        // Check if there is a connecting way
        if (connecting_way == NULL) {
            printf("error: there are no roads between node %d and node %d.\n", node_id, next_node_id);
            return -1.0;
        }

        // Check if the nodes are connected in the correct order
        bool adjacency_found = false;
        for (int j = 0; j < connecting_way->num_nodes - 1; ++j) {
            if (connecting_way->node_ids[j] == node_id && connecting_way->node_ids[j + 1] == next_node_id) {
                adjacency_found = true;
                break;
            }
            if (connecting_way->node_ids[j] == next_node_id && connecting_way->node_ids[j + 1] == node_id) {
                adjacency_found = true;
                break;
            }
        }
        if (!adjacency_found) {
            printf("error: cannot go directly from node %d to node %d.\n", node_id, next_node_id);
            return -1.0;
        }

        // Check for one-way violation
        if (connecting_way->oneway) {
            bool oneway_violation = true;
            for (int j = 0; j < connecting_way->num_nodes - 1; j++) {
                if (connecting_way->node_ids[j] == node_id && connecting_way->node_ids[j + 1] == next_node_id) {
                    oneway_violation = false;
                    break;
                }
            }
            if (oneway_violation) {
                printf("error: cannot go in reverse from node %d to node %d.\n", node_id, next_node_id);
                return -1.0;
            }
        }

        if (size > 0 && node_visited[node_ids[size - 1]]) {
            printf("error: node %d appeared more than once.\n", node_ids[size - 1]);
            return -1.0;
        }

        double segment_distance = distance_between_nodes(&current_node, &next_node); // Distance in kilometers
        double segment_time = segment_distance / (connecting_way->maxspeed / 60.0); // Time = Distance / Speed, converted from hours to minutes
        total_time += segment_time;
    }

    return total_time; // Return the total travel time
}

// Priority Queue Function Implementations
#include <float.h> // For DBL_MAX

// Define a structure for the priority queue nodes
typedef struct {
    int node_id;
    double priority; // Lower values have higher priority
} PriorityQueueNode;

// Define a priority queue (min-heap) structure
typedef struct {
    PriorityQueueNode *nodes;
    int capacity;
    int size;
} PriorityQueue;

PriorityQueue *priority_queue_create(int capacity) {
    PriorityQueue *pq = (PriorityQueue *)malloc(sizeof(PriorityQueue));
    if (!pq) {
        printf("error: out of memory.\n");
        return NULL;
    }
    pq->nodes = (PriorityQueueNode *)malloc(sizeof(PriorityQueueNode) * capacity);
    if (!pq->nodes) {
        printf("error: out of memory.\n");
        free(pq);
        return NULL;
    }
    pq->capacity = capacity;
    pq->size = 0;
    return pq;
}

void priority_queue_destroy(PriorityQueue *pq) {
    free(pq->nodes);
    free(pq);
}

void swap_nodes(PriorityQueueNode *a, PriorityQueueNode *b) {
    PriorityQueueNode temp = *a;
    *a = *b;
    *b = temp;
}

void priority_queue_push(PriorityQueue *pq, int node_id, double priority) {
    if (pq->size >= pq->capacity) {
        printf("error: priority queue full.\n");
        return;
    }
    pq->nodes[pq->size] = (PriorityQueueNode){node_id, priority};
    int index = pq->size;
    pq->size++;

    while (index != 0) {
        int parent_index = (index - 1) / 2;
        if (pq->nodes[parent_index].priority <= pq->nodes[index].priority) {
            break;
        }
        swap_nodes(&pq->nodes[parent_index], &pq->nodes[index]);
        index = parent_index;
    }
}

int priority_queue_pop(PriorityQueue *pq) {
    if (pq->size == 0) {
        printf("error: priority queue empty.\n");
        return -1;
    }
    int pop_id = pq->nodes[0].node_id;
    pq->size--;
    pq->nodes[0] = pq->nodes[pq->size];
    int index = 0;

    while (true) {
        int left_child = 2 * index + 1;
        int right_child = 2 * index + 2;
        int smallest = index;

        if (left_child < pq->size && pq->nodes[left_child].priority < pq->nodes[smallest].priority) {
            smallest = left_child;
        }
        if (right_child < pq->size && pq->nodes[right_child].priority < pq->nodes[smallest].priority) {
            smallest = right_child;
        }
        if (smallest == index) {
            break;
        }
        swap_nodes(&pq->nodes[index], &pq->nodes[smallest]);
        index = smallest;
    }

    return pop_id;
}

bool priority_queue_is_empty(const PriorityQueue *pq) {
    return pq->size == 0;
}

// Standard Dijkstra's Algorithm from https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm in the assignment instruction.
void ssmap_path_create(const struct ssmap * m, int start_id, int end_id) {
    // Check for valid node IDs
    if (start_id < 0 || start_id >= m->nr_nodes || end_id < 0 || end_id >= m->nr_nodes) {
        printf("error: node %d does not exist.\n", start_id < 0 || start_id >= m->nr_nodes ? start_id : end_id);
        return;
    }

    // Initialize data structures for Dijkstra's algorithm
    double distances[m->nr_nodes];
    int previous[m->nr_nodes];
    bool visited[m->nr_nodes];

    for (int i = 0; i < m->nr_nodes; i++) {
        distances[i] = DBL_MAX;  // Set all distances to "infinity"
        previous[i] = -1;        // No predecessors initially
        visited[i] = false;      // All nodes are initially unvisited
    }

    distances[start_id] = 0.0;   // Distance to the start node is zero

    // Priority queue for the nodes to visit
    PriorityQueue *pq = priority_queue_create(m->nr_nodes);
    if (!pq) {
        // If we run out of memory when creating the priority queue
        printf("error: out of memory.\n");
        return;
    }

    priority_queue_push(pq, start_id, 0.0);  // Start with the start_id

    while (!priority_queue_is_empty(pq)) {
        // Pop the node with the smallest distance
        int current_id = priority_queue_pop(pq);
        if (visited[current_id]) continue; // Skip nodes that have been visited

        visited[current_id] = true; // Mark the node as visited

        if (current_id == end_id) break; // We've reached our destination

        struct node current_node = m->nodes[current_id];

        // Iterate over all adjacent nodes
        for (int i = 0; i < current_node.num_ways; i++) {
            struct way *w = &m->ways[current_node.way_ids[i]];
            for (int j = 0; j < w->num_nodes - 1; j++) {
                // Check both directions if the way is not one-way
                if (w->node_ids[j] == current_id || (!w->oneway && w->node_ids[j + 1] == current_id)) {
                    int neighbour_id = w->node_ids[j] == current_id ? w->node_ids[j + 1] : w->node_ids[j];
                    if (visited[neighbour_id]) continue; // Skip visited neighbors

                    // Calculate the time to reach the neighbour via current node
                    double travel_time = distance_between_nodes(&current_node, &m->nodes[neighbour_id]) / w->maxspeed;
                    double new_dist = distances[current_id] + travel_time;

                    // If the new distance is smaller, update it
                    if (new_dist < distances[neighbour_id]) {
                        distances[neighbour_id] = new_dist;
                        previous[neighbour_id] = current_id;
                        priority_queue_push(pq, neighbour_id, new_dist);
                    }
                }
            }
        }
    }

    if (distances[end_id] == DBL_MAX) {
        printf("error: could not find a path from node %d to node %d.\n", start_id, end_id);
        priority_queue_destroy(pq);
        return;
    }

    // Reconstruct the path
    int path[m->nr_nodes];
    int path_size = 0;
    for (int at = end_id; at != -1; at = previous[at]) {
        path[path_size++] = at;
    }

    // Print the path in the correct order
    for (int i = path_size - 1; i >= 0; i--) {
        if (i < path_size - 1) printf(" ");
        printf("%d", path[i]);
    }
    printf("\n");

    priority_queue_destroy(pq);
}
/*
 * node_arrays.c
 * Global per-field arrays for the structure-peeled node layout.
 * Definitions live here; extern declarations are in defines.h.
 */
#include "defines.h"

cost_t  *node_potentials    = NULL;
int     *node_orientations  = NULL;
node_p  *node_childs        = NULL;
node_p  *node_preds         = NULL;
node_p  *node_siblings      = NULL;
node_p  *node_sibling_prevs = NULL;
arc_p   *node_basic_arcs    = NULL;
arc_p   *node_firstouts     = NULL;
arc_p   *node_firstins      = NULL;
flow_t  *node_flows         = NULL;
LONG    *node_depths        = NULL;
int     *node_numbers       = NULL;
int     *node_times         = NULL;

LONG     g_num_nodes        = 0;

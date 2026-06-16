/**************************************************************************
PSTART.C of ZIB optimizer MCF, SPEC version
Modified: Structure Peeling -- node_p is int32_t index.
          root node is net->n (last index), non-root nodes are 1..n-1.
**************************************************************************/

#include "pstart.h"


#ifdef _PROTO_
LONG primal_start_artificial( network_t *net )
#else
LONG primal_start_artificial( net )
    network_t *net;
#endif
{
    node_p root;   /* index of root node */
    node_p node;   /* loop index         */
    arc_t *arc;
    int i;

    /* Root is at index net->n  (same layout as original &net->nodes[net->n]) */
    root = (node_p)net->n;

    NODE_BASIC_ARC(root)   = NULL;
    NODE_PRED(root)        = INVALID_NODE;
    NODE_CHILD(root)       = (node_p)1;   /* first non-root node */
    NODE_SIBLING(root)     = INVALID_NODE;
    NODE_SIBLING_PREV(root)= INVALID_NODE;
    NODE_DEPTH(root)       = (LONG)(net->n) + 1;
    NODE_ORIENTATION(root) = 0;
    NODE_POTENTIAL(root)   = (cost_t)-MAX_ART_COST;
    NODE_FLOW(root)        = ZERO;

    arc = net->arcs;
#if (defined(_OPENMP) || defined(SPEC_OPENMP)) && !defined(SPEC_SUPPRESS_OPENMP) && !defined(SPEC_AUTO_SUPPRESS_OPENMP)
#pragma omp parallel for
#endif
    for( i = 0; i < net->m; i++ )
        if( arc[i].ident != FIXED )
            arc[i].ident = AT_LOWER;

    arc = net->dummy_arcs;
#if (defined(_OPENMP) || defined(SPEC_OPENMP)) && !defined(SPEC_SUPPRESS_OPENMP) && !defined(SPEC_AUTO_SUPPRESS_OPENMP)
#pragma omp parallel for
#endif
    for( i = 0; i < net->n; i++ )
    {
        node = (node_p)i;   /* node index == i */

        NODE_BASIC_ARC(node)   = &arc[i];
        NODE_PRED(node)        = root;
        NODE_CHILD(node)       = INVALID_NODE;
        /* sibling linked list: node i -> node i+1 */
        NODE_SIBLING(node)     = (node_p)(i + 1);
        NODE_SIBLING_PREV(node)= (node_p)(i - 1);
        NODE_DEPTH(node)       = 1;
        NODE_ORIENTATION(node) = UP;
        NODE_POTENTIAL(node)   = ZERO;
        NODE_FLOW(node)        = (flow_t)0;

        arc[i].cost  = (cost_t)MAX_ART_COST;
        arc[i].ident = BASIC;
        arc[i].tail  = node;
        arc[i].head  = root;
        arc[i].id    = DUMMY_ARC;
    }

    /* Fix up endpoints of the sibling list */
    NODE_SIBLING((node_p)(net->n - 1)) = INVALID_NODE;
    NODE_SIBLING_PREV((node_p)0)       = INVALID_NODE;

    return 0;
}

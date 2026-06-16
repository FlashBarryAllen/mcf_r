/**************************************************************************
READMIN.C of ZIB optimizer MCF, SPEC version
Modified: Structure Peeling -- allocate per-field node arrays instead of
          a monolithic node_t array. node_p is now int32_t (index).
**************************************************************************/

#include "readmin.h"

/* Helper: allocate all per-field node arrays for (count) nodes */
static int alloc_node_arrays(LONG count)
{
    node_potentials    = (cost_t *) calloc(count, sizeof(cost_t));
    node_orientations  = (int *)    calloc(count, sizeof(int));
    node_childs        = (node_p *) calloc(count, sizeof(node_p));
    node_preds         = (node_p *) calloc(count, sizeof(node_p));
    node_siblings      = (node_p *) calloc(count, sizeof(node_p));
    node_sibling_prevs = (node_p *) calloc(count, sizeof(node_p));
    node_basic_arcs    = (arc_p *)  calloc(count, sizeof(arc_p));
    node_firstouts     = (arc_p *)  calloc(count, sizeof(arc_p));
    node_firstins      = (arc_p *)  calloc(count, sizeof(arc_p));
    node_flows         = (flow_t *) calloc(count, sizeof(flow_t));
    node_depths        = (LONG *)   calloc(count, sizeof(LONG));
    node_numbers       = (int *)    calloc(count, sizeof(int));
    node_times         = (int *)    calloc(count, sizeof(int));

    if (!node_potentials || !node_orientations || !node_childs ||
        !node_preds      || !node_siblings     || !node_sibling_prevs ||
        !node_basic_arcs || !node_firstouts    || !node_firstins ||
        !node_flows      || !node_depths       || !node_numbers  ||
        !node_times)
        return 0;

    /* Initialise all node_p fields to INVALID_NODE */
    LONG i;
    for (i = 0; i < count; i++) {
        node_childs[i]        = INVALID_NODE;
        node_preds[i]         = INVALID_NODE;
        node_siblings[i]      = INVALID_NODE;
        node_sibling_prevs[i] = INVALID_NODE;
    }
    return 1;
}

#ifdef _PROTO_
LONG read_min( network_t *net )
#else
LONG read_min( net )
     network_t *net;
#endif
{                                       
    FILE *in = NULL;
    char instring[201];
    LONG t, h, c;
    LONG i, actArc = 0;
    arc_t *arc;
    /* node indices */
    node_p node_root;   /* index 0 */

    if(( in = fopen( net->inputfile, "r")) == NULL )
        return -1;

    fgets( instring, 200, in );
#ifdef SPEC
    if( sscanf( instring, "%" PRId64 " %" PRId64 , &t, &h ) != 2 )
#else
    if( sscanf( instring, "%ld %ld", &t, &h ) != 2 )
#endif
        return -1;

    net->n_trips = t;
    net->m_org   = h;
    net->n       = (t + t + 1);
    net->m       = (t + t + t + h);

    net->max_elems = K;
    net->nr_group  = ((net->m - 1) / K) + 1;
    if (net->m % K != 0)
        net->full_groups = net->nr_group - (K - (net->m % K));
    else
        net->full_groups = net->nr_group;
    while (net->full_groups < 0) {
        net->full_groups = net->nr_group + net->full_groups;
        net->max_elems--;
    }

    if( net->n_trips <= MAX_NB_TRIPS_FOR_SMALL_NET )
    {
        net->max_m              = net->m;
        net->max_new_m          = MAX_NEW_ARCS_SMALL_NET;
        net->max_residual_new_m = net->max_m - net->m;
    }
    else
    {
        net->max_m     = MAX(net->m + MAX_NEW_ARCS, STRECHT(STRECHT(net->m)));
        net->max_new_m = MAX_NEW_ARCS_LARGE_NET;
    }

    assert( net->max_new_m >= 3 );

    /* --- Allocate per-field node arrays (net->n + 1 entries, 0-based) --- */
    g_num_nodes = net->n + 1;
    if (!alloc_node_arrays(g_num_nodes))
    {
        printf("read_min(): not enough memory for node arrays\n");
        getfree(net);
        return -1;
    }

    /* net->nodes is now the index of the first node (0).
     * net->stop_nodes is the one-past-end index. */
    net->nodes      = (node_p)0;
    net->stop_nodes = (node_p)(net->n + 1);

    /* --- Allocate arc pools (unchanged) --- */
    net->dummy_arcs  = (arc_t *) calloc(net->n,       sizeof(arc_t));
    net->sorted_arcs = (arc_t *) calloc(net->max_m,   sizeof(arc_t));
    net->arcs        = (arc_t *) calloc(net->max_m,   sizeof(arc_t));

    if( !(net->arcs && net->dummy_arcs && net->sorted_arcs) )
    {
        printf("read_min(): not enough memory for arc arrays\n");
        getfree(net);
        return -1;
    }

    net->stop_arcs  = net->arcs + net->m;
    net->stop_dummy = net->dummy_arcs + net->n;

    /* Root node is index net->n (last slot, same as original &node[net->n]) */
    node_root = (node_p)net->n;

    arc = net->arcs;

    for( i = 1; i <= net->n_trips; i++ )
    {
        fgets( instring, 200, in );
#ifdef SPEC
        if( sscanf( instring, "%" PRId64 " %" PRId64 , &t, &h ) != 2 || t > h )
#else
        if( sscanf( instring, "%ld %ld", &t, &h ) != 2 || t > h )
#endif
            return -1;

        /* node index i */
        NODE_NUMBER(i)              = (int)-i;
        NODE_FLOW(i)                = (flow_t)-1;
        NODE_TIME(i)                = (int)t;

        /* node index i + net->n_trips */
        NODE_NUMBER(i + net->n_trips) = (int)i;
        NODE_FLOW(i + net->n_trips)   = (flow_t)1;
        NODE_TIME(i + net->n_trips)   = (int)h;

        /* Arc: root -> node[i] */
        arc->id       = (int)actArc;
        arc->tail     = node_root;
        arc->head     = (node_p)i;
        arc->org_cost = arc->cost = (cost_t)(net->bigM + 15);
        arc->nextout  = NODE_FIRSTOUT(node_root);
        NODE_FIRSTOUT(node_root) = arc;
        arc->nextin   = NODE_FIRSTIN(i);
        NODE_FIRSTIN(i) = arc;
        arc = net->arcs + getArcPosition(net, ++actArc);

        /* Arc: node[i+n_trips] -> root */
        arc->id       = (int)actArc;
        arc->tail     = (node_p)(i + net->n_trips);
        arc->head     = node_root;
        arc->org_cost = arc->cost = (cost_t)15;
        arc->nextout  = NODE_FIRSTOUT(i + net->n_trips);
        NODE_FIRSTOUT(i + net->n_trips) = arc;
        arc->nextin   = NODE_FIRSTIN(node_root);
        NODE_FIRSTIN(node_root) = arc;
        arc = net->arcs + getArcPosition(net, ++actArc);

        /* Arc: node[i] -> node[i+n_trips] */
        arc->id       = (int)actArc;
        arc->tail     = (node_p)i;
        arc->head     = (node_p)(i + net->n_trips);
        arc->org_cost = arc->cost = (cost_t)(2 * MAX(net->bigM, (LONG)BIGM));
        arc->nextout  = NODE_FIRSTOUT(i);
        NODE_FIRSTOUT(i) = arc;
        arc->nextin   = NODE_FIRSTIN(i + net->n_trips);
        NODE_FIRSTIN(i + net->n_trips) = arc;
        arc = net->arcs + getArcPosition(net, ++actArc);
    }

    if( i != net->n_trips + 1 )
        return -1;

    for( i = 0; i < net->m_org; i++, arc = net->arcs + getArcPosition(net, ++actArc))
    {
        fgets( instring, 200, in );
#ifdef SPEC
        if( sscanf( instring, "%" PRId64 " %" PRId64 " %" PRId64 , &t, &h, &c ) != 3 )
#else
        if( sscanf( instring, "%ld %ld %ld", &t, &h, &c ) != 3 )
#endif
            return -1;

        arc->id       = (int)actArc;
        arc->tail     = (node_p)(t + net->n_trips);
        arc->head     = (node_p)h;
        arc->org_cost = arc->cost = (cost_t)c;
        arc->nextout  = NODE_FIRSTOUT(t + net->n_trips);
        NODE_FIRSTOUT(t + net->n_trips) = arc;
        arc->nextin   = NODE_FIRSTIN(h);
        NODE_FIRSTIN(h) = arc;
    }
    arc = net->stop_arcs;

    if( net->stop_arcs != arc )
    {
        net->stop_arcs = arc;
        arc = net->arcs;
        for( net->m = 0; arc < net->stop_arcs; arc++ )
            (net->m)++;
        net->m_org = net->m;
    }

    fclose( in );

    net->clustfile[0] = (char)0;
    for( i = 1; i <= net->n_trips; i++ )
    {
        arc = net->arcs + getArcPosition(net, 3 * i - 1);
        arc->cost     = (cost_t)((-2) * MAX(net->bigM, (LONG)BIGM));
        arc->org_cost = (cost_t)((-2) * (MAX(net->bigM, (LONG)BIGM)));
    }

    return 0;
}

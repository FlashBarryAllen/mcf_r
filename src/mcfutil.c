/**************************************************************************
MCFUTIL.C of ZIB optimizer MCF, SPEC version
Modified: Structure Peeling -- getfree frees per-field node arrays.
          refresh_neighbour_lists uses NODE_FIRSTIN/FIRSTOUT macros.
          flow_cost / dual_feasible / primal_feasible use node_p index.
**************************************************************************/

#include "mcfutil.h"

#ifdef _PROTO_
void refresh_neighbour_lists( network_t *net, LONG (*getPos)(network_t*, LONG) )
#else
void refresh_neighbour_lists( net )
    network_t *net;
#endif
{
    node_p node;
    arc_t *arc;
    int i;

    for( node = net->nodes; node < net->stop_nodes; node++ )
    {
        NODE_FIRSTIN(node)  = (arc_p)NULL;
        NODE_FIRSTOUT(node) = (arc_p)NULL;
    }

    arc = net->arcs;
    for( i = 0; i < net->m; i++, arc = &net->arcs[getPos(net, i)] )
    {
        arc->nextout            = NODE_FIRSTOUT(arc->tail);
        NODE_FIRSTOUT(arc->tail) = arc;
        arc->nextin             = NODE_FIRSTIN(arc->head);
        NODE_FIRSTIN(arc->head)  = arc;
    }

    return;
}


#ifdef _PROTO_
double flow_cost( network_t *net )
#else
double flow_cost( net )
    network_t *net;
#endif
{
    arc_t *arc;
    node_p node;

    LONG fleet = 0;
    int i;
    cost_t operational_cost = 0;

    arc = net->arcs;
#if (defined(_OPENMP) || defined(SPEC_OPENMP)) && !defined(SPEC_SUPPRESS_OPENMP) && !defined(SPEC_AUTO_SUPPRESS_OPENMP)
#pragma omp parallel for
#endif
    for( i = 0; i < net->m; i++ )
    {
        if( arc[i].ident == AT_UPPER )
            arc[i].flow = (flow_t)1;
        else
            arc[i].flow = (flow_t)0;
    }

#if (defined(_OPENMP) || defined(SPEC_OPENMP)) && !defined(SPEC_SUPPRESS_OPENMP) && !defined(SPEC_AUTO_SUPPRESS_OPENMP)
#pragma omp parallel for
#endif
    /* Root node is at index net->n and has no basic_arc; iterate non-root nodes only */
    for( i = 1; i < net->n; i++ )
        NODE_BASIC_ARC(i)->flow = NODE_FLOW(i);

    arc = net->arcs;
#if (defined(_OPENMP) || defined(SPEC_OPENMP)) && !defined(SPEC_SUPPRESS_OPENMP) && !defined(SPEC_AUTO_SUPPRESS_OPENMP)
#pragma omp parallel for reduction(+ : fleet, operational_cost)
#endif
    for( i = 0; i < net->m; i++ )
    {
        if( arc[i].flow )
        {
            node_p t = arc[i].tail;
            node_p h = arc[i].head;
            if( !(NODE_NUMBER(t) < 0 && NODE_NUMBER(h) > 0) )
            {
                if( !NODE_NUMBER(t) )
                {
                    operational_cost += (arc[i].cost - net->bigM);
                    fleet++;
                }
                else
                    operational_cost += arc[i].cost;
            }
        }
    }

    return (double)fleet * (double)net->bigM + (double)operational_cost;
}

static LONG old_group = 0;
static LONG old_Arc   = 0;

#ifdef _PROTO_
LONG start()
#else
start()
#endif
{
    old_group = 0;
    old_Arc   = 0;
    return 0;
}

#ifdef _PROTO_
LONG getArcPosition(network_t *net, LONG actArc)
#else
int getArcPosition(net, actArc)
network_t *net;
LONG actArc;
#endif
{
    LONG result, akt_group;
    akt_group = actArc % net->nr_group;
    if (akt_group > net->full_groups) {
        result = (actArc / net->nr_group) + (net->full_groups * net->max_elems + (akt_group - net->full_groups) * (net->max_elems - 1));
    } else {
        result = (actArc / net->nr_group) + (akt_group * net->max_elems);
    }
    return (int)result;
}

#ifdef _PROTO_
LONG getOriginalArcPosition(network_t *net, LONG actArc)
#else
getOriginalArcPosition(net, actArc)
network_t *net;
LONG actArc;
#endif
{
    return actArc;
}


#ifdef _PROTO_
double flow_org_cost( network_t *net )
#else
double flow_org_cost( net )
    network_t *net;
#endif
{
    arc_t *arc;
    int i;

    LONG fleet = 0;
    cost_t operational_cost = 0;

    arc = net->arcs;
#if (defined(_OPENMP) || defined(SPEC_OPENMP)) && !defined(SPEC_SUPPRESS_OPENMP) && !defined(SPEC_AUTO_SUPPRESS_OPENMP)
#pragma omp parallel for
#endif
    for( i = 0; i < net->m; i++ )
    {
        if( arc[i].ident == AT_UPPER )
            arc[i].flow = (flow_t)1;
        else
            arc[i].flow = (flow_t)0;
    }

#if (defined(_OPENMP) || defined(SPEC_OPENMP)) && !defined(SPEC_SUPPRESS_OPENMP) && !defined(SPEC_AUTO_SUPPRESS_OPENMP)
#pragma omp parallel for reduction(+: fleet, operational_cost)
#endif
    for( i = 0; i < net->n; i++ )
        NODE_BASIC_ARC(i)->flow = NODE_FLOW(i);

    arc = net->arcs;
    for( i = 0; i < net->m; i++ )
    {
        if( arc[i].flow )
        {
            node_p t = arc[i].tail;
            node_p h = arc[i].head;
            if( !(NODE_NUMBER(t) < 0 && NODE_NUMBER(h) > 0) )
            {
                if( !NODE_NUMBER(t) )
                {
                    operational_cost += (arc[i].org_cost - net->bigM);
                    fleet++;
                }
                else
                    operational_cost += arc[i].org_cost;
            }
        }
    }

    printf("ORG_COST: %f", (double)fleet * (double)net->bigM + (double)operational_cost);
    return (double)fleet * (double)net->bigM + (double)operational_cost;
}


#ifdef _PROTO_
LONG primal_feasible( network_t *net )
#else
LONG primal_feasible( net )
    network_t *net;
#endif
{
    node_p node;
    arc_t *dummy     = net->dummy_arcs;
    arc_t *stop_dummy = net->stop_dummy;
    arc_t *arc;
    flow_t flow;

    for( node = net->nodes + 1; node < net->stop_nodes; node++ )
    {
        arc  = NODE_BASIC_ARC(node);
        flow = NODE_FLOW(node);
        if( arc >= dummy && arc < stop_dummy )
        {
            if( ABS(flow) > (flow_t)net->feas_tol )
            {
                printf("PRIMAL NETWORK SIMPLEX: ");
                printf("artificial arc with nonzero flow, node %d (%" PRId64 ")\n",
                       NODE_NUMBER(node), flow);
            }
        }
        else
        {
            if( flow < (flow_t)(-net->feas_tol)
               || flow - (flow_t)1 > (flow_t)net->feas_tol )
            {
                printf("PRIMAL NETWORK SIMPLEX: ");
                printf("basis primal infeasible (%" PRId64 ")\n", flow);
                net->feasible = 0;
                return 1;
            }
        }
    }

    net->feasible = 1;
    return 0;
}


#ifdef _PROTO_
LONG dual_feasible( network_t *net )
#else
LONG dual_feasible( net )
    network_t *net;
#endif
{
    arc_t  *arc;
    arc_t  *stop    = net->stop_arcs;
    cost_t  red_cost;
    LONG i = 0;

    for( i = 0, arc = net->arcs; arc < stop; arc++, i++ )
    {
        red_cost = arc->cost - NODE_POTENTIAL(arc->tail)
                             + NODE_POTENTIAL(arc->head);
        switch( arc->ident )
        {
        case BASIC:
        case AT_LOWER:
            if( red_cost < (cost_t)-net->feas_tol )
                goto DUAL_INFEAS;
            break;
        case AT_UPPER:
            if( red_cost > (cost_t)net->feas_tol )
                goto DUAL_INFEAS;
            break;
        case FIXED:
        default:
            break;
        }
    }
    return 0;

DUAL_INFEAS:
    fprintf(stderr, "DUAL NETWORK SIMPLEX: ");
    fprintf(stderr, "basis dual infeasible\n");
    return 1;
}


#ifdef _PROTO_
LONG getfree( network_t *net )
#else
LONG getfree( net )
     network_t *net;
#endif
{
    /* Free per-field node arrays */
    FREE(node_potentials);
    FREE(node_orientations);
    FREE(node_childs);
    FREE(node_preds);
    FREE(node_siblings);
    FREE(node_sibling_prevs);
    FREE(node_basic_arcs);
    FREE(node_firstouts);
    FREE(node_firstins);
    FREE(node_flows);
    FREE(node_depths);
    FREE(node_numbers);
    FREE(node_times);
    node_potentials = node_flows = NULL;
    node_depths = NULL;
    node_childs = node_preds = node_siblings = node_sibling_prevs = NULL;
    node_basic_arcs = node_firstouts = node_firstins = NULL;
    node_orientations = node_numbers = node_times = NULL;
    g_num_nodes = 0;

    /* Free arc arrays */
    FREE(net->arcs);
    FREE(net->dummy_arcs);
    FREE(net->sorted_arcs);
    net->nodes = net->stop_nodes = INVALID_NODE;
    net->arcs = net->stop_arcs = NULL;
    net->dummy_arcs = net->stop_dummy = NULL;
    net->sorted_arcs = NULL;

    return 0;
}

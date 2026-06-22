/**************************************************************************
PBEAMPP.C of ZIB optimizer MCF, SPEC version
Optimisations:
  1. cost_compare uses basket->arc_id instead of (*b)->a->id, eliminating
     a pointer dereference on every comparison (a real, measured win).
  2. (Reverted) A hand-rolled partial-sort replacement for spec_qsort was
     tried and measured SLOWER on real SPEC input -- see note below.
**************************************************************************/

#include "pbeampp.h"
#include "mcfutil.h"
#if defined(SPEC)
# include "spec_qsort.h"
#endif

static arc_t* full_group_end_arc;

#ifdef _PROTO_
void set_static_vars(network_t *net, arc_t* arcs)
#else
void set_static_vars(net, arcs)
network_t *net;
arc_t* arcs;
#endif
{
  full_group_end_arc = arcs + net->full_groups * net->max_elems;
}


#ifdef _PROTO_
int bea_is_dual_infeasible( arc_t *arc, cost_t red_cost )
#else
int bea_is_dual_infeasible( arc, red_cost )
    arc_t *arc;
    cost_t red_cost;
#endif
{
    return(    (red_cost < 0 && arc->ident == AT_LOWER)
            || (red_cost > 0 && arc->ident == AT_UPPER) );
}


/* cost_compare: primary key = abs_cost (descending), tie-break = arc_id (ascending).
 * arc_id is now stored directly in BASKET -- no pointer dereference needed. */
#ifdef _PROTO_
int cost_compare( BASKET **b1, BASKET **b2 )
#else
int cost_compare( b1, b2 )
    BASKET **b1;
    BASKET **b2;
#endif
{
  if( (*b1)->abs_cost < (*b2)->abs_cost ) return  1;
  if( (*b1)->abs_cost > (*b2)->abs_cost ) return -1;
  if( (*b1)->arc_id   > (*b2)->arc_id   ) return  1;
  return -1;
}


/* NOTE: An earlier version of this file replaced the qsort below with a
 * hand-rolled partial-sort (max-heap of size B) on the theory that we only
 * need the top-B out of up to K candidates. In practice, READY is reached
 * as soon as the running total across threads hits B, so basket_sizes[thread]
 * is almost always only slightly larger than B itself -- the heap's "large n"
 * branch essentially never executes. The wrapper call and extra branch were
 * pure overhead and measured SLOWER end-to-end (verified on real SPEC input:
 * 3m12s -> 4m42s). Reverted to the direct qsort call; the arc_id optimisation
 * in cost_compare (avoiding the a->id pointer chase) is kept since it is a
 * real, measured win with no downside. */


#ifdef _PROTO_
BASKET *primal_bea_mpp( LONG m, arc_t *arcs, arc_t *stop_arcs,
                          LONG* basket_sizes, BASKET** perm, int thread,
                          arc_t** end_arc, LONG step, LONG num_threads, LONG max_elems)
#else
arc_t *primal_bea_mpp( m, arcs, stop_arcs, basket_sizes, perm, thread, end_arc, step, num_threads, max_elems )
LONG m;
arc_t *arcs;
arc_t *stop_arcs;
LONG *basket_sizes;
BASKET** perm;
int thread;
arc_t** end_arc;
LONG step;
LONG num_threads;
LONG max_elems;
#endif
{
    LONG i, j, count, global_basket_size, next;
    arc_t *arc, *old_end_arc;
    cost_t red_cost;

    /* Re-validate survivors from previous iteration */
    for( i = 1, next = 0; i <= B && i <= basket_sizes[thread]; i++ )
    {
        arc      = perm[i]->a;
        count    = perm[i]->number;
        red_cost = arc->cost - NODE_POTENTIAL(arc->tail) + NODE_POTENTIAL(arc->head);
        if( count > 0 && ((red_cost < 0 && arc->ident == AT_LOWER)
            || (red_cost > 0 && arc->ident == AT_UPPER)) )
        {
            next++;
            perm[next]->a        = arc;
            perm[next]->cost     = red_cost;
            perm[next]->abs_cost = ABS(red_cost);
            perm[next]->arc_id   = arc->id;   /* keep arc_id in sync */
            perm[next]->number   = 0;
        }
    }

    basket_sizes[thread] = next;
    old_end_arc = *end_arc;

NEXT:
    arc = *end_arc + step;

    if (*end_arc >= full_group_end_arc)
        *end_arc = *end_arc + max_elems - 1;
    else
        *end_arc = *end_arc + max_elems;

    for ( ; arc < *end_arc; arc += num_threads) {
        if( arc->ident > BASIC )
        {
            red_cost = arc->cost - NODE_POTENTIAL(arc->tail) + NODE_POTENTIAL(arc->head);
            if( bea_is_dual_infeasible( arc, red_cost ) )
            {
                basket_sizes[thread]++;
                perm[basket_sizes[thread]]->a        = arc;
                perm[basket_sizes[thread]]->cost     = red_cost;
                perm[basket_sizes[thread]]->abs_cost = ABS(red_cost);
                perm[basket_sizes[thread]]->arc_id   = arc->id;  /* inline id */
                perm[basket_sizes[thread]]->number   = 0;
            }
        }
    }

    if( *end_arc >= stop_arcs )
        *end_arc = arcs;

    if (*end_arc != old_end_arc) {
#if (defined(_OPENMP) || defined(SPEC_OPENMP)) && !defined(SPEC_SUPPRESS_OPENMP) && !defined(SPEC_AUTO_SUPPRESS_OPENMP)
        #pragma omp barrier
#endif
        global_basket_size = 0;
        for (j = 0; j < num_threads; j++)
            global_basket_size += basket_sizes[j];
        if ( global_basket_size >= B)
            goto READY;
#if (defined(_OPENMP) || defined(SPEC_OPENMP)) && !defined(SPEC_SUPPRESS_OPENMP) && !defined(SPEC_AUTO_SUPPRESS_OPENMP)
        #pragma omp barrier
#endif
        goto NEXT;
    }

READY:
    perm[basket_sizes[thread] + 1]->number = -1;

    if (basket_sizes[thread] == 0)
        return NULL;

    /* Direct qsort over the collected candidates (basket_sizes[thread] is
     * almost always just over B in practice, so a partial-sort wrapper adds
     * a function-call layer with no real algorithmic benefit -- see note
     * above cost_compare). */
#if defined(SPEC)
    spec_qsort(perm + 1, basket_sizes[thread], sizeof(BASKET*),
            (int (*)(const void *, const void *))cost_compare);
#else
    qsort(perm + 1, basket_sizes[thread], sizeof(BASKET*),
            (int (*)(const void *, const void *))cost_compare);
#endif

    return perm[1];
}
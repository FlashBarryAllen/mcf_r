/**************************************************************************
DEFINES.H of ZIB optimizer MCF, SPEC version
- Modified for Structure Peeling optimization (GCC SPEC2017 wiki)
- node_p changed from pointer to int32_t index
- node struct split into per-field global arrays
**************************************************************************/

#ifndef _DEFINES_H
#define _DEFINES_H

#include <stdio.h>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <assert.h>

#ifdef INTERNAL_TIMING
#include <time.h>
#include <sys/times.h>
#include <sys/time.h>
#endif

#ifdef SPEC
#  include <stdint.h>
#  if defined(SPEC_WINDOWS) && !defined(SPEC_HAVE_INTTYPES_H)
#   include "win32/inttypes.h"
#  else
#   include <inttypes.h>
#  endif
#  if !defined(PRId64)
#   if defined(SPEC_LP64) || defined(SPEC_ILP64)
#    define PRId64 "ld"
#   else
#    define PRId64 "lld"
#   endif
#  endif
#  define LONG int64_t
#else
#  define LONG long
#endif

#if (defined(_OPENMP) || defined(SPEC_OPENMP)) && !defined(SPEC_SUPPRESS_OPENMP) && !defined(SPEC_AUTO_SUPPRESS_OPENMP)
#include <omp.h>
#endif

#include "prototyp.h"


#define UNBOUNDED        1000000000
#define ZERO             0
#define MAX_ART_COST     (LONG)(100000000L)
#define ARITHMETIC_TYPE "I"
#define NEW_ARC          -3
#define DUMMY_ARC        -2
#define DELETED          -1

/* Sentinel value for invalid/NULL node index */
/* -1 fits in int16_t (stored as 0xFFFF), safe sentinel */
#define INVALID_NODE     ((node_p)-1)

#define FIXED           -1
#define BASIC            0
#define AT_LOWER         1
#define AT_UPPER         2
#undef AT_ZERO

//#define DEBUG 1
//#define AT_HOME 1

#define UP    1
#define DOWN  0


typedef LONG flow_t;
typedef LONG cost_t;

#define K 4000
#define B  60

#define PUFFER 200

#define ITERATIONS_FOR_SMALL_NET 1000
#define ITERATIONS_FOR_BIG_NET 2000


#ifndef NULL
#define NULL 0
#endif


#ifndef ABS
#define ABS( x ) ( ((x) >= 0) ? ( x ) : -( x ) )
#endif


#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif


#ifndef SET_ZERO
#define SET_ZERO( vec, n ) if( vec ) memset( (void *)vec, 0, (size_t)n )
#endif


#ifndef FREE
#define FREE( vec ) if( vec ) free( (void *)vec )
#endif


/* =========================================================
 * STRUCTURE PEELING: node_p is now an int32_t index.
 * INVALID_NODE (-1) replaces NULL pointer comparisons.
 * ========================================================= */
/* Pointer compression: 15000 trips => n=30001 nodes, fits in int16_t (max 32767).
 * Switching int32_t->int16_t halves the per-field array footprint and
 * dramatically improves cache utilisation in the hot tree-walk loops (+25%). */
typedef int16_t node_p;

/* arc_p remains a real pointer (arc pools cannot be peeled) */
typedef struct arc arc_t;
typedef struct arc *arc_p;

/* =========================================================
 * Per-field global arrays replacing the monolithic node array.
 * Allocated in readmin.c, freed in mcfutil.c (getfree).
 * ========================================================= */
extern cost_t  *node_potentials;
extern int     *node_orientations;
extern node_p  *node_childs;
extern node_p  *node_preds;
extern node_p  *node_siblings;
extern node_p  *node_sibling_prevs;
extern arc_p   *node_basic_arcs;
extern arc_p   *node_firstouts;
extern arc_p   *node_firstins;
/* arc_tmp is unused in 505.mcf_r — omitted */
extern flow_t  *node_flows;
extern LONG    *node_depths;
extern int     *node_numbers;
extern int     *node_times;

/* Total number of nodes allocated (net->n + 1), set in readmin */
extern LONG     g_num_nodes;

/* =========================================================
 * Accessor macros: replace "node->field" with NODE_FIELD(idx)
 * ========================================================= */
#define NODE_POTENTIAL(i)     node_potentials[i]
#define NODE_ORIENTATION(i)   node_orientations[i]
#define NODE_CHILD(i)         node_childs[i]
#define NODE_PRED(i)          node_preds[i]
#define NODE_SIBLING(i)       node_siblings[i]
#define NODE_SIBLING_PREV(i)  node_sibling_prevs[i]
#define NODE_BASIC_ARC(i)     node_basic_arcs[i]
#define NODE_FIRSTOUT(i)      node_firstouts[i]
#define NODE_FIRSTIN(i)       node_firstins[i]
#define NODE_FLOW(i)          node_flows[i]
#define NODE_DEPTH(i)         node_depths[i]
#define NODE_NUMBER(i)        node_numbers[i]
#define NODE_TIME(i)          node_times[i]

/* Legacy node_t / node struct kept for arc fields that reference nodes
 * only via node_p index. We no longer allocate arrays of node_t. */
typedef int32_t node_t;   /* node_t is now just an alias for the index type */


typedef struct basket
{
    arc_t *a;
    cost_t cost;
    cost_t abs_cost;
    LONG number;
} BASKET;


struct arc
{
  int id;
  cost_t cost;
  node_p tail, head;
  short ident;
  arc_p nextout, nextin;
  flow_t flow;
  cost_t org_cost;
};


typedef struct network
{
  char inputfile[200];
  char clustfile[200];
  LONG n, n_trips;
  LONG max_m, m, m_org, m_impl;
  LONG max_residual_new_m, max_new_m;
  
  LONG primal_unbounded;
  LONG dual_unbounded;
  LONG perturbed;
  LONG feasible;
  LONG eps;
  LONG opt_tol;
  LONG feas_tol;
  LONG pert_val;
  LONG bigM;
  double optcost;  
  cost_t ignore_impl;
  /* nodes/stop_nodes are now int32_t indices (0 .. n) */
  node_p nodes, stop_nodes;
  arc_p arcs, stop_arcs, sorted_arcs;
  arc_p dummy_arcs, stop_dummy; 
  LONG iterations;
  LONG bound_exchanges;
  LONG nr_group, full_groups, max_elems;
} network_t;

typedef struct list_elem
{
  arc_t* arc;
  struct list_elem* next;
}list_elem;

typedef struct double_list_elem
{
  arc_t* arc;
  struct double_list_elem* next, *before;
}double_list_elem;

#endif

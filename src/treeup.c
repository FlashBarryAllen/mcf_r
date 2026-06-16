/**************************************************************************
TREEUP.C of ZIB optimizer MCF, SPEC version
Modified: Structure Peeling -- node pointers replaced by node_p (int32_t)
          indices; field accesses use NODE_xxx() macros.
**************************************************************************/

#include "treeup.h"


#ifdef _PROTO_
void update_tree( 
                 LONG cycle_ori,
                 LONG new_orientation,
                 flow_t delta,
                 flow_t new_flow,
                 node_p iplus,
                 node_p jplus,
                 node_p iminus,
                 node_p jminus,
                 node_p w,
                 arc_t *bea,
                 cost_t sigma,
                 flow_t feas_tol
                )
#else
void update_tree( cycle_ori, new_orientation, delta, new_flow,
                 iplus, jplus, iminus, jminus, w, bea, sigma, feas_tol )
     LONG cycle_ori;
     LONG new_orientation;
     flow_t delta;
     flow_t new_flow;
     node_p iplus, jplus;
     node_p iminus, jminus;
     node_p w;
     arc_t *bea;
     cost_t sigma;
     flow_t feas_tol;
#endif
{
    arc_t  *basic_arc_temp;
    arc_t  *new_basic_arc;
    node_p  father;
    node_p  temp;
    node_p  new_pred;
    LONG    orientation_temp;
    LONG    depth_temp;
    LONG    depth_iminus;
    LONG    new_depth;
    flow_t  flow_temp;

    if( (bea->tail == jplus && sigma < 0) ||
        (bea->tail == iplus && sigma > 0) )
        sigma = ABS(sigma);
    else
        sigma = -(ABS(sigma));

    father = iminus;
    NODE_POTENTIAL(father) += sigma;

 RECURSION:
    temp = NODE_CHILD(father);
    if( temp != INVALID_NODE )
    {
    ITERATION:
        NODE_POTENTIAL(temp) += sigma;
        father = temp;
        goto RECURSION;
    }
 TEST:
    if( father == iminus )
        goto CONTINUE;
    temp = NODE_SIBLING(father);
    if( temp != INVALID_NODE )
        goto ITERATION;
    father = NODE_PRED(father);
    goto TEST;

 CONTINUE:

    temp      = iplus;
    father    = NODE_PRED(temp);
    new_depth = depth_iminus = NODE_DEPTH(iminus);
    new_pred  = jplus;
    new_basic_arc = bea;

    while( temp != jminus )
    {
        node_p sib      = NODE_SIBLING(temp);
        node_p sib_prev = NODE_SIBLING_PREV(temp);

        if( sib != INVALID_NODE )
            NODE_SIBLING_PREV(sib) = sib_prev;
        if( sib_prev != INVALID_NODE )
            NODE_SIBLING(sib_prev) = sib;
        else
            NODE_CHILD(father) = sib;

        NODE_PRED(temp)    = new_pred;
        NODE_SIBLING(temp) = NODE_CHILD(new_pred);
        if( NODE_SIBLING(temp) != INVALID_NODE )
            NODE_SIBLING_PREV(NODE_SIBLING(temp)) = temp;
        NODE_CHILD(new_pred)    = temp;
        NODE_SIBLING_PREV(temp) = INVALID_NODE;

        orientation_temp = !(NODE_ORIENTATION(temp));
        if( orientation_temp == cycle_ori )
            flow_temp = NODE_FLOW(temp) + delta;
        else
            flow_temp = NODE_FLOW(temp) - delta;

        basic_arc_temp = NODE_BASIC_ARC(temp);
        depth_temp     = NODE_DEPTH(temp);

        NODE_ORIENTATION(temp) = (int)new_orientation;
        NODE_FLOW(temp)        = new_flow;
        NODE_BASIC_ARC(temp)   = new_basic_arc;
        NODE_DEPTH(temp)       = new_depth;

        new_pred        = temp;
        new_orientation = orientation_temp;
        new_flow        = flow_temp;
        new_basic_arc   = basic_arc_temp;
        new_depth       = depth_iminus - depth_temp;
        temp   = father;
        father = NODE_PRED(temp);
    }

    if( delta > feas_tol )
    {
        for( temp = jminus; temp != w; temp = NODE_PRED(temp) )
        {
            NODE_DEPTH(temp) -= depth_iminus;
            if( NODE_ORIENTATION(temp) != cycle_ori )
                NODE_FLOW(temp) += delta;
            else
                NODE_FLOW(temp) -= delta;
        }
        for( temp = jplus; temp != w; temp = NODE_PRED(temp) )
        {
            NODE_DEPTH(temp) += depth_iminus;
            if( NODE_ORIENTATION(temp) == cycle_ori )
                NODE_FLOW(temp) += delta;
            else
                NODE_FLOW(temp) -= delta;
        }
    }
    else
    {
        for( temp = jminus; temp != w; temp = NODE_PRED(temp) )
            NODE_DEPTH(temp) -= depth_iminus;
        for( temp = jplus; temp != w; temp = NODE_PRED(temp) )
            NODE_DEPTH(temp) += depth_iminus;
    }
}

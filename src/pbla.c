/**************************************************************************
PBLA.C of ZIB optimizer MCF, SPEC version
Modified: Structure Peeling -- node_p is int32_t index.
**************************************************************************/

#include "pbla.h"


#define TEST_MIN( nod, ex, comp, op ) \
{ \
      if( *delta op (comp) ) \
      { \
            iminus = nod; \
            *delta = (comp); \
            *xchange = ex; \
      } \
}


#ifdef _PROTO_
node_p primal_iminus( 
                      flow_t *delta,
                      LONG *xchange,
                      node_p iplus,
                      node_p jplus,
                      node_p *w
                    )
#else
node_p primal_iminus( delta, xchange, iplus, jplus, w )
    flow_t *delta;
    LONG *xchange;
    node_p iplus, jplus;
    node_p *w;
#endif
{
    node_p iminus = INVALID_NODE;

    while( iplus != jplus )
    {
        if( NODE_DEPTH(iplus) < NODE_DEPTH(jplus) )
        {
            if( NODE_ORIENTATION(iplus) )
                TEST_MIN( iplus, 0, NODE_FLOW(iplus), > )
            else if( NODE_PRED(NODE_PRED(iplus)) != INVALID_NODE )
                TEST_MIN( iplus, 0, (flow_t)1 - NODE_FLOW(iplus), > )
            iplus = NODE_PRED(iplus);
        }
        else
        {
            if( !NODE_ORIENTATION(jplus) )
                TEST_MIN( jplus, 1, NODE_FLOW(jplus), >= )
            else if( NODE_PRED(NODE_PRED(jplus)) != INVALID_NODE )
                TEST_MIN( jplus, 1, (flow_t)1 - NODE_FLOW(jplus), >= )
            jplus = NODE_PRED(jplus);
        }
    }

    *w = iplus;

    return iminus;
}

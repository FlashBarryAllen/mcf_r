/**************************************************************************
PFLOWUP.C of ZIB optimizer MCF, SPEC version
Modified: Structure Peeling -- node_p is int32_t index.
**************************************************************************/

#include "pflowup.h"


#ifdef _PROTO_
void primal_update_flow( 
                 node_p iplus,
                 node_p jplus,
                 node_p w
                 )
#else
void primal_update_flow( iplus, jplus, w )
    node_p iplus, jplus;
    node_p w;
#endif
{
    for( ; iplus != w; iplus = NODE_PRED(iplus) )
    {
        if( NODE_ORIENTATION(iplus) )
            NODE_FLOW(iplus) = (flow_t)0;
        else
            NODE_FLOW(iplus) = (flow_t)1;
    }

    for( ; jplus != w; jplus = NODE_PRED(jplus) )
    {
        if( NODE_ORIENTATION(jplus) )
            NODE_FLOW(jplus) = (flow_t)1;
        else
            NODE_FLOW(jplus) = (flow_t)0;
    }
}

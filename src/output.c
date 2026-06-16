/**************************************************************************
OUTPUT.C of ZIB optimizer MCF, SPEC version
Modified: Structure Peeling -- node field accesses use NODE_xxx() macros.
**************************************************************************/

#include "output.h"


#ifdef _PROTO_
LONG write_circulations(
                   char *outfile,
                   network_t *net
                   )
#else
LONG write_circulations( outfile, net )
     char *outfile;
     network_t *net;
#endif
{
  FILE *out = NULL;
  arc_t *block;
  arc_t *arc;
  arc_t *arc2;
  arc_t *first_impl = net->stop_arcs - net->m_impl;
  node_p root_node = (node_p)net->n;  /* root index */

  if(( out = fopen( outfile, "w" )) == NULL )
    return -1;

  refresh_neighbour_lists( net, &getArcPosition );

  for( block = NODE_FIRSTOUT(root_node); block; block = block->nextout )
  {
    if( block->flow )
    {
      fprintf( out, "()\n" );

      arc = block;
      while( arc )
      {
        if( arc >= first_impl )
          fprintf( out, "***\n" );

        node_p h = arc->head;
        fprintf( out, "%d\n", -NODE_NUMBER(h) );

        /* arc->head[net->n_trips].firstout  =>  NODE_FIRSTOUT(h + net->n_trips) */
        arc2 = NODE_FIRSTOUT(h + net->n_trips);
        for( ; arc2; arc2 = arc2->nextout )
          if( arc2->flow )
            break;
        if( !arc2 )
        {
          fclose( out );
          return -1;
        }

        if( NODE_NUMBER(arc2->head) )
          arc = arc2;
        else
          arc = NULL;
      }
    }
  }

  fclose(out);
  return 0;
}


#ifdef _PROTO_
LONG write_objective_value(
                   char *outfile,
                   network_t *net
                   )
#else
LONG write_objective_value( outfile, net )
     char *outfile;
     network_t *net;
#endif
{
  FILE *out = NULL;

  if(( out = fopen( outfile, "w" )) == NULL )
    return -1;

  fprintf( out, "%.0f\n", flow_cost(net) );

  fclose(out);
  return 0;
}

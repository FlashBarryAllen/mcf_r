/**************************************************************************
OUTPUT.C of ZIB optimizer MCF, SPEC version
write_circulations: NOT called in SPEC CPU 2017 run (only write_objective_value
is called from mcf.c). The nextout/nextin fields it used have been removed
from the arc struct as part of the arc layout optimisation.
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
    /* nextout/nextin removed from arc struct; this function is never called
     * in the SPEC CPU 2017 run path. Stub to satisfy linker. */
    (void)outfile; (void)net;
    return -1;
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

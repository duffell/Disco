
#include "../paul.h"

void restart( struct domain * theDomain ){

   //This file is included so that the code can easily
   //compile without hdf5.  Just don't try to restart
   //from file because that's impossible without hdf5.

   int rank = theDomain->rank;
   if( rank==0 ) printf("**********\nDanger! You should not be here.\nYou are trying to restart from file but you didn't compile with HDF5 input\n**********\n");

}


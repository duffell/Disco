
#include "paul.h"

int mpiSetup( struct domain * theDomain , int argc, char * argv[] ){

   MPI_Comm_size(MPI_COMM_WORLD,&(theDomain->size));
   int size = theDomain->size;
   int * dim_rank   = theDomain->dim_rank;
   int * dim_size   = theDomain->dim_size;
   int * left_rank  = theDomain->left_rank;
   int * right_rank = theDomain->right_rank;

   int NUM_Z = theDomain->theParList.Num_Z;

   //A bunch of crap I need to create the cart.
   int wraparound[2] = {1,1};
   dim_size[0] = 0;
   dim_size[1] = 0;
   dim_rank[0] = 0;
   dim_rank[1] = 0;
   int reorder = 1, ndims = 2;

   if( NUM_Z == 1 ){ ndims = 1; dim_size[1] = 1; }

   //This gives me the most efficient factorization of
   //the number of processes being used into a j-by-k grid.
   MPI_Dims_create(size,ndims,dim_size);

   //Just in case the dimensions don't make sense; I don't know how
   //the MPI_Dims_create funcition works, so you never know.
   if( dim_size[0]*dim_size[1] != size ){
      printf("Error: dimensions don't jive.\n");
      return(1);
   }

   //Creating the cart!  Woohoo.
   MPI_Cart_create(MPI_COMM_WORLD,ndims,dim_size,wraparound,reorder,&(theDomain->theComm));
   MPI_Comm_rank(theDomain->theComm,&(theDomain->rank));
   int rank = theDomain->rank;
   MPI_Cart_coords(theDomain->theComm,rank,ndims,dim_rank);

   if(rank==0){printf("dim1 = %d, dim2 = %d\n",dim_size[0],dim_size[1]);}

   //Determine ranks to my left and right in all 2 dimensions.
   int next_rank[2];
   int i;
   for( i=0 ; i<2 ; ++i ){
      next_rank[i] = dim_rank[i];
   }
   for( i=0 ; i<2 ; ++i ){
      next_rank[i] += 1;
      MPI_Cart_rank(theDomain->theComm,next_rank,&(right_rank[i]));
      next_rank[i] -= 2;
      MPI_Cart_rank(theDomain->theComm,next_rank,&(left_rank[i]));
      next_rank[i] += 1;
   }

   return(0);
}


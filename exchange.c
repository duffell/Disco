
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "paul.h"

struct cell_lite{
   double prim[NUM_Q];
   double cons[NUM_Q];
   double RKcons[NUM_Q];
   double piph;
   double RKpiph;
   double wiph;
};

void generate_mpi_cell( MPI_Datatype * cell_mpi ){

   struct cell_lite test;
   int count = 6;
   int blocksize[]      = {NUM_Q,NUM_Q,NUM_Q,1,1,1};
   MPI_Datatype types[] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
   MPI_Aint offsets[6];

   offsets[0] = (char *)&(test.prim)   - (char *)(&test);
   offsets[1] = (char *)&(test.cons)   - (char *)(&test);
   offsets[2] = (char *)&(test.RKcons) - (char *)(&test);
   offsets[3] = (char *)&(test.piph)   - (char *)(&test);
   offsets[4] = (char *)&(test.RKpiph) - (char *)(&test);
   offsets[5] = (char *)&(test.wiph)   - (char *)(&test);

   MPI_Type_create_struct( count , blocksize , offsets , types , cell_mpi );
   MPI_Type_commit( cell_mpi );

}

void copy_cell_to_lite( struct cell * c , struct cell_lite * cl ){
  
   memcpy( cl->prim   , c->prim   , NUM_Q*sizeof(double) ); 
   memcpy( cl->cons   , c->cons   , NUM_Q*sizeof(double) ); 
   memcpy( cl->RKcons , c->RKcons , NUM_Q*sizeof(double) );
   cl->piph   = c->piph;
   cl->RKpiph = c->RKpiph;
   cl->wiph   = c->wiph; 

}

void copy_lite_to_cell( struct cell_lite * cl , struct cell * c ){ 

   memcpy( c->prim   , cl->prim   , NUM_Q*sizeof(double) ); 
   memcpy( c->cons   , cl->cons   , NUM_Q*sizeof(double) ); 
   memcpy( c->RKcons , cl->RKcons , NUM_Q*sizeof(double) );
   c->piph   = cl->piph;
   c->RKpiph = cl->RKpiph;
   c->wiph   = cl->wiph;

}

void generate_sendbuffer( struct domain * theDomain , int rnum , int znum , int dim , int * nijk , int * indexL , int * indexR , struct cell_lite * pl , struct cell_lite * pr , int dn1 , int dn2 , int mode ){

   struct cell ** theCells = theDomain->theCells;
   int Periodic = theDomain->theParList.Z_Periodic;
   if( dim == 0 ) Periodic = 0;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   int Nr = theDomain->Nr;
   int * Np = theDomain->Np;
   int i,j,k;

   int iL = 0;
   int iR = 0;
   for( j=0 ; j<rnum ; ++j ){
      nijk[0]=j;
      for( k=0 ; k<znum ; ++k ){
         nijk[1]=k;
         nijk[dim] += dn1;

         int jk = nijk[0]+Nr*nijk[1];
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            if( mode==1 ){
               copy_cell_to_lite( c , pl+iL );
            }else if( mode==2 && ( dim_rank[dim] != 0 || Periodic ) ){
               copy_lite_to_cell( pl+iL , c );
            }
            ++iL;
         }

         nijk[dim] += dn2-dn1;

         jk = nijk[0]+Nr*nijk[1];
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            if( mode==1 ){
               copy_cell_to_lite( c , pr+iR );
            }else if( mode==2 && ( dim_rank[dim] != dim_size[dim]-1 || Periodic ) ){
               copy_lite_to_cell( pr+iR , c );
            }
            ++iR;
         }

         nijk[dim] -= dn2;
      }
   }

   *indexL = iL;
   *indexR = iR;
}

void generate_intbuffer( struct domain * theDomain , int rnum , int znum , int dim , int * nijk , int * indexL , int * indexR , int * Npl , int * Npr , int dn1 , int dn2 , int mode ){

   struct cell ** theCells = theDomain->theCells;
   int Periodic = theDomain->theParList.Z_Periodic;
   if( dim == 0 ) Periodic = 0;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   int Nr = theDomain->Nr;
   int * Np = theDomain->Np;
   int j,k; 

   int iL = 0;
   int iR = 0; 
   for( j=0 ; j<rnum ; ++j ){
      nijk[0]=j;
      for( k=0 ; k<znum ; ++k ){
         nijk[1]=k;
         nijk[dim] += dn1; 

         int jk = nijk[0]+Nr*nijk[1];
         if( mode==1 ){
            Npl[ iL ] = Np[jk];
         }else if( mode==2 && ( dim_rank[dim] != 0 || Periodic ) ){ 
            Np[jk] = Npl[ iL ];
            theCells[jk] = (struct cell *) realloc( theCells[jk] , Np[jk]*sizeof(struct cell) );
         }    
         ++iL;

         nijk[dim] += dn2-dn1;

         jk = nijk[0]+Nr*nijk[1];
         if( mode==1 ){
            Npr[ iR ] = Np[jk];
         }else if( mode==2 && ( dim_rank[dim] != dim_size[dim]-1 || Periodic ) ){
            Np[jk] = Npr[ iR ];
            theCells[jk] = (struct cell *) realloc( theCells[jk] , Np[jk]*sizeof(struct cell) );
         }    
         ++iR;

         nijk[dim] -= dn2; 

      }    
   }
   *indexL = iL; 
   *indexR = iR;

}

void exchangeData( struct domain * theDomain , int dim ){

   MPI_Datatype cell_mpi = {0}; 
   generate_mpi_cell( &cell_mpi );

   MPI_Comm grid_comm = theDomain->theComm;
   int * left_rank = theDomain->left_rank;
   int * right_rank = theDomain->right_rank;
   int Ng = theDomain->Ng;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;

   int tag = 0;
   MPI_Status status;
   int nijk[2];
   int rnum = Nr;
   int znum = Nz;
   int NN;

   if( dim == 0 ){
      NN = Nr;
      rnum = Ng;
   }else{
      NN = Nz;
      znum = Ng;
   }

   int indexL,indexR;
   int send_sizeL = 0;
   int send_sizeR = 0;
   int recv_sizeL = 0;
   int recv_sizeR = 0;

////////////////
//Send Np[jk]...
////////////////

//Count the number of Np's to send...
   generate_intbuffer( theDomain , rnum , znum , dim , nijk , &indexL , &indexR , NULL    , NULL    , Ng , NN-2*Ng , 0 );
   send_sizeL = indexL;
   send_sizeR = indexR;
//Tell your neighbor how many to expect...
   MPI_Sendrecv( &send_sizeL , 1 , MPI_INT ,  left_rank[dim] , tag   , 
                 &recv_sizeR , 1 , MPI_INT , right_rank[dim] , tag  , grid_comm , &status);
   MPI_Sendrecv( &send_sizeR , 1 , MPI_INT , right_rank[dim] , tag+1 , 
                 &recv_sizeL , 1 , MPI_INT ,  left_rank[dim] , tag+1, grid_comm , &status);

   int Nl_send[send_sizeL];
   int Nr_send[send_sizeR];
   int Nl_recv[recv_sizeL];
   int Nr_recv[recv_sizeR];
//Build up list of ints to send...
   generate_intbuffer( theDomain , rnum , znum , dim , nijk , &indexL , &indexR , Nl_send , Nr_send , Ng , NN-2*Ng , 1 );
//Send!
   MPI_Sendrecv( Nl_send , send_sizeL , MPI_INT ,  left_rank[dim] , tag+2 ,
                 Nr_recv , recv_sizeR , MPI_INT , right_rank[dim] , tag+2, grid_comm , &status);
   MPI_Sendrecv( Nr_send , send_sizeR , MPI_INT , right_rank[dim] , tag+3 ,
                 Nl_recv , recv_sizeL , MPI_INT ,  left_rank[dim] , tag+3, grid_comm , &status);
//Now take the list of ints and put them where they belong...
   generate_intbuffer( theDomain , rnum , znum , dim , nijk , &indexL , &indexR , Nl_recv , Nr_recv , 0  , NN-Ng   , 2 );

////////////
//Send Cells
////////////

//Count the number of cell_lites to send...
   generate_sendbuffer( theDomain , rnum , znum , dim , nijk , &indexL , &indexR , NULL    , NULL    , Ng , NN-2*Ng , 0 );
   send_sizeL = indexL;
   send_sizeR = indexR;
//Tell you neighbor how many cells to expect...
   MPI_Sendrecv( &send_sizeL , 1 , MPI_INT ,  left_rank[dim] , tag   ,
                 &recv_sizeR , 1 , MPI_INT , right_rank[dim] , tag  , grid_comm , &status);
   MPI_Sendrecv( &send_sizeR , 1 , MPI_INT , right_rank[dim] , tag+1 ,
                 &recv_sizeL , 1 , MPI_INT ,  left_rank[dim] , tag+1, grid_comm , &status);

   struct cell_lite pl_send[send_sizeL];
   struct cell_lite pr_send[send_sizeR];
   struct cell_lite pl_recv[recv_sizeL];
   struct cell_lite pr_recv[recv_sizeR];
//Build up list of cells to send...
   generate_sendbuffer( theDomain , rnum , znum , dim , nijk , &indexL , &indexR , pl_send , pr_send , Ng , NN-2*Ng , 1 );
//Send!
   MPI_Sendrecv( pl_send , send_sizeL , cell_mpi ,  left_rank[dim] , tag+2 ,
                 pr_recv , recv_sizeR , cell_mpi , right_rank[dim] , tag+2, grid_comm , &status);
   MPI_Sendrecv( pr_send , send_sizeR , cell_mpi , right_rank[dim] , tag+3 ,
                 pl_recv , recv_sizeL , cell_mpi ,  left_rank[dim] , tag+3, grid_comm , &status);
//Now take the list of cells and put them into the appropriate locations...
   generate_sendbuffer( theDomain , rnum , znum , dim , nijk , &indexL , &indexR , pl_recv , pr_recv , 0  , NN-Ng   , 2 );

   MPI_Type_free( &cell_mpi );
}


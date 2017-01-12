
#include "../paul.h"
#include <hdf5.h>
#include <string.h>

void getH5dims( char * file , char * group , char * dset , hsize_t * dims ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );
   hid_t h5spc = H5Dget_space( h5dst );

   H5Sget_simple_extent_dims( h5spc , dims , NULL);

   H5Sclose( h5spc );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void readSimple( char * file , char * group , char * dset , void * data , hid_t type ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   H5Dread( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void readPatch( char * file , char * group , char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   hsize_t mdims[dim];
   hsize_t fdims[dim];

   hsize_t fstart[dim];
   hsize_t fstride[dim];
   hsize_t fcount[dim];
   hsize_t fblock[dim];

   int d;
   for( d=0 ; d<dim ; ++d ){
      mdims[d] = loc_size[d];
      fdims[d] = glo_size[d];

      fstart[d]  = start[d];
      fstride[d] = 1;
      fcount[d]  = loc_size[d];
      fblock[d]  = 1;
   }
   hid_t mspace = H5Screate_simple(dim,mdims,NULL);
   hid_t fspace = H5Screate_simple(dim,fdims,NULL);

   H5Sselect_hyperslab( fspace , H5S_SELECT_SET , fstart , fstride , fcount , fblock );

   H5Dread( h5dst , type , mspace , fspace , H5P_DEFAULT , data );

   H5Sclose( mspace );
   H5Sclose( fspace );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}


void Doub2Cell( double * Q , struct cell * c ){
   int q;
   for( q=0 ; q<NUM_Q ; ++q ) c->prim[q] = Q[q];
   for( q=0 ; q<NUM_FACES ; ++q ) c->Phi[q] = Q[NUM_Q+q];
   c->piph = Q[NUM_Q+NUM_FACES];
}

int getN0( int , int , int );
void freeDomain( struct domain * );

void setPlanetParams( struct domain * );
void initializePlanets( struct planet * );
int num_diagnostics( void );
int get_num_rzFaces( int , int , int );

void restart( struct domain * theDomain ){

   //This code has not been bug-tested in 3D.
   //The ordering of indices in particular hasn't been checked
   //i.e. should I use jk = j + Nt*k or jk = j*Np + k?

   freeDomain( theDomain );

   int Ng = theDomain->Ng;
   int Z_Periodic = theDomain->theParList.Z_Periodic;
   int rank = theDomain->rank;
   int size = theDomain->size;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;

   int i,j,k;

   char filename[256];
   strcpy(filename,"input.h5");
   char group1[256];
   strcpy(group1,"Grid");
   char group2[256];
   strcpy(group2,"Data");

   hsize_t dims[3];

   if( rank==0 ) printf("Restarting from file...\n");

   int NUM_R,NUM_Z;
   double tstart;
   //Read the time from "T" and get Nt and Np from 
   //the dimensions of "Index".  Broadcast to all ranks.
   if(rank==0){
      readSimple( filename , group1 ,"T", &tstart , H5T_NATIVE_DOUBLE );
      getH5dims( filename , group1 ,"Index", dims );
      NUM_Z = dims[0];
      NUM_R = dims[1];
   }
   MPI_Bcast( &NUM_R  , 1 , MPI_INT    , 0 , theDomain->theComm );
   MPI_Bcast( &NUM_Z  , 1 , MPI_INT    , 0 , theDomain->theComm );
   MPI_Bcast( &tstart , 1 , MPI_DOUBLE , 0 , theDomain->theComm );
   theDomain->theParList.Num_R = NUM_R;
   theDomain->theParList.Num_Z = NUM_Z;
   theDomain->t = tstart;
  
   //The following is very similar to equivalent code in gridsetup.c
   //Now you're just doing the process over because you're restarting
   //from file. 
   int N0r = getN0( dim_rank[0]   , dim_size[0] , NUM_R );
   int N1r = getN0( dim_rank[0]+1 , dim_size[0] , NUM_R );
   if( dim_rank[0] != 0 ) N0r -= Ng;
   if( dim_rank[0] != dim_size[0]-1 ) N1r += Ng;
   int Nr = N1r-N0r;

   int N0z = getN0( dim_rank[1]   , dim_size[1] , NUM_Z );
   int N1z = getN0( dim_rank[1]+1 , dim_size[1] , NUM_Z );
   if( NUM_Z > 1 ){
      if( dim_rank[1] != 0 || Z_Periodic ) N0z -= Ng;
      if( dim_rank[1] != dim_size[1]-1 || Z_Periodic ) N1z += Ng;
   }
   int Nz = N1z-N0z;

   theDomain->Nr = Nr;
   theDomain->Nz = Nz;

   theDomain->Np    = (int *)    malloc( Nr*Nz*sizeof(int) );
   theDomain->r_jph = (double *) malloc( (Nr+1)*sizeof(double) );
   theDomain->z_kph = (double *) malloc( (Nz+1)*sizeof(double) );

   theDomain->theCells = (struct cell **) malloc( Nr*Nz*sizeof( struct cell * ) );

   struct cell ** theCells = theDomain->theCells;

   //The following must happen in serial because different processors
   //will try to read from the same file.  In principle we can write
   //this using parallel HDF5, but I'd rather cross that bridge 
   //when I come to it.
   int nrk;
   int Nq=0;
   for( nrk=0 ; nrk<size ; ++nrk ){
   if( rank==nrk ){

      //Read the R values of the grid...
      int start1[1] = {N0r};
      int loc_size1[1] = {Nr+1};
      int glo_size1[1] = {NUM_R+1};
      double r_jph[Nr+1];
      readPatch( filename , group1 ,"r_jph", r_jph , H5T_NATIVE_DOUBLE , 1 , start1 , loc_size1 , glo_size1 ); 
      memcpy( theDomain->r_jph , r_jph , (Nr+1)*sizeof(double) );
 
      //Read the Z values of the grid...
      start1[0]    = N0z;
      loc_size1[0] = Nz+1;
      glo_size1[0] = NUM_Z+1;
      double z_kph[Nz+1];
      readPatch( filename , group1 ,"z_kph", z_kph , H5T_NATIVE_DOUBLE , 1 , start1 , loc_size1 , glo_size1 );
      memcpy( theDomain->z_kph , z_kph , (Nz+1)*sizeof(double) );

      ++(theDomain->r_jph);
      ++(theDomain->z_kph);

      //Read the indexing information so you know how to read in
      //The radial tracks of data which are coming up...
      int start2[2]   = {N0z,N0r};
      int loc_size2[2] = {Nz,Nr};
      int glo_size2[2] = {NUM_Z,NUM_R};
      int Np[Nr*Nz];
      int Index[Nr*Nz];

      printf("%d %d\n", NUM_Z, NUM_R);
      printf("%d %d\n", Nz, Nr);
      printf("%d %d\n", N0z, N0r);

      readPatch( filename , group1 ,"Np"   , Np    , H5T_NATIVE_INT , 2 , start2 , loc_size2 , glo_size2 );
      readPatch( filename , group1 ,"Index", Index , H5T_NATIVE_INT , 2 , start2 , loc_size2 , glo_size2 );
      memcpy( theDomain->Np , Np , Nr*Nz*sizeof(int) );

      getH5dims( filename , group2 ,"Cells", dims );
      int Nc = dims[0];
      Nq = dims[1];

      //Read in each radial track one at a time, because
      //you don't really know where the different radial
      //tracks belong in memory, they might not be 
      //contiguous, and they definitely are not in a rectangular block.
      for( k=0 ; k<Nz ; ++k ){
         for( j=0 ; j<Nr ; ++j ){
            int jk = j+Nr*k;
            start2[0] = Index[jk];
            start2[1] = 0;
            loc_size2[0] = Np[jk];
            loc_size2[1] = Nq;
            glo_size2[0] = Nc;
            glo_size2[1] = Nq;
            double TrackData[Np[jk]*Nq];
            readPatch( filename , group2 ,"Cells", TrackData , H5T_NATIVE_DOUBLE , 2 , start2 , loc_size2 , glo_size2 );
            theDomain->theCells[jk] = (struct cell *) malloc( Np[jk]*sizeof(struct cell) );
            for( i=0 ; i<Np[jk] ; ++i ){
               struct cell * c = theCells[jk]+i;
               Doub2Cell( TrackData + i*Nq , c );
            }
         }
      }

   }
   MPI_Barrier(theDomain->theComm);
   }
   if( Nq != NUM_Q+NUM_FACES+1 ){ if(rank==0)printf("Ummm, I got an hdf5 read error. Check NUM_Q.\n"); exit(1); }

   setPlanetParams( theDomain );
   int Npl = theDomain->Npl;
   theDomain->thePlanets = (struct planet *) malloc( Npl*sizeof(struct planet) );
   initializePlanets( theDomain->thePlanets );

   double num_tools = num_diagnostics();
   theDomain->num_tools = num_tools;
   theDomain->theTools.t_avg = 0.0; 
   theDomain->theTools.Qr = (double *) calloc( Nr*num_tools , sizeof(double) );

   theDomain->N_ftracks_r = get_num_rzFaces( theDomain->Nr , theDomain->Nz , 1 ); 
   theDomain->N_ftracks_z = get_num_rzFaces( theDomain->Nr , theDomain->Nz , 2 ); 
   theDomain->fIndex_r = (int *) malloc( (theDomain->N_ftracks_r+1)*sizeof(int) );
   theDomain->fIndex_z = (int *) malloc( (theDomain->N_ftracks_z+1)*sizeof(int) );

}


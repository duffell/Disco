
#include "../paul.h"
#include <hdf5.h>

void createFile( char * fname ){
   hid_t h5file = H5Fcreate( fname , H5F_ACC_TRUNC , H5P_DEFAULT , H5P_DEFAULT );
   H5Fclose( h5file );
}

void createGroup( char * fname , char * gname ){
   hid_t h5file = H5Fopen( fname , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5group = H5Gcreate1( h5file , gname , 0 );
   H5Gclose( h5group );
   H5Fclose( h5file );
}

void createDataset( char * fname , char * gname , char * dname , int dim , hsize_t * fdims , hid_t type ){
   hid_t h5file  = H5Fopen( fname , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5group = H5Gopen1( h5file , gname );
   hid_t fspace  = H5Screate_simple(dim,fdims,NULL);
   hid_t h5dset  = H5Dcreate1( h5group , dname , type , fspace , H5P_DEFAULT );
   H5Sclose( fspace );
   H5Dclose( h5dset );
   H5Gclose( h5group );
   H5Fclose( h5file );
}

void writeSimple( char * file , char * group , char * dset , void * data , hid_t type ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   H5Dwrite( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void writePatch( char * file , char * group , char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size){
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

   H5Dwrite( h5dst , type , mspace , fspace , H5P_DEFAULT , data );

   H5Sclose( mspace );
   H5Sclose( fspace );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

int Cell2Doub( struct cell * c , double * Q , int mode ){
   if( mode==0 ) return(NUM_Q+NUM_FACES+1); else{
      int q;
      for( q=0 ; q<NUM_Q ; ++q ) Q[q] = c->prim[q];
      for( q=0 ; q<NUM_FACES ; ++q ) Q[NUM_Q+q] = c->Phi[q];
      Q[NUM_Q+NUM_FACES] = c->piph;
      return(0);
   }
}

double get_dV( double * , double * );
void prim2cons( double * , double * , double *, double );
void cons2prim( double * , double * , double *, double );

void zero_diagnostics( struct domain * );
void avg_diagnostics( struct domain * );

void output( struct domain * theDomain , char * filestart ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Npl = theDomain->Npl;
   int Ng = theDomain->Ng;
   int Nr_Tot = theDomain->theParList.Num_R;
   int Nz_Tot = theDomain->theParList.Num_Z;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int rank = theDomain->rank;
   int size = theDomain->size;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   int Z_Periodic = theDomain->theParList.Z_Periodic;

   int NpDat = 6;
   int Ntools = theDomain->num_tools;
   avg_diagnostics( theDomain );

   char filename[256];
   sprintf(filename,"%s.h5",filestart);
   int jmin = 0;
   if( dim_rank[0] != 0 ) jmin = Ng;
   int jmax = Nr;
   if( dim_rank[0] != dim_size[0]-1 ) jmax = Nr-Ng;
   int kmin = 0;
   if( dim_rank[1] != 0 || Z_Periodic ) kmin = Ng;
   int kmax = Nz;
   if( dim_rank[1] != dim_size[1]-1 || Z_Periodic ) kmax = Nz-Ng;

   int Ntot = 0;
   int j,k;
   for( j=jmin ; j<jmax ; ++j ){
      for( k=kmin ; k<kmax ; ++k ){
         int jk = j+Nr*k;
         Ntot += Np[jk];
      }
   }
   int myNtot = Ntot;
   MPI_Allreduce( MPI_IN_PLACE , &Ntot  , 1 , MPI_INT , MPI_SUM , theDomain->theComm );

   int Ndoub = Cell2Doub(NULL,NULL,0);

   hsize_t fdims1[1];
   hsize_t fdims2[2];
   if( rank==0 ){
      printf("Writing Checkpoint...\n");
      
      createFile(filename);
      createGroup(filename,"Grid");

      fdims1[0] = 1;
      createDataset(filename,"Grid","T",1,fdims1,H5T_NATIVE_DOUBLE);
      fdims1[0] = Nr_Tot+1;
      createDataset(filename,"Grid","r_jph",1,fdims1,H5T_NATIVE_DOUBLE);
      fdims1[0] = Nz_Tot+1;
      createDataset(filename,"Grid","z_kph",1,fdims1,H5T_NATIVE_DOUBLE);
      fdims2[0] = Nz_Tot;
      fdims2[1] = Nr_Tot;
      createDataset(filename,"Grid","Index",2,fdims2,H5T_NATIVE_INT);
      createDataset(filename,"Grid","Np",2,fdims2,H5T_NATIVE_INT);
      createDataset(filename,"Grid","Id_phi0",2,fdims2,H5T_NATIVE_INT);

      createGroup(filename,"Data");

      fdims2[0] = Ntot;
      fdims2[1] = Ndoub;
      createDataset(filename,"Data","Cells",2,fdims2,H5T_NATIVE_DOUBLE);
      fdims2[0] = Npl;
      fdims2[1] = NpDat;
      createDataset(filename,"Data","Planets",2,fdims2,H5T_NATIVE_DOUBLE);
      fdims2[0] = Nr_Tot;
      fdims2[1] = Ntools;
      createDataset(filename,"Data","Radial_Diagnostics",2,fdims2,H5T_NATIVE_DOUBLE);
   }
   MPI_Barrier( theDomain->theComm );
   if( rank==0 ){
      writeSimple(filename,"Grid","T",&(theDomain->t),H5T_NATIVE_DOUBLE);
      double PlanetData[Npl*NpDat];
      int p;
      for( p=0 ; p<Npl ; ++p ){
         struct planet * pl = theDomain->thePlanets+p;
         PlanetData[NpDat*p + 0] = pl->M;
         PlanetData[NpDat*p + 1] = pl->vr;
         PlanetData[NpDat*p + 2] = pl->omega;
         PlanetData[NpDat*p + 3] = pl->r;
         PlanetData[NpDat*p + 4] = pl->phi;
         PlanetData[NpDat*p + 5] = pl->eps;
      }
      writeSimple(filename,"Data","Planets",PlanetData,H5T_NATIVE_DOUBLE);
   }

   int jSize = jmax-jmin;
   int kSize = kmax-kmin;
   int nrk;
   int j0,k0;
   int jSum = 0;
   for( nrk=0 ; nrk < dim_size[0] ; ++nrk ){
      if( nrk == dim_rank[0] ){
         j0 = jSum;
         if( dim_rank[1] == 0 ) jSum += jSize;
      }
      MPI_Allreduce( MPI_IN_PLACE , &jSum , 1 , MPI_INT , MPI_MAX , theDomain->theComm );
   }
   int kSum = 0;
   for( nrk=0 ; nrk < dim_size[1] ; ++nrk ){
      if( nrk == dim_rank[1] ){
         k0 = kSum;
         if( dim_rank[0] == 0 ) kSum += kSize;
      }
      MPI_Allreduce( MPI_IN_PLACE , &kSum , 1 , MPI_INT , MPI_MAX , theDomain->theComm );
   }


   if( Nr_Tot == 1 ){ j0 = 0; jSum = 1; }
   if( Nz_Tot == 1 ){ k0 = 0; kSum = 1; }


   int * Index   = (int *) malloc( jSize*kSize*sizeof(int) );
   int * Size    = (int *) malloc( jSize*kSize*sizeof(int) );
   int * Id_phi0 = (int *) malloc( jSize*kSize*sizeof(int) );
   double * Qwrite = (double *) malloc( myNtot*Ndoub*sizeof(double) );

   int index = 0;
   for( k=kmin ; k<kmax ; ++k ){
      for( j=jmin ; j<jmax ; ++j ){
         int jk = (k-kmin)*jSize + (j-jmin);
         Index[jk] = index;
         Size[jk] = Np[j+Nr*k];
 
         double phi0 = M_PI;
         int Id = 0;
         int i;
         for( i=0 ; i<Np[j+Nr*k] ; ++i ){
            struct cell * c = &(theCells[j+Nr*k][i]);
            Cell2Doub( c , Qwrite+index*Ndoub , 1 );
            double phi = c->piph-.5*c->dphi;
            if( cos(phi0) < cos(phi) ){ phi0 = phi; Id = index; }
            ++index;
         }
         Id_phi0[jk] = Id;
      }
   }

   int runningTot = 0;
   for( nrk=0 ; nrk < size ; ++nrk ){
      int thisTot = myNtot;
      MPI_Bcast( &thisTot , 1 , MPI_INT , nrk , theDomain->theComm );
      if( rank > nrk ) runningTot += thisTot;
   }

   int jk;
   for( jk=0 ; jk<jSize*kSize ; ++jk ){
      Index[jk] += runningTot;
      Id_phi0[jk] += runningTot;
   }

   for( nrk=0 ; nrk < size ; ++nrk ){
      if( rank==nrk ){      
         //Write Cell Data
         int start2[2]    = {runningTot,0};
         int loc_size2[2] = {myNtot,Ndoub};
         int glo_size2[2] = {Ntot,Ndoub};
         writePatch( filename , "Data" , "Cells" , Qwrite , H5T_NATIVE_DOUBLE , 2 , start2 , loc_size2 , glo_size2 );
         //Write Indices and Sizes for each radial track
         start2[0] = k0;
         start2[1] = j0;
         loc_size2[0] = kSize;
         loc_size2[1] = jSize;
         glo_size2[0] = Nz_Tot;
         glo_size2[1] = Nr_Tot;
         writePatch( filename , "Grid" , "Index"   , Index   , H5T_NATIVE_INT , 2 , start2 , loc_size2 , glo_size2 );
         writePatch( filename , "Grid" , "Np"      , Size    , H5T_NATIVE_INT , 2 , start2 , loc_size2 , glo_size2 );
         writePatch( filename , "Grid" , "Id_phi0" , Id_phi0 , H5T_NATIVE_INT , 2 , start2 , loc_size2 , glo_size2 );
         //Write 1D Radial Data
         if( dim_rank[1] == 0 ){
            int offset = Ng;
            if( dim_rank[0] == 0 ) offset = 0;
            int start1[1]    = {j0};
            int loc_size1[1] = {jSize};
            if( dim_rank[0] == dim_size[0]-1 ) loc_size1[0]++;
            int glo_size1[1] = {Nr_Tot+1};
            writePatch( filename , "Grid" , "r_jph" , r_jph-1+offset , H5T_NATIVE_DOUBLE , 1 , start1 , loc_size1 , glo_size1 );
            int start2[2] = {j0,0};
            int loc_size2[2] = {jSize,Ntools};
            int glo_size2[2] = {Nr_Tot,Ntools};
            double * Q = theDomain->theTools.Qr;
            writePatch( filename , "Data" , "Radial_Diagnostics" , Q + offset*Ntools , H5T_NATIVE_DOUBLE , 2 , start2 , loc_size2 , glo_size2 );
         }
         //Write 1D Vertical Data
         if( dim_rank[0] == 0 ){
            int offset = Ng;
            if( dim_rank[1] == 0 ) offset = 0;
            if( Z_Periodic && dim_rank[1] == 0 ) offset += Ng;
            int start1[1]    = {k0};
            int loc_size1[1] = {kSize};
            if( dim_rank[1] == dim_size[1]-1 ) loc_size1[0]++;
            int glo_size1[1] = {Nz_Tot+1};
            writePatch( filename , "Grid" , "z_kph" , z_kph-1+offset , H5T_NATIVE_DOUBLE , 1 , start1 , loc_size1 , glo_size1 );
         }
      }
      MPI_Barrier( theDomain->theComm );
   }
   zero_diagnostics( theDomain );

   free(Index);
   free(Size);
   free(Id_phi0);
   free(Qwrite);
   MPI_Barrier(theDomain->theComm);
}



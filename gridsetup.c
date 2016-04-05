
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "paul.h"

int getN0( int drank , int dsize , int dnum ){
   int N0 = (dnum*drank)/dsize;
   return(N0);
}

void setupGrid( struct domain * theDomain ){

   int Ng = NUM_G;
   theDomain->Ng = Ng;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   int Num_R = theDomain->theParList.Num_R;
   int Num_Z = theDomain->theParList.Num_Z;
   int LogZoning = theDomain->theParList.LogZoning;
   int Z_Periodic = theDomain->theParList.Z_Periodic;
   double aspect = theDomain->theParList.aspect;

   double Rmin = theDomain->theParList.rmin;
   double Rmax = theDomain->theParList.rmax;
   double Zmin = theDomain->theParList.zmin;
   double Zmax = theDomain->theParList.zmax;
   double Pmax = theDomain->theParList.phimax;

   int N0r = getN0( dim_rank[0]   , dim_size[0] , Num_R );
   int N1r = getN0( dim_rank[0]+1 , dim_size[0] , Num_R );
   if( dim_rank[0] != 0 ) N0r -= Ng;
   if( dim_rank[0] != dim_size[0]-1 ) N1r += Ng;
   int Nr = N1r-N0r;

   int N0z = getN0( dim_rank[1]   , dim_size[1] , Num_Z );
   int N1z = getN0( dim_rank[1]+1 , dim_size[1] , Num_Z );
   if( Num_Z > 1 ){
      if( dim_rank[1] != 0 || Z_Periodic ) N0z -= Ng;
      if( dim_rank[1] != dim_size[1]-1 || Z_Periodic ) N1z += Ng;
   }
   int Nz = N1z-N0z;

   theDomain->Nr = Nr;
   theDomain->Nz = Nz;
   printf("Rank = %d, Nr = %d, Nz = %d\n",theDomain->rank,Nr,Nz);

   theDomain->Np    = (int *)    malloc( Nr*Nz*sizeof(int) );
   theDomain->r_jph = (double *) malloc( (Nr+1)*sizeof(double) );
   theDomain->z_kph = (double *) malloc( (Nz+1)*sizeof(double) );

   ++(theDomain->r_jph);
   ++(theDomain->z_kph);

   int j,k;

   double dx = 1./(double)Num_R;
   double x0 = (double)N0r/(double)Num_R;
   double R0 = theDomain->theParList.LogRadius;
   for( j=-1 ; j<Nr ; ++j ){
      double x = x0 + ((double)j+1.)*dx;
      if( LogZoning == 0 ){
         theDomain->r_jph[j] = Rmin + x*(Rmax-Rmin);
      }else if( LogZoning == 1 ){
         theDomain->r_jph[j] = Rmin*pow(Rmax/Rmin,x);
      }else{
         theDomain->r_jph[j] = R0*pow(Rmax/R0,x) + Rmin-R0 + (R0-Rmin)*x;
      }
   }
   double dz = (Zmax-Zmin)/(double)Num_Z;
   double z0 = Zmin + (double)N0z*dz;
   for( k=-1 ; k<Nz ; ++k ){
      theDomain->z_kph[k] = z0 + ((double)k+1.)*dz;
   }

   for( j=0 ; j<Nr ; ++j ){
      double rp = theDomain->r_jph[j];
      double rm = theDomain->r_jph[j-1];
      double dr = rp-rm;
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         double dp = dr/rp*aspect;
         int Np = (int)(Pmax/dp);
         if( Np<4 ) Np=4;
//         if( Np<40 ) Np=40;
         theDomain->Np[jk] = Np;
      }
   }
}



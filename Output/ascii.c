
#include "../paul.h"

double get_dV( double * , double * );
void prim2cons( double * , double * , double , double );
void cons2prim( double * , double * , double , double );

void output( struct domain * theDomain , char * filestart ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Ng = theDomain->Ng;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int rank = theDomain->rank;
   int size = theDomain->size;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   int Z_Periodic = theDomain->theParList.Z_Periodic;

   char filename[256];
   sprintf(filename,"%s.dat",filestart);
/*
   char filename_1d[256];
   sprintf(filename_1d,"%s_1d.dat",filestart);
*/
   if( rank==0 ){
      FILE * pFile = fopen( filename , "w" );
      fclose(pFile);
   }
   MPI_Barrier( theDomain->theComm );

   int j_min = 0;
   int j_max = Nr;
   int k_min = 0;
   int k_max = Nz;

   if( dim_rank[0] != 0 ) j_min = Ng;
   if( dim_rank[0] != dim_size[0]-1 ) j_max = Nr-Ng;
   if( dim_rank[1] != 0 || Z_Periodic ) k_min = Ng;
   if( dim_rank[1] != dim_size[1]-1 || Z_Periodic ) k_max = Nz-Ng;

   int rk;
   for( rk=0 ; rk<size ; ++rk ){
      if( rank==rk ){
         FILE * pFile = fopen( filename , "a" );
         int i,j,k,q;
         for( j=j_min ; j<j_max ; ++j ){
            double r  = .5*(r_jph[j]+r_jph[j-1]);
            double dr = r_jph[j]-r_jph[j-1];
            for( k=k_min ; k<k_max ; ++k ){
               double z = .5*(z_kph[k]+z_kph[k-1]);
               double dz = z_kph[k]-z_kph[k-1];
               int jk = j+Nr*k;
               for( i=0 ; i<Np[jk] ; ++i ){
                  struct cell * c = &(theCells[jk][i]);
                  double phi  = c->piph - .5*c->dphi;
                  double dphi = c->dphi;
                  fprintf(pFile,"%e %e %e ",r,phi,z);
                  fprintf(pFile,"%e %e %e ",dr,dphi,dz);
                  for( q=0 ; q<NUM_Q ; ++q ){
                     fprintf(pFile,"%e ",c->prim[q]);
                  }
                  fprintf(pFile,"\n");
               }
            }
         }
         fclose( pFile );
      }
      MPI_Barrier( theDomain->theComm );
   }
/*
   double Rmin = theCells[0][0].riph;
   double Rmax = theCells[0][Nr[0]-1].riph;
 
   int NUM_R = theDomain->theParList.Num_R;
   int i,j,k,q;
   double cons_1d_avg[NUM_R*NUM_Q];
   double prim_1d_avg[NUM_R*NUM_Q];
   for( i=0 ; i<NUM_R*NUM_Q ; ++i ){
      cons_1d_avg[i] = 0.0;
      prim_1d_avg[i] = 0.0;
   }
   
   for( j=j_min ; j<j_max ; ++j ){
      double thp = t_jph[j];
      double thm = t_jph[j-1];
      for( k=k_min ; k<k_max ; ++k ){
         double php = p_kph[k];
         double phm = p_kph[k-1];
         int jk = j+Nt*k;
         for( i=1 ; i<Nr[jk] ; ++i ){
            double rp = theCells[jk][i].riph;
            double rm = rp-theCells[jk][i].dr;
            double ip = (NUM_R)*(rp-Rmin)/(Rmax-Rmin);//(NUM_R)*log(rp/Rmin)/log(Rmax/Rmin);
            double im = (NUM_R)*(rm-Rmin)/(Rmax-Rmin);//(NUM_R)*log(rm/Rmin)/log(Rmax/Rmin);
            double xp[3] = {rp,thp,php};
            double xm[3] = {rm,thm,phm};
            double dV = get_dV( xp , xm );
          
            double delta_i = ip-im; 
            int iip = (int) ip;
            int iim = ((int) im);
            double dip = ip-(double)iip;
            double dim = (double)(iim+1)-im;

            if( iip == iim ) dip = delta_i;
        
            int iCons;
            for( iCons = iim ; iCons <= iip ; ++iCons ){
               double frac = 1./delta_i;
               if( iCons == iim ) frac = dim/delta_i;
               if( iCons == iip ) frac = dip/delta_i;

               int i2 = iCons;
               if( i2 < 0 ) i2 = 0;
               if( i2 > NUM_R-1 ) i2 = NUM_R-1;

               for( q=0 ; q<NUM_Q ; ++q ){
                  cons_1d_avg[i2*NUM_Q + q] += frac*theCells[jk][i].cons[q];
                  prim_1d_avg[i2*NUM_Q + q] += frac*theCells[jk][i].prim[q]*dV;
               }
            }
         }
      }
   }

   MPI_Allreduce( MPI_IN_PLACE , cons_1d_avg , NUM_R*NUM_Q , MPI_DOUBLE , MPI_SUM , theDomain->theComm );
   MPI_Allreduce( MPI_IN_PLACE , prim_1d_avg , NUM_R*NUM_Q , MPI_DOUBLE , MPI_SUM , theDomain->theComm );

   double THETA_MIN = theDomain->theParList.thmin;
   double THETA_MAX = theDomain->theParList.thmax;
   double PHI_MAX   = theDomain->theParList.phimax;

   if( rank==0 ){
      FILE * pFile_1d = fopen(filename_1d,"w");
      for( i=0 ; i<NUM_R ; ++i ){
         double P_out[NUM_Q];
         double xxm = (double)i/(double)NUM_R;
         double xxp = (double)(i+1)/(double)NUM_R;
         double rm = xxm*(Rmax-Rmin)+Rmin;//pow(Rmax/Rmin,xxm)*Rmin;
         double rp = xxp*(Rmax-Rmin)+Rmin;//pow(Rmax/Rmin,xxp)*Rmin;
         double xp[3] = {rp,THETA_MAX,PHI_MAX};
         double xm[3] = {rm,THETA_MIN,0.0};
         double r = (2./3.)*(rp*rp*rp-rm*rm*rm)/(rp*rp-rm*rm);
         double dV = get_dV( xp , xm );
         cons2prim( &(cons_1d_avg[i*NUM_Q]) , P_out , r , dV );
         fprintf(pFile_1d,"%e ",r);
         for( q=0 ; q<NUM_Q ; ++q ){
            fprintf(pFile_1d,"%e ",P_out[q]);
            fprintf(pFile_1d,"%e ",prim_1d_avg[i*NUM_Q+q]/dV);
         }
         fprintf(pFile_1d,"\n");
      }
      fclose( pFile_1d );
   }
*/
}

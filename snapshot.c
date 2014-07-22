
#include "paul.h"

double get_dV( double * , double * );
void prim2cons( double * , double * , double , double );
void cons2prim( double * , double * , double , double );

void snapshot( struct domain * theDomain , char * filestart ){

   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   int Ng = theDomain->Ng;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   int rank = theDomain->rank;
   int size = theDomain->size;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;

   if(rank==0) printf("Generating Snapshot...\n");

   int j_min = 0; 
   int j_max = Nt;
   int k_min = 0; 
   int k_max = Np;

   if( dim_rank[0] != 0 ) j_min = NUM_G;
   if( dim_rank[0] != dim_size[0]-1 ) j_max = Nt-NUM_G;
   if( dim_rank[1] != 0 ) k_min = NUM_G;
   if( dim_rank[1] != dim_size[1]-1 ) k_max = Np-NUM_G;

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
            struct cell * c = theCells[jk]+i;
            double rp = c->riph;
            double rm = rp - c->dr;
            double ip = (NUM_R)*(rp-Rmin)/(Rmax-Rmin);
            double im = (NUM_R)*(rm-Rmin)/(Rmax-Rmin);
            double xp[3] = {rp,thp,php};
            double xm[3] = {rm,thm,phm};
            double dV = get_dV( xp , xm );
          
            double delta_i = ip-im; 
            int iip = (int) ip;
            int iim = (int) im;
            if(iim<0) iim=0;
            if(iip<0) iip=0;
            if(iim>NUM_R-1) iim=NUM_R-1;
            if(iip>NUM_R-1) iip=NUM_R-1;
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
                  cons_1d_avg[i2*NUM_Q + q] += frac*c->cons[q];
                  prim_1d_avg[i2*NUM_Q + q] += frac*c->prim[q]*dV;
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

   char fname_rad[256];
   strcpy( fname_rad , filestart );
   strcat( fname_rad , "_radial.dat" );
   char fname_th0[256];
   strcpy( fname_th0 , filestart );
   strcat( fname_th0 , "_rad0.dat" );
   char fname_th1[256];
   strcpy( fname_th1 , filestart );
   strcat( fname_th1 , "_rad1.dat" );

   if( rank==0 ){
      FILE * pFile_1d = fopen(fname_rad,"w");
      for( i=0 ; i<NUM_R ; ++i ){
         double P_out[NUM_Q];
         double xxm = (double)i/(double)NUM_R;
         double xxp = (double)(i+1)/(double)NUM_R;
         double rm = xxm*(Rmax-Rmin)+Rmin;
         double rp = xxp*(Rmax-Rmin)+Rmin;
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

      pFile_1d = fopen(fname_th0,"w");
      for( i=0 ; i<Nr[0] ; ++i ){
         struct cell * c = theCells[0]+i;
         fprintf(pFile_1d,"%e ",c->riph-.5*c->dr);
         for( q=0 ; q<NUM_Q ; ++q ) fprintf(pFile_1d,"%e ",c->prim[q]);
         fprintf(pFile_1d,"\n");
      }
      fclose( pFile_1d );
   }

   double thetaSlice = 0.05;
   double thP = t_jph[Nt-1];
   double thM = t_jph[-1];
   double Dth = thP-thM;
   double Mth = .5*(thP+thM);
   double delta = fabs(thetaSlice-Mth)/Dth;

   struct { double value ; int index ; } minbuf;
   minbuf.value = delta;
   minbuf.index = rank;
   MPI_Allreduce( MPI_IN_PLACE , &minbuf , 1 , MPI_DOUBLE_INT , MPI_MINLOC , theDomain->theComm );

   int thRank = minbuf.index;
   if( rank==thRank ){
      j=0;
      while( j<Nt && t_jph[j-1]<thetaSlice ) ++j;
      if( j>0 ) --j;
      FILE * pFile_1d = fopen(fname_th1,"w");
      for( i=0 ; i<Nr[j] ; ++i ){
         struct cell * c = theCells[j]+i;
         fprintf(pFile_1d,"%e ",c->riph-.5*c->dr);
         for( q=0 ; q<NUM_Q ; ++q ) fprintf(pFile_1d,"%e ",c->prim[q]);
         fprintf(pFile_1d,"\n");         
      }
      fclose( pFile_1d );
   }

//////////////////////////Stupid Ebins Bullshit
   double maxgam = 0.0;
   for( j=0 ; j<Nt ; ++j ){
      for( i=0 ; i<Nr[j] ; ++i ){
         double ur = theCells[j][i].prim[UU1];
         double up = theCells[j][i].prim[UU2];
         if( maxgam < ur ) maxgam = sqrt(ur*ur+up*up);
      }
   }

   MPI_Allreduce( MPI_IN_PLACE , &maxgam , 1 , MPI_DOUBLE , MPI_MAX , theDomain->theComm );
   int Ngam = 1000;
   double Ebin[Ngam];
   double Mtot = 0.0;
   int n;
   for( n=0 ; n<Ngam ; ++n ) Ebin[n] = 0.0;
   int jmin = Ng;
   int jmax = Nt-Ng;
   if( dim_rank[0]==0 ) jmin = 0;
   if( dim_rank[0]==dim_size[0]-1 ) jmax = Nt;
   
   for( j=jmin ; j<jmax ; ++j ){
      for( i=0 ; i<Nr[j] ; ++i ){
         double ur = theCells[j][i].prim[UU1];
         double up = theCells[j][i].prim[UU2];
         double gam = sqrt(ur*ur+up*up);
         double ndoub = (double)Ngam*gam/maxgam;
         int nn = (int) ndoub;
         if( nn > Ngam ) nn = Ngam;
         if( nn < 0 ) nn = 0;
         Ebin[nn] += theCells[j][i].cons[TAU];
         Mtot += theCells[j][i].cons[DEN];
      }
   }
   MPI_Allreduce( MPI_IN_PLACE , Ebin , Ngam , MPI_DOUBLE , MPI_SUM , theDomain->theComm );
   MPI_Allreduce( MPI_IN_PLACE , &Mtot ,  1  , MPI_DOUBLE , MPI_SUM , theDomain->theComm );

   char fname_ebins[256];
   strcpy( fname_ebins , filestart );
   strcat( fname_ebins , "_ebins.dat" );

   if( rank==0 ){
      FILE * gamFile = fopen( fname_ebins ,"w");
      for( n=Ngam-2 ; n>=0 ; --n ){
         Ebin[n] += Ebin[n+1];
      }
      fprintf(gamFile,"# E/M = %e\n",Ebin[0]/Mtot);
      for( n=0 ; n<Ngam ; ++n ){
         double gam = ((double)n/(double)Ngam)*maxgam;
         fprintf(gamFile,"%e %e\n",gam,Ebin[n]/Ebin[0]);
      }
      fclose(gamFile);
   }

///////////////////////////////////////////////

   char fname_theta[256];
   strcpy( fname_theta , filestart );
   strcat( fname_theta , "_theta.dat" );
   FILE * thFile;

   if(rank==0){
      thFile = fopen( fname_theta ,"w");
      fclose(thFile);
   }
   int nrk;
   double Eth[Nt];
   double XEth[Nt];
   double g_max[Nt];
   double r_max[Nt];

   for( j=j_min ; j<j_max ; ++j ){
      double maxgam = 0.0;
      double maxr   = 0.0;
      Eth[j]  = 0.0;
      XEth[j] = 0.0;
      for( i=0 ; i<Nr[j] ; ++i ){
         struct cell * c = theCells[j]+i;
         Eth[j]  += c->cons[TAU];
         if(NUM_N>0) XEth[j] += c->cons[TAU]*c->prim[NUM_C];
         double ur = c->prim[UU1];
         double up = c->prim[UU2];
         double gam = sqrt(ur*ur+up*up);
         if( maxgam < gam ){
            maxgam = gam;
         }
         double r = theCells[j][i].riph;
         if( maxr<r && gam>.1*maxgam ){
            maxr = r;
         }
      }
      g_max[j] = maxgam;
      r_max[j] = maxr;
   }
   for( nrk=0 ; nrk<size ; ++nrk ){
   if( rank==nrk ){
      FILE * thFile = fopen( fname_theta ,"a");
      double xp[3] = {1.0,0.0,p_kph[0]};
      double xm[3] = {0.0,0.0,p_kph[-1]};
      for( j=jmin ; j<j_max ; ++j ){
         xp[1] = t_jph[j];
         xm[1] = t_jph[j-1];
         double dOmega = get_dV(xp,xm);
         fprintf(thFile,"%e ",.5*(t_jph[j]+t_jph[j-1]));
         fprintf(thFile,"%e %e %e %e",r_max[j],g_max[j],Eth[j]/dOmega,XEth[j]/dOmega);
         fprintf(thFile,"\n");
      }
      fclose(thFile);
   }
   MPI_Barrier(theDomain->theComm);
   }

}



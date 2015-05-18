#include "paul.h"
#include <string.h>

double get_dA( double * , double * , int );
double get_dV( double * , double * );
double get_moment_arm( double * , double * );

void clean_pi( struct domain * theDomain ){
   
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Npl = theDomain->Npl;
   double phi_max = theDomain->phi_max;

   int i,jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         double phi = c->piph;
         while( phi > phi_max ){ phi -= phi_max; c->RKpiph -= phi_max; }
         while( phi < 0.0     ){ phi += phi_max; c->RKpiph += phi_max; }
         c->piph = phi; 
      }    
   }
   int p;
   for( p=0 ; p<Npl ; ++p ){
      struct planet * pl = theDomain->thePlanets + p;
      double phi = pl->phi;
      while( phi > phi_max ){ phi -= phi_max; pl->RK_phi -= phi_max; }
      while( phi < 0.0     ){ phi += phi_max; pl->RK_phi += phi_max; }
      pl->phi = phi; 
   }

}

double mindt( double * , double , double * , double * );

double getmindt( struct domain * theDomain ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   double dt = 1e100;
   int i,j,k;
   for( j=1 ; j<Nr-1 ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = phip - c->dphi;
            double xp[3] = {r_jph[j  ] , phip , z_kph[k  ]};
            double xm[3] = {r_jph[j-1] , phim , z_kph[k-1]};
            int im = i-1;
            if( i==0 ) im = Np[jk]-1; 
            double wm = theCells[jk][im].wiph;
            double wp = c->wiph;
            double w = .5*(wm+wp);
            double dt_temp = mindt( c->prim , w , xp , xm );
            if( dt > dt_temp ) dt = dt_temp;
         }
      }
   }
   dt *= theDomain->theParList.CFL; 
   MPI_Allreduce( MPI_IN_PLACE , &dt , 1 , MPI_DOUBLE , MPI_MIN , theDomain->theComm );

   return( dt );
}

void initial( double * , double * );
void prim2cons( double * , double * , double , double );
void cons2prim( double * , double * , double , double );
void restart( struct domain * );
/*
void clear_w( struct domain * theDomain ){
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int i,jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         theDomain->theCells[jk][i].wiph = 0.0;
      }
   }
}*/

double get_omega( double * );
double mesh_om( double );

void set_wcell( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int mesh_motion = theDomain->theParList.Mesh_Motion;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;

   int i,j,k;
   for( j=0 ; j<Nr ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * cL = &(theCells[jk][i ]);  
            double w = 0.0;
            if( mesh_motion ){
               int ip = (i+1)%Np[jk];
               struct cell * cR = &(theCells[jk][ip]);
               double wL = get_omega( cL->prim );
               double wR = get_omega( cR->prim );
               w = .5*(wL + wR); 
            }
            cL->wiph = w;
         }
      }    
   }
   if( mesh_motion == 3 ){
      for( j=0 ; j<Nr ; ++j ){
         for( k=0 ; k<Nz ; ++k ){
            int jk = j+Nr*k;
            double w = 0.0;
            for( i=0 ; i<Np[jk] ; ++i ){
               w += theCells[jk][i].wiph; 
            }    
            w /= (double)Np[jk];
            for( i=0 ; i<Np[jk] ; ++i ){
               theCells[jk][i].wiph = w; 
            }    
         }    
      } 
   }
   if( mesh_motion == 4 ){
      for( j=0 ; j<Nr ; ++j ){
         double r = .5*(r_jph[j]+r_jph[j-1]);
         for( k=0 ; k<Nz ; ++k ){
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i ){
               theCells[jk][i].wiph = r*mesh_om(r); 
            }    
         }    
      } 
   }

}

void initial( double * , double * );
void clear_cell( struct cell * );

void regrid( struct domain * theDomain ){
/*
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   double jthresh = 0.1;
   int i,j,k;
   for( j=0 ; j<Nt ; ++j ){
   for( k=0 ; k<Np ; ++k ){
      int jk = j+Nt*k;
      for( i=2 ; i<Nr[jk]-1 ; ++i ){
         double jump = log( theCells[jk][i-1].prim[RHO]/theCells[jk][i+1].prim[RHO] );
         if( fabs(jump) > jthresh ){
            struct cell * c = theCells[jk]+i;
            double rp = c->riph;
            double rm = (c-1)->riph;
            int Nsplit = (int)(fabs(jump/jthresh));
            if(Nsplit>5) Nsplit=5;
            if(Nsplit>3) printf("r=%.2e jump=%.2e split = %d (No Cause For Alarm)\n",theCells[jk][i].riph,jump,Nsplit);
            int blocksize = (Nr[jk]-1) - i;
            Nr[jk] += Nsplit;
            theCells[jk] = (struct cell *) realloc( theCells[jk] , Nr[jk]*sizeof(struct cell) );
            memmove( theCells[jk]+i+1+Nsplit , theCells[jk]+i+1 , blocksize*sizeof(struct cell) );
            int l;
            for( l=0 ; l<Nsplit+1 ; ++l ){
               int m = l+i;
               struct cell * cnew = theCells[jk]+m;
               clear_cell( cnew );
               double rlp = rm + (rp-rm)*( (double)(l+1)/(double)(Nsplit+1) );
               double rlm = rm + (rp-rm)*( (double)(l  )/(double)(Nsplit+1) );
               double xp[3] = {rlp,t_jph[j]  ,p_kph[k]  };
               double xm[3] = {rlm,t_jph[j-1],p_kph[k-1]};
               double x[3];
               int d;
               for( d=0 ; d<3 ; ++d ) x[d] = .5*(xp[d]+xm[d]);
               initial( cnew->prim , x );
               cnew->riph = rlp;
               cnew->dr = rlp-rlm;
               double dV = get_dV(xp,xm);
               double rr = (2./3.)*(rlp*rlp*rlp-rlm*rlm*rlm)/(rlp*rlp-rlm*rlm);
               prim2cons( cnew->prim , cnew->cons , rr , dV );
            }
            i += Nsplit;
         }
      }
   }
   }
*/
}


void adjust_RK_cons( struct domain * theDomain , double RK ){
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;

   int i,jk,q;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         for( q=0 ; q<NUM_Q ; ++q ){
            c->cons[q] = (1.-RK)*c->cons[q] + RK*c->RKcons[q];
         }
      }
   }
}

void planet_RK_adjust( struct planet * , double );

void adjust_RK_planets( struct domain * theDomain , double RK ){
   int Npl = theDomain->Npl;
   int p;
   for( p=0 ; p<Npl ; ++p ){
      planet_RK_adjust( theDomain->thePlanets+p , RK );
   }
}

void move_cells( struct domain * theDomain , double RK , double dt){
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int i,jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         c->piph = (1.-RK)*c->piph + RK*c->RKpiph + c->wiph*dt;
      }
   }
}

double get_dp( double , double );

void calc_dp( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;

   int i,jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         int im = i-1;
         if( i == 0 ) im = Np[jk]-1;
         double phim = theCells[jk][im].piph;
         double phip = theCells[jk][i ].piph;
         double dphi = get_dp(phip,phim);
         theCells[jk][i].dphi = dphi;
      }
   }
}

void calc_prim( struct domain * theDomain ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   int i,j,k;
   for( j=0 ; j<Nr ; ++j ){
      double rm = r_jph[j-1];
      double rp = r_jph[j];
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = phip-c->dphi;
            double xp[3] = {rp,phip,z_kph[k]  };
            double xm[3] = {rm,phim,z_kph[k-1]};
            double r = get_moment_arm( xp , xm );
            double dV = get_dV( xp , xm );
            cons2prim( c->cons , c->prim , r , dV );
         }
      }
   }
}

void plm_phi( struct domain * );
void riemann_phi( struct cell * , struct cell * , double , double );

void phi_flux( struct domain * theDomain , double dt ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int i,j,k;
   plm_phi( theDomain );
   for( j=0 ; j<Nr ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            int ip = (i+1)%Np[jk];
            struct cell * cp = theCells[jk];
            double phi = cp[i].piph;
            double xp[3] = {r_jph[j]  ,phi,z_kph[k]  };
            double xm[3] = {r_jph[j-1],phi,z_kph[k-1]};
            double r = get_moment_arm(xp,xm);
            double dA = get_dA(xp,xm,0); 
            riemann_phi( &(cp[i]) , &(cp[ip]) , r , dA*dt );
         }
      }
   }

}

void buildfaces( struct domain * , struct face * , int * , int , int );
void plm_trans( struct domain * , struct face * , int , int );
void riemann_trans( struct face * , double , int );

void trans_flux( struct domain * theDomain , struct face * theFaces , int Nf , double dt , int dim ){

   plm_trans( theDomain , theFaces , Nf , dim );
   int f;
   for( f=0 ; f<Nf ; ++f ){
      riemann_trans( &(theFaces[f]) , dt , dim );
   }

}

int get_num_rzFaces( int , int , int );

void setup_faces( struct domain * theDomain , struct face ** theFaces , int * nn , int dim ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int NN = get_num_rzFaces( Nr , Nz , dim );

   buildfaces( theDomain , NULL      , nn , dim , 0 );
   int Nf = nn[NN];
   *theFaces = (struct face *) malloc( Nf*sizeof(struct face) );
   buildfaces( theDomain , *theFaces , nn , dim , 1 );

}

void source( double * , double * , double * , double * , double );
void planet_src( struct planet * , double * , double * , double * , double * , double );

void add_source( struct domain * theDomain , double dt ){

   struct cell ** theCells = theDomain->theCells;
   struct planet * thePlanets = theDomain->thePlanets;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Npl = theDomain->Npl;

   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   int i,j,k,p;
   for( j=0 ; j<Nr ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = c->piph - c->dphi;
            double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            double dV = get_dV(xp,xm);
            source( c->prim , c->cons , xp , xm , dV*dt );
            for( p=0 ; p<Npl ; ++p ){
               planet_src( thePlanets+p , c->prim , c->cons , xp , xm , dV*dt );
            }
         }    
      }    
   }   

}

void longandshort( struct domain * theDomain , double * L , double * S , int * iL , int * iS , struct cell * sweep , int j , int k ){

   int Nr = theDomain->Nr;
   int jk = j + Nr*k;
   int Np = theDomain->Np[jk];
   double * r_jph = theDomain->r_jph;
   double dr = r_jph[j]-r_jph[j-1];
   double r  = r_jph[j];
   double Long = 0.0;
   double Short = 0.0;
   int iLong = 0;
   int iShort = 0;
   int i;
   for( i=0 ; i<Np ; ++i ){
      double dx = dr;
      double dy = r*sweep[i].dphi;
      double l = dy/dx;
      double s = dx/dy;
      if( Long  < l ){ Long  = l; iLong  = i; } 
      if( Short < s ){ Short = s; iShort = i; } 
   }
   *L  = Long;
   *iL = iLong;
   *S  = Short;
   *iS = iShort;

}

void AMRsweep( struct domain * theDomain , struct cell ** swptr , int jk ){

   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int Nr = theDomain->Nr;
   int * Np = theDomain->Np;
   int j = jk%Nr;
   int k = jk/Nr;

   double MaxShort = theDomain->theParList.MaxShort;
   double MaxLong  = theDomain->theParList.MaxLong;

   struct cell * sweep = *swptr;
   double L,S;
   int iL=0;
   int iS=0;
   longandshort( theDomain , &L , &S , &iL , &iS , sweep , j , k );
   //printf("Long = %e, Short = %e\n",L,S);

   if( S>MaxShort ){
      int iSp = (iS+1)%Np[jk];

      //Possibly shift iS backwards by 1
      int iSm = iS-1;
      if( iSm == -1 ) iSm = Np[jk]-1;
      double dpL = sweep[iSm].dphi;
      double dpR = sweep[iSp].dphi;
      if( dpL < dpR ){
         --iS;
         --iSm;
         --iSp;
         if( iS  == -1 ) iS  = Np[jk]-1;
         if( iSm == -1 ) iSm = Np[jk]-1;
         if( iSp == -1 ) iSp = Np[jk]-1;
      }

      //Remove Zone at iS+1
      sweep[iS].dphi  += sweep[iSp].dphi;
      sweep[iS].piph   = sweep[iSp].piph;
      sweep[iS].RKpiph = sweep[iSp].RKpiph;
      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         sweep[iS].cons[q]   += sweep[iSp].cons[q];
         sweep[iS].RKcons[q] += sweep[iSp].RKcons[q];
      }
      double phip = sweep[iS].piph;
      double phim = phip - sweep[iS].dphi;
      double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
      double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
      double r  = get_moment_arm( xp , xm );
      double dV = get_dV( xp , xm );
      cons2prim( sweep[iS].cons , sweep[iS].prim , r , dV );
      //Shift Memory
      int blocksize = Np[jk]-iSp-1;
      if( iSp != Np[jk]-1 ) memmove( sweep+iSp , sweep+iSp+1 , blocksize*sizeof(struct cell) );
      Np[jk] -= 1;
      *swptr = (struct cell *) realloc( sweep , Np[jk]*sizeof(struct cell) );
      sweep = *swptr;
      if( iS < iL ) iL--;

   }

   if( L>MaxLong ){
      Np[jk] += 1;
      *swptr = (struct cell *) realloc( sweep , Np[jk]*sizeof(struct cell) );
      sweep = *swptr;
      int blocksize = Np[jk]-iL-1;
      memmove( sweep+iL+1 , sweep+iL , blocksize*sizeof(struct cell) );

      double dphi = sweep[iL].dphi;
      double phip = sweep[iL].piph;
      double phim = phip - dphi;
      double phi0 = .5*(phip+phim);

      sweep[iL].piph   = phi0;
      sweep[iL].RKpiph = phi0;
      sweep[iL].dphi   = .5*dphi;
      sweep[iL+1].dphi = .5*dphi;

      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         sweep[iL].cons[q]     *= .5;
         sweep[iL].RKcons[q]   *= .5;
         sweep[iL+1].cons[q]   *= .5;
         sweep[iL+1].RKcons[q] *= .5;
      }

      double xp[3] = {r_jph[j]  ,phi0,z_kph[k]  };
      double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
      double dV = get_dV( xp , xm );
      double r  = get_moment_arm( xp , xm );
      cons2prim( sweep[iL].cons , sweep[iL].prim , r , dV );

      xp[1] = phip;
      xm[1] = phi0;
      dV = get_dV( xp , xm );
      r  = get_moment_arm( xp , xm );
      cons2prim( sweep[iL+1].cons , sweep[iL+1].prim , r , dV );

   }
}

void AMR( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      AMRsweep( theDomain , theCells+jk , jk );
   }
}


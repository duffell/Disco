
#include "paul.h"

void initial( double * , double * ); 
double get_moment_arm( double * , double * );
double get_dA( double * , double * , int );
double get_dV( double * , double * );
void setup_faces( struct domain * , int );
int get_num_rzFaces( int , int , int );
void B_faces_to_cells( struct domain * , int );
 
void set_B_fields( struct domain * theDomain ){

   int i,j,k;
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   for( j=0 ; j<Nr ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = phip-c->dphi;
            double xp[3] = { r_jph[j]   , phip , z_kph[k]   };
            double xm[3] = { r_jph[j-1] , phim , z_kph[k-1] };
            double r = get_moment_arm( xp , xm );
            double x[3] = { r , phip , .5*(z_kph[k]+z_kph[k-1])};
            double prim[NUM_Q];
            initial( prim , x ); 
            double dA = get_dA( xp , xm , 0 );
            double Phi = 0.0;
            if( NUM_Q > BPP ) Phi = prim[BPP]*dA;
            c->Phi[0] = Phi;
         }    
      }    
   }

   int NRZ1 = theDomain->N_ftracks_r;
   setup_faces( theDomain , 1 ); 

   int n;
   for( n=0 ; n<theDomain->fIndex_r[NRZ1] ; ++n ){
      struct face * f = theDomain->theFaces_1 + n;
      double prim[NUM_Q];
      initial( prim , f->cm );
      double Phi = 0.0;
      if( NUM_Q > BRR ) Phi = prim[BRR]*f->dA;
      if( f->LRtype == 0 ){
         f->L->Phi[2] = Phi;
      }else{
         f->R->Phi[1] = Phi;
      }
   }

   if( NUM_FACES==5 && theDomain->Nz > 1 ){
      int NRZ2 = theDomain->N_ftracks_z;
      setup_faces( theDomain , 2 );
      for( n=0 ; n<theDomain->fIndex_z[NRZ2] ; ++n ){
         struct face * f = theDomain->theFaces_2 + n; 
         double prim[NUM_Q];
         initial( prim , f->cm );
         double Phi = 0.0; 
         if( NUM_Q > BZZ ) Phi = prim[BZZ]*f->dA;

         if( f->LRtype == 0 ){ 
            f->L->Phi[4] = Phi; 
         }else{
            f->R->Phi[3] = Phi; 
         }    
      }
   }

   B_faces_to_cells( theDomain , 0 );
   B_faces_to_cells( theDomain , 1 );
 
   free( theDomain->theFaces_1 );
   if( theDomain->theFaces_2 ) free( theDomain->theFaces_2 );

}

void B_faces_to_cells( struct domain * theDomain , int type ){

   if( NUM_Q > BZZ ){
      struct cell ** theCells = theDomain->theCells;
      struct face * theFaces_1 = theDomain->theFaces_1;
      struct face * theFaces_2 = theDomain->theFaces_2;
   
      int Nr = theDomain->Nr;
      int Nz = theDomain->Nz;
      int * Np = theDomain->Np;
      double * r_jph = theDomain->r_jph;
      double * z_kph = theDomain->z_kph;

      int Nf = theDomain->fIndex_r[theDomain->N_ftracks_r];  
 
      int i,j,k;

      for( j=0 ; j<Nr ; ++j ){
         for( k=0 ; k<Nz ; ++k ){
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i ){
               struct cell * c = &(theCells[jk][i]);
               c->tempDoub = 0.0;
               if( type==0 ){
                  c->prim[BRR] = 0.0;
                  c->prim[BPP] = 0.0;
                  if( NUM_FACES==5 ) c->prim[BZZ] = 0.0;
               }else{
                  c->cons[BRR] = 0.0;
                  c->cons[BPP] = 0.0;
                  if( NUM_FACES==5 ) c->cons[BZZ] = 0.0;
               }
            }
         }
      }   

      int n;
      for( n=0 ; n<Nf ; ++n ){
         struct face * f = theFaces_1 + n;
         struct cell * cL = f->L;
         struct cell * cR = f->R;
         double Phi;
         if( f->LRtype==0 ){
            Phi = cL->Phi[2];
         }else{
            Phi = cR->Phi[1];
         }
         cL->tempDoub += f->dA;
         cR->tempDoub += f->dA;

         if( type==0 ){
            cL->prim[BRR] += Phi;
            cR->prim[BRR] += Phi;
         }else{
            cL->cons[BRR] += Phi;
            cR->cons[BRR] += Phi;
         }
      }

      for( j=0 ; j<Nr ; ++j ){
         for( k=0 ; k<Nz ; ++k ){
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i ){
               int im = i-1;
               if( im==-1 ) im = Np[jk]-1;

               struct cell * c = &(theCells[jk][i]);
               struct cell * cm = &(theCells[jk][im]);

               double xp[3] = { r_jph[j]   , c->piph  , z_kph[k]   };
               double xm[3] = { r_jph[j-1] , cm->piph , z_kph[k-1] };
               double r = get_moment_arm( xp , xm );
               double dA = get_dA( xp , xm , 0 );
               double dV = get_dV( xp , xm );

               if( type==0 ){
                  c->prim[BRR] /= c->tempDoub;
                  c->prim[BPP] = .5*(c->Phi[0]+cm->Phi[0])/dA;
               }else{
                  c->cons[BRR] *= dV/r/c->tempDoub;
                  c->cons[BPP] = .5*(c->Phi[0]+cm->Phi[0])*dV/dA/r;
               }
            }
         }
      }

      if( NUM_FACES == 5 ){
         for( j=0 ; j<Nr ; ++j ){
            for( k=0 ; k<Nz ; ++k ){
               int jk = j+Nr*k;
               for( i=0 ; i<Np[jk] ; ++i ){
                  theCells[jk][i].tempDoub = 0.0;
               }
            }
         }
         int Nfz = theDomain->fIndex_z[theDomain->N_ftracks_z];  
         for( n=0 ; n<Nfz ; ++n ){
            struct face * f = theFaces_2 + n;
            struct cell * cL = f->L;
            struct cell * cR = f->R;
            double Phi;
            if( f->LRtype==0 ){
               Phi = cL->Phi[4];
            }else{
               Phi = cR->Phi[3];
            }
            cL->tempDoub += f->dA;
            cR->tempDoub += f->dA;

            if( type==0 ){
               cL->prim[BZZ] += Phi;
               cR->prim[BZZ] += Phi;
            }else{
               cL->cons[BZZ] += Phi;
               cR->cons[BZZ] += Phi;
            }
         }

         for( j=0 ; j<Nr ; ++j ){
            for( k=0 ; k<Nz ; ++k ){
               int jk = j+Nr*k;
               for( i=0 ; i<Np[jk] ; ++i ){
                  int im = i-1; 
                  if( im==-1 ) im = Np[jk]-1;

                  struct cell * c = &(theCells[jk][i]);
                  struct cell * cm = &(theCells[jk][im]);

                  double xp[3] = { r_jph[j]   , c->piph  , z_kph[k]   };   
                  double xm[3] = { r_jph[j-1] , cm->piph , z_kph[k-1] };
                  double dV = get_dV( xp , xm );

                  if( type==0 ){
                     c->prim[BZZ] /= c->tempDoub;
                  }else{
                     c->cons[BZZ] *= dV/c->tempDoub;
                  }    
               }    
            }    
         }  
      }
   }
}

void make_edge_adjust( struct domain * , double );

void update_B_fluxes( struct domain * theDomain , double dt ){

   struct face * theFaces = theDomain->theFaces_1;
   int * Nf = theDomain->fIndex_r;
   int Njk = theDomain->N_ftracks_r;
   int n;
   int jk;
   int n0 = 0;
   for( jk=0 ; jk<Njk ; ++jk ){
      for( n=0 ; n<Nf[jk]-n0 ; ++n ){
         struct face * f  = theFaces + n0 + n;
         int np = n+1;
         if( np==Nf[jk]-n0 ) np=0; 
         struct face * fp = theFaces + n0 + np;

         double E;
         double dl = f->dl;
         if( f->LRtype == 0 ){
            E = f->L->E[1];
            f->L->Phi[2] -= E*dl*dt;
            f->L->Phi[0] += E*dl*dt;
         }else{
            E = f->R->E[0];
            f->R->Phi[1] -= E*dl*dt;
            f->R->Phi[0] -= E*dl*dt;
         }
         if( fp->LRtype == 0 ){
            fp->L->Phi[2] += E*dl*dt;
         }else{
            fp->R->Phi[1] += E*dl*dt;
         }
      }
      n0 = Nf[jk];
   }
   if( NUM_FACES == 5 && NUM_EDGES == 8 ){
      theFaces = theDomain->theFaces_2;
      Nf = theDomain->fIndex_z;
      Njk = theDomain->N_ftracks_z;
      n0 = 0; 
      for( jk=0 ; jk<Njk ; ++jk ){
         for( n=0 ; n<Nf[jk]-n0 ; ++n ){
            struct face * f  = theFaces + n0 + n; 
            int np = n+1; 
            if( np==Nf[jk]-n0 ) np=0; 
            struct face * fp = theFaces + n0 + np;

            double E;
            double dl = f->dl;
            if( f->LRtype == 0 ){ 
               E = f->L->E[5];
               f->L->Phi[4] -= E*dl*dt;
               f->L->Phi[0] += E*dl*dt;
            }else{
               E = f->R->E[4];
               f->R->Phi[3] -= E*dl*dt;
               f->R->Phi[0] -= E*dl*dt;
            }    
            if( fp->LRtype == 0 ){ 
               fp->L->Phi[4] += E*dl*dt;
            }else{
               fp->R->Phi[3] += E*dl*dt;
            }    
         }    
         n0 = Nf[jk];
      }

      if( NUM_AZ_EDGES == 4 ) make_edge_adjust( theDomain , dt );
   } 
}

void add_E_phi( double * phiL , double * phiR , double * phiD , double * phiU , double Edldt ){
   *phiL -= Edldt;
   *phiR += Edldt;
   *phiU -= Edldt;
   *phiD += Edldt;
}

double get_dp( double , double );

void avg_Efields( struct domain * theDomain ){

   int i,j,k;
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;

   for( j=0 ; j<Nr ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c  = theCells[jk]+i;
            int ip = (i+1)%Np[jk];
            struct cell * cp = theCells[jk]+ip;

            double El_avg = .5*( c->E[0] + cp->E[2] );
            double Er_avg = .5*( c->E[1] + cp->E[3] );

             c->E[0] = El_avg;
             c->E[1] = Er_avg;
            cp->E[2] = El_avg;
            cp->E[3] = Er_avg;

            double Bl_avg = .5*( c->B[0] + cp->B[2] );
            double Br_avg = .5*( c->B[1] + cp->B[3] );

             c->B[0] = Bl_avg;
             c->B[1] = Br_avg;
            cp->B[2] = Bl_avg;
            cp->B[3] = Br_avg;

            if( NUM_EDGES == 8 ){
               double El_avg = .5*( c->E[4] + cp->E[6] );
               double Er_avg = .5*( c->E[5] + cp->E[7] );

                c->E[4] = El_avg;
                c->E[5] = Er_avg;
               cp->E[6] = El_avg;
               cp->E[7] = Er_avg;

               double Bl_avg = .5*( c->B[4] + cp->B[6] );
               double Br_avg = .5*( c->B[5] + cp->B[7] );

                c->B[4] = Bl_avg;
                c->B[5] = Br_avg;
               cp->B[6] = Bl_avg;
               cp->B[7] = Br_avg;
            }

         }
      }
   }

   int Nf = theDomain->fIndex_r[theDomain->N_ftracks_r];
   struct face * theFaces = theDomain->theFaces_1;
   int n;
   for( n=0 ; n<Nf ; ++n ){
      struct face * f = theFaces+n;
      struct cell * c1;
      struct cell * c2;
      if( f->LRtype == 0 ){
         c1 = f->L;
         c2 = f->R;
      }else{
         c1 = f->R;
         c2 = f->L;
      }
      double p1 = c1->piph;
      double p2 = c2->piph;
      double dp1 = get_dp(p2,p1);
      double dp2 = c2->dphi - dp1;
      if( f->LRtype == 0 ){
         double Eavg = ( dp2*c2->E[0] + dp1*c2->E[2] )/(dp1+dp2);
         double Bavg = ( dp2*c2->B[0] + dp1*c2->B[2] )/(dp1+dp2);
         f->E = .5*(f->L->E[1] + Eavg);
         f->B = .5*(f->L->B[1] + Bavg);
      }else{
         double Eavg = ( dp2*c2->E[1] + dp1*c2->E[3] )/(dp1+dp2);
         double Bavg = ( dp2*c2->B[1] + dp1*c2->B[3] )/(dp1+dp2);
         f->E = .5*(f->R->E[0] + Eavg);
         f->B = .5*(f->R->B[0] + Bavg);
      }
   }

   for( n=0 ; n<Nf ; ++n ){
      struct face * f = theFaces+n;
      if( f->LRtype==0 ){
         f->L->E[1] = f->E;
         f->L->B[1] = f->B;
      }else{
         f->R->E[0] = f->E;
         f->R->B[0] = f->B;
      }
      f->E = 0.0;
      f->B = 0.0;
   }

   if( NUM_EDGES == 8 ){
//REPEAT THE ABOVE FOR VERTICALLY-ORIENTED FACES & RADIAL EDGES
      Nf = theDomain->fIndex_z[theDomain->N_ftracks_z];
      theFaces = theDomain->theFaces_2;
      int n;
      for( n=0 ; n<Nf ; ++n ){
         struct face * f = theFaces+n;
         struct cell * c1;
         struct cell * c2;
         if( f->LRtype == 0 ){ 
            c1 = f->L;
            c2 = f->R;
         }else{
            c1 = f->R;
            c2 = f->L;
         }    
         double p1 = c1->piph;
         double p2 = c2->piph;
         double dp1 = get_dp(p2,p1);
         double dp2 = c2->dphi - dp1; 
         if( f->LRtype == 0 ){ 
            double Eavg = ( dp2*c2->E[4] + dp1*c2->E[6] )/(dp1+dp2);
            double Bavg = ( dp2*c2->B[4] + dp1*c2->B[6] )/(dp1+dp2);
            f->E = .5*(f->L->E[5] + Eavg);
            f->B = .5*(f->L->B[5] + Bavg);
         }else{
            double Eavg = ( dp2*c2->E[5] + dp1*c2->E[7] )/(dp1+dp2);
            double Bavg = ( dp2*c2->B[5] + dp1*c2->B[7] )/(dp1+dp2);
            f->E = .5*(f->R->E[4] + Eavg);
            f->B = .5*(f->R->B[4] + Bavg);
         }    
      }

      for( n=0 ; n<Nf ; ++n ){
         struct face * f = theFaces+n;
         if( f->LRtype==0 ){
            f->L->E[5] = f->E;
            f->L->B[5] = f->B;
         }else{
            f->R->E[4] = f->E;
            f->R->B[4] = f->B;
         }    
         f->E = 0.0; 
         f->B = 0.0; 
      }

   }

//E_Z ALONG THE POLE...
   j=0;
   for( k=0 ; k<Nz ; ++k ){
      int jk = j+Nr*k;
      double E = 0.0;
      double B = 0.0;
      for( i=0 ; i<Np[jk] ; ++i ){
         struct cell * c = theCells[jk]+i;
         E += c->E[0]/(double)Np[jk];
      //   B += c->B[0]/(double)Np[jk];
      }
      for( i=0 ; i<Np[jk] ; ++i ){
         struct cell * c = theCells[jk]+i;
         c->E[0] = E;
         c->B[0] = B;
      }
   }

   for( j=0 ; j<Nr ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c  = theCells[jk]+i;
            int ip = (i+1)%Np[jk];
            struct cell * cp = theCells[jk]+ip;

            cp->E[2] = c->E[0];   
            cp->E[3] = c->E[1];   
            cp->B[2] = c->B[0];   
            cp->B[3] = c->B[1];  

            if( NUM_EDGES == 8 ){
 
               cp->E[6] = c->E[4];   
               cp->E[7] = c->E[5];   
               cp->B[6] = c->B[4];   
               cp->B[7] = c->B[5];  

            }
         }
      }
   }

   if( NUM_AZ_EDGES == 4 ){
//AVERAGE E_PHI SOMEHOW
   }

}

void subtract_advective_B_fluxes( struct domain * theDomain ){

   int i,j,k;
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;

   for( j=0 ; j<Nr ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c  = theCells[jk]+i;

            double wm = c->wiph;
            double wp = c->wiph;

            c->E[0] -= wm*c->B[0];
            c->E[1] -= wp*c->B[1];

            if( NUM_EDGES == 8 ){
               c->E[4] -= wm*c->B[4];
               c->E[5] -= wp*c->B[5];
            }
         }
      }
   }

}

double get_signed_dp( double , double );

void check_flipped( struct domain * theDomain ){

   struct face * theFaces = theDomain->theFaces_1;
   int * fI = theDomain->fIndex_r;
   int Nf_t = theDomain->N_ftracks_r;

   int n;
   int n0 = 0; 
   int jk;
   for( jk=0 ; jk<Nf_t ; ++jk ){
      for( n=0 ; n<fI[jk]-n0 ; ++n ){
         struct face * f  = theFaces + n0 + n; 
         int np = n+1; 
         if( np==fI[jk]-n0 ) np=0; 
         struct face * fp = theFaces + n0 + np;

         double pp,pm;
         if( f->LRtype==0 ){
            pm = f->L->piph;
         }else{
            pm = f->R->piph;
         }
         if( fp->LRtype==0 ){
            pp = fp->L->piph;
         }else{
            pp = fp->R->piph;
         }
         double dp = get_signed_dp( pp , pm );
         if( dp<0. ){ fp->flip_flag=1; }//printf("FLIPPED YO\n");}
      }    
      n0 = fI[jk];
   }

}

void get_phi_pointer( struct face * f , double ** P , double ** RK_P ){

   if( f->LRtype == 0 ){ 
      *P    = &(f->L->Phi[2]);
      *RK_P = &(f->L->RK_Phi[2]);
   }else{
      *P    = &(f->R->Phi[1]);
      *RK_P = &(f->R->RK_Phi[1]);
   }

}

void flip_fluxes( struct domain * theDomain ){

   struct face * theFaces = theDomain->theFaces_1;
   int * fI = theDomain->fIndex_r;
   int Nf_t = theDomain->N_ftracks_r;

   int n;
   int n0 = 0;
   int jk;
   for( jk=0 ; jk<Nf_t ; ++jk ){
      for( n=0 ; n<fI[jk]-n0 ; ++n ){
         struct face * f  = theFaces + n0 + n;
         if( f->flip_flag ){

            int np = n+1;
            if( np==fI[jk]-n0 ) np=0;
            struct face * fp = theFaces + n0 + np;
            int nm = n-1;
            if( nm==-1 ) nm = fI[jk]-n0-1;
            struct face * fm = theFaces + n0 + nm;

            double *P, *Pm, *Pp, *RK_P, *RK_Pm, *RK_Pp;

            get_phi_pointer( f  , &P  , &RK_P  );
            get_phi_pointer( fm , &Pm , &RK_Pm );
            get_phi_pointer( fp , &Pp , &RK_Pp );

            double Phi = *P;
            double RK_Phi = *RK_P;

            *P      = *Pm    + Phi;
            *RK_P   = *RK_Pm + RK_Phi;
            *Pm     = -Phi;
            *RK_Pm  = -RK_Phi;
            *Pp     = *Pp    + Phi;
            *RK_Pp  = *RK_Pp + RK_Phi;

            f->LRtype  = !f->LRtype;
            fm->LRtype = !fm->LRtype;

         }
      }    
      n0 = fI[jk];
   }

}


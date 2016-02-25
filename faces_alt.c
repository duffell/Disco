
#include "paul.h"

double get_dA( double * , double * , int );
double get_dL( double * , double * , int );
double get_dp( double , double );

int get_num_rzFaces( int Nr , int Nz , int dim ){
   if( dim==1 ) return( (Nr-1)*Nz );
   else return( (Nz-1)*Nr );
}

void addFace( struct face * theFaces , int n , struct cell * cL , struct cell * cR , double dxL , double dxR , double * xp , double * xm , int dim , int LRtype ){
   int d;
   for( d=0 ; d<3 ; ++d ) theFaces[n].cm[d] = .5*(xp[d]+xm[d]); //Consider calculating center of mass in geometry.c
   double dp = get_dp(xp[1],xm[1]);
   double phic = get_dp(xp[1],.5*dp);
   double rp = xp[0];
   double rm = xm[0];
   double r2 = (rp*rp+rm*rm+rp*rm)/3.;
   theFaces[n].cm[0] = r2/(.5*(rp+rm));
   theFaces[n].cm[1] = phic;
   theFaces[n].L   = cL;
   theFaces[n].R   = cR;
   theFaces[n].dxL = dxL;
   theFaces[n].dxR = dxR;
   theFaces[n].dphi= dp;//get_dp(xp[1],xm[1]);
   theFaces[n].dA  = get_dA(xp,xm,dim);

   int dim_trans = 3-dim;//4-2*dim;//3-dim;
   theFaces[n].dl = get_dL( xp , xm , dim_trans );//xp[dim_trans] - xm[dim_trans];

   theFaces[n].E = 0.0;
   theFaces[n].B = 0.0;
   theFaces[n].LRtype = LRtype;
   theFaces[n].flip_flag = 0;

}

void buildfaces( struct domain * theDomain , int dim , int mode ){
  
   struct cell ** theCells = theDomain->theCells;
   struct face * theFaces;
   int * ntj;
   if( dim==1 ){
      theFaces = theDomain->theFaces_1;
      ntj = theDomain->fIndex_r;
   }else{
      theFaces = theDomain->theFaces_2;
      ntj = theDomain->fIndex_z;
   }

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   double Pmax = theDomain->phi_max;
   int i,j,k; 

   int I0[Nr*Nz];

   for( k=0 ; k<Nz ; ++k ){
      for( j=0 ; j<Nr ; ++j ){
         int jk = j+Nr*k;
         int found=0;
         int quad_prev=0;
         for( i=0 ; i<Np[jk] && !found ; ++i ){
            struct cell * c = theCells[jk]+i;
            double convert = 2.*M_PI/Pmax;
            double sn = sin(c->piph*convert);
            double cs = cos(c->piph*convert);
            if( sn>0. && cs>0. && quad_prev ){
               quad_prev = 0;
               found = 1;
               I0[jk] = i; 
            }else if( sn<0. && cs>0. ){
               quad_prev=1;
            }
         }
         if( !found ) I0[jk]=0;
      }
   }

   int n=0;
   int Nrmax = Nr-1;
   if( dim==2 ) Nrmax = Nr;
   int Nzmax = Nz-1;
   if( dim==1 ) Nzmax = Nz;

   if( mode==0 ){
      int n=0;
      for( k=0 ; k<Nzmax ; ++k ){
         for( j=0 ; j<Nrmax ; ++j ){
            int JK  = j+Nrmax*k;
            int jk  = j+Nr*k;
            int jkp = j+1 + Nr*k;
            if( dim==2 ) jkp = j+Nr*(k+1);

            ntj[JK] = n;
            n += Np[jk]+Np[jkp];
         }
      }
      ntj[Nrmax*Nzmax] = n;
   }else{
      for( k=0 ; k<Nzmax ; ++k ){
         for( j=0 ; j<Nrmax ; ++j ){
            int jp = j+1;
            int kp = k+1;

            int jk  = j  + Nr*k;
            int jkp = jp + Nr*k;
            if( dim==2 ) jkp = j + Nr*kp;

            double dxL,dxR;
            if( dim==1 ){
               dxL = .5*(r_jph[j]  - r_jph[j-1]);
               dxR = .5*(r_jph[jp] - r_jph[j]  );
            }else{
               dxL = .5*(z_kph[k]  - z_kph[k-1]);
               dxR = .5*(z_kph[kp] - z_kph[k]  );
            }
            double xp[3] = {r_jph[j],0.0,z_kph[k  ]};
            double xm[3] = {r_jph[j],0.0,z_kph[k-1]};
            if( dim==2 ){ xm[0] = r_jph[j-1] ; xm[2] = z_kph[k] ; }

            int i  = I0[jk];
            int ip = I0[jkp];

            int f;
            for( f=0 ; f<Np[jk]+Np[jkp] ; ++f ){
 
               struct cell * cL = &(theCells[jk ][i] );
               struct cell * cR = &(theCells[jkp][ip]);
 
               double dphi = cL->piph - cR->piph;
               while( dphi > .5*Pmax ) dphi -= Pmax;
               while( dphi <-.5*Pmax ) dphi += Pmax;
 
               int LR = 0;
               if( dphi > 0. ) LR = 1;

               dphi = cL->piph-cL->dphi - cR->piph+cR->dphi;
               while( dphi > .5*Pmax ) dphi -= Pmax;
               while( dphi <-.5*Pmax ) dphi += Pmax;

               int LR_back = 0;
               if( dphi < 0. ) LR_back = 1;

               if( LR == 0 ){
                  xp[1] = cL->piph;
               }else{
                  xp[1] = cR->piph;
               }
               if( LR_back == 0 ){
                  xm[1] = cL->piph-cL->dphi;
               }else{
                  xm[1] = cR->piph-cR->dphi;
               }
               addFace( theFaces , n , cL , cR , dxL , dxR , xp , xm , dim , LR );
               ++n;
            
               if( LR == 0 ){
                  ++i;
                  if( i == Np[jk] ) i=0;
               }else{
                  ++ip;
                  if( ip== Np[jkp] ) ip=0;
               }
            }
         }
      }
   }

}

double get_dL( double * , double * , int );
//*phiL,*phiR,*phiD,*phiU,Edldt
void add_E_phi( double * , double * , double * , double * , double );

int phi_switch( double dphi , double Pmax , int mode ){

   while( dphi > .5*Pmax ) dphi -= Pmax;
   while( dphi <-.5*Pmax ) dphi += Pmax;
   if( mode == 1 ) dphi = -dphi;

   int LR = 0;
   if( dphi > 0.) LR = 1;

   return( LR );

}

int get_which4( double phi , double phiR , double phiU , double phiUR , int * LR_alt , int * UD_alt , int mode , double Pmax ){

   int which4;
   double dphi;
   
   dphi = phi - phiR;
   int LR_D = phi_switch( dphi , Pmax , mode );

   dphi = phiU - phiUR;
   int LR_U = phi_switch( dphi , Pmax , mode );

   double phi1 = phi;
   if( LR_D ) phi1 = phiR;
   double phi2 = phiU;
   if( LR_U ) phi2 = phiUR;

   dphi = phi1-phi2;
   int UD = phi_switch( dphi , Pmax , mode );

   if( UD==0 ){
      if( LR_D==0 ) which4 = 0;
      else which4 = 1;
   }else{
      if( LR_U==0 ) which4 = 2;
      else which4 = 3;
   }

   if( mode == 0 ){
      if( which4 == 0 || which4 == 1 ) *LR_alt = LR_U;
      else                             *LR_alt = LR_D;

      if( which4 == 0 || which4 == 2 ){
         dphi = phiR - phiUR;
         *UD_alt = phi_switch( dphi , Pmax , mode );
      }else{
         dphi = phi - phiU;
         *UD_alt = phi_switch( dphi , Pmax , mode );
      }
   }

   return( which4 );
}

void make_edge_adjust( struct domain * theDomain , double dt ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   double Pmax = theDomain->phi_max;
   int i,j,k;

   int I0[Nr*Nz];
   for( k=0 ; k<Nz ; ++k ){
      for( j=0 ; j<Nr ; ++j ){
         int jk = j+Nr*k;
         int found=0;
         int quad_prev=0;
         for( i=0 ; i<Np[jk] && !found ; ++i ){
            struct cell * c = theCells[jk]+i;
            double convert = 2.*M_PI/Pmax;
            double sn = sin(c->piph*convert);
            double cs = cos(c->piph*convert);
            if( sn>0. && cs>0. && quad_prev ){
               quad_prev = 0; 
               found = 1; 
               I0[jk] = i; 
            }else if( sn<0. && cs>0. ){
               quad_prev=1;
            }    
         }    
         if( !found ) I0[jk]=0;
      }    
   }

   for( k=0 ; k<Nz-1 ; ++k ){
      for( j=0 ; j<Nr-1 ; ++j ){
         int jk   = j  +Nr*k;
         int jkR  = j+1+Nr*k;
         int jkU  = j  +Nr*(k+1);
         int jkUR = j+1+Nr*(k+1);

         double xp[3] = {r_jph[j],0.0,z_kph[k]};
         double xm[3] = {r_jph[j],0.0,z_kph[k]};

         int i   = I0[jk  ];
         int iR  = I0[jkR ];
         int iU  = I0[jkU ];
         int iUR = I0[jkUR];

         int Ne = Np[jk] + Np[jkR] + Np[jkU] + Np[jkUR];
         int e;
         for( e=0 ; e<Ne ; ++e ){

            struct cell * c   = &(theCells[jk  ][i  ]);
            struct cell * cR  = &(theCells[jkR ][iR ]);
            struct cell * cU  = &(theCells[jkU ][iU ]);
            struct cell * cUR = &(theCells[jkUR][iUR]);

            int LR_alt;
            int UD_alt;
            int buffer;
            int which4      = get_which4( c->piph         , cR->piph          , cU->piph          , cUR->piph           , &LR_alt , &UD_alt , 0 , Pmax );
            int which4_back = get_which4( c->piph-c->dphi , cR->piph-cR->dphi , cU->piph-cU->dphi , cUR->piph-cUR->dphi , &buffer , &buffer , 1 , Pmax );

            double * PhiL;
            double * PhiR;
            double * PhiU;
            double * PhiD;

            double E;

            if( which4 == 0 ){
               PhiL = c->Phi+4;
               PhiD = c->Phi+2;
               E = .25*( c->E_phi[3] + c->E_phi[1] );
               if( LR_alt==0 ){
                  PhiU = cU->Phi+2;
                  E += .25*cU->E_phi[1];
               }else{
                  PhiU = cUR->Phi+1;
                  E += .25*cUR->E_phi[0];
               }
               if( UD_alt==0 ){
                  PhiR = cR->Phi+4;
                  E += .25*cR->E_phi[3];
               }else{
                  PhiR = cUR->Phi+3;
                  E += .25*cUR->E_phi[2];
               }
            }else if( which4 == 1 ){
               PhiR = cR->Phi+4;
               PhiD = cR->Phi+1;
               E = .25*( cR->E_phi[3] + cR->E_phi[0] );
               if( LR_alt==0 ){
                  PhiU = cU->Phi+2;
                  E += .25*cU->E_phi[1];
               }else{
                  PhiU = cUR->Phi+1;
                  E += .25*cUR->E_phi[0];
               }
               if( UD_alt==0 ){
                  PhiL = c->Phi+4;
                  E += .25*c->E_phi[3];
               }else{
                  PhiL = cU->Phi+3;
                  E += .25*cU->E_phi[2];
               }
            }else if( which4 == 2 ){
               PhiL = cU->Phi+3;
               PhiU = cU->Phi+2;
               E = .25*( cU->E_phi[2] + cU->E_phi[1] );
               if( LR_alt==0 ){
                  PhiD = c->Phi+2;
                  E += .25*c->E_phi[1];
               }else{
                  PhiD = cR->Phi+1;
                  E += .25*cR->E_phi[0];
               }
               if( UD_alt==0 ){
                  PhiR = cR->Phi+4;
                  E += .25*cR->E_phi[3];
               }else{
                  PhiR = cUR->Phi+3;
                  E += .25*cUR->E_phi[2];
               }
            }else{
               PhiR = cUR->Phi+3;
               PhiU = cUR->Phi+1;
               E = .25*( cU->E_phi[2] + cU->E_phi[0] );
               if( LR_alt==0 ){
                  PhiD = c->Phi+2;
                  E += .25*c->E_phi[1];
               }else{
                  PhiD = cR->Phi+1;
                  E += .25*cR->E_phi[0];
               }
               if( UD_alt==0 ){
                  PhiL = c->Phi+4;
                  E += .25*cR->E_phi[3];
               }else{
                  PhiL = cU->Phi+3;
                  E += .25*cUR->E_phi[2];
               }
            } 

            if( which4 == 0 ){
               xp[1] = c->piph;
            }else if( which4 == 1 ){
               xp[1] = cR->piph;
            }else if( which4 == 2 ){
               xp[1] = cU->piph;
            }else{
               xp[1] = cUR->piph;
            }

            if( which4_back == 0 ){
               xm[1] = c->piph  - c->dphi;
            }else if( which4_back == 1 ){
               xm[1] = cR->piph - cR->dphi;
            }else if( which4_back == 2 ){
               xm[1] = cU->piph - cU->dphi;
            }else{
               xm[1] = cUR->piph- cUR->dphi;
            }

            double dl = get_dL( xp , xm , 0 );
//if( e==0 ) printf("dl = %e which4 = %d, which4_back = %d, phip = %e phim = %e dphi=%e \n",dl,which4,which4_back,xp[1],xm[1],xp[1]-xm[1]);
            add_E_phi( PhiL , PhiR , PhiD , PhiU , E*dl*dt );

            if( which4 == 0 ){
               ++i;
               if( i   == Np[jk  ] ) i  =0;
            }else if( which4 == 1 ){
               ++iR;
               if( iR  == Np[jkR ] ) iR =0;
            }else if( which4 == 2 ){
               ++iU;
               if( iU  == Np[jkU ] ) iU =0;
            }else{
               ++iUR;
               if( iUR == Np[jkUR] ) iUR=0;
            }
         }

      }
   }
}


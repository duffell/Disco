
enum{RHX,PPX,UXX,UYY,UUZ,BXX,BYY,BBZ};
enum{DEN,TAX,SXX,SYY,SSZ};

#include <stdio.h>
#include <math.h>
#include "paul.h"

//#define GAMMA_LAW (5./3.)
#define TOL 1e-10

static double GAMMA_LAW = 0.0;

void setHlldParams( struct domain * theDomain ){
   GAMMA_LAW = theDomain->theParList.Adiabatic_Index;
}

void get_velocities( double * pL , double * pR , double * n , double * vel , double * Pt_star ){
  
   double gam = GAMMA_LAW;

   double rho_L = pL[RHO];
   double P_L   = pL[PPP];
   double u_L   = pL[UXX]*n[0]    + pL[UYY]*n[1]    + pL[UZZ]*n[2];
   double Bn_L  = pL[BXX]*n[0]    + pL[BYY]*n[1]    + pL[BZZ]*n[2];
   double B2_L  = pL[BXX]*pL[BXX] + pL[BYY]*pL[BYY] + pL[BZZ]*pL[BZZ];
   double Pt_L  = P_L + .5*B2_L;

   double rho_R = pR[RHO];
   double P_R   = pR[PPP];
   double u_R   = pR[UXX]*n[0]    + pR[UYY]*n[1]    + pR[UZZ]*n[2];
   double Bn_R  = pR[BXX]*n[0]    + pR[BYY]*n[1]    + pR[BZZ]*n[2];
   double B2_R  = pR[BXX]*pR[BXX] + pR[BYY]*pR[BYY] + pR[BZZ]*pR[BZZ];
   double Pt_R  = P_R + .5*B2_R;

   double Bn = .5*(Bn_L+Bn_R);
 
   double cf2_L = .5*(gam*P_L + B2_L + sqrt( fabs( pow( gam*P_L + B2_L , 2. ) - 4.*gam*P_L*Bn_L*Bn_L ) ) )/rho_L;
   double cf2_R = .5*(gam*P_R + B2_R + sqrt( fabs( pow( gam*P_R + B2_R , 2. ) - 4.*gam*P_R*Bn_R*Bn_R ) ) )/rho_R;

   double cf2;
   if( cf2_L > cf2_R ) cf2 = cf2_L; else cf2 = cf2_R;
   double umin,umax;
   if( u_L < u_R ){ umin = u_L; umax = u_R; }else{ umin = u_R; umax = u_L; }

   double S_R = umax + sqrt(fabs(cf2));
   double S_L = umin - sqrt(fabs(cf2));

   double S_M = ( ( S_R - u_R )*rho_R*u_R - ( S_L - u_L )*rho_L*u_L - Pt_R + Pt_L )/( ( S_R - u_R )*rho_R - (S_L-u_L)*rho_L );

   double rhostar_L = rho_L*( S_L - u_L )/( S_L - S_M );
   double rhostar_R = rho_R*( S_R - u_R )/( S_R - S_M );

   double Ss_L = S_M - fabs(Bn)/sqrt(rhostar_L);
   double Ss_R = S_M + fabs(Bn)/sqrt(rhostar_R);

   vel[0] = S_L;
   vel[1] = Ss_L;
   vel[2] = S_M;
   vel[3] = Ss_R;
   vel[4] = S_R;
   
   double num   = ( S_R - u_R )*rho_R*Pt_L - ( S_L - u_L )*rho_L*Pt_R + rho_L*rho_R*( S_R - u_R )*( S_L - u_L )*( u_R - u_L );
   double denom = ( S_R - u_R )*rho_R - ( S_L - u_L )*rho_L;

   *Pt_star = num/denom;
}

void prim_to_cons( double * prim , double * cons ){

   double rho  = prim[RHO];
   double Pp   = prim[PPP];
   double vx   = prim[UXX];
   double vy   = prim[UYY];
   double vz   = prim[UZZ];

   double Bx  = prim[BXX];
   double By  = prim[BYY];
   double Bz  = prim[BZZ];

   double v2 = vx*vx + vy*vy + vz*vz;
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);

   double B2 = Bx*Bx + By*By + Bz*Bz;

   cons[DEN] = rho;
   cons[SXX] = rho*vx;
   cons[SYY] = rho*vy;
   cons[SZZ] = rho*vz;
   cons[TAU] = ( .5*rho*v2 + rhoe + .5*B2 );

   cons[BXX] = Bx;
   cons[BYY] = By;
   cons[BZZ] = Bz;

}

void get_flux( double * prim , double * flux , double * n ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vx  = prim[UXX];
   double vy  = prim[UYY];
   double vz  = prim[UZZ];

   double Bx   = prim[BXX];
   double By   = prim[BYY];
   double Bz   = prim[BZZ];

   double v2 = vx*vx+vy*vy+vz*vz;
   double vB = vx*Bx + vy*By + vz*Bz;
   double B2 = Bx*Bx + By*By + Bz*Bz;
   double vn = vx*n[0] + vy*n[1] + vz*n[2];
   double Bn = Bx*n[0] + By*n[1] + Bz*n[2];

   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);

   flux[DEN] = rho*vn;
   flux[SXX] = rho*vx*vn + ( Pp + .5*B2 )*n[0] - Bx*Bn;
   flux[SYY] = rho*vy*vn + ( Pp + .5*B2 )*n[1] - By*Bn;
   flux[SZZ] = rho*vz*vn + ( Pp + .5*B2 )*n[2] - Bz*Bn;
   flux[TAU] = (.5*rho*v2 + rhoe + Pp + B2 )*vn - vB*Bn;

   flux[BXX] = Bx*vn - vx*Bn;
   flux[BYY] = By*vn - vy*Bn;
   flux[BZZ] = Bz*vn - vz*Bn;

}

void get_single_star( double * prim , double * Ustar , double S_K , double S_M , double Pt_star , double * n ){

   double rho = prim[RHO];
   double vx  = prim[UXX];
   double vy  = prim[UYY];
   double vz  = prim[UZZ];
   double Pp  = prim[PPP];
   double Bx  = prim[BXX];
   double By  = prim[BYY];
   double Bz  = prim[BZZ];

   double u_K   = vx*n[0] + vy*n[1] + vz*n[2];
   double Bn    = Bx*n[0] + By*n[1] + Bz*n[2];
   double vdotB = vx*Bx + vy*By + vz*Bz;
   double v2    = vx*vx + vy*vy + vz*vz;
   double B2    = Bx*Bx + By*By + Bz*Bz;

   double Pt = Pp + .5*B2;
 
   double gam = GAMMA_LAW;
   double e = .5*rho*v2 + Pp/(gam-1.) + .5*B2;

   double rhostar = rho*( S_K - u_K )/( S_K - S_M );

   double denom = rho*(S_K-u_K)*(S_K-S_M) - Bn*Bn ;
   if( fabs(denom) < TOL ) denom = TOL;

   double vx_star = vx - Bn*Bx*( S_M - u_K )/denom;
   double vy_star = vy - Bn*By*( S_M - u_K )/denom;
   double vz_star = vz - Bn*Bz*( S_M - u_K )/denom;
 
   double Bx_star = Bx*( rho*pow(S_K-u_K,2.) - Bn*Bn )/denom;
   double By_star = By*( rho*pow(S_K-u_K,2.) - Bn*Bn )/denom;
   double Bz_star = Bz*( rho*pow(S_K-u_K,2.) - Bn*Bn )/denom;

   double vn_star = vx_star*n[0] + vy_star*n[1] + vz_star*n[2];
   double Bn_star = Bx_star*n[0] + By_star*n[1] + Bz_star*n[2];

   vx_star += ( S_M - vn_star )*n[0];
   vy_star += ( S_M - vn_star )*n[1];
   vz_star += ( S_M - vn_star )*n[2];

   Bx_star += ( Bn  - Bn_star )*n[0];
   By_star += ( Bn  - Bn_star )*n[1];
   Bz_star += ( Bn  - Bn_star )*n[2];

   double vBstar = vx_star*Bx_star + vy_star*By_star + vz_star*Bz_star;
   double e_star = ( (S_K-u_K)*e - Pt*u_K + Pt_star*S_M + Bn*( vdotB - vBstar ) )/( S_K - S_M );

   Ustar[DEN] = rhostar;
   Ustar[SXX] = rhostar*vx_star;
   Ustar[SYY] = rhostar*vy_star;
   Ustar[SZZ] = rhostar*vz_star;
   Ustar[TAU] = e_star;
   Ustar[BXX] = Bx_star;
   Ustar[BYY] = By_star;
   Ustar[BZZ] = Bz_star;

}

void get_double_star( double * UsL , double * UsR , double * Uss , double * vel , double * n , int LR ){

   double rho_L = UsL[RHO];
   double rho_R = UsR[RHO];
   double rhostar_K;
   if( LR==0 ) rhostar_K = rho_L; else rhostar_K = rho_R;
   double rrh_L = sqrt( rho_L );
   double rrh_R = sqrt( rho_R );
   double rrh_K = sqrt( rhostar_K );
  
   double Sx_L = UsL[SXX];
   double Sy_L = UsL[SYY];
   double Sz_L = UsL[SZZ];
   double Bx_L = UsL[BXX];
   double By_L = UsL[BYY];
   double Bz_L = UsL[BZZ];

   double vx_L = Sx_L/rho_L;
   double vy_L = Sy_L/rho_L;
   double vz_L = Sz_L/rho_L;
 
   double Sx_R = UsR[SXX];
   double Sy_R = UsR[SYY];
   double Sz_R = UsR[SZZ];
   double Bx_R = UsR[BXX];
   double By_R = UsR[BYY];
   double Bz_R = UsR[BZZ];

   double vx_R = Sx_R/rho_R;
   double vy_R = Sy_R/rho_R;
   double vz_R = Sz_R/rho_R;
 
   double Bn_L = Bx_L*n[0] + By_L*n[1] + Bz_L*n[2];
   double Bn_R = Bx_R*n[0] + By_R*n[1] + Bz_R*n[2];
   double Bn = .5*(Bn_L+Bn_R);
   double signBn;
   if( Bn>0.0 ) signBn = 1.0; else signBn = -1.0;
 
   double denom = rrh_L+rrh_R;

   double vx_ss = ( rrh_L*vx_L + rrh_R*vx_R + ( Bx_R - Bx_L )*signBn )/denom;
   double vy_ss = ( rrh_L*vy_L + rrh_R*vy_R + ( By_R - By_L )*signBn )/denom;
   double vz_ss = ( rrh_L*vz_L + rrh_R*vz_R + ( Bz_R - Bz_L )*signBn )/denom;

   double vn_ss = vx_ss*n[0]+vy_ss*n[1]+vz_ss*n[2];
   vx_ss += (vel[2]-vn_ss)*n[0];
   vy_ss += (vel[2]-vn_ss)*n[1];
   vz_ss += (vel[2]-vn_ss)*n[2];

   double Bx_ss = ( rrh_L*Bx_R + rrh_R*Bx_L + rrh_L*rrh_R*( vx_R - vx_L )*signBn )/denom;
   double By_ss = ( rrh_L*By_R + rrh_R*By_L + rrh_L*rrh_R*( vy_R - vy_L )*signBn )/denom;
   double Bz_ss = ( rrh_L*Bz_R + rrh_R*Bz_L + rrh_L*rrh_R*( vz_R - vz_L )*signBn )/denom;

   double Bn_ss = Bx_ss*n[0]+By_ss*n[1]+Bz_ss*n[2];
   Bx_ss += (Bn-Bn_ss)*n[0];
   By_ss += (Bn-Bn_ss)*n[1];
   Bz_ss += (Bn-Bn_ss)*n[2];

   double vB_ss;
   vB_ss = vx_ss*Bx_ss + vy_ss*By_ss + vz_ss*Bz_ss;

   double vBstar;
   double e_star;
   if( LR==0 ){
      vBstar = vx_L*Bx_L + vy_L*By_L + vz_L*Bz_L;
      e_star = UsL[TAU];
   }else{
      vBstar = vx_R*Bx_R + vy_R*By_R + vz_R*Bz_R;
      e_star = UsR[TAU];
   }
   double plusminus = -1.0;
   if( LR==1 ) plusminus = 1.0;

   double e_ss = e_star + plusminus*rrh_K*( vBstar - vB_ss )*signBn;

   Uss[DEN] = rhostar_K;
   Uss[SXX] = rhostar_K*vx_ss;
   Uss[SYY] = rhostar_K*vy_ss;
   Uss[SZZ] = rhostar_K*vz_ss;
   Uss[TAU] = e_ss;
   Uss[BXX] = Bx_ss;
   Uss[BYY] = By_ss;
   Uss[BZZ] = Bz_ss;

}

void get_Ustar_HLLD( double w , double * pL , double * pR , double * F , double * U , double r , double * n ){

   //Get Velocities
   double vel[5];
   double Pt_star;
   pL[UYY] *= r;
   pR[UYY] *= r;
   get_velocities( pL , pR , n , vel , &Pt_star );

   int q;

   if( w < vel[2] ){
      //Left Side
      double UL[8];
      get_flux( pL , F , n );
      prim_to_cons( pL , UL );
      if( w < vel[0] ){
         //Upwind Left
         prim_to_cons( pL , U );
      }else if( w < vel[1] ){
         //Single Star Left
         get_single_star( pL , U , vel[0] , vel[2] , Pt_star , n );
         for( q=0 ; q<8 ; ++q ){
            F[q] += vel[0]*( U[q] - UL[q] );
         }
      }else{
         //Double Star Left
         double UsL[8];
         double UsR[8];
         get_single_star( pL , UsL , vel[0] , vel[2] , Pt_star , n );
         get_single_star( pR , UsR , vel[4] , vel[2] , Pt_star , n );
         get_double_star( UsL , UsR , U , vel , n , 0 );
         for( q=0 ; q<8 ; ++q ){
            F[q] += vel[1]*( U[q] - UsL[q] ) + vel[0]*( UsL[q] - UL[q] );
         }
      }
   }else{
      //Right Side
      double UR[8];
      get_flux( pR , F , n );
      prim_to_cons( pR , UR );
      if( w > vel[4] ){
         //Upwind Right
         prim_to_cons( pR , U );
      }else if( w > vel[3] ){
         //Single Star Right
         get_single_star( pR , U , vel[4] , vel[2] , Pt_star , n );
         for( q=0 ; q<8 ; ++q ){
            F[q] += vel[4]*( U[q] - UR[q] );
         }
      }else{
         //Double Star Right
         double UsL[8];
         double UsR[8];
         get_single_star( pL , UsL , vel[0] , vel[2] , Pt_star , n );
         get_single_star( pR , UsR , vel[4] , vel[2] , Pt_star , n );
         get_double_star( UsL , UsR , U , vel , n , 1 );
         for( q=0 ; q<8 ; ++q ){
            F[q] += vel[3]*( U[q] - UsR[q] ) + vel[4]*( UsR[q] - UR[q] );
         }
      }
   }
   F[BXX] /= r;
   F[BYY] /= r;
   F[LLL] *= r;
   U[BXX] /= r;
   U[BYY] /= r;
   U[LLL] *= r;
}



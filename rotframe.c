#include "paul.h"

static int om_flag = 0;
static double Omega0 = 0.0;
static double d = 0.0;

void setOmegaParams( struct domain * theDomain ){
   om_flag = 0;
   Omega0 = 1.0;
   d = 0.5;
}

double get_dp( double , double );

void subtract_omega( double * prim ){
   if( om_flag ) prim[UPP] -= Omega0;
}

void omegaForce( double r , double phi , double vr , double omega , double * fr , double * fp ){

   //Omega0^2 vec R + 2 Omega0 x v

   double Rr = r-d*cos(phi);
   double Rp = d*sin(phi);

   *fr = Omega0*Omega0*Rr + 2.*Omega0*r*omega;
   *fp = Omega0*Omega0*Rp - 2.*Omega0*vr;

}

void omega_src( double * prim , double * cons , double * xp , double * xm , double dVdt ){

   if( om_flag ){

      double rp = xp[0];
      double rm = xm[0];
      double rho = prim[RHO];
      double vr  = prim[URR];
      double omega = prim[UPP];
   
      double r = 0.5*(rp+rm);
      double vp  = r*omega;
      double dphi = get_dp(xp[1],xm[1]);
      double phi = xm[1] + 0.5*dphi;

      double Fr,Fp;
      omegaForce( r , phi , vr , omega , &Fr , &Fp );

      cons[SRR] += rho*Fr*dVdt;
      cons[LLL] += rho*Fp*r*dVdt;
      cons[TAU] += rho*( Fr*vr + Fp*vp )*dVdt;

   }

}


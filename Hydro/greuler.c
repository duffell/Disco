
#include "../paul.h"
#include "metric.h"
#include "frame.h"

double get_cs2( double );
double get_dp( double , double );
double get_dL( double * , double * , int );

static double gamma_law = 0.0; 
static double RHO_FLOOR = 0.0; 
static double PRE_FLOOR = 0.0; 
static int isothermal = 0;

void setHydroParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
   isothermal = theDomain->theParList.isothermal_flag;
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
}

int set_B_flag(void){
   return(0);
}

double get_omega(double *prim, double *x)
{
    int i,j;
    double l[3] = {prim[URR], prim[UPP], prim[UZZ]};
    
    double lapse;
    double shift[3];
    double igam[9];

    lapse = metric_lapse(x);
    metric_shift(x, shift);
    metric_igam(x, igam);

    double u2 = 0.0;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            u2 += igam[3*i+j]*l[i]*l[j];
    double w = sqrt(1.0 + u2);
    double vp = lapse*(igam[3]*l[0]+igam[4]*l[1]+igam[5]*l[2])/w - shift[1];
    
    return vp;
}

void prim2cons( double *prim, double *cons, double *x, double dV)
{
    double r = x[0];
    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double l[3] = {prim[URR], prim[UPP], prim[UZZ]};

    double lapse;
    double shift[3];
    double igam[9];
    double jac;
    double U[4];

    lapse = metric_lapse(x);
    metric_shift(x, shift);
    metric_igam(x, igam);
    jac = metric_jacobian(x) / r;
    frame_U(x, U);

    double w, u0, u2;
    double u[3];
    double igaml[3];
    int i,j;
    for(i=0; i<3; i++)
    {
        igaml[i] = 0.0;
        for(j=0; j<3; j++)
            igaml[i] += igam[3*i+j]*l[j];
    }
    u2 = l[0]*igaml[0] + l[1]*igaml[1] + l[2]*igaml[2];
    w = sqrt(1.0 + u2);
    u0 = w/lapse;
    for(i=0; i<3; i++)
        u[i] = igaml[i] - shift[i]*u0;

    double l0 = -lapse*w + shift[0]*l[0] + shift[1]*l[1] + shift[2]*l[2];
    double uU = U[0]*l0 + U[1]*l[0] + U[2]*l[1] + U[3]*l[2];
    
    double rhoh = rho + gamma_law/(gamma_law-1.0)*Pp;
    double rhoe = Pp / (gamma_law-1.0);

    cons[DDD] = jac * rho*u0 * dV;
    cons[SRR] = jac * rhoh*u0*l[0] * dV;
    cons[LLL] = jac * rhoh*u0*l[1] * dV;
    cons[SZZ] = jac * rhoh*u0*l[2] * dV;
    cons[TAU] = jac * (-rhoe*uU*u0 - Pp*(uU*u0+U[0]) - rho*(uU+1)*u0) * dV;

    int q;
    for(q = NUM_C; q < NUM_Q; q++)
        cons[q] = prim[q]*cons[DDD];
}

void getUstar(double *prim, double *Ustar, double *x, double Sk, double Ss, 
                double *n, double *Bpack)
{
    Ustar[DDD] = 0.0;
    Ustar[SRR] = 0.0;
    Ustar[LLL] = 0.0;
    Ustar[SZZ] = 0.0;
    Ustar[TAU] = 0.0;

    int q;
    for(q = NUM_C; q < NUM_Q; q++)
        Ustar[q] = prim[q]*Ustar[DDD];
}

void cons2prim( double * cons , double * prim , double r , double dV ){
   
   double rho = cons[DDD]/dV;
   if( rho < RHO_FLOOR )   rho = RHO_FLOOR;
   double Sr  = cons[SRR]/dV;
   double Sp  = cons[LLL]/dV/r;
   double Sz  = cons[SZZ]/dV;
   double E   = cons[TAU]/dV;
   double om  = get_om( r );
   
   double vr = Sr/rho;
   double vp = Sp/rho;
   double vp_off = vp - om*r;
   double vz = Sz/rho;

   double KE = .5*( Sr*vr + rho*vp_off*vp_off + Sz*vz );
   double rhoe = E-KE;
   double Pp = (gamma_law - 1.)*rhoe;

   if( Pp  < PRE_FLOOR*rho ) Pp = PRE_FLOOR*rho;
   if( isothermal ){
      double cs2 = get_cs2( r );
      Pp = cs2*rho/gamma_law;
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = vp/r;
   prim[UZZ] = vz;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      prim[q] = cons[q]/cons[DDD];
   }

}

void flux(double *prim, double *flux, double *x, double *n)
{
    double r = x[0];
    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double l[3] = {prim[URR], prim[UPP], prim[UZZ]};

    double lapse;
    double shift[3];
    double igam[9];
    double jac;
    double U[4];

    lapse = metric_lapse(x);
    metric_shift(x, shift);
    metric_igam(x, igam);
    jac = metric_jacobian(x) / r;
    frame_U(x, U);

    double w, u0, u2;
    double u[3];
    double igaml[3];
    int i,j;
    for(i=0; i<3; i++)
    {
        igaml[i] = 0.0;
        for(j=0; j<3; j++)
            igaml[i] += igam[3*i+j]*l[j];
    }
    u2 = l[0]*igaml[0] + l[1]*igaml[1] + l[2]*igaml[2];
    w = sqrt(1.0 + u2);
    u0 = w/lapse;
    for(i=0; i<3; i++)
        u[i] = igaml[i] - shift[i]*u0;

    double l0 = -lapse*w + shift[0]*l[0] + shift[1]*l[1] + shift[2]*l[2];
    double uU = U[0]*l0 + U[1]*l[0] + U[2]*l[1] + U[3]*l[2];
    double un = u[0]*n[0] + r*u[1]*n[1] + u[2]*n[2];
    double Un = U[1]*n[0] + r*u[2]*n[1] + u[3]*n[2];
    
    double rhoh = rho + gamma_law/(gamma_law-1.0)*Pp;
    double rhoe = Pp / (gamma_law-1.0);

    flux[DDD] = jac * rho*un;
    flux[SRR] = jac * rhoh*un*l[0];
    flux[LLL] = jac * rhoh*un*l[1];
    flux[SZZ] = jac * rhoh*un*l[2];
    flux[TAU] = jac * (-rhoe*uU*un - Pp*(uU*un+Un) - rho*(uU+1)*un);

    int q;
    for(q = NUM_C; q < NUM_Q; q++)
        flux[q] = prim[q]*flux[DDD];
}

void source(double *prim, double *cons, double *xp, double *xm, double dVdt)
{
    double x[3] = {0.5*(xm[0]+xp[0]), 0.5*(xm[1]+xp[1]), 0.5*(xm[2]+xp[2])};
    double r = x[0];
    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double l[4] = {0.0, prim[URR], prim[UPP], prim[UZZ]};

    int i,mu,nu;
    double lapse;
    double shift[3];
    double igam[9];
    double ig[16];
    double jac;
    double U[4], dU[16];

    lapse = metric_lapse(x);
    metric_shift(x, shift);
    metric_igam(x, igam);
    jac = metric_jacobian(x) / r;
    frame_U(x, U);
    frame_der_U(x, dU);

    double ia2 = 1.0/(lapse*lapse);
    ig[0] = -ia2;
    for(mu=0; mu<3; mu++)
    {
        ig[mu+1] = shift[i]*ia2;
        ig[4*(mu+1)] = ig[mu+1];
        for(nu=0; nu<3; nu++)
            ig[4*(mu+1)+nu+1] = igam[3*mu+nu]-shift[mu]*shift[nu]*ia2;
    }

    double w, u[4], u2;
    double igaml[3];
    igaml[0] = igam[0]*l[1] + igam[1]*l[2] + igam[2]*l[3];
    igaml[1] = igam[3]*l[1] + igam[4]*l[2] + igam[5]*l[3];
    igaml[2] = igam[6]*l[1] + igam[7]*l[2] + igam[8]*l[3];
    u2 = l[1]*igaml[0] + l[2]*igaml[1] + l[3]*igaml[2];
    w = sqrt(1.0 + u2);
    
    u[0] = w/lapse;
    for(i=0; i<3; i++)
        u[i+1] = igaml[i] - shift[i]*u[0];
    l[0] = -lapse*w + shift[0]*l[1] + shift[1]*l[2] + shift[2]*l[3];
    
    double rhoh = rho + gamma_law/(gamma_law-1.0)*Pp;

    double S0, Sk[3];
    for(i=0; i<3; i++)
    {
        Sk[i] = 0.0;
        if(metric_killing(i+1))
            continue;

        double dg[16];
        metric_der_g(x, i+1, dg);
        for(mu=0; mu<4; mu++)
            for(nu=0; nu<4; nu++)
                Sk[i] += (rhoh*u[mu]*u[nu]+ig[4*mu+nu]*Pp)*dg[4*mu+nu];
        Sk[i] *= 0.5;
    }
    S0 = -U[1]*Sk[0] - U[2]*Sk[1] - U[3]*Sk[2];

    for(mu=1; mu<4; mu++)
        for(nu=0; nu<4; nu++)
        {
            if(mu == nu)
                S0 += -(rhoh*u[mu]*l[nu] + Pp) * dU[4*mu+nu];
            else
                S0 += -(rhoh*u[mu]*l[nu])*dU[4*mu+nu];
        }

    cons[SRR] += jac * Sk[0] * dVdt;
    cons[LLL] += jac * Sk[1] * dVdt;
    cons[SZZ] += jac * Sk[2] * dVdt;
    cons[TAU] += jac * S0 * dVdt;
}

void visc_flux(double *prim, double *gprim, double *flux, double *x, 
                double *n){}

void vel(double *prim1, double *prim2, double *Sl, double *Sr, double *Ss, 
            double *n, double *x, double *Bpack)
{

    double r = x[0];
   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];
   double vn1  = prim1[URR]*n[0] + prim1[UPP]*n[1]*r + prim1[UZZ]*n[2];

   double cs1 = sqrt(gamma_law*P1/rho1);

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];
   double vn2  = prim2[URR]*n[0] + prim2[UPP]*n[1]*r + prim2[UZZ]*n[2];

   double cs2 = sqrt(gamma_law*P2/rho2);

   *Ss = ( P2 - P1 + rho1*vn1*(-cs1) - rho2*vn2*cs2 )/( rho1*(-cs1) - rho2*cs2 );

   *Sr =  cs1 + vn1;
   *Sl = -cs1 + vn1;

   if( *Sr <  cs2 + vn2 ) *Sr =  cs2 + vn2;
   if( *Sl > -cs2 + vn2 ) *Sl = -cs2 + vn2;

}

double mindt(double *prim, double w, double *xp, double *xm)
{
   double r = .5*(xp[0]+xm[0]);
   double Pp  = prim[PPP];
   double rho = prim[RHO];
   double vp  = (prim[UPP]-w)*r;
   double vr  = prim[URR];
   double vz  = prim[UZZ];
   double cs  = sqrt(gamma_law*Pp/rho);

   double maxvr = cs + fabs(vr);
   double maxvp = cs + fabs(vp);
   double maxvz = cs + fabs(vz);

   double dtr = get_dL(xp,xm,0)/maxvr;
   double dtp = get_dL(xp,xm,1)/maxvp;
   double dtz = get_dL(xp,xm,2)/maxvz;
   
   double dt = dtr;
   if( dt > dtp ) dt = dtp;
   if( dt > dtz ) dt = dtz;

   return( dt );
}

double getReynolds(double *prim, double w, double *x, double dx)
{
    return 0.0;
}


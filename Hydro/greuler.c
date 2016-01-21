
#include "../paul.h"
#include "metric.h"
#include "frame.h"

//Global Functions
double get_cs2( double );
double get_dp( double , double );
double get_dL( double * , double * , int );

//Local Functions
void cons2prim_prep(double *cons, double *x);
void cons2prim_solve_isothermal(double *cons, double *prim, double *x);
void cons2prim_solve_adiabatic(double *cons, double *prim, double *x);
void cons2prim_finalize(double *prim, double *x);

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

void cons2prim(double *cons, double *prim, double *x, double dV)
{
    int q;
    double cons1[NUM_Q];
    for(q=0; q<NUM_Q; q++)
        cons1[NUM_Q] = cons[NUM_Q]/dV;

    cons2prim_prep(cons1, x);
    if(isothermal)
        cons2prim_solve_isothermal(cons1, prim, x);
    else
        cons2prim_solve_adiabatic(cons1, prim, x);
    cons2prim_finalize(prim, x);


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
        ig[mu+1] = shift[mu]*ia2;
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
    double rho1 = prim1[RHO];
    double P1   = prim1[PPP];
    double l1[3]  = {prim1[URR], prim1[UPP], prim1[UZZ]};

    double cs21 = gamma_law*P1/(rho1+gamma_law/(gamma_law-1.0)*P1);

    double rho2 = prim2[RHO];
    double P2   = prim2[PPP];
    double l2[3]  = {prim2[URR], prim2[UPP], prim2[UZZ]};

    double cs22 = gamma_law*P2/(rho2+gamma_law/(gamma_law-1.0)*P2);

    //*Ss = ( P2 - P1 + rho1*vn1*(-cs1) - rho2*vn2*cs2 )/( rho1*(-cs1) - rho2*cs2 );
    
    double a, b[3], igam[9];
    a = metric_lapse(x);
    metric_shift(x, b);
    metric_igam(x, igam);
    double bn = b[0]*n[0] + r*b[1]*n[1] + b[2]*n[2];

    int i,j;
    double u21 = 0.0;
    double u22 = 0.0;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
        {
            u21 += igam[3*i+j]*l1[i]*l1[j];
            u22 += igam[3*i+j]*l2[i]*l2[j];
        }
    double w1 = sqrt(1.0+u21);
    double w2 = sqrt(1.0+u22);
    double v1[3], v2[3], vn1, vn2;
    for(i=0; i<3; i++)
    {
        v1[i] = -b[i];
        v2[i] = -b[i];
        for(j=0; j<3; j++)
        {
            v1[i] += igam[3*i+j]*l1[j]*a/w1;
            v2[i] += igam[3*i+j]*l2[j]*a/w2;
        }
    }
    vn1 = n[0]*v1[0] + r*n[1]*v1[1] + n[2]*v1[2];
    vn2 = n[0]*v2[0] + r*n[1]*v2[1] + n[2]*v2[2];

    double sig1 = cs21/(w1*w1*(1.0-cs21));
    double sig2 = cs22/(w2*w2*(1.0-cs22));
    double dv1 = sqrt(sig1*(1.0+sig1)*a*a*igam[

    double vn1 = 

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

void cons2prim_prep(double *cons, double *x)
{

}

void cons2prim_solve_isothermal(double *cons, double *prim, double *x)
{

    

    int q;
    for( q=NUM_C ; q<NUM_Q ; ++q )
        prim[q] = cons[q]/cons[DDD];
}

void cons2prim_solve_adiabatic(double *cons, double *prim, double *x)
{
    double prec = 1.0e6;
    double max_iter = 100;

    double r = x[0];

    double D = cons[DDD];
    double S[3] = {cons[SRR], cons[LLL], cons[SZZ]};
    double tau = cons[TAU];

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

    double s2 = 0.0;
    double Us = 0.0;
    double e = (tau/D + Us + 1.0) / (lapse*U[0]);
    double n = (gamma_law-1.0)/gamma_law;

    int i,j;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            s2 += igam[3*i+j]*S[i]*S[j];
    s2 /= D*D;

    for(i=0; i<3; i++)
        Us += S[i]*(shift[i]*U[0] + U[i+1]);
    Us /= D;

    double wmo;
    if(s2 == 0.0)
        wmo = 0.0;
    else
    {
        //Newton-Raphson on a quartic polynomial to solve for w-1
        //TODO: Minimize truncation error.
        double c[5];
        c[0] = -s2*(n-1)*(n-1);
        c[1] = 2*((e-n)*(e-n)+2*(n-1)*s2);
        c[2] = (5*e-n)*(e-n)-2*(3-n)*s2;
        c[3] = 4*(e*e-s2) - 2*n*e;
        c[4] = e*e-s2;

        //Bounds
        double wmomin = 0.0; //u=0
        double wmomax = sqrt(1.0+s2)-1.0; //eps = P = 0

        //Initial guess: previous w
        double u2 = 0.0;
        double l[3] = {prim[URR], prim[UPP], prim[UZZ]};
        for(i=0; i<3; i++)
            for(j=0; j<3; j++)
                u2 += igam[3*i+j]*l[i]*l[j];
        double wmo0 = u2 / (1.0+sqrt(1.0+u2)); // sqrt(1+u2)-1

        //Run Newton-Raphson
        double wmo1 = wmo0;
        i = 0;
        do
        {
            //TODO: Telescoping evaluation.
            wmo = wmo1;
            double f = c[0] + c[1]*wmo + c[2]*wmo*wmo + c[3]*wmo*wmo*wmo
                        + c[4]*wmo*wmo*wmo*wmo;
            double df = c[1] + 2*c[2]*wmo + 3*c[3]*wmo*wmo
                        + 4*c[4]*wmo*wmo*wmo;
            wmo1 = wmo - f/df;

            if(f > 0.0 && wmo<wmomax)
                wmomax = wmo;
            else if(f < 0.0 && wmo > wmomin)
                wmomin = wmo;
            if(wmo1 < wmomin || wmo1 > wmomax)
                wmo1 = 0.5*(wmomin+wmomax);
            i++;
        }
        while(fabs(wmo-wmo1) > prec && i < max_iter);

        if(i == max_iter)
            printf("ERROR: NR failed to converge\n");
    }

    //Prim recovery
    double w = wmo + 1.0;
    double u0 = w/lapse;
    double hmo = w*(e-w) / (w*w-n);

    double rho = D / (jac*u0);
    if(rho < RHO_FLOOR)
        rho = RHO_FLOOR;
    double Pp = n * rho * hmo;
    if(Pp < PRE_FLOOR*rho)
        Pp = PRE_FLOOR*rho;
    
    double h = 1.0 + gamma_law/(gamma_law-1.0) * Pp/rho;
    double l[3] = {S[0]/(D*h), S[1]/(D*h), S[2]/(D*h)};

    prim[RHO] = rho;
    prim[URR] = l[0];
    prim[UPP] = l[1];
    prim[UZZ] = l[2];
    prim[PPP] = Pp;

    int q;
    for( q=NUM_C ; q<NUM_Q ; ++q )
        prim[q] = cons[q]/cons[DDD];
}

void cons2prim_finalize(double *prim, double *x)
{

}

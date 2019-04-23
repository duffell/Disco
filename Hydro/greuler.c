
#include "../paul.h"
#include "metric.h"
#include "frame.h"

#define DEBUG 0

//Global Functions
double get_cs2( double * );
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
    double r = x[0];
    double rho = prim[RHO];
    double Pp = prim[PPP];
    double l[3] = {prim[URR], prim[UPP], prim[UZZ]};

    double a, b[3], igam[9], jac, U[4];
    a = metric_lapse(x);
    metric_shift(x, b);
    metric_igam(x, igam);
    jac = metric_jacobian(x) / r;
    frame_U(x, U);

    double bn = n[0]*b[0] + n[1]*b[1] + n[2]*b[2];
    double ign = n[0]*igam[0] + n[1]*igam[4] + n[2]*igam[8];
    double Un = n[0]*U[1] + n[1]*U[2] + n[2]*U[3];
    double hn = n[0] + r*n[1] + n[2];

    double ss = Ss/hn;
    double sk = Sk/hn;

    double uS[3];
    double u2 = 0.0;
    int i,j;
    for(i=0; i<3; i++)
    {
        uS[i] = 0.0;;
        for(j=0; j<3; j++)
            uS[i] += igam[3*i+j]*l[j];
        u2 += l[i]*uS[i];
    }

    double w = sqrt(1+u2);
    double u0 = w/a;
    double l0 = -a*w + b[0]*l[0] + b[1]*l[1] + b[2]*l[2];
    double vn = a * (n[0]*uS[0] + n[1]*uS[1] + n[2]*uS[2]) / w - bn;
    double uU = l0*U[0] + l[0]*U[1] + l[1]*U[2] + l[2]*U[3];

    // q == F - s * U
    double rhoh = rho + gamma_law/(gamma_law-1.0) * Pp;
    double qE = rhoh*w*w*(vn-sk) + Pp*(sk+bn);

    double mn = rhoh*w * (n[0]*uS[0] + n[1]*uS[1] + n[2]*uS[2]);
    double qMn = mn*(vn-sk) + a*Pp*ign;

    double ssS = (ss+bn)/a;
    double skS = (sk+bn)/a;

    // P star!
    double Pstar = (qMn - ssS*qE) / (a*(ign - ssS*skS));

    double kappa = (vn - sk) / (ss - sk);
    double alpha1 = (Pp - Pstar) / (ss - sk);
    double alpha2 = (vn*Pp - ss*Pstar) / (ss - sk);

    double rhoe = Pp / (gamma_law - 1.0);
    double tau = -rhoe*uU*u0 - Pp*(u0*uU+U[0]) - rho*(uU+1.0)*u0;

    Ustar[DDD] = jac * rho*u0 * kappa;
    Ustar[SRR] = jac * (rhoh*u0*l[0] * kappa + alpha1 * n[0]);
    Ustar[LLL] = jac * (rhoh*u0*l[1] * kappa + alpha1 * n[1]);
    Ustar[SZZ] = jac * (rhoh*u0*l[2] * kappa + alpha1 * n[2]);
    Ustar[TAU] = jac * (tau * kappa - Un*alpha1 + U[0]*alpha2);

    int q;
    for(q = NUM_C; q < NUM_Q; q++)
        Ustar[q] = prim[q]*Ustar[DDD];
}

void cons2prim(double *cons, double *prim, double *x, double dV)
{
    int q;
    double cons1[NUM_Q];
    for(q=0; q<NUM_Q; q++)
        cons1[q] = cons[q]/dV;

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
    //double un = u[0]*n[0] + r*u[1]*n[1] + u[2]*n[2];
    //double Un = U[1]*n[0] + r*U[2]*n[1] + U[3]*n[2];
    double un = u[0]*n[0] + u[1]*n[1] + u[2]*n[2];
    double Un = U[1]*n[0] + U[2]*n[1] + U[3]*n[2];
    double hn = n[0] + r*n[1] + n[2];
    
    double rhoh = rho + gamma_law/(gamma_law-1.0)*Pp;
    double rhoe = Pp / (gamma_law-1.0);

    //flux[DDD] = jac * rho*un;
    //flux[SRR] = jac * (rhoh*un*l[0] + Pp*n[0]);
    //flux[LLL] = jac * (rhoh*un*l[1] + Pp*n[1]*r);
    //flux[LLL] = jac * (rhoh*un*l[1] + Pp*n[1]);
    //flux[SZZ] = jac * (rhoh*un*l[2] + Pp*n[2]);
    //flux[TAU] = jac * (-rhoe*uU*un - Pp*(uU*un+Un) - rho*(uU+1)*un);

    flux[DDD] = jac * hn * rho*un;
    flux[SRR] = jac * hn * (rhoh*un*l[0] + Pp*n[0]);
    flux[LLL] = jac * hn * (rhoh*un*l[1] + Pp*n[1]);
    flux[SZZ] = jac * hn * (rhoh*un*l[2] + Pp*n[2]);
    flux[TAU] = jac * hn * (-rhoe*uU*un - Pp*(uU*un+Un) - rho*(uU+1)*un);
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

    double a, b[3], gam[9], igam[9];
    a = metric_lapse(x);
    metric_shift(x, b);
    metric_gam(x, gam);
    metric_igam(x, igam);

    int i,j;
    double u21 = 0.0;
    double u22 = 0.0;
    double uS1[3], uS2[3];
    for(i=0; i<3; i++)
    {
        uS1[i] = 0.0;
        uS2[i] = 0.0;
        for(j=0; j<3; j++)
        {
            uS1[i] += igam[3*i+j]*l1[j];
            uS2[i] += igam[3*i+j]*l2[j];
        }
        u21 += uS1[i]*l1[i];
        u22 += uS2[i]*l2[i];
    }

    double w1 = sqrt(1.0+u21);
    double w2 = sqrt(1.0+u22);
    double v21 = u21/(w1*w1);
    double v22 = u22/(w2*w2);

    //TODO: Use n[] PROPERLY.  This only works if n = (1,0,0) or some 
    //      permutation.
    double vSn1 = (uS1[0]*n[0]+uS1[1]*n[1]+uS1[2]*n[2]) / w1;
    double vSn2 = (uS2[0]*n[0]+uS2[1]*n[1]+uS2[2]*n[2]) / w2;
    double bn = (b[0]*n[0]+b[1]*n[1]+b[2]*n[2]);
    double ign = igam[3*0+0]*n[0] + igam[3*1+1]*n[1] + igam[3*2+2]*n[2];

    double dv1 = sqrt(cs21*(ign - vSn1*vSn1 - cs21*(ign*v21-vSn1*vSn1))) / w1;
    double dv2 = sqrt(cs22*(ign - vSn2*vSn2 - cs22*(ign*v22-vSn2*vSn2))) / w2;
    double hn = n[0] + r*n[1] + n[2];

    double sl1 = hn * (a * (vSn1*(1.0-cs21) - dv1) / (1.0-v21*cs21) - bn);
    double sr1 = hn * (a * (vSn1*(1.0-cs21) + dv1) / (1.0-v21*cs21) - bn);
    double sl2 = hn * (a * (vSn2*(1.0-cs22) - dv2) / (1.0-v22*cs22) - bn);
    double sr2 = hn * (a * (vSn2*(1.0-cs22) + dv2) / (1.0-v22*cs22) - bn);

/*
    printf("cs2(L/R): %.12lg %.12lg\n", cs21, cs22);
    printf("vSn(L/R): %.12lg %.12lg\n", vSn1, vSn2);
    printf("w  (L/R): %.12lg %.12lg\n", w1, w2);
    printf("u2 (L/R): %.12lg %.12lg\n", u21, u22);
    printf("v2 (L/R): %.12lg %.12lg\n", v21, v22);
    printf("sl (L/R): %.12lg %.12lg\n", sl1, sl2);
    printf("sr (L/R): %.12lg %.12lg\n", sr1, sr2);
*/

    *Sr = sr1 > sr2 ? sr1 : sr2;
    *Sl = sl1 < sl2 ? sl1 : sl2;

    //double maxv = fabs(*Sr) > fabs(*Sl) ? fabs(*Sr) : fabs(*Sl);
    //*Sr = maxv;
    //*Sl = -maxv;

    //Now for the contact wave speed.
    double sL = *Sl / hn;
    double sR = *Sr / hn;

    double rhohL = rho1 + gamma_law/(gamma_law-1)*P1;
    double rhohR = rho2 + gamma_law/(gamma_law-1)*P2;
    double vnL = a*vSn1 - bn;
    double vnR = a*vSn2 - bn;

    double ML[3] = {rhohL*w1*l1[0], rhohL*w1*l1[1], rhohL*w1*l1[2]};
    double MR[3] = {rhohR*w2*l2[0], rhohR*w2*l2[1], rhohR*w2*l2[2]};
    double EL = rhohL*w1*w1-P1;
    double ER = rhohR*w2*w2-P2;

    double FML[3] = {ML[0]*vnL+n[0]*a*P1, ML[1]*vnL+n[1]*a*P1, 
                        ML[2]*vnL+n[2]*a*P1};
    double FMR[3] = {MR[0]*vnR+n[0]*a*P2, MR[1]*vnR+n[1]*a*P2, 
                        MR[2]*vnR+n[2]*a*P2};
    double FEL = EL*vnL + P1*(vnL+bn);
    double FER = ER*vnR + P2*(vnR+bn);

    double UE = (sR*ER - sL*EL - FER + FEL) / (sR - sL);
    double FE = ((sR*FEL - sL*FER + sL*sR*(ER-EL)) / (sR - sL) - bn*UE) / a;

    double UM_hll[3], FM_hll[3];
    for(i=0; i<3; i++)
    {
        UM_hll[i] = (sR*MR[i]-sL*ML[i]-FMR[i]+FML[i]) / (sR-sL);
        FM_hll[i] = ((sR*FML[i]-sL*FMR[i]+sL*sR*(MR[i]-ML[i])) / (sR-sL)
                        - bn*UM_hll[i]) / a;
    }
    double UM = 0.0;
    double FM = 0.0;
    for(i=0; i<3; i++)
    {
        double igi = n[0]*igam[3*i] + n[1]*igam[3*i+1] + n[2]*igam[3*i+2];
        UM += igi * UM_hll[i];
        FM += igi * FM_hll[i];
    }

    double A = FE;
    double B = -FM-ign*UE;
    double C = ign*UM;

    double sS;
    if(fabs(4*A*C/(B*B)) < 1.0e-7)
        sS = -C/B * (1.0 + A*C/(B*B) + 2*A*A*C*C/(B*B*B*B));
    else
        sS = (-B - sqrt(B*B-4*A*C)) / (2*A);

    *Ss = hn * (a*sS - bn);
}

double mindt(double *prim, double wc, double *xp, double *xm)
{
    double x[3] = {0.5*(xm[0]+xp[0]), 0.5*(xm[1]+xp[1]), 0.5*(xm[2]+xp[2])};
    double r = x[0];
    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double l[3] = {prim[URR], prim[UPP], prim[UZZ]};
    double cs  = sqrt(gamma_law*Pp/(rho+gamma_law/(gamma_law-1)*Pp));

    double a, b[3], gam[9], igam[9];
    a = metric_lapse(x);
    metric_shift(x, b);
    metric_gam(x, gam);
    metric_igam(x, igam);

    int i,j;
    double uS[3];
    double u2 = 0.0;

    for(i=0; i<3; i++)
    {
        uS[i] = 0.0;
        for(j=0; j<3; j++)
            uS[i] += igam[3*i+j]*l[j];
        u2 += uS[i]*l[i];
    }
    double w = sqrt(1.0+u2);

    double v2, vS[3];
    for(i=0; i<3; i++)
        vS[i] = uS[i]/w;
    v2 = u2/(w*w);

    double sig = 1-cs*cs;

    double dvr = cs * sqrt(igam[0]*(1-cs*cs*v2) - sig*vS[0]*vS[0]) / w;
    double vrl = fabs(a * (vS[0]*sig - dvr) / (1-v2*cs*cs) - b[0]);
    double vrr = fabs(a * (vS[0]*sig + dvr) / (1-v2*cs*cs) - b[0]);

    double dvp = cs * sqrt(igam[4]*(1-cs*cs*v2) - sig*vS[1]*vS[1]) / w;
    double vpl = fabs(r * (a * (vS[1]*sig - dvp) / (1-v2*cs*cs) - b[1] - wc));
    double vpr = fabs(r * (a * (vS[1]*sig + dvp) / (1-v2*cs*cs) - b[1] - wc));
    
    double dvz = cs * sqrt(igam[8]*(1-cs*cs*v2) - sig*vS[2]*vS[2]) / w;
    double vzl = fabs(a * (vS[2]*sig - dvz) / (1-v2*cs*cs) - b[2]);
    double vzr = fabs(a * (vS[2]*sig + dvz) / (1-v2*cs*cs) - b[2]);

    double maxvr = vrr > vrl ? vrr : vrl;
    double maxvp = vpr > vpl ? vpr : vpl;
    double maxvz = vzr > vzl ? vzr : vzl;

    double dtr = get_dL(xp,xm,1)/maxvr;
    double dtp = get_dL(xp,xm,0)/maxvp;
    double dtz = get_dL(xp,xm,2)/maxvz;

    double dt = dtr;
    dt = dt < dtp ? dt : dtp;
    dt = dt < dtz ? dt : dtz;

    return dt;
}

double getReynolds(double *prim, double w, double *x, double dx)
{
    return 0.0;
}

void cons2prim_prep(double *cons, double *x)
{
    //TODO: complete this.
}

void cons2prim_solve_isothermal(double *cons, double *prim, double *x)
{
    //TODO: complete this.
    int q;
    for( q=NUM_C ; q<NUM_Q ; ++q )
        prim[q] = cons[q]/cons[DDD];
}

void cons2prim_solve_adiabatic(double *cons, double *prim, double *x)
{
    double prec = 1.0e-15;
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

    int i,j;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            s2 += igam[3*i+j]*S[i]*S[j];
    s2 /= D*D;

    for(i=0; i<3; i++)
        Us += S[i]*(shift[i]*U[0] + U[i+1]);
    Us /= D;
    
    double e = (tau/D + Us + 1.0) / (lapse*U[0]);
    double n = (gamma_law-1.0)/gamma_law;

    if(e*e < s2 && DEBUG)
    {
        printf("Not enough thermal energy (r=%.12lg, e2=%.12lg, s2=%.12lg)\n",
                r, e*e, s2);

        double cons0[NUM_Q];
        prim2cons(prim, cons0, x, 1.0);

        printf("prim: %.16lg %.16lg %.16lg %.16lg %.16lg\n",
                prim[RHO], prim[PPP], prim[URR], prim[UPP], prim[UZZ]);
        printf("cons0: %.16lg %.16lg %.16lg %.16lg %.16lg\n",
                cons0[DDD], cons0[TAU], cons0[SRR], cons0[LLL], cons0[SZZ]);
        printf("cons: %.16lg %.16lg %.16lg %.16lg %.16lg\n",
                cons[DDD], cons[TAU], cons[SRR], cons[LLL], cons[SZZ]);
    }

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
        while(fabs((wmo-wmo1)/(wmo+1.0)) > prec && i < max_iter);

        if(i == max_iter)
            printf("ERROR: NR failed to converge: err = %.12lg\n", fabs((wmo-wmo1)/(wmo+1.0)));
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
    
    if(e*e < s2 && DEBUG)
    {
        double cons1[NUM_Q];
        prim2cons(prim, cons1, x, 1.0);

        printf("prim1: %.16lg %.16lg %.16lg %.16lg %.16lg\n",
                prim[RHO], prim[PPP], prim[URR], prim[UPP], prim[UZZ]);
        printf("cons1: %.16lg %.16lg %.16lg %.16lg %.16lg\n",
                cons1[DDD], cons1[TAU], cons1[SRR], cons1[LLL], cons1[SZZ]);
    }
}

void cons2prim_finalize(double *prim, double *x)
{

}

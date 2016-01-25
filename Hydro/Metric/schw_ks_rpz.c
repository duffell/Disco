#include "../../paul.h"
#include "../metric.h"

/*
 * This is the Schwarzschild metric in Kerr-Schild cylindrical coordinates 
 * (t, r, phi, z). These are related to the usual spherical Schwarzschild 
 * coordinates (t_sc, R, theta, phi) by: t = t_sc + 2M log|R/2M-1|, 
 * r = R sin(theta), phi = phi, and z = R cos(theta).
 */

static double om = 0.0; 
static double M = 1.0; 

void setMetricParams(struct domain *theDomain)
{
   //om = theDomain->theParList.MetricPar1;
   //M = theDomain->theParList.MetricPar2;
}

double metric_lapse(double x[3])
{
    double R = sqrt(x[0]*x[0] + x[2]*x[2]);
    return 1.0/sqrt(1.0 + 2*M/R);
}

void metric_shift(double x[3], double b[3])
{
    double r = x[0];
    double z = x[2];
    double R = sqrt(r*r+z*z);
    double a2 = 1.0 / (1.0 + 2*M/R);

    b[0] = 2*M*r/(R*R) * a2;
    b[1] = 0.0;
    b[2] = 2*M*z/(R*R) * a2;
}

void metric_gam(double x[3], double gam[9])
{
    double r = x[0];
    double z = x[2];
    double R = sqrt(r*r+z*z);

    double H = 2*M/R;
    double st = r/R;
    double ct = z/R;

    gam[0] = 1.0 + H*st*st;
    gam[1] = 0.0;
    gam[2] = H*st*ct;
    gam[3] = 0.0;
    gam[4] = r*r;
    gam[5] = 0.0;
    gam[6] = H*st*ct;
    gam[7] = 0.0;
    gam[8] = 1.0 + H*ct*ct;
}

void metric_igam(double x[3], double igam[9])
{
    double r = x[0];
    double z = x[2];
    double R = sqrt(r*r+z*z);

    double H = 2*M/R;
    double st = r/R;
    double ct = z/R;
    double a2 = 1.0/(1.0 + H);

    igam[0] = (1.0 + H*ct*ct) * a2;
    igam[1] = 0.0;
    igam[2] = -H*st*ct*a2;
    igam[3] = 0.0;
    igam[4] = 1.0/(r*r);
    igam[5] = 0.0;
    igam[6] = -H*st*ct*a2;
    igam[7] = 0.0;
    igam[8] = (1.0 + H*st*st) * a2;
}

double metric_jacobian(double x[3])
{
    return x[0];
}

void metric_der_g(double x[3], int i, double dg[16])
{
    int mu;
    for(mu=0; mu<16; mu++)
        dg[mu] = 0.0;

    if(i != 1 && i != 3)
        return;

    double r = x[0];
    double z = x[2];
    double R = sqrt(r*r+z*z);

    double st = r/R;
    double ct = z/R;

    double MoR2 = M/(R*R);
    double MoR3 = M/(R*R*R);
    double MoR4 = M/(R*R*R*R);

    if(i == 1)
    {
        dg[0]  = -2*MoR2 * st;               // 00
        dg[1]  =  2*MoR2 * (1.0-2*st*st);    // 01
        dg[3]  = -4*MoR3*z * st;          // 03
        dg[5]  =  2*MoR3*r * (2.0-3.0*st*st);// 11
        dg[7]  =  2*MoR3*z * (1.0-3.0*st*st);// 13
        dg[10] =  2*r;                       // 22
        dg[15] = -6*MoR4*z*z * st;           // 33
        dg[4]  = dg[1];
        dg[12] = dg[3];
        dg[13] = dg[7];
    }
    else if(i == 3)
    {
        dg[0]  = -2*MoR2 * ct;               // 00
        dg[1]  = -4*MoR3*r * ct*ct;          // 01
        dg[3]  =  2*MoR2 * (1.0-2*ct*ct);    // 03
        dg[5]  = -6*MoR4*r*r * ct;           // 11
        dg[7]  =  2*MoR3*r * (1.0-3.0*ct*ct);// 13
        dg[15] =  2*MoR3*z * (2.0-3.0*ct*ct);// 33
        dg[4]  = dg[1];
        dg[12] = dg[3];
        dg[13] = dg[7];
    }
}

void metric_der_lapse(double x[3], double da[4])
{
    double r = x[0];
    double z = x[2];
    double R = sqrt(r*r+z*z);

    double st = r/R;
    double ct = z/R;

    double MoR2 = M/(R*R);
    double a = 1.0/sqrt(1.0+2*M/R);

    da[0] = 0.0;
    da[1] = a*a*a * MoR2 * st;
    da[2] = 0.0;
    da[3] = a*a*a * MoR2 * ct;
}

void metric_der_shift(double x[3], double db[12])
{
    double r = x[0];
    double z = x[2];
    double R = sqrt(r*r+z*z);

    double st = r/R;
    double ct = z/R;

    double denom = 1.0/(R*R + 2*M*R);
    double d_denom = -denom*denom * 2*(R+M);

    int i;

    //dt
    for(i=0; i<3; i++)
        db[i] = 0;

    //dr
    db[3] = 2*M * denom + 2*M*r * d_denom*st;
    db[4] = 0.0;
    db[5] = 2*M*z * d_denom*st;

    //dphi
    for(i=0; i<3; i++)
        db[6+i] = 0;

    //dz
    db[9] = 2*M*r * d_denom*ct;
    db[10] = 0.0;
    db[11] = 2*M * denom + 2*M*z * d_denom*ct;
}

int metric_killing(int mu)
{
    if(mu == 1 || mu == 3)
        return 0;
    return 1;
}

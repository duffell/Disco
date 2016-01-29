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
    return sqrt(1.0 - 2*M/R);
}

void metric_shift(double x[3], double b[3])
{
    b[0] = 0.0;
    b[1] = 0.0;
    b[2] = 0.0;
}

void metric_gam(double x[3], double gam[9])
{
    double r = x[0];
    double z = x[2];
    double R = sqrt(r*r+z*z);

    gam[0] = r*r/(R*R-2*M*R) + z*z/(R*R);
    gam[1] = 0.0;
    gam[2] = 2*M*r*z/(R*R*(R-2*M));
    gam[3] = 0.0;
    gam[4] = r*r;
    gam[5] = 0.0;
    gam[6] = 2*M*r*z/(R*R*(R-2*M));
    gam[7] = 0.0;
    gam[8] = r*r/(R*R) + z*z/(R*R-2*M*R);
}

void metric_igam(double x[3], double igam[9])
{
    double r = x[0];
    double z = x[2];
    double R = sqrt(r*r+z*z);

    igam[0] = (r*r*(R-2*M) + z*z*R) / (R*R*R);
    igam[1] = 0.0;
    igam[2] = -2*M*r*z/(R*R*R);
    igam[3] = 0.0;
    igam[4] = 1.0/(r*r);
    igam[5] = 0.0;
    igam[6] = -2*M*r*z/(R*R*R);
    igam[7] = 0.0;
    igam[8] = (r*r*R + z*z*(R-2*M)) / (R*R*R);
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

    if(i == 1)
    {
        dg[0]  = -2*M*r/(R*R*R);             // 00
        dg[5]  = -2*M*r*(r*r*R-2*z*z*(R-2*M)) / (R*R*R*R*(R-2*M)*(R-2*M));//11
        dg[7] = 2*M*z*(R*R*(R-2*M)-r*r*(3*R-4*M)) / (R*R*R*R*(R-2*M)*(R-2*M));
        dg[10] =  2*r;                       // 22
        dg[15]  = -2*M*r*z*z*(3*R-4*M) / (R*R*R*R*(R-2*M)*(R-2*M)); //33
        dg[13] = dg[7];
    }
    else if(i == 3)
    {
        dg[0]  = -2*M*z/(R*R*R);             // 00
        dg[5]  = -2*M*r*r*z*(3*R-4*M) / (R*R*R*R*(R-2*M)*(R-2*M)); //33
        dg[7] = 2*M*r*(R*R*(R-2*M)-z*z*(3*R-4*M)) / (R*R*R*R*(R-2*M)*(R-2*M));
        dg[15]  = -2*M*z*(z*z*R-2*r*r*(R-2*M)) / (R*R*R*R*(R-2*M)*(R-2*M));//11
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
    double ia = 1.0/sqrt(1.0-2*M/R);

    da[0] = 0.0;
    da[1] = ia * MoR2 * st;
    da[2] = 0.0;
    da[3] = ia * MoR2 * ct;
}

void metric_der_shift(double x[3], double db[12])
{
    int i;
    for(i=0; i<12; i++)
        db[i] = 0;
}

int metric_killing(int mu)
{
    if(mu == 1 || mu == 3)
        return 0;
    return 1;
}

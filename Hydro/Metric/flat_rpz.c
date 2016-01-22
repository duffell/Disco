#include "../metric.h"

static double om = 0.0; 

void setMetricParams(struct domain *theDomain)
{
   //om = theDomain->theParList.MetricPar1;
}

double metric_lapse(double x[3])
{
    return 1.0;
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
    gam[0] = 1.0;
    gam[1] = 0.0;
    gam[2] = 0.0;
    gam[3] = 0.0;
    gam[4] = r*r;
    gam[5] = 0.0;
    gam[6] = 0.0;
    gam[7] = 0.0;
    gam[8] = 1.0;
}

void metric_igam(double x[3], double igam[9])
{
    double r = x[0];
    igam[0] = 1.0;
    igam[1] = 0.0;
    igam[2] = 0.0;
    igam[3] = 0.0;
    igam[4] = 1.0/(r*r);
    igam[5] = 0.0;
    igam[6] = 0.0;
    igam[7] = 0.0;
    igam[8] = 1.0;
}

double metric_jacobian(double x[3])
{
    return x[0];
}

void metric_der_g(double x[3], int i, double dg[16])
{
    double r = x[0];
    int mu;
    for(mu=0; mu<16; mu++)
        dg[mu] = 0.0;

    if(i == 1)
        dg[4*2+2] = 2*r;
}

void metric_der_lapse(double x[3], double da[4])
{
    int i;
    for(i=0; i<4; i++)
        da[i] = 0.0;
}

void metric_der_shift(double x[3], double db[12])
{
    int i;
    for(i=0; i<12; i++)
        db[i] = 0.0;
}

int metric_killing(int mu)
{
    if(mu == 1)
        return 0;
    return 1;
}

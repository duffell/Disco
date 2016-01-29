#include "../metric.h"
#include "../frame.h"

static double M = 1.0;

void frame_U(double x[3], double U[4])
{
    double r = x[0];
    double Risco = 6.0*M;

    double U0, Ur, Up;

    if(r > 6.0*M)
    {
        U0 = sqrt(r/(r-3*M));
        Ur = 0.0;
        Up = sqrt(M/(r*r*(r-3*M)));
    }
    else
    {
        double x = Risco/r - 1.0;
        U0 = 2*(sqrt(2.0)*r - M*sqrt(x*x*x)) / (3*(r-2*M));
        //U0 = -2*sqrt(2.0)/(3*(-1+2*M/r));
        Ur = -sqrt(x*x*x) / 3.0;
        Up = 2.0*sqrt(3.0) / (r*r);
    }

    U[0] = U0;
    U[1] = Ur;
    U[2] = Up;
    U[3] = 0.0;
}

void frame_der_U(double x[3], double dU[16])
{
    int i;
    for(i=0; i<16; i++)
        dU[i] = 0.0;

    double r = x[0];
    double Risco = 6*M;

    double dU0dr, dUrdr, dUpdr;

    if(r > Risco)
    {
        dU0dr = -1.5*M / sqrt(r*(r-3*M)*(r-3*M)*(r-3*M));
        dUrdr = 0.0;
        dUpdr = -1.5 * (r-2*M) * sqrt(M/((r-3*M)*(r-3*M)*(r-3*M))) / (r*r);
    }
    else
    {
        double x = sqrt(Risco/r - 1.0);
        double y = M/r;

        dU0dr = -2*M * (x*(18*y*y-15*y+1) + 2*sqrt(2.0)) / (3*(r-2*M)*(r-2*M));
        //dU0dr = -2*sqrt(2.0) / (3.0*(-1+2*M/r)*(-1+2*M/r)) * (-2*M/(r*r));
        dUrdr = 0.5 * x * Risco/(r*r);
        dUpdr = -4*sqrt(3.0) * M / (r*r*r);
    }

    dU[4] = dU0dr;
    dU[5] = dUrdr;
    dU[6] = dUpdr;
}


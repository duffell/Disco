#include "../metric.h"
#include "../frame.h"

void frame_U(double x[3], double U[4])
{
    double ia, b[3];
    ia = 1.0/metric_lapse(x);
    metric_shift(x, b);

    U[0] = ia;
    U[1] = -ia*b[0];
    U[2] = -ia*b[1];
    U[3] = -ia*b[2];
}

void frame_der_U(double x[3], double dU[16])
{
    int i,j;
    double ia, da[4], b[3], db[12];
    ia = 1.0/metric_lapse(x);
    metric_der_lapse(x, da);
    metric_shift(x, b);
    metric_der_shift(x, db);

    for(i=0; i<16; i++)
        dU[i] = 0.0;

    for(i=0; i<4; i++)
    {
        if(metric_killing(i))
            continue;

        dU[4*i] = -ia*ia*da[i];
        for(j=0; j<3; j++)
            dU[4*i+j+1] = b[j]*ia*ia*da[i] - db[3*i+j]*ia;
    }
}


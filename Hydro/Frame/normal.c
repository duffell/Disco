#include "../frame.h"

void frame_U(double x[3], double U[4])
{
    U[0] = 1.0;
    U[1] = 0.0;
    U[2] = 0.0;
    U[3] = 0.0;
}

void frame_der_U(double x[3], double dU[16])
{
    int i;
    for(i=0; i<16; i++)
        dU[i] = 0.0;
}


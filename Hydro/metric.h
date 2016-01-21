#ifndef DISCO_METRIC
#define DISCO_METRIC

double metric_lapse(double x[3]);
void metric_shift(double x[3], double b[3]);
void metric_gam(double x[3], double gam[9]);
void metric_igam(double x[3], double igam[9]);
double metric_jacobian(double x[3]);
void metric_der_g(double x[3], int i, double dg[16]);
void metric_der_lapse(double x[3], double da[4]);
void metric_der_shift(double x[3], double db[12]);
int metric_killing(int mu);

#endif

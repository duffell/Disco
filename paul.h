enum{RHO,PPP,URR,UPP,UZZ,BRR,BPP,BZZ};
enum{DDD,TAU,SRR,LLL,SZZ};
enum{C_FIXED,C_WCELL,C_WRIEMANN};

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define MOVE_CELLS C_WCELL

#define NUM_C 8
#define NUM_N 1
#define NUM_Q (NUM_C+NUM_N)
#define NUM_G 2

//Magnetic field tracking things.  Can be set to zero if there is no MHD.
#define NUM_EDGES 8    //0, 4 or 8 
#define NUM_FACES 5    //0, 3 or 5
#define NUM_AZ_EDGES 4 //0, 0 or 4

struct param_list{

   double t_min, t_max;
   int Num_R, Num_Z;
   double aspect;
   int NumRepts, NumSnaps, NumChecks;
   int Out_LogTime;

   double rmin, rmax;
   double zmin, zmax;
   double phimax;

   int LogZoning, Z_Periodic;
   double LogRadius;
   double MaxShort, MaxLong;
   int Mesh_Motion, Riemann_Solver;
   int Absorb_BC, Initial_Regrid, visc_flag, include_atmos;

   double CFL, PLM;
   double Density_Floor, Pressure_Floor;

   double Adiabatic_Index;
   double viscosity;
   int isothermal_flag;

   double Disk_Mach;
   double Mass_Ratio;
   double Eccentricity;
   double Drift_Rate,Drift_Exp;
   int alpha_flag;

   int restart_flag;
   int CT;

};

struct diagnostic_avg{
   double * Qr;
   double t_avg;
};

struct domain{

   struct cell ** theCells;
   struct face * theFaces_1;
   struct face * theFaces_2;
   struct planet * thePlanets;
   int * Np;
   int Nr,Nz,Ng;
   int N_ftracks_r;
   int N_ftracks_z;
   int Npl;
   double * r_jph;
   double * z_kph;
   double phi_max;
   int * fIndex_r;
   int * fIndex_z;

   time_t Wallt_init;
   int rank,size;
   int dim_rank[2];
   int dim_size[2];
   int left_rank[2];
   int right_rank[2];
   MPI_Comm theComm;

   struct param_list theParList;
   int num_tools;
   struct diagnostic_avg theTools;

   double t;
   int count_steps;
   double t_init, t_fin;
   int nrpt;
   int N_rpt;
   int nsnp;
   int N_snp;
   int nchk;
   int N_chk;

   int final_step;
   int check_plz;

};

struct cell{

   double prim[NUM_Q];
   double cons[NUM_Q];
   double RKcons[NUM_Q];
   double grad[NUM_Q];
   double gradp[NUM_Q];
   double piph;
   double dphi;
   double wiph;

   double E[NUM_EDGES];
   double B[NUM_EDGES];
   double E_phi[NUM_AZ_EDGES];
   double    Phi[NUM_FACES];
   double RK_Phi[NUM_FACES];
   double tempDoub;

};

struct edge{
   struct cell * LU;
   struct cell * RU;
   struct cell * LD;
   struct cell * RD;

   int Prim_Zone;
   int Alt_LR;
   int Alt_UD;

   double E_dl;
};

struct face{
   struct cell * L;
   struct cell * R;
   double dxL;
   double dxR;
   double cm[3];
   double dphi;
   double dl;
   double dA;

   double E,B;
   int LRtype;
   int flip_flag;
};

struct planet{
   double r;
   double phi; 
   double M;
   double omega;
   double vr;
   double RK_r;
   double RK_phi;
   double RK_M;
   double RK_omega;
   double RK_vr;

   double eps;
   double Fr;
   double Fp;
};


#include "paul.h"

int between( double phi , double phip , double phim , double phi_max ){
   double dp1 = phi-phim;
   while( dp1 > phi_max ) dp1 -= phi_max;
   while( dp1 < 0.0     ) dp1 += phi_max;
   double dp2 = phip-phim;
   while( dp2 > phi_max ) dp2 -= phi_max;
   while( dp2 < 0.0 ) dp2 += phi_max;

   int between = 0;
   if( dp1 < dp2 ) between = 1;
   return(between);
}

double get_dA( double * , double * , int );
double get_dp( double , double );

int get_num_rzFaces( int Nr , int Nz , int dim ){
   if( dim==1 ) return( (Nr-1)*Nz );
   else return( (Nz-1)*Nr );
}

void addFace( struct face * theFaces , int n , struct cell * cL , struct cell * cR , double dxL , double dxR , double * xp , double * xm , int dim ){
   int d;
   for( d=0 ; d<3 ; ++d ) theFaces[n].cm[d] = .5*(xp[d]+xm[d]); //Consider calculating center of mass in geometry.c
   double dp = get_dp(xp[1],xm[1]);
   double phic = get_dp(xp[1],.5*dp);
   double rp = xp[0];
   double rm = xm[0];
   double r2 = (rp*rp+rm*rm+rp*rm)/3.;
   theFaces[n].cm[0] = r2/(.5*(rp+rm));
   theFaces[n].cm[1] = phic;
   theFaces[n].L   = cL;
   theFaces[n].R   = cR;
   theFaces[n].dxL = dxL;
   theFaces[n].dxR = dxR;
   theFaces[n].dphi= dp;//get_dp(xp[1],xm[1]);
   theFaces[n].dA  = get_dA(xp,xm,dim);
}

void buildfaces( struct domain * theDomain , struct face * theFaces , int * ntj , int dim , int mode ){
 
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   double Pmax = theDomain->phi_max;
   int i,j,k; 
   int n=0;

   int Nrmax = Nr-1;
   if( dim==2 ) Nrmax = Nr;
   int Nzmax = Nz-1;
   if( dim==1 ) Nzmax = Nz;

   for( j=0 ; j<Nrmax ; ++j ){
      int jp = j+1;
      for( k=0 ; k<Nzmax ; ++k ){
         int JK = j+Nrmax*k;
         if( mode == 0 ) ntj[JK] = n;
         int kp = k+1;

         int jk  = j  + Nr*k;
         int jkp = jp + Nr*k;
         if( dim==2 ) jkp = j + Nr*kp;

         double dxL = .5*(r_jph[j]  - r_jph[j-1]);
         double dxR = .5*(r_jph[jp] - r_jph[j]  );
         if( dim==2 ){
            dxL = .5*(z_kph[k]  - z_kph[k-1]);
            dxR = .5*(z_kph[kp] - z_kph[k]  );
         }
         double xp[3] = {r_jph[j],0.0,z_kph[k  ]};
         double xm[3] = {r_jph[j],0.0,z_kph[k-1]};
         if( dim==2 ){ xm[0] = r_jph[j-1] ; xm[2] = z_kph[k] ; }

         int ip=0;
         int found = 0;
         double phi0 = theCells[jk][0].piph-theCells[jk][0].dphi;
         while( !found && ip<Np[jkp] ){
            //Find zone with smallest piph which overlaps with zone 0
            int im = ip-1;
            if( im == -1 ) im = Np[jkp]-1;
            double phip = theCells[jkp][ip].piph;
            double phim = theCells[jkp][im].piph;
            found = between( phi0 , phip , phim , Pmax );
            if( !found ) ++ip;
         }
         if( !found ) ip = 0;
         for( i=0 ; i<Np[jk] ; ++i ){

            struct cell * cL = &(theCells[jk ][i] );
            struct cell * cR = &(theCells[jkp][ip]);
            double dphi = cL->piph - cR->piph;
            while( dphi > .5*Pmax ) dphi -= Pmax;
            while( dphi <-.5*Pmax ) dphi += Pmax;
            //First figure out if cell+ covers all of cell-, if so create one face out of cell-, and move on to next i.
            if( dphi < 0.0 ){
               xp[1] = cL->piph;
               xm[1] = cL->piph-cL->dphi;
               if( mode==1 ) addFace( theFaces , n , cL , cR , dxL , dxR , xp , xm , dim );
               ++n;
            }else{
            //Otherwise, three steps:
               //Step A: face formed out of beginning of cell- and end of cell+. ++ip;
               xp[1] = cR->piph;
               xm[1] = cL->piph-cL->dphi;
               if( mode==1 ) addFace( theFaces , n , cL , cR , dxL , dxR , xp , xm , dim );
               ++n;

               ++ip;
               if( ip == Np[jkp] ) ip = 0;
               cR = &(theCells[jkp][ip]);
               double dphi = cL->piph - cR->piph;
               while( dphi > .5*Pmax ) dphi -= Pmax;
               while( dphi <-.5*Pmax ) dphi += Pmax;
               while( dphi > 0.0 ){
                  //Step B: (optional) all faces formed out of part of cell- and all of cell+. ++ip;
                  xp[1] = cR->piph;
                  xm[1] = cR->piph-cR->dphi;
                  if( mode==1 ) addFace( theFaces , n , cL , cR , dxL , dxR , xp , xm , dim );
                  ++n;
               
                  ++ip;
                  if( ip == Np[jkp] ) ip = 0;
                  cR = &(theCells[jkp][ip]);
                  dphi = cL->piph - cR->piph;
                  while( dphi > .5*Pmax ) dphi -= Pmax;
                  while( dphi <-.5*Pmax ) dphi += Pmax;
               }
               //Step C: face formed out of end of cell- and beginning of cell+.
               xp[1] = cL->piph;
               xm[1] = cR->piph-cR->dphi;
               if( mode==1 ) addFace( theFaces , n , cL , cR , dxL , dxR , xp , xm , dim );
               ++n;
            }
         }
      }
   }

   if( mode==0 ){
      int NN = get_num_rzFaces( Nr , Nz , dim );
      ntj[NN] = n;
   }

}



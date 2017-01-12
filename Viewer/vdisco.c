
//
// This code was created by Jeff Molofee '99 (ported to Linux/GLUT by Richard Campbell '99)
// If you've found this code useful, please let me know.
//
// Visit me at www.demonews.com/hosted/nehe 
// (email Richard Campbell at ulmont@bellsouth.net)
//
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <fstream>

#include <hdf5.h>

#include <GL/glut.h>    // Header File For The GLUT Library 
#include <GL/gl.h>	// Header File For The OpenGL32 Library
#include <GL/glu.h>	// Header File For The GLu32 Library
#include <unistd.h>     // needed to sleep

/* ASCII code for the escape key. */
#define ESCAPE 27

#define VAL_FLOOR -1   //0.95//0 //-8e-3//-3e-2 //(-HUGE_VAL)  //.96
#define VAL_CEIL  1 //4.5e-3 //1.05//5.25e-21 //5.25e-9 //8e-3//3e-2 //5.24e-5 //(HUGE_VAL)  //1.04
#define FIXMAXMIN 1
#define COLORMAX 7
#define CAM_BACKUP  1.5
#define ZRORDER 1 // 1: checkpoints organized with faster index r (default,new)
                  // 0: checkpoints organized with slower index r (old)

static int WindowWidth  = 600;
static int WindowHeight = 600;
 
int CommandMode;
int FullScreenMode=0;

int dim3d = 0;
int t_off = 0;
int p_off = 0;
int cmap = 4;
int draw_1d = 0;
int draw_bar = 0;
int draw_t   = 0;
int draw_spiral = 0;
int draw_jet = 0;
int draw_planet = 0;
int draw_scale = 0;
int reflect  = 0;
int valq=0;
int draw_border = 0;
int logscale = 0;
int floors=0;
int help_screen=0;
int print_vals=0;
int fix_zero=0;

double rotate_angle = M_PI/2.;

double offx, offy, rescale, maxval, minval;

double t;
int Nr,Nz,Nq,Npl,N1d,midz;
int * Np;
double * r_jph;
double * z_kph;
double ** p_iph;
double *** theZones;
double ** rzZones;
double ** thePlanets;
double ** theRadialData;

double max_1d = 1.0;

void get_rgb( double , float * , float * , float * , int );

double getval( double * thisZone , int q ){
   if( q!=-1 ) return( thisZone[q] );
   double rho = thisZone[0];
//   double X   = thisZone[5];
   double P   = thisZone[1];
   double Br = thisZone[5];
   double Bp = thisZone[6];
   double Bz = thisZone[7];
//   double gam = sqrt(1.+ur*ur+up*up);
//   double e = (rho+4.*P)*gam*gam-P - rho*gam;
   return( .5*(Br*Br+Bp*Bp+Bz*Bz) ); //fabs(P/pow(rho,5./3.)-1.) );// fabs(thisZone[1]/pow(thisZone[0],5./3.)-1.) );
//   return( rho*P );
}

void getMaxMin(void){
   
   maxval = -HUGE_VAL;
   minval = HUGE_VAL;
   int q = valq;
   double val;
   int i,j,k;
   for( j=0 ; j<Nr ; ++j ){
      for( i=0 ; i<Np[j] ; ++i ){
         val = getval( theZones[j][i], q );
         if( logscale ) val = log(val)/log(10.);
         if( maxval < val ) maxval = val;
         if( minval > val ) minval = val;
      }
      if( dim3d ){
         for( k=0 ; k<Nz ; ++k ){
            int jk;
            if(ZRORDER)
                jk = k*Nr+j;
            else
                jk = j*Nz+k;
            val = getval( rzZones[jk], q );
            if( logscale ) val = log(val)/log(10.);
            if( maxval < val ) maxval = val;
            if( minval > val ) minval = val;
         }
      }
   }
   if( floors ){
      if( maxval > VAL_CEIL  ) maxval = VAL_CEIL;
      if( minval < VAL_FLOOR ) minval = VAL_FLOOR;
   }
   if( FIXMAXMIN && floors ){
      maxval = VAL_CEIL;
      minval = VAL_FLOOR;
   }
   if( fix_zero ){
      maxval = .5*(maxval-minval);
      minval = -maxval;
   }
   //if( floors ) minval = maxval-5.0;
   //if( floors ) minval = log( getval( theZones[0][Np[0]-1] , q ) )/log(10.);
}

void getH5dims( char * file , char * group , char * dset , hsize_t * dims ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );
   hid_t h5spc = H5Dget_space( h5dst );

   H5Sget_simple_extent_dims( h5spc , dims , NULL);

   H5Sclose( h5spc );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void readSimple( char * file , char * group , char * dset , void * data , hid_t type ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   H5Dread( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void readPatch( char * file , char * group , char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   hsize_t mdims[dim];
   hsize_t fdims[dim];

   hsize_t fstart[dim];
   hsize_t fstride[dim];
   hsize_t fcount[dim];
   hsize_t fblock[dim];

   int d;
   for( d=0 ; d<dim ; ++d ){
      mdims[d] = loc_size[d];
      fdims[d] = glo_size[d];

      fstart[d]  = start[d];
      fstride[d] = 1;
      fcount[d]  = loc_size[d];
      fblock[d]  = 1;
   }
   hid_t mspace = H5Screate_simple(dim,mdims,NULL);
   hid_t fspace = H5Screate_simple(dim,fdims,NULL);

   H5Sselect_hyperslab( fspace , H5S_SELECT_SET , fstart , fstride , fcount , fblock );

   H5Dread( h5dst , type , mspace , fspace , H5P_DEFAULT , data );

   H5Sclose( mspace );
   H5Sclose( fspace );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

/* The number of our GLUT window */
int window; 

// Here are the fonts: 
void ** glutFonts[7] = { 
    GLUT_BITMAP_9_BY_15, 
    GLUT_BITMAP_8_BY_13, 
    GLUT_BITMAP_TIMES_ROMAN_10, 
    GLUT_BITMAP_TIMES_ROMAN_24, 
    GLUT_BITMAP_HELVETICA_10, 
    GLUT_BITMAP_HELVETICA_12, 
    GLUT_BITMAP_HELVETICA_18 
}; 

// This function prints some text wherever you want it. 
void glutPrint(float x, float y, float z , void ** font, char* text, float r, float g, float b, float a) 
{ 
    if(!text || !strlen(text)) return; 
    int blending = 0; 
    if(glIsEnabled(GL_BLEND)) blending = 1; 
    glEnable(GL_BLEND); 
    glColor4f(r,g,b,a); 
    glRasterPos3f(x,y,z); 
    while (*text) { 
        glutBitmapCharacter(font, *text); 
        text++; 
    } 
    if(!blending) glDisable(GL_BLEND); 
}  
 
// A general OpenGL initialization function.  Sets all of the initial parameters. 
void InitGL(int Width, int Height)	        // We call this right after our OpenGL window is created.
{
  glClearColor(0.8f, 0.8f, 0.8f, 0.0f);		// This Will Clear The Background Color To Black
  glClearDepth(1.0);				// Enables Clearing Of The Depth Buffer
  glDepthFunc(GL_LESS);			        // The Type Of Depth Test To Do
  glEnable(GL_DEPTH_TEST);		        // Enables Depth Testing
  glShadeModel(GL_SMOOTH);			// Enables Smooth Color Shading
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();				// Reset The Projection Matrix
  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);	// Calculate The Aspect Ratio Of The Window
  glMatrixMode(GL_MODELVIEW);
}

/* The function called when our window is resized (which shouldn't happen, because we're fullscreen) */
void ReSizeGLScene(int Width, int Height)
{
  if (Height==0)				// Prevent A Divide By Zero If The Window Is Too Small
    Height=1;

  glViewport(0, 0, Width, Height);		// Reset The Current Viewport And Perspective Transformation

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);
  glMatrixMode(GL_MODELVIEW);
}

void TakeScreenshot(const char *useless){
   int dimx = WindowWidth;
   int dimy = WindowHeight;
   printf("Writing out.ppm... ");
   float *pixels = new float[3*dimx*dimy];
   glReadBuffer(GL_BACK);
   glPixelStorei(GL_PACK_ALIGNMENT,1);
   glReadPixels(0, 0, dimx, dimy, GL_RGB, GL_FLOAT, pixels);

   char *filename = (char *)"out.ppm";
   std::ofstream fp(filename);
   fp << "P3\n" << dimx << " " << dimy << std::endl << "255\n";
   int i,j;
   for( i=dimy-1; i>=0; i--){    
      for( j=0; j<dimx; ++j ){    
         int pixelPos = (i*dimx+j)*3;
         int r = (int)(pixels[pixelPos+0]*255.0);
         int g = (int)(pixels[pixelPos+1]*255.0);
         int b = (int)(pixels[pixelPos+2]*255.0);
         fp << r << " " << g << " " << b << std::endl;
      }    
   }    
   fp.close();
   printf("done!\n");
}

/* The main drawing function. */
void DrawGLScene(){

   glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
   glLoadIdentity();

   double RotationAngleX = 0.0;
   double RotationAngleY = 0.0;
   double RotationAngleZ = 0.0;

   if( dim3d==1 ) RotationAngleX = -60.0;
   if( dim3d==2 ) RotationAngleX = -90.0;
 
   glTranslatef( 0.0 , 0.0 , -CAM_BACKUP );

   glRotatef(RotationAngleX, 1, 0, 0);
   glRotatef(RotationAngleY, 0, 1, 0);
   glRotatef(RotationAngleZ, 0, 0, 1);

   int q = valq;

   getMaxMin();
   if( print_vals==0 ){
      printf("Max = %e Min = %e\n",maxval,minval);
      print_vals=1;
   }

   double camdist = 0.0;//-CAM_BACKUP;
   double xoff    = offx;
   double yoff    = offy;
   double zoff    = 0.0;

   if( draw_1d ){
      glColor3f(1.0,1.0,1.0);
      glBegin(GL_QUADS);
      glVertex3f(-1./rescale-xoff,-1./rescale-yoff, camdist + zoff );
      glVertex3f( 1./rescale-xoff,-1./rescale-yoff, camdist + zoff );
      glVertex3f( 1./rescale-xoff, 1./rescale-yoff, camdist + zoff );
      glVertex3f(-1./rescale-xoff, 1./rescale-yoff, camdist + zoff );
      glEnd();
      glLineWidth( 3.0f );
      glColor3f(0.0,0.0,0.0);
      glBegin( GL_LINE_STRIP );
      int j;
//      for( j=0 ; j<Nr ; ++j ){
      for( j=0 ; j<Nr ; ++j ){
         //if(logscale) val = log(getval(theZones[i],q))/log(10.);
         double rm = r_jph[j-1];
         double rp = r_jph[j]; 
         double val = theRadialData[j][q]/max_1d;//*pow(.5*(rp+rm),1.5);//2.*(theRadialData[j][q]-minval)/(maxval-minval) - 1.;//getval(theZones[i],q);
         if( q==1 ) val *= 1e2;
         val = 2.*val - 1.;
         double c0 = 2.*((.5*(rp+rm)-r_jph[-1])/(r_jph[Nr-1]-r_jph[-1]) - .5)/rescale;
         double c1 = val/rescale;
         glVertex3f( c0-xoff, c1-yoff, camdist-zoff+.001 );
      }    
      glEnd();
/*
      glColor3f(0.5,0.5,0.5);
      glBegin( GL_LINE_STRIP );

      for( j=0 ; j<Nr ; ++j ){
         //if(logscale) val = log(getval(theZones[i],q))/log(10.);
         double rm = r_jph[j-1];
         double rp = r_jph[j];
         double r = .5*(rp+rm); 
         double val = pow(r,-1.5)/max_1d;// *pow(.5*(rp+rm),1.5);//2.*(theRadialData[j][q]-minval)/(maxval-minval) - 1.;//getval(theZones[i],q);         if( q==1 ) val *= 1e2;
         val = 2.*val - 1.;
         double c0 = 2.*((.5*(rp+rm)-r_jph[-1])/(r_jph[Nr-1]-r_jph[-1]) - .5)/rescale;                 double c1 = val/rescale;
         glVertex3f( c0-xoff, c1-yoff, camdist-zoff+.001 );
      }       
      glEnd();
*/
      if( draw_planet ){
         double eps = 0.01;
         int p;
         for( p=0 ; p<Npl ; ++p ){
            double r   = thePlanets[p][0];
            double x = 2.*(( r - r_jph[-1] )/( r_jph[Nr-1]-r_jph[-1] ) - .5)/rescale;
            glColor3f(0.0,0.0,0.0);
            glLineWidth(2.0);
            glBegin(GL_LINE_LOOP);
            glVertex3f( x-xoff+eps , 1./rescale-yoff+eps , camdist+.001 );
            glVertex3f( x-xoff-eps , 1./rescale-yoff+eps , camdist+.001 );
            glVertex3f( x-xoff-eps , 1./rescale-yoff-eps , camdist+.001 );
            glVertex3f( x-xoff+eps , 1./rescale-yoff-eps , camdist+.001 );
            glEnd();
         }
      }

   }else{
      int Num_Draws = 1+reflect+draw_border;
      int count;
      int i,j;
      for( count = 0 ; count < Num_Draws ; ++count ){
         int draw_border_now = 0;
         int draw_reflected_now = 0;
         if( count==reflect+1 ) draw_border_now = 1;
         if( reflect && count==1 ) draw_reflected_now = 1;
         if( draw_border_now ) zoff -= .0001;
      for( j=0 ; j<Nr ; ++j ){
         double rm = r_jph[j-1]/rescale;
         double rp = r_jph[j]/rescale;
         double rc = .5*(rp+rm);
         double dr = rp-rm;
      
         if( !draw_border ){
            rm = rc-.55*dr;
            rp = rc+.55*dr;
         }
         for( i=0 ; i<Np[j] ; ++i ){

            double phip = p_iph[j][i];
            double phim;
            if( i==0 ) phim = p_iph[j][Np[j]-1]; else phim = p_iph[j][i-1];
            if( t_off && !p_off ){ phip -= t; phim -= t; }
            if( p_off ){ phip -= thePlanets[1][1]; phim -= thePlanets[1][1]; }
            double dp   = phip-phim;
            while( dp < 0.0 ) dp += 2.*M_PI;
            while( dp > 2.*M_PI ) dp -= 2.*M_PI;
            double phi = phim + 0.5*dp;
            if( !draw_border ){
               phip = phi+.55*dp;
               phim = phi-.55*dp;
            }

            double val = (getval(theZones[j][i],q)-minval)/(maxval-minval);
            if(logscale) val = (log(getval(theZones[j][i],q))/log(10.)-minval)/(maxval-minval);
            if( val > 1.0 ) val = 1.0;
            if( val < 0.0 ) val = 0.0;
            //double u = getval( theZones[j][i] , 2 );
            //if( uMax < u ) uMax = u;

            float rrr,ggg,bbb;
            get_rgb( val , &rrr , &ggg , &bbb , cmap );

            //if( (!dim3d || (sin(phi)>0 || cos(phi+.25)<0.0)) && dim3d !=2 ){
            if( dp < 1.5 && (!dim3d || (sin(phi)>0 && cos(phi)>0.0)) && dim3d !=2 ){
               if( !draw_border_now ){ 
                  glColor3f( rrr , ggg , bbb );
                  glBegin(GL_POLYGON);
               }else{
                  glLineWidth(3.0f);
                  glColor3f(0.0,0.0,0.0);
                  glBegin(GL_LINE_LOOP);
               }

               double z0 = 0.0;
               if( dim3d ) z0 = z_kph[Nz-1]/rescale;
     
               double c0 = rm*cos(phim);
               double c1 = rm*sin(phim);
               glVertex3f( c0-xoff, c1-yoff, camdist-zoff+z0 );

               c0 = rp*cos(phim);
               c1 = rp*sin(phim);
               glVertex3f( c0-xoff, c1-yoff, camdist-zoff+z0 );

               c0 = rp*cos(phi)/cos(.5*dp);
               c1 = rp*sin(phi)/cos(.5*dp);
               glVertex3f( c0-xoff, c1-yoff, camdist-zoff+z0 );

               c0 = rp*cos(phip);
               c1 = rp*sin(phip);
               glVertex3f( c0-xoff, c1-yoff, camdist-zoff+z0 );

               c0 = rm*cos(phip);
               c1 = rm*sin(phip);
               glVertex3f( c0-xoff, c1-yoff, camdist-zoff+z0 );

               glEnd();
            }
         }
         if( dim3d ){
            int k;
            for( k=0 ; k<Nz ; ++k ){
               int jk;
               if(ZRORDER)
                  jk = k*Nr+j;
               else
                  jk = j*Nz+k;
               double rp = r_jph[j]/rescale;
               double rm = r_jph[j-1]/rescale;
               double zp = z_kph[k]/rescale;
               double zm = z_kph[k-1]/rescale;

               double phi = 0.0;//rzZones[jk][Nq];

            double val = (getval(rzZones[jk],q)-minval)/(maxval-minval);
            if(logscale) val = (log(getval(rzZones[jk],q))/log(10.)-minval)/(maxval-minval);
            if( val > 1.0 ) val = 1.0;
            if( val < 0.0 ) val = 0.0;
            float rrr,ggg,bbb;
            get_rgb( val , &rrr , &ggg , &bbb , cmap );
               if( !draw_border_now ){ 
                  glColor3f( rrr , ggg , bbb );
                  glBegin(GL_POLYGON);
               }else{
                  glLineWidth(2.0f);
                  glColor3f(0.0,0.0,0.0);
                  glBegin(GL_LINE_LOOP);
               } 
                  
                  glVertex3f( rp*cos(phi) - xoff , rp*sin(phi) - yoff + zoff, zp );
                  glVertex3f( rm*cos(phi) - xoff , rm*sin(phi) - yoff + zoff, zp );
                  glVertex3f( rm*cos(phi) - xoff , rm*sin(phi) - yoff + zoff, zm );
                  glVertex3f( rp*cos(phi) - xoff , rp*sin(phi) - yoff + zoff, zm );
                  glEnd();
            }
            if( draw_border_now || dim3d == 1 ){
               glLineWidth(2.0f);
               glColor3f(0.0,0.0,0.0);
               glBegin(GL_LINE_LOOP);

               double rp = r_jph[Nr-1]/rescale;
               //double rm =-r_jph[Nr-1]/rescale;
               double rm = r_jph[-1]/rescale;
               double zp = z_kph[Nz-1]/rescale;
               double zm = z_kph[  -1]/rescale;

               glVertex3f( rp-xoff , -yoff+zoff , zp );
               glVertex3f( rm-xoff , -yoff+zoff , zp );
               glVertex3f( rm-xoff , -yoff+zoff , zm );
               glVertex3f( rp-xoff , -yoff+zoff , zm );
               glEnd();
            }
         }
      }
   }


   if( draw_bar ){
      double xb = 0.6;
      double hb = 1.0;
      double wb = 0.02;
      int Nb = 1000;
      double dy = hb/(double)Nb;
      int k;
      for( k=0 ; k<Nb ; ++k ){
         double y = (double)k*dy - .5*hb;
         double val = (double)k/(double)Nb;
         float rrr,ggg,bbb;
         get_rgb( val , &rrr , &ggg , &bbb , cmap );
         glLineWidth(0.0f);
         glColor3f( rrr , ggg , bbb );
         glBegin(GL_POLYGON);
         glVertex3f( xb    , y    , camdist + .001 ); 
         glVertex3f( xb+wb , y    , camdist + .001 ); 
         glVertex3f( xb+wb , y+dy , camdist + .001 ); 
         glVertex3f( xb    , y+dy , camdist + .001 ); 
         glEnd();
      }
      int Nv = 8;
      for( k=0 ; k<Nv ; ++k ){
         double y = (double)k*hb/(double)(Nv-1) - .5*hb;
         double val = (double)k/(double)(Nv-1)*(maxval-minval) + minval;
         char valname[256];
         sprintf(valname,"%+.2e",val);
         glLineWidth(1.0f);
         glColor3f(0.0,0.0,0.0);
         glBegin(GL_LINE_LOOP);
         glVertex3f( xb    , y , camdist + .0011 );
         glVertex3f( xb+wb , y , camdist + .0011 );
         glEnd();
         glutPrint( xb+1.5*wb , y-.007 , camdist + .001 ,glutFonts[6] , valname , 0.0f, 0.0f , 0.0f , 0.5f );
      }
   } 

   if( draw_planet ){
      double eps = 0.01;
      int p;
      for( p=0 ; p<Npl ; ++p ){
         double r   = thePlanets[p][0]/rescale;
         double phi = thePlanets[p][1];
         if( t_off && !p_off ) phi -= t;
         if( p_off ) phi -= thePlanets[1][1];
         glColor3f(0.0,0.0,0.0);
         glLineWidth(2.0);
         glBegin(GL_LINE_LOOP);
         glVertex3f( r*cos(phi)-xoff+eps , r*sin(phi)-yoff+eps , camdist+.001 );
         glVertex3f( r*cos(phi)-xoff-eps , r*sin(phi)-yoff+eps , camdist+.001 );
         glVertex3f( r*cos(phi)-xoff-eps , r*sin(phi)-yoff-eps , camdist+.001 );
         glVertex3f( r*cos(phi)-xoff+eps , r*sin(phi)-yoff-eps , camdist+.001 );
         glEnd();
      }
   }

   if( draw_spiral ){

      double rp   = 5.0;//thePlanets[1][0];
      double Mach = 3.65;
      double p0 = 1.2;
      double dr = .07;
      int Nr = 200;
      double Rmin = 0.5;
      double Rmax = 1.1;
      int k;
      glLineWidth(3.0f);
      glColor3f(0.0,0.0,0.0);
      //glColor3f(1.0,1.0,1.0);
      glBegin(GL_LINE_LOOP);
      for( k=0 ; k<Nr ; ++k ){
         //double phi0 = ((double)k+.5)/(double)Nr*2.*M_PI;//(3.-2.*sqrt(1./r)-r)*20.;
         double x = ((double)k+.5)/(double)Nr;
         double r = .001*pow(.5/.001,x);
         //double r = .5*pow(1.5/.5,x);
         double phi0 = p0-log(r/0.001)*Mach;//(3.-2.*sqrt(1./r)-r)*5.;
         //double phi0 = (3.-2.*sqrt(1./r)-r)*20.;
         //if( r<1. ) phi0 = -phi0;

//         double x0 = rp + dr*cos(phi0);
//         double y0 = dr*sin(phi0);
  
         double phi = phi0;
//         double r = rp;       

//         double phi = atan2(y0,x0);
//         double phi = ((double)k+.5)/(double)Nr*2.*M_PI*9.32 + 4.*M_PI;
//         double r   = sqrt(x0*x0+y0*y0);
//         double r = pow(phi/10.,-2.);
//         double r   = 2./(1.+sin(phi));//1.0;//((double)k+0.5)/(double)Nr*(Rmax-Rmin) + Rmin;
//         if( r<1. ) phi = -phi;
         r /= rescale;
         
         glVertex3f( r*cos(phi)-xoff , r*sin(phi)-yoff , camdist + .0011 );
         if( k%2==1 ){
            glEnd();
            glBegin(GL_LINE_LOOP);
         }
      }
      glEnd();


/*
      e += 0.01;
      glBegin(GL_LINE_LOOP);
      for( k=0 ; k<Nr ; ++k ){
         double phi0 = ((double)k-.5)/(double)Nr*2.*M_PI;//(3.-2.*sqrt(1./r)-r)*20.;
         double x0 = 1.+e*cos(phi0);
         double y0 = e*sin(phi0);
         
         double phi = atan2(y0,x0);
         double r   = sqrt(x0*x0+y0*y0);
//         double phi = (double)k/(double)Nr*2.*M_PI;//(3.-2.*sqrt(1./r)-r)*20.;
//         double r   = 2./(1.+sin(phi));//1.0;//((double)k+0.5)/(double)Nr*(Rmax-Rmin) + Rmin;
//         if( r<1. ) phi = -phi;
         r /= rescale;
         
         glVertex3f( r*cos(phi)-xoff , r*sin(phi)-yoff , camdist + .0011 );
         if( k%2==1 ){
            glEnd();
            glBegin(GL_LINE_LOOP);
         }
      }
      glEnd();
*/
   }
   }

   if( draw_t ){
      char tprint[256];
      sprintf(tprint,"t = %.2e",t/2./M_PI);
      glutPrint( -.6 , .5 , camdist + .001 , glutFonts[6] , tprint , 0.0f, 0.0f , 0.0f , 0.5f );
  //    sprintf(tprint,"uMax = %.1f",uMax);
  //    glutPrint( -.8 , .4 , camdist + .001 , glutFonts[6] , tprint , 0.0f, 0.0f , 0.0f , 0.5f );
   }
   if( help_screen ){
      int NLines = 9;
      double dy = -.04;
      char help[NLines][256];
      sprintf(help[0],"Help Display:");
      sprintf(help[1],"b - Toggle Colorbar");
      sprintf(help[2],"c - Change Colormap");
      sprintf(help[3],"f - Toggle Max/Min Floors");
      sprintf(help[4],"p - Toggle Planet Data");
      sprintf(help[5],"1-9 - Choose Primitive Variable to Display");
      sprintf(help[6],"wasd - Move Camera");
      sprintf(help[7],"z/x - Zoom in/out");
      sprintf(help[8],"h - Toggle Help Screen");
      int i;
      for( i=0 ; i<NLines ; ++i ){
         glutPrint( -.8 , i*dy , camdist + .001 , glutFonts[6] , help[i] , 0.0f, 0.0f , 0.0f , 0.5f );
      }
   } 

   if( CommandMode ){
      TakeScreenshot("out.ppm");
      exit(1);
   }

   glutSwapBuffers();

}

/* The function called whenever a key is pressed. */
void keyPressed(unsigned char key, int x, int y) 
{
   usleep(100);
   if (key == ESCAPE){
      glutDestroyWindow(window); 
      exit(0);                   
   }
   if (key == '[' ){
       --cmap;
       cmap = (cmap+COLORMAX) % COLORMAX;
   }  
   if (key == ']' ){
       ++cmap;
       cmap = cmap % COLORMAX;
   }
   if( key >= (int)'0' && key < (int)'1'+Nq ) valq = (int)key-(int)'1';
   if( key == 'b' ) draw_bar = !draw_bar;
   if( key == 'd' ) {++dim3d; if(dim3d==3) dim3d=0;}
   if( key == 'f' ) floors = !floors;
   if( key == 'g' ) draw_border = !draw_border;
   if( key == 'h' ) help_screen = !help_screen;
   if( key == 'l' ) logscale = !logscale;
   if( key == 'r' ) reflect = !reflect;
   if( key == 's' ) draw_scale = !draw_scale;
   if( key == 't' ) draw_t = !draw_t;
   if( key == 'u' ) draw_1d = !draw_1d;
   if( key == 'x' ){
      rescale *= 1.3;
      offx /= 1.3;
      offy /= 1.3;
   }
   if( key == 'z' ){
      rescale /= 1.3;
      offx *= 1.3;
      offy *= 1.3;
   }
   if( key == 'j' ){
      max_1d *= 1.3;
   }
   if( key == 'J' ){
      max_1d /= 1.3;
   }
   if (key == 'F') {
      if (FullScreenMode) {
         glutReshapeWindow(WindowWidth, WindowHeight);
         glutPositionWindow(0,0);
         FullScreenMode = 0; 
      }else{
         glutFullScreen();
         FullScreenMode = 1; 
      }
   }
   if( key == 'p' ) draw_planet = !draw_planet;
   if( key == 's' ) draw_spiral = !draw_spiral;
   if( key == 'Z' ) fix_zero = !fix_zero;
   glutPostRedisplay();
}

void specialKeyPressed(int key, int x, int y){
   usleep(100);
   if (key == GLUT_KEY_LEFT ) offx -= .1;
   if (key == GLUT_KEY_RIGHT) offx += .1;

   if (key == GLUT_KEY_UP   ) offy += .1;
   if (key == GLUT_KEY_DOWN ) offy -= .1;
}

int main(int argc, char **argv) 
{
   if( argc < 2 ){
      printf("Please specify the input file.\n");
      exit(1);
   }
   char filename[256];
   if( argv[1] ){
      strcpy( filename , argv[1] );
   }
   CommandMode=0;
   if( argc>2 ){ CommandMode=1; FullScreenMode=1; }
   char group1[256];
   char group2[256];
   strcpy( group1 , "Grid" );
   strcpy( group2 , "Data" );

   hsize_t dims[3];
   
   readSimple( filename , group1 , (char *)"T" , &t , H5T_NATIVE_DOUBLE );
   getH5dims( filename , group1 , (char *)"r_jph" , dims );
   Nr = dims[0]-1;
   getH5dims( filename , group1 , (char *)"z_kph" , dims );
   Nz = dims[0]-1;

   Np = (int *) malloc( Nr*sizeof(int) );
   r_jph = (double *) malloc( (Nr+1)*sizeof(double) );
   z_kph = (double *) malloc( (Nz+1)*sizeof(double) );
   int Tindex[Nr];
   
   printf("t = %.2f, Nr = %d Nz = %d\n",t,Nr,Nz);

   readSimple( filename , group1 , (char *)"r_jph" , r_jph , H5T_NATIVE_DOUBLE );
   readSimple( filename , group1 , (char *)"z_kph" , z_kph , H5T_NATIVE_DOUBLE );

   //midz = Nz/2;
   midz = Nz-1;

   int start[2]    = {0,0};
   int loc_size[2] = {Nz,Nr};
   int glo_size[2] = {Nz,Nr};
   if(!ZRORDER)
   {
       loc_size[0] = Nr; loc_size[1] = Nz;
       glo_size[0] = Nr; glo_size[1] = Nz;
   }

   int Np_All[Nr*Nz];
   int Tindex_All[Nr*Nz];
   int Id_phi0[Nr*Nz];

   readPatch( filename , group1 , (char *)"Np"      , Np_All     , H5T_NATIVE_INT , 2 , start , loc_size , glo_size);
   readPatch( filename , group1 , (char *)"Index"   , Tindex_All , H5T_NATIVE_INT , 2 , start , loc_size , glo_size);
   readPatch( filename , group1 , (char *)"Id_phi0" , Id_phi0    , H5T_NATIVE_INT , 2 , start , loc_size , glo_size);

   int i,j,k;
   for( j=0 ; j<Nr ; ++j ){
      k = midz;
      int jk = k*Nr+j;
      if(!ZRORDER)
          jk = j*Nz+k;
      Np[j]     = Np_All[jk];
      Tindex[j] = Tindex_All[jk];
   }

   getH5dims( filename , group2 , (char *)"Cells" , dims );
   int Nc = dims[0];
   Nq = dims[1]-1;

   getH5dims( filename , group2 , (char *)"Planets" , dims );
   Npl = dims[0];
   int NpDat = dims[1];
   printf("Nc = %d Nr = %d Nq=%d Npl=%d NpDat=%d\n",Nc,Nr,Nq,Npl,NpDat);

   theZones = (double ***) malloc( Nr*sizeof(double **) );
   int q;
   for( j=0 ; j<Nr ; ++j ){
      theZones[j] = (double **) malloc( Np[j]*sizeof( double * ) );
      for( i=0 ; i<Np[j] ; ++i ){
         theZones[j][i] = (double *) malloc( Nq*sizeof( double ) );
      }
   }

   rzZones = (double **) malloc( Nr*Nz*sizeof(double *) );
   for( j=0 ; j<Nr*Nz ; ++j ){
      rzZones[j] = (double *) malloc( (Nq+1)*sizeof( double ) );
   }

   p_iph = (double **) malloc( Nr*sizeof( double * ) );
   for( j=0 ; j<Nr ; ++j ){
      p_iph[j] = (double *) malloc( Np[j]*sizeof( double ) );
   }
   thePlanets = (double **) malloc( Npl*sizeof(double *) );
   int p;
   for( p=0 ; p<Npl ; ++p ){ 
      thePlanets[p] = (double *) malloc( 2*sizeof(double) );
   }

   getH5dims( filename , group2 , (char *)"Radial_Diagnostics" , dims );
   N1d = dims[1];
   theRadialData = (double **) malloc( Nr*sizeof(double *) );
   for( j=0 ; j<Nr ; ++j ){
      theRadialData[j] = (double *) malloc( N1d*sizeof(double) );
   }

   printf("Zones Allocated\n");
   loc_size[1] = Nq+1;
   glo_size[0] = Nc;
   glo_size[1] = Nq+1;

   for( j=0 ; j<Nr ; ++j ){
      loc_size[0] = Np[j];
      start[0] = Tindex[j];

      double TrackData[Np[j]*(Nq+1)];
      readPatch( filename , group2 , (char *)"Cells" , TrackData , H5T_NATIVE_DOUBLE , 2 , start , loc_size , glo_size);
      for( i=0 ; i<Np[j] ; ++i ){
         p_iph[j][i] = TrackData[i*(Nq+1) + Nq];
         for( q=0 ; q<Nq ; ++q ){
            theZones[j][i][q] = TrackData[i*(Nq+1) + q];
         }
      }
   }
   printf("theZones built\n");
/**/
   loc_size[0] = 1;
   if(ZRORDER)
   {
       for( k=0 ; k<Nz ; ++k ){
          for( j=0 ; j<Nr ; ++j ){
             int jk = k*Nr+j;
             start[0] = Id_phi0[jk];
             readPatch( filename , group2 , (char *)"Cells" , rzZones[jk] , H5T_NATIVE_DOUBLE , 2 , start , loc_size , glo_size );
          }
       }
   }
   else
   {
       for( j=0 ; j<Nr ; ++j ){
          for( k=0 ; k<Nz ; ++k ){
             int jk = j*Nz+k;
             start[0] = Id_phi0[jk];
             readPatch( filename , group2 , (char *)"Cells" , rzZones[jk] , H5T_NATIVE_DOUBLE , 2 , start , loc_size , glo_size );
          }
       }
   }
/**/
   printf("rzZones built\n");

   start[1] = 0;
   loc_size[0] = 1;
   glo_size[0] = Npl;
   loc_size[1] = NpDat;
   glo_size[1] = NpDat;
   for( p=0 ; p<Npl ; ++p ){
      start[0] = p;
      double thisPlanet[6];
      readPatch( filename , group2 , (char *)"Planets" , thisPlanet , H5T_NATIVE_DOUBLE , 2 , start , loc_size , glo_size );
      int dp;
      printf("Planet = ");
      for( dp=0 ; dp<NpDat ; ++dp ) printf("%e ",thisPlanet[dp]);
      printf("\n");
      thePlanets[p][0] = thisPlanet[3];
      thePlanets[p][1] = thisPlanet[4];
   }

   loc_size[0] = 1;
   glo_size[0] = Nr;
   loc_size[1] = N1d;
   glo_size[1] = N1d;
   start[1] = 0;
   for( j=0 ; j<Nr ; ++j ){
      start[0] = j;
      double thisQ[N1d];
      readPatch( filename , group2 , (char *)"Radial_Diagnostics" , thisQ , H5T_NATIVE_DOUBLE , 2 , start , loc_size , glo_size );
      int nq;
      for( nq = 0 ; nq < N1d ; ++nq ) theRadialData[j][nq] = thisQ[nq];
   }

   r_jph++;
   z_kph++;

   //double RR = r_iph[0][Nr[0]-1];
   //double r_max = .98*RR;
   //double r_min = 0.0;//r_iph[0][0];
   //double r_min = .82*RR;

   //printf("Rmin = %.2e Rmax = %.2e\n",r_min,r_max);
   //double thalf = .5*t_jph[Nt-1];
   rescale = 1.7*r_jph[Nr-1];  //(r_max-r_min);
   //offx = .5*(r_min+r_max)/rescale;
   //offy = 0.5;//.5*(r_min+r_max)*sin(thalf)/rescale;
   //offy = .5*(r_min+r_max)/rescale;
   //rescale *= .5;
   //rescale *= .4;
   rescale *= 1.2;
   //offy = 0.75/rescale;
   //offx = 0.6/rescale;//0.0;//1.0/rescale;
   offy = 0.0;
   offx = 0.0;//2.5/rescale;//0.0;//1.0/rescale;

//////////////////////////////
   glutInit(&argc, argv);  
   glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);  
   glutInitWindowSize(WindowWidth, WindowHeight);
   glutInitWindowPosition(0, 0);
   window = glutCreateWindow("DISCO Viewer");
   glutDisplayFunc(&DrawGLScene);  
   if(FullScreenMode) glutFullScreen();
   glutIdleFunc(&DrawGLScene);
   glutReshapeFunc(&ReSizeGLScene);
   glutKeyboardFunc(&keyPressed);
   glutSpecialFunc(&specialKeyPressed);
   InitGL(WindowWidth, WindowHeight);
   glutMainLoop();  

   for( p=0 ; p<Npl ; ++p ){
      free(thePlanets[p]);
   }
   free(thePlanets);

   for( j=0 ; j<Nr ; ++j ){
      for( i=0 ; i<Np[j] ; ++i ){
         free( theZones[j][i] );
      }
      free( theZones[j] );
   }
   for( j=0 ; j<Nr*Nz ; ++j ){
      free( rzZones[j] );
   }
   free( rzZones );

   free( theZones );
   for( j=0 ; j<Nr ; ++j ){
      free( p_iph[j] );
      free( theRadialData[j] );
   }
   free( p_iph );
   free( theRadialData );
   free( Np );
   r_jph--;
   free( r_jph );
   z_kph--;
   free( z_kph );

   return (0);

}

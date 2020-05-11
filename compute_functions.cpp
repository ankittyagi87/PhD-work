#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cmath>
//#include <math.h>
#include <complex>
//#include <sys/types.h>
//#include <unistd.h>
#include <iostream>

#include "channelflow/flowfield.h"
#include "channelflow/vector.h"
#include "channelflow/dns.h"
#include "channelflow/periodicfunc.h"
#include "channelflow/def"

using namespace std;
using namespace channelflow;

extern int  list[n_atom], head[Nof_lattice], map[map_size];
extern int  backward_index,Nof_backward_step,overlaping_a, overlaping_b;
extern int  cell_of_part[n_atom+2];
extern int npair[n_atom],cell_of_part[n_atom+2];
extern int atom_a, atom_b,time_step;
extern int previous_atom_a, previous_atom_b;
extern int  N_diff;
extern int count_Rep;

extern bool event_noise_injection;



extern double sigma,sigma_sq,radius,Reynolds;
extern double p_x[n_atom+2],p_y[n_atom+2],p_z[n_atom+2];//position
extern double p_vx[n_atom+2],p_vy[n_atom+2],p_vz[n_atom+2],p_v[n_atom+2];//velocity
extern double  p_angvelx[n_atom+2],p_angvely[n_atom+2],p_angvelz[n_atom+2],p_angvel[n_atom+2];//angular velocity
extern double vx_max_air,vy_air,vz_air,b_ly;
extern double fdrag_x[n_atom],fdrag_y[n_atom],fdrag_z[n_atom];
extern double  Total_Rep_x, Total_Rep_y, Total_Rep_z;

extern double den_gas, vis_gas,part_mass,den_part,moment_inertia;

extern double g_x,g_y,g_z;

extern double const_diff_x,const_diff_y,const_diff_z,const_diff_xy;
extern double length_scale,time_scale,acc_scale,vel_scale;

extern double acc_x[n_atom+2],acc_y[n_atom+2],acc_z[n_atom+2];
extern double angaccx[n_atom+2],angaccy[n_atom+2],angaccz[n_atom+2];
extern double coll_time,ctime_wall;
extern bool overlap;

extern double p_x_previous[n_atom+2],p_y_previous[n_atom+2],p_z_previous[n_atom+2],
  p_vx_previous[n_atom+2],p_vy_previous[n_atom+2],p_vz_previous[n_atom+2];
double p_angvelx_previous[n_atom+2],p_angvely_previous[n_atom+2],p_angvelz_previous[n_atom+2];

extern double t_next,run_time;
extern double t_noise,del_tb,del_tc;
extern double b_lx,b_ly,b_lz;
extern double et,en,ep;
extern double tau_vp;

extern double tot_KE,tc[n_atom];double tot_angularKE;
extern double  p_xi, p_yi, p_zi, p_vxi, p_vyi, p_vzi;


extern double u_noise,v_noise,w_noise;
extern double diff_x[diff_nmax],diff_y[diff_nmax],diff_z[diff_nmax],
  diff_xy[diff_nmax],y_ch[diff_nmax];


extern  void    maping_neighbour();
extern  void    links ( ); 
extern  int     index_cell(int, int, int);
extern  void    drag(int i , const FlowField& u,const FlowField& omega,const Vector& x_grid, const Vector&  z_grid);
extern  double  collision_time(int i,int j);
extern  double  wall_coll_time(int i,int j);
extern  void    backward_step();
extern  void	advance_particle();
extern  void    post_coll_vel();
extern  void	update_coll_table(const FlowField& u ,const FlowField& omega,const Vector& x_grid , const Vector& z_grid  );
extern  void    adding_noise();

extern int n_order; // --5 th order lagrangian interpolation
extern void interpolation(const FlowField& , const Vector& ,const Vector& ,double ,double , double , int, double[]);
void interpolation_tensor(FlowField& , const Vector& , const Vector& ,double ,double , double , int , double [3][3] );
//extern void reverseforce(int i, const FlowField& force);
double freverse_x[n_atom],freverse_y[n_atom],freverse_z[n_atom];
double torque_x[n_atom],torque_y[n_atom],torque_z[n_atom];  //Torque of particles
double treverse_x[n_atom],treverse_y[n_atom],treverse_z[n_atom];
//double ***freverse_x2; double ***freverse_y2; double ***freverse_z2;

double ***freverse_x2, ***freverse_y2, ***freverse_z2;
extern double **frev_cheby_x, **frev_cheby_y, **frev_cheby_z; 
extern double **frev_cheby_diff2_x, **frev_cheby_diff2_y, **frev_cheby_diff2_z;
double **frev_fourier_x,**frev_fourier_y,**frev_fourier_z;
double ***freverse_x2_cheby_diff2, ***freverse_y2_cheby_diff2, ***freverse_z2_cheby_diff2;
extern double **tmp3, **tmp_diff, **tmp_diff2;

//*************Functions for calc of stress moments by interpolation****************************************************************//
extern void reverseforce(int i,FlowField& force);

extern void third_order_stress_moments(int i,FlowField& force, FlowField& Dxx1, FlowField& Dyy1, FlowField& Dzz1);
extern void rotational_stress_tensor(int i, FlowField& force, FlowField& dx1, FlowField& dy1, FlowField& dz1);
extern void symmetric_stress_tensor(int i,FlowField& u,FlowField& utmp,FlowField& force,FlowField& sx1, const Vector& x_grid, const Vector&  z_grid);

extern void chebyshev_transform(const FlowField& force);
extern void chebyshev_differentiation(const FlowField& force);
void reverse_quantity (int i, const FlowField& force, double **in, double ***out);
extern void reverseforce_spectral(int i, FlowField& force);

extern double symmetric_air_p_xx[n_atom],symmetric_air_p_xy[n_atom],symmetric_air_p_xz[n_atom],
                  symmetric_air_p_yx[n_atom],symmetric_air_p_yy[n_atom],symmetric_air_p_yz[n_atom],
                  symmetric_air_p_zx[n_atom],symmetric_air_p_zy[n_atom],symmetric_air_p_zz[n_atom];

//**************End of Functions for calc of stress moments by interpolation**********************************************************//

//**************Functions for calc of stress moments by spectral method***************************************************************//
//extern void reverseforce_spectral_method(int i, FlowField& force);
//extern void third_order_stress_moments_spectral(int i,FlowField& Dxx,FlowField& Dyy,FlowField& Dzz,FlowField& Dxx1,FlowField& Dyy1,FlowField& Dzz1);
//extern void rotational_stress_tensor_spectral(int i,FlowField& dx, FlowField& dy, FlowField& dz, FlowField& dx1, FlowField& dy1, FlowField& dz1);
//extern void symmetric_stress_tensor_spectral_method(int i,FlowField& u,FlowField& sx,FlowField& sx1,const Vector& x_grid, const Vector&  z_grid);
//void reversetorque_spectral_method(int i, FlowField& torque);
//**************End Functions for calc of stress moments by spectral method***************************************************************//

extern int Max_dns_grid;
//extern channelflow::Complex ***tmp4;
//extern void Fourier_chebhyshev_transform(int i, const FlowField& u);



extern int Max_dns_grid;

void  maping_neighbour()
{

	int ix,iy,iz;
	int imap,imap_init,icell_0;
  if(dimension==3)
{ 
  for(iy=0;iy<Z-1;iy++)
  
  for(iz = 0;iz<Z;iz++)
    
       for(ix=0;ix<Z;ix++)
	   {
             
	     imap = ( index_cell(ix, iy, iz )  ) * 13;
	     //cout<<" in maps"<<imap<<endl;
	     map[imap      ] = index_cell( ix + 1, iy    , iz     );
	     map[imap + 1  ] = index_cell( ix + 1, iy + 1, iz     );
	     map[imap + 2  ] = index_cell( ix    , iy + 1, iz     );
	     map[imap + 3  ] = index_cell( ix - 1, iy + 1, iz     );
	     map[imap + 4  ] = index_cell( ix + 1, iy    , iz - 1 );
	     map[imap + 5  ] = index_cell( ix + 1, iy + 1, iz - 1 );
	     map[imap + 6  ] = index_cell( ix    , iy + 1, iz - 1 );
	     map[imap + 7  ] = index_cell( ix - 1, iy + 1, iz - 1 );
	     map[imap + 8  ] = index_cell( ix + 1, iy    , iz + 1 );
	     map[imap + 9  ] = index_cell( ix + 1, iy + 1, iz + 1 );
	     map[imap + 10 ] = index_cell( ix    , iy + 1, iz + 1 );
	     map[imap + 11 ] = index_cell( ix - 1, iy + 1, iz + 1 );
	     map[imap + 12 ] = index_cell( ix    , iy    , iz + 1 );
	     
	  }

    imap_init= imap+13;
    icell_0=index_cell(0,Z-1, 0);
    
  for(iy=Z-1;iy<Z;iy++)  
     
   for(iz =0;iz<Z;iz++)    
               
         for(ix=0;ix<Z;ix++)
	   {
	     imap = imap_init+( index_cell (ix, iy, iz )- icell_0 )*4;
             map[imap      ] = index_cell( ix + 1, iy    , iz     );
	     map[imap + 1  ] = index_cell( ix + 1, iy    , iz - 1 );
	     map[imap + 2  ] = index_cell( ix + 1, iy    , iz + 1 );
	     map[imap + 3  ] = index_cell( ix    , iy    , iz + 1 );
	    
           //cout<<ix<<" "<<iy<<" "<<iz<<endl;
	   //  for(i=0;i<4;i++)
	  // cout<<map[imap+i   ] <<endl;

	   }
  /*
for(i=0;i<13;i++)
    cout<<"neighbour=  "<<map[1196+i]<<endl;
    cin.get();
  */
}

  

  else if(dimension ==2){
    iz=0;
     for(iy=0;iy<Z-1;iy++)
          
         for(ix=0;ix<Z;ix++)
	   {
             
	     imap = ( index_cell (ix, iy, iz )  ) *4;
	     //cout<<" in maps"<<imap<<endl;
	     map[imap      ] = index_cell( ix + 1, iy    , iz     );
	     map[imap + 1  ] = index_cell( ix + 1, iy + 1, iz     );
	     map[imap + 2  ] = index_cell( ix    , iy + 1, iz     );
	     map[imap + 3  ] = index_cell( ix - 1, iy + 1, iz     );
	     
	     
		   }
     imap_init= imap+4;
     icell_0=index_cell (0,Z-1, iz );
	
     for(iy=Z-1;iy<Z;iy++)
         for(ix=0;ix<Z;ix++)
	   {          
	     imap = imap_init+( index_cell (ix, iy, iz )- icell_0 ) *1;
	     
	     map[imap      ] = index_cell( ix + 1, iy    , iz     );
	    
	   }  

  }
  /*for(i=0;i<13;i++)
    {
    cout<<"neighbour=  "<<map[1196+i]<<endl;
    cin.get();
    }
  */
  cout<<" returning from map"<<endl;
return;
}

//------------------------------------ index cell -------------------

int index_cell (int IX1, int IY1, int IZ1 )
{
  int p;
  p= (int)(fmod ((double) (IX1 + Z),(double) Z )+ fmod ( (double) (IZ1 + Z),(double) Z ) * Z
       + fmod( (double)(IY1 + Z),(double) Z ) * Z * Z);

 return (p);
}



/*------------------------------------------------------------
     function to calculate te drag on the particle
 --------------------------------------------------------------*/

void drag(int i,const FlowField& u, const FlowField& omega, const Vector& x_grid, const Vector&  z_grid)
	       
         {
	  
	  double vx_air,  vy_air, vz_air,rel_vx,rel_angvelx,Rep_x,Ap;
           double angvelx_air,angvely_air,angvelz_air;
	   double rel_vy, rel_angvely,Rep_y;
	   double rel_vz, rel_angvelz,Rep_z;
           double interpolated_u[3];
           double interpolated_omega[3];
	   int n_order=5;
	   double Cd_x,Cd_y,Cd_z;

	   double a_=u.a();
	   double b_=u.b();
	   double ubase_air;
	   
		//----------------------------------------	
		 //Calculation of drag in the X-direction
		 //---------------------------------------

		 sigma=sqrt(sigma_sq);
		 double xp=p_x[i];
		 double yp=p_y[i];
		 double zp=p_z[i];

		 count_Rep=count_Rep+1;  //to calculate average Rep
   
  interpolation(u, x_grid,z_grid,xp,yp,zp, n_order, interpolated_u);
  interpolation(omega, x_grid, z_grid,xp,yp,zp, n_order, interpolated_omega); 

  ubase_air = (1.0 - square(abs(p_y[i]-(b_+a_)/2.0)/((b_-a_)/2.0))); //for channel flow

  //cout<<" in compute "<< p_y[i]<<"   " <<ubase_air<<endl;
 

  //ubase_air =p_y[i]-(a_+b_)/2.0;    // for couette flow 

  vx_air=ubase_air+interpolated_u[0];
  // vx_air=interpolated_u[0];
  // cout<<"vx_air  "<<vx_air<<endl;
  // cin.get();
  vy_air=interpolated_u[1];  // nondimensionalized with centerline air velocity of prabolic profile
  vz_air=interpolated_u[2];

  // cout<<vx_air<<" "<<vy_air<<" "<<vz_air<<endl;
  //cin.get();
  angvelx_air=0.5*interpolated_omega[0];
  angvely_air=0.5*interpolated_omega[1];
  angvelz_air=0.5*interpolated_omega[2];
//cout<<angvelx_air<<endl;
  // cout<<vx_air<<" "<<vy_air<<" "<<vz_air<<endl;
  //cin.get();
   
   rel_vx= fabs(p_vx[i]-vx_air);//rel_vx=Relative velocity of the particle in x-direction
       
   Rep_x=(sigma* rel_vx*den_gas)/vis_gas; 
   Rep_x=Rep_x*(length_scale*vel_scale);//particle Reynolds number 

   Total_Rep_x=Total_Rep_x+Rep_x;

   rel_angvelx=fabs(p_angvelx[i]-angvelx_air);
   rel_angvely=fabs(p_angvely[i]-angvely_air);
   rel_angvelz=fabs(p_angvelz[i]-angvelz_air); 

              

   Cd_x=0;   //just to initiallize

  if(Rep_x<=1.0e-16) 
    Cd_x=0.0;
      	         
  else if(Rep_x>=1.0e-16 && Rep_x <= 0.20)
      Cd_x=24.0/Rep_x;

    else if(Rep_x>0.20 && Rep_x <=1000.0)
      Cd_x=24.0*(1.0+0.15*(pow(Rep_x,0.67)))/Rep_x;

		 else if(Rep_x>1000.0 && Rep_x <=2.e5) 
		 Cd_x=0.440;

		    else if(Rep_x>2.e5)
		      Cd_x=0.050;
		 
    Ap=((22.0/7.0)*(sigma*sigma))/4.0;

    // fdrag_x[i]=Cd_x*(rel_vx*rel_vx)*den_gas*Ap/2.0;//in this expression, rel_vx and Ap is nondimensional
    // fdrag_x[i]=fdrag_x[i]/part_mass; // if it is non-stokes flow comment the line below and activate this line and the line above

   fdrag_x[i]=rel_vx/tau_vp;   //here fdrag_x is not drag it is acceleration , use above two lines for Rep>1
     
     // fdrag_x[i] = 3*pi*vis_gas*rel_vx*sigma/part_mass;
          
               if(p_vx[i]>vx_air)
		 fdrag_x[i]=-fdrag_x[i];
// ---------------------------------------------------------------//
//--------------------Reverse force-------------------------------//   
               freverse_x[i] = -6*pi*radius*rel_vx*b_ly*b_ly/Reynolds;
                if(p_vx[i]>vx_air)
		 freverse_x[i]=-freverse_x[i];
               //  freverse_x[i] = 1.0;
                 
//--------------------------------------------------------------//
//---------------------------------------------------------------//
//------------------Particle Torque-----------------------//
                     torque_x[i] = 10*rel_angvelx/(3*tau_vp);
                   if(p_angvelx[i]>angvelx_air)
                      torque_x[i] = -torque_x[i];
                
 //----------------------------------------------------------------//
//--------------------Angular Acceleration------------------------//
                angaccx[i]=torque_x[i];
//----------------------------------------------------------------//	

//-------------------------------------------------------------------//
//--------------------Reverse Torque---------------------------------//
                treverse_x[i] = 4*pi*radius*radius*radius*rel_angvelx*b_ly*b_ly/Reynolds;
                    if(p_angvelx[i]>angvelx_air)
                     treverse_x[i] = -treverse_x[i];
	 
	       //---------------------------------------
	       //Calculation of drag in the y-direction
	       //---------------------------------------

    rel_vy=fabs(p_vy[i]-vy_air);   //rel_vy=Relative velocity of the particle in p_y-direction
               
    Rep_y=(sigma* rel_vy*den_gas)/vis_gas; //particle Reynolds number 

    Rep_y=Rep_y*(length_scale*vel_scale);

 Total_Rep_y=Total_Rep_y+Rep_y;

  Cd_y=0;   //just to initiallize

        if(Rep_y<=1.0e-16) 
	  Cd_y=0.0;        
	else  if(Rep_y>=1.0e-16 && Rep_y <= 0.20)
		 Cd_y=24.0/Rep_y;

		 else if(Rep_y>0.20 && Rep_y <=1000.0)
		 Cd_y=24.0*(1.0+0.150*(pow(Rep_y,0.67)))/Rep_y;

		 else if(Rep_y>1000.0 && Rep_y <=2.e5) 
		 Cd_y=0.44;

		 else if(Rep_y>2.e5)
		 Cd_y=0.05;
         
		 Ap=((22.0/7.0)*(sigma*sigma))/4.0;

		 // fdrag_y[i]=Cd_y*(rel_vy*rel_vy)*den_gas*Ap/2.0; //in this expression, rel_vy and Ap is nondimensional
		 // fdrag_y[i]=fdrag_y[i]/part_mass;//here fdrag_y is not drag it is acceleration

		fdrag_y[i]=rel_vy/tau_vp; //here fdrag_y is not drag it is acceleration, use above two lines for Rep>1

                // fdrag_y[i] = 3*pi*vis_gas*rel_vy*sigma/part_mass;
 
		 if(p_vy[i]>vy_air)
		 fdrag_y[i]=- fdrag_y[i];
//--------------------------------------------------------//
//----------------------Reverse force---------------------//
                 freverse_y[i] = -6*pi*radius*rel_vy*b_ly*b_ly/Reynolds;
                    if(p_vy[i]>vy_air)
		       freverse_y[i]=-freverse_y[i];
//-------------------------------------------------------//
	//---------------------------------------------------------------//
//---------------------Torque on particles-----------------------//
                 torque_y[i] = 10*rel_angvely/(3*tau_vp);
                  if(p_angvely[i]>angvely_air)
                     torque_y[i] = -torque_y[i];
 
//----------------------------------------------------------------//
//----------------------------------------------------------------//
//--------------------Angular Acceleration------------------------//
               angaccy[i]=torque_y[i];
//----------------------------------------------------------------//

//-------------------------------------------------------------------//
//--------------------Reverse Torque---------------------------------//
                treverse_y[i] = 4*pi*radius*radius*radius*rel_angvely*b_ly*b_ly/Reynolds;
                    if(p_angvely[i]>angvely_air)
                      treverse_y[i] = -treverse_y[i];
	
		
              //--------------------------------------
                //Calculation of drag in the z-direction
                //---------------------------------------

  rel_vz=fabs(p_vz[i]-vz_air);    //V_R_Z=Relative velocity of the particle in z-direction
               
   Rep_z=(sigma*rel_vz*den_gas)/vis_gas; //particle Reynolds number 
    Rep_z=Rep_z*(length_scale*vel_scale);

    Total_Rep_z=Total_Rep_z+Rep_z;

    Cd_z=0;   //just to initiallize

       if(Rep_z<=1.0e-16) 
	  Cd_z=0.0;          
      if(Rep_z>=1.0e-16 && Rep_z <= 0.20)
		 Cd_z=24.0/Rep_z;

		 else if(Rep_z>0.20 && Rep_z <=1000.0)
		 Cd_z=24.0*(1.0+0.15*pow(Rep_z,0.67))/Rep_z;

		 else if(Rep_z>1000.0 && Rep_z <=2.e5) 
		 Cd_z=0.44;

		 else if(Rep_z>2.e5)
		 Cd_z=0.05;
		 
               Ap=((22.0/7.0)*(sigma*sigma))/4.0;

	       // fdrag_z[i]=Cd_z*(rel_vz*rel_vz)*den_gas*Ap/2.0;//in this expression, rel_vz and Ap is nondimensional
	       // fdrag_z[i]=fdrag_z[i]/part_mass;//here fdrag_y is not drag it is acceleration

	       fdrag_z[i]=rel_vz/tau_vp; //here fdrag_y is not drag it is acceleration

               // fdrag_z[i] = 3*pi*vis_gas*rel_vz*sigma/part_mass;

               if(p_vz[i]>vz_air)
		 fdrag_z[i]=- fdrag_z[i];
//---------------------------------------------------------------------//
//------------------------------Reverse force-------------------------//
                    freverse_z[i] = -6*pi*radius*rel_vz*b_ly*b_ly/Reynolds;
                    if(p_vz[i]>vz_air)
		       freverse_z[i]=-freverse_z[i];
//----------------------------------------------------------------------//
//---------------------------------------------------------------//
//---------------------Torque on particles-----------------------//
               torque_z[i] = 10*rel_angvelz/(3*tau_vp);
                     if(p_angvelz[i]>angvelz_air)
                     torque_z[i] = -torque_z[i];
//----------------------------------------------------------------//	
//----------------------------------------------------------------//
//--------------------Angular Acceleration------------------------//
                angaccz[i]=torque_z[i];
//----------------------------------------------------------------//

//-------------------------------------------------------------------//
//--------------------Reverse Torque---------------------------------//
                treverse_z[i] = 4*pi*radius*radius*radius*rel_angvelz*b_ly*b_ly/Reynolds;
               if(p_angvelz[i]>angvelz_air)
                     treverse_z[i] = -treverse_z[i];
          
		 return;
		 
	 }   
 

 void reverseforce_spectral(int i, FlowField& force)
{
              int Ny = force.Ny();             
                  /*     for (int ny=0; ny<Ny; ++ny)
                           {
                              frev_cheby_x[ny]  = 0.0;
                              frev_cheby_x[ny]  = 0.0; 
                              frev_cheby_x[ny]  = 0.0;
 
                              frev_cheby_diff2_x[ny] = 0.0;
                              frev_cheby_diff2_x[ny] = 0.0;
                              frev_cheby_diff2_x[ny] = 0.0;
                            }  */ 

         for(i=0;i<n_atom;++i)            
                 for (int ny=0; ny<Ny; ++ny) {

                     frev_cheby_x[i][ny] = freverse_x[i]*tmp3[i][ny];
                     frev_cheby_y[i][ny] = freverse_y[i]*tmp3[i][ny];
                     frev_cheby_z[i][ny] = freverse_z[i]*tmp3[i][ny];

                     /*frev_cheby_diff2_x[i][ny] =freverse_x[i]*tmp3[i][ny]*tmp_diff2[i][ny];
                     frev_cheby_diff2_y[i][ny] =freverse_y[i]*tmp3[i][ny]*tmp_diff2[i][ny];
                     frev_cheby_diff2_z[i][ny] =freverse_z[i]*tmp3[i][ny]*tmp_diff2[i][ny]; */    
       
                     frev_cheby_diff2_x[i][ny] =freverse_x[i]*tmp_diff2[i][ny];
                     frev_cheby_diff2_y[i][ny] =freverse_y[i]*tmp_diff2[i][ny];
                     frev_cheby_diff2_z[i][ny] =freverse_z[i]*tmp_diff2[i][ny];                    
                         }
                
            /*   for(int ny=0;ny<Ny;++ny)
                        {
                          frev_cheby_x[ny] /= n_atom;
                          frev_cheby_y[ny] /= n_atom;
                          frev_cheby_z[ny] /= n_atom;
                         
                          frev_cheby_diff2_x[ny]  /= n_atom;
                          frev_cheby_diff2_y[ny]  /= n_atom;
                          frev_cheby_diff2_z[ny]  /= n_atom;                          
                        }
                 */
             return;
 }


  

void reverseforce(int i,FlowField& force)

{

     int Nx = force.Nx();
     int Ny = force.Ny();
     int Nz = force.Nz();     
      
          freverse_x2=new double **[Nx]; freverse_y2=new double **[Nx]; freverse_z2=new double **[Nx];          
              for(int nx=0;nx<Nx;++nx){
          freverse_x2[nx]=new double *[Ny]; freverse_y2[nx]=new double *[Ny]; freverse_z2[nx]=new double *[Ny];    
                for(int ny=0;ny<Ny;++ny){
          freverse_x2[nx][ny]=new double [Nz]; freverse_y2[nx][ny]=new double [Nz]; freverse_z2[nx][ny]=new double [Nz];
                    }
                 }

       double factor = (Nx*Nz)/(force.Lx()*force.Lz());      

       reverse_quantity(i,force,frev_cheby_x, freverse_x2);
       reverse_quantity(i,force,frev_cheby_y, freverse_y2);
       reverse_quantity(i,force,frev_cheby_z, freverse_z2);


         force.setState(Physical,Spectral);

         for (int nx=0; nx<Nx; ++nx)
                for (int nz=0; nz<Nz; ++nz)
                   for (int ny=0; ny<Ny; ++ny) {                 
            force(nx,ny,nz,0) = freverse_x2[nx][ny][nz];
            force(nx,ny,nz,1) = freverse_y2[nx][ny][nz];
            force(nx,ny,nz,2) = freverse_z2[nx][ny][nz];
             
                        }
                 
     force.makeSpectral();

      for(int j=0;j<force.Nd();++j)
            for(int mx=0;mx<force.Mx();++mx)
               for(int mz=0;mz<force.Mz();++mz)
                  for(int ny=0;ny<Ny;++ny)
                     force.cmplx(mx,ny,mz,j) = factor* force.cmplx(mx,ny,mz,j);
  

return;

}


void third_order_stress_moments(int i,FlowField& force, FlowField& Dxx1, FlowField& Dyy1, FlowField& Dzz1)

{

     int Nx = force.Nx();
     int Ny = force.Ny();
     int Nz = force.Nz(); 

     int Mx = force.Mx();
     int Mz = force.Mz(); 

     

      freverse_x2_cheby_diff2=new double **[Nx]; freverse_y2_cheby_diff2=new double **[Nx]; freverse_z2_cheby_diff2=new double **[Nx];          
              for(int nx=0;nx<Nx;++nx){
          freverse_x2_cheby_diff2[nx]=new double *[Ny]; freverse_y2_cheby_diff2[nx]=new double *[Ny]; freverse_z2_cheby_diff2[nx]=new double *[Ny];    
                for(int ny=0;ny<Ny;++ny){
          freverse_x2_cheby_diff2[nx][ny]=new double [Nz]; freverse_y2_cheby_diff2[nx][ny]=new double [Nz]; freverse_z2_cheby_diff2[nx][ny]=new double [Nz];
                    }
                 }

              for (int nx=0; nx<Nx; ++nx)
                for (int nz=0; nz<Nz; ++nz)
                   for (int ny=0; ny<Ny; ++ny) {
                        freverse_x2_cheby_diff2[nx][ny][nz] = 0.0;
                        freverse_y2_cheby_diff2[nx][ny][nz] = 0.0;
                        freverse_z2_cheby_diff2[nx][ny][nz] = 0.0;
                             }

            reverse_quantity(i,force,frev_cheby_diff2_x,freverse_x2_cheby_diff2);
            reverse_quantity(i,force,frev_cheby_diff2_y,freverse_y2_cheby_diff2);
            reverse_quantity(i,force,frev_cheby_diff2_z,freverse_z2_cheby_diff2);

             double factor = (Nx*Nz)/(force.Lx()*force.Lz());

      Dxx1.setState(Physical,Spectral);
         
         for (int nx=0; nx<Nx; ++nx){
                 
                for (int nz=0; nz<Nz; ++nz){
                 
                   for (int ny=0; ny<Ny; ++ny) {
                   
            Dxx1(nx,ny,nz,0) = radius*radius*freverse_x2[nx][ny][nz]/3.0;
            Dxx1(nx,ny,nz,1) = radius*radius*freverse_y2[nx][ny][nz]/3.0;
            Dxx1(nx,ny,nz,2) = radius*radius*freverse_z2[nx][ny][nz]/3.0;
              
                        }
                  }
         }
                                                        
     Dxx1.makeSpectral();

     for(int j=0;j<force.Nd();++j)
            for(int mx=0;mx<Mx;++mx)
               for(int mz=0;mz<Mz;++mz)
                  for(int ny=0;ny<Ny;++ny)
                     Dxx1.cmplx(mx,ny,mz,j) = factor* Dxx1.cmplx(mx,ny,mz,j);


            for(int mx=0;mx<Mx;++mx){
                for(int mz=0;mz<Mz;++mz){
                     for(int ny=0;ny<Ny;++ny){
                         const channelflow::Complex Dx = Dxx1.Dx(mx,2);
         Dxx1.cmplx(mx,ny,mz,0)=Dx*Dxx1.cmplx(mx,ny,mz,0);
         Dxx1.cmplx(mx,ny,mz,1)=Dx*Dxx1.cmplx(mx,ny,mz,1);
         Dxx1.cmplx(mx,ny,mz,2)=Dx*Dxx1.cmplx(mx,ny,mz,2);
                             }
                          }
                       }   

          Dyy1.setState(Physical,Spectral);
         for (int nx=0; nx<Nx; ++nx)
                for (int nz=0; nz<Nz; ++nz)
                   for (int ny=0; ny<Ny; ++ny) {                  
            Dyy1(nx,ny,nz,0) = radius*radius*freverse_x2_cheby_diff2[nx][ny][nz]/3.0;
            Dyy1(nx,ny,nz,1) = radius*radius*freverse_y2_cheby_diff2[nx][ny][nz]/3.0;     
            Dyy1(nx,ny,nz,2) = radius*radius*freverse_z2_cheby_diff2[nx][ny][nz]/3.0;
              
                        }             
     Dyy1.makeSpectral();

      for(int j=0;j<force.Nd();++j)
            for(int mx=0;mx<Mx;++mx)
               for(int mz=0;mz<Mz;++mz)
                  for(int ny=0;ny<Ny;++ny)
                     Dyy1.cmplx(mx,ny,mz,j) = factor* Dyy1.cmplx(mx,ny,mz,j);

      

  Dzz1.setState(Physical,Spectral);  

         for (int nx=0; nx<Nx; ++nx)
                for (int nz=0; nz<Nz; ++nz)
                   for (int ny=0; ny<Ny; ++ny) {             
            Dzz1(nx,ny,nz,0) = radius*radius*freverse_x2[nx][ny][nz]/3.0;;
            Dzz1(nx,ny,nz,1) = radius*radius*freverse_y2[nx][ny][nz]/3.0;
            Dzz1(nx,ny,nz,2) = radius*radius*freverse_z2[nx][ny][nz]/3.0;
                        }        

     Dzz1.makeSpectral();

       for(int j=0;j<force.Nd();++j)
            for(int mx=0;mx<Mx;++mx)
               for(int mz=0;mz<Mz;++mz)
                  for(int ny=0;ny<Ny;++ny)
                     Dzz1.cmplx(mx,ny,mz,j) = factor* Dzz1.cmplx(mx,ny,mz,j);
                
                for(int mx=0;mx<Mx;++mx){
                  for(int mz=0;mz<Mz;++mz){
                     for(int ny=0;ny<Ny;++ny){
                         const channelflow::Complex Dz = Dzz1.Dz(mz,2);
         Dzz1.cmplx(mx,ny,mz,0)=Dz*Dzz1.cmplx(mx,ny,mz,0);
         Dzz1.cmplx(mx,ny,mz,1)=Dz*Dzz1.cmplx(mx,ny,mz,1);
         Dzz1.cmplx(mx,ny,mz,2)=Dz*Dzz1.cmplx(mx,ny,mz,2);
                             }
                          }
                       }

      
      for(int nx=0;nx<Nx;++nx){
           for(int ny=0;ny<Ny;++ny){
        delete [] freverse_x2[nx][ny];delete [] freverse_y2[nx][ny];delete [] freverse_z2[nx][ny];
                  }
        delete [] freverse_x2[nx]; delete [] freverse_y2[nx]; delete  [] freverse_z2[nx];
             }
        delete [] freverse_x2; delete [] freverse_y2; delete  [] freverse_z2; 
         
                  for(int nx=0;nx<Nx;++nx){
           for(int ny=0;ny<Ny;++ny){
        delete [] freverse_x2_cheby_diff2[nx][ny];delete [] freverse_y2_cheby_diff2[nx][ny];delete [] freverse_z2_cheby_diff2[nx][ny];
                  }
        delete [] freverse_x2_cheby_diff2[nx]; delete [] freverse_y2_cheby_diff2[nx]; delete  [] freverse_z2_cheby_diff2[nx];
             }
        delete [] freverse_x2_cheby_diff2; delete [] freverse_y2_cheby_diff2; delete  [] freverse_z2_cheby_diff2; 
     
return;

}


void rotational_stress_tensor(int i, FlowField& force, FlowField& dx1, FlowField& dy1, FlowField& dz1)
{
    
              int nx=0;
              int ny=0;
              
              int Nx = force.Nx();
              int Ny = force.Ny();
              int Nz = force.Nz();
              int Mx = force.Mx();
              int Mz = force.Mz();

            double **trev_cheby_x,**trev_cheby_y,**trev_cheby_z;
            double **trev_cheby_diff_x,**trev_cheby_diff_y,**trev_cheby_diff_z;
    
          trev_cheby_x = new double *[n_atom];trev_cheby_y = new double *[n_atom];trev_cheby_z = new double *[n_atom];
            trev_cheby_diff_x = new double *[n_atom];trev_cheby_diff_y = new double *[n_atom];trev_cheby_diff_z = new double *[n_atom];
            for(i=0;i<n_atom;++i)
             {
            trev_cheby_x[i] = new double [Ny];trev_cheby_y[i] = new double [Ny];trev_cheby_z[i] = new double [Ny];
            trev_cheby_diff_x[i] = new double [Ny];trev_cheby_diff_y[i] = new double [Ny];trev_cheby_diff_z[i] = new double [Ny];
             }
              

         double ***treverse_x2, ***treverse_y2, ***treverse_z2;
         double ***treverse_x2_cheby_diff, ***treverse_y2_cheby_diff, ***treverse_z2_cheby_diff;
 
         double factor = (Nx*Nz)/(force.Lx()*force.Lz());
          
          treverse_x2=new double **[Nx]; treverse_y2=new double **[Nx]; treverse_z2=new double **[Nx];          
              for(nx=0;nx<Nx;++nx){
          treverse_x2[nx]=new double *[Ny]; treverse_y2[nx]=new double *[Ny]; treverse_z2[nx]=new double *[Ny];    
                for(ny=0;ny<Ny;++ny){
          treverse_x2[nx][ny]=new double [Nz]; treverse_y2[nx][ny]=new double [Nz]; treverse_z2[nx][ny]=new double [Nz];
                    }
                 }

          treverse_x2_cheby_diff=new double **[Nx]; treverse_y2_cheby_diff=new double **[Nx]; treverse_z2_cheby_diff=new double **[Nx];          
              for(nx=0;nx<Nx;++nx){
          treverse_x2_cheby_diff[nx]=new double *[Ny]; treverse_y2_cheby_diff[nx]=new double *[Ny]; treverse_z2_cheby_diff[nx]=new double *[Ny];    
                for(ny=0;ny<Ny;++ny){
          treverse_x2_cheby_diff[nx][ny]=new double [Nz]; treverse_y2_cheby_diff[nx][ny]=new double [Nz]; treverse_z2_cheby_diff[nx][ny]=new double [Nz];
                    }
                 }

         for(i=0;i<n_atom;++i)            
              for(int ny=0;ny<Ny;++ny)
               {

                     trev_cheby_x[i][ny] = treverse_x[i]*tmp3[i][ny];
                     trev_cheby_y[i][ny] = treverse_y[i]*tmp3[i][ny];
                     trev_cheby_z[i][ny] = treverse_z[i]*tmp3[i][ny];

                     trev_cheby_diff_x[i][ny] =-treverse_x[i]*tmp_diff[i][ny];
                     trev_cheby_diff_y[i][ny] =-treverse_y[i]*tmp_diff[i][ny];
                     trev_cheby_diff_z[i][ny] =-treverse_z[i]*tmp_diff[i][ny];
                         
             }
                
             

      reverse_quantity(i,force,trev_cheby_x,treverse_x2);
      reverse_quantity(i,force,trev_cheby_y,treverse_y2);
      reverse_quantity(i,force,trev_cheby_z,treverse_z2);

      reverse_quantity(i,force,trev_cheby_diff_x,treverse_x2_cheby_diff);
      reverse_quantity(i,force,trev_cheby_diff_y,treverse_y2_cheby_diff);
      reverse_quantity(i,force,trev_cheby_diff_z,treverse_z2_cheby_diff);

       
           dx1.setState(Physical,Spectral);  

         for(int nx=0;nx<Nx;++nx)
                 for(int ny=0;ny<Ny;++ny)
                    for(int nz=0;nz<Nz;++nz){
                        dx1(nx,ny,nz,0) = 0.0;                       
                        dx1(nx,ny,nz,1) = treverse_z2_cheby_diff[nx][ny][nz];
                        dx1(nx,ny,nz,2) = -treverse_y2[nx][ny][nz];
            }

      dx1.makeSpectral();      

       for(int j=0;j<dx1.Nd();++j)
         for(int mx=0;mx<Mx;++mx)
               for(int mz=0;mz<Mz;++mz)
                  for(int ny=0;ny<Ny;++ny)
                     dx1.cmplx(mx,ny,mz,j) = factor* dx1.cmplx(mx,ny,mz,j); 

            for(int mx=0;mx<Mx;++mx)
               for(int mz=0;mz<Mz;++mz)
                  for(int ny=0;ny<Ny;++ny){
                       const channelflow::Complex Dz=dx1.Dz(mz);                        
                         dx1.cmplx(mx,ny,mz,2)=-Dz*dx1.cmplx(mx,ny,mz,2);
                             }

        
           dy1.setState(Physical,Spectral);  


           for(int nx=0;nx<Nx;++nx){
                 for(int ny=0;ny<Ny;++ny){
                    for(int nz=0;nz<Nz;++nz){                        
                        dy1(nx,ny,nz,0)=-treverse_z2[nx][ny][nz];
                        dy1(nx,ny,nz,1) = 0.0;
                        dy1(nx,ny,nz,2)=treverse_x2[nx][ny][nz];               
                }
              }
            }

         dy1.makeSpectral();    

          for(int j=0;j<dx1.Nd();++j)
             for(int mx=0;mx<Mx;++mx)
               for(int mz=0;mz<Mz;++mz)
                  for(int ny=0;ny<Ny;++ny)
                     dy1.cmplx(mx,ny,mz,j) = factor* dy1.cmplx(mx,ny,mz,j);                            
                              
                for(int mx=0;mx<Mx;++mx)
                   for(int mz=0;mz<Mz;++mz)
                     for(int ny=0;ny<Ny;++ny){
                         const channelflow::Complex Dx = dy1.Dx(mx);
                         const channelflow::Complex Dz = dy1.Dz(mz);
         dy1.cmplx(mx,ny,mz,0)=-Dx*dy1.cmplx(mx,ny,mz,0);
         dy1.cmplx(mx,ny,mz,2)=-Dz*dy1.cmplx(mx,ny,mz,2);
                             }
                          
       dz1.setState(Physical,Spectral); 
 
        for(int nx=0;nx<Nx;++nx)
                 for(int ny=0;ny<Ny;++ny)
                    for(int nz=0;nz<Nz;++nz){                        
                        dz1(nx,ny,nz,0)=treverse_y2[nx][ny][nz];
                        dz1(nx,ny,nz,1)=-treverse_x2_cheby_diff[nx][ny][nz];
                        dz1(nx,ny,nz,2) = 0.0;
                }
            
              dz1.makeSpectral();

            for(int j=0;j<dx1.Nd();++j)
              for(int mx=0;mx<Mx;++mx)
                for(int mz=0;mz<Mz;++mz)
                  for(int ny=0;ny<Ny;++ny)
                     dz1.cmplx(mx,ny,mz,j) = factor* dz1.cmplx(mx,ny,mz,j);
            
           for(int mx=0;mx<Mx;++mx)
                  for(int mz=0;mz<Mz;++mz)
                     for(int ny=0;ny<Ny;++ny){
                         const channelflow::Complex Dx = dz1.Dx(mx);                         
                         dz1.cmplx(mx,ny,mz,0)=-Dx*dz1.cmplx(mx,ny,mz,0);                                                  
                          }
                              
                   for(i=0;i<n_atom;++i){  
                  delete [] trev_cheby_x[i]; delete [] trev_cheby_y[i];  delete [] trev_cheby_z[i];
                  delete [] trev_cheby_diff_x[i]; delete [] trev_cheby_diff_y[i];  delete [] trev_cheby_diff_z[i];
                   }
	          delete [] trev_cheby_x; delete [] trev_cheby_y;  delete [] trev_cheby_z;
                  delete [] trev_cheby_diff_x; delete [] trev_cheby_diff_y;  delete [] trev_cheby_diff_z;           

                  for(int nx=0;nx<Nx;++nx){
           for(int ny=0;ny<Ny;++ny){
        delete [] treverse_x2[nx][ny];delete [] treverse_y2[nx][ny];delete [] treverse_z2[nx][ny];
                  }
        delete [] treverse_x2[nx]; delete [] treverse_y2[nx]; delete  [] treverse_z2[nx];
             }
        delete [] treverse_x2; delete [] treverse_y2; delete  [] treverse_z2;  
       
        for(int nx=0;nx<Nx;++nx){
           for(int ny=0;ny<Ny;++ny){
        delete [] treverse_x2_cheby_diff[nx][ny];delete [] treverse_y2_cheby_diff[nx][ny];delete [] treverse_z2_cheby_diff[nx][ny];
                  }
        delete [] treverse_x2_cheby_diff[nx]; delete [] treverse_y2_cheby_diff[nx]; delete  [] treverse_z2_cheby_diff[nx];
             }
        delete [] treverse_x2_cheby_diff; delete [] treverse_y2_cheby_diff; delete  [] treverse_z2_cheby_diff;  

          return;
        }


//**************************************************************************************************************************************//
//****************************************Symmetric stress tensor calc by Interpolation*************************************************//
//**************************************************************************************************************************************//

void symmetric_stress_tensor(int i,FlowField& u,FlowField& utmp,FlowField& force, FlowField& sx1, const Vector& x_grid, const Vector&  z_grid)
{
  
              int nx=0;
              int ny=0;
              int nz=0;
             
              int Nx = u.Nx();
              int Ny = u.Ny();
              int Nz = u.Nz();
             
              
              Real a = u.a();   
              Real b = u.b();   
              Real Lx = u.Lx();
              Real Lz = u.Lz(); 
              Vector x = periodicpoints(Nx, Lx);
              Vector y = chebypoints(Ny,a,b);
              Vector z = periodicpoints(Nz, Lz);

             
           double interpolated_sx[3][3];
           
 
           double ***symmetric_air_p_grid_xx,***symmetric_air_p_grid_xy,***symmetric_air_p_grid_xz,
                     ***symmetric_air_p_grid_yx,***symmetric_air_p_grid_yy,***symmetric_air_p_grid_yz,
                     ***symmetric_air_p_grid_zx,***symmetric_air_p_grid_zy,***symmetric_air_p_grid_zz;

           double **symmetric_cheby_diff_xx,**symmetric_cheby_diff_xy,**symmetric_cheby_diff_xz,
                     **symmetric_cheby_diff_yx,**symmetric_cheby_diff_yy,**symmetric_cheby_diff_yz,   
                     **symmetric_cheby_diff_zx,**symmetric_cheby_diff_zy,**symmetric_cheby_diff_zz;   

          double factor = (Nx*Nz)/(force.Lx()*force.Lz());  

          symmetric_air_p_grid_xx=new double **[Nx]; symmetric_air_p_grid_xy=new double **[Nx]; symmetric_air_p_grid_xz=new double **[Nx];
          symmetric_air_p_grid_yx=new double **[Nx]; symmetric_air_p_grid_yy=new double **[Nx]; symmetric_air_p_grid_yz=new double **[Nx];
          symmetric_air_p_grid_zx=new double **[Nx]; symmetric_air_p_grid_zy=new double **[Nx]; symmetric_air_p_grid_zz=new double **[Nx];
          
          
              for(nx=0;nx<Nx;++nx){
          symmetric_air_p_grid_xx[nx]=new double *[Ny]; symmetric_air_p_grid_xy[nx]=new double *[Ny]; symmetric_air_p_grid_xz[nx]=new double *[Ny];
          symmetric_air_p_grid_yx[nx]=new double *[Ny]; symmetric_air_p_grid_yy[nx]=new double *[Ny]; symmetric_air_p_grid_yz[nx]=new double *[Ny];
          symmetric_air_p_grid_zx[nx]=new double *[Ny]; symmetric_air_p_grid_zy[nx]=new double *[Ny]; symmetric_air_p_grid_zz[nx]=new double *[Ny];

          
                for(ny=0;ny<Ny;++ny){
          symmetric_air_p_grid_xx[nx][ny]=new double [Nz]; symmetric_air_p_grid_xy[nx][ny]=new double [Nz]; symmetric_air_p_grid_xz[nx][ny]=new double [Nz];
          symmetric_air_p_grid_yx[nx][ny]=new double [Nz]; symmetric_air_p_grid_yy[nx][ny]=new double [Nz]; symmetric_air_p_grid_yz[nx][ny]=new double [Nz];
          symmetric_air_p_grid_zx[nx][ny]=new double [Nz]; symmetric_air_p_grid_zy[nx][ny]=new double [Nz]; symmetric_air_p_grid_zz[nx][ny]=new double [Nz];
       
                    }
                 }

            
          symmetric_cheby_diff_xx = new double *[n_atom]; symmetric_cheby_diff_xy = new double *[n_atom]; symmetric_cheby_diff_xz = new double *[n_atom];
          symmetric_cheby_diff_yx = new double *[n_atom]; symmetric_cheby_diff_yy = new double *[n_atom]; symmetric_cheby_diff_yz = new double *[n_atom];   
          symmetric_cheby_diff_zx = new double *[n_atom]; symmetric_cheby_diff_zy = new double *[n_atom]; symmetric_cheby_diff_zz = new double *[n_atom];
            for(i=0;i<n_atom;++i){
          symmetric_cheby_diff_xx[i] = new double [Ny]; symmetric_cheby_diff_xy[i] = new double [Ny]; symmetric_cheby_diff_xz[i] = new double [Ny];
          symmetric_cheby_diff_yx[i] = new double [Ny]; symmetric_cheby_diff_yy[i] = new double [Ny]; symmetric_cheby_diff_yz[i] = new double [Ny];   
          symmetric_cheby_diff_zx[i] = new double [Ny]; symmetric_cheby_diff_zy[i] = new double [Ny]; symmetric_cheby_diff_zz[i] = new double [Ny];          
               }

    
      utmp = u;
     utmp.makeSpectral();

      ChebyCoeff U(Ny,a,b,Physical);
       

  for (int ny=0; ny<Ny; ++ny) 
    U[ny] = (1.0 - square(abs(y[ny]-(b+a)/2.0)/((b-a)/2.0)));
   
    U.makeSpectral();
   
      utmp += U;    

     FlowField sx(u.Nx(), u.Ny(), u.Nz(), 9, u.Lx(),
                     u.Lz(), u.a(), u.b(), Physical, Spectral);

     FlowField ucurl(u.Nx(), u.Ny(), u.Nz(), 9, u.Lx(),
                     u.Lz(), u.a(), u.b(), Physical, Physical);

     sx.setState(Spectral,Spectral);
       sx.setToZero();  

          ComplexChebyCoeff g(u.Ny(),u.a(),u.b(),Spectral);
          ComplexChebyCoeff gy(u.Ny(),u.a(),u.b(),Spectral);  
             
           for(int mx=0;mx<u.Mx();++mx){
              for(int mz=0;mz<u.Mz();++mz){
                for(int i=0;i<3;++i){ 
                 for(int ny=0;ny<u.Ny();++ny){
                       g.set(ny,utmp.cmplx(mx,ny,mz,i));                        
                        }
                           diff(g,gy);
                          for(int ny=0;ny<u.Ny();++ny){
                          sx.cmplx(mx,ny,mz,i3j(i,1))=gy[ny];  
                         
                       } 
                      }
                      
                     for(int ny=0;ny<u.Ny();++ny){ 
                       const channelflow::Complex Dx = u.Dx(mx);
                       const channelflow::Complex Dz = u.Dz(mz);
               sx.cmplx(mx,ny,mz,i3j(0,0)) = Dx*utmp.cmplx(mx,ny,mz,0);      
               sx.cmplx(mx,ny,mz,i3j(0,2)) = Dz*utmp.cmplx(mx,ny,mz,0);
               sx.cmplx(mx,ny,mz,i3j(1,0)) = Dx*utmp.cmplx(mx,ny,mz,1);                   //(i,j)--> i is the dir of velocity and j is the dir of derivatives              
               sx.cmplx(mx,ny,mz,i3j(1,2)) = Dz*utmp.cmplx(mx,ny,mz,1);
               sx.cmplx(mx,ny,mz,i3j(2,0)) = Dx*utmp.cmplx(mx,ny,mz,2);             
               sx.cmplx(mx,ny,mz,i3j(2,2)) = Dz*utmp.cmplx(mx,ny,mz,2);          
                     }
                 }
               }

         sx.makePhysical_xz();

         int n_order=5;
        
               for(i=0;i<n_atom;++i){
		 double xp=p_x[i];
		 double yp=p_y[i];
		 double zp=p_z[i];
            
       interpolation_tensor(sx, x_grid, z_grid,xp,yp,zp, n_order, interpolated_sx);

       symmetric_air_p_xx[i] = interpolated_sx[0][0];    
       symmetric_air_p_xy[i] = 0.5*(interpolated_sx[0][1]+interpolated_sx[1][0]);   //strain = E_ij = du_i/dx_j
       symmetric_air_p_xz[i] = 0.5*(interpolated_sx[0][2]+interpolated_sx[2][0]);     
       symmetric_air_p_yx[i] = 0.5*(interpolated_sx[1][0]+interpolated_sx[0][1]);
       symmetric_air_p_yy[i] = interpolated_sx[1][1];                                 //symmmetric part S_ij = 0.5(du_i/dx_j + du_j/dx_i)
       symmetric_air_p_yz[i] = 0.5*(interpolated_sx[1][2]+interpolated_sx[2][1]); 
       symmetric_air_p_zx[i] = 0.5*(interpolated_sx[2][0]+interpolated_sx[0][2]);     //antisymmetric part A_ij = 0.5(du_i/dx_j - du_j/dx_i)
       symmetric_air_p_zy[i] = 0.5*(interpolated_sx[2][1]+interpolated_sx[1][2]);    
       symmetric_air_p_zz[i] = interpolated_sx[2][2]; 
       
        }  

 
           
                    
                     for(nx=0;nx<Nx;++nx)
                        for(ny=0;ny<Ny;++ny)
                          for(nz=0;nz<Nz;++nz){
                               symmetric_air_p_grid_xx[nx][ny][nz]=0;
                               symmetric_air_p_grid_xy[nx][ny][nz]=0;
                               symmetric_air_p_grid_xz[nx][ny][nz]=0;
                               symmetric_air_p_grid_yx[nx][ny][nz]=0;
                               symmetric_air_p_grid_yy[nx][ny][nz]=0;
                               symmetric_air_p_grid_yz[nx][ny][nz]=0;
                               symmetric_air_p_grid_zx[nx][ny][nz]=0;
                               symmetric_air_p_grid_zy[nx][ny][nz]=0;
                               symmetric_air_p_grid_zz[nx][ny][nz]=0;
                                        }
                                      
                               
                    for(i=0;i<n_atom;++i)            
                       for(int ny=0;ny<Ny;++ny)
                        {
                     
                     symmetric_cheby_diff_xx[i][ny] = symmetric_air_p_xx[i]*tmp3[i][ny];
                     symmetric_cheby_diff_xy[i][ny] = symmetric_air_p_xy[i]*tmp_diff[i][ny];
                     symmetric_cheby_diff_xz[i][ny] = symmetric_air_p_xz[i]*tmp3[i][ny];
                     symmetric_cheby_diff_yx[i][ny] = symmetric_air_p_yx[i]*tmp3[i][ny];
                     symmetric_cheby_diff_yy[i][ny] = symmetric_air_p_yy[i]*tmp_diff[i][ny];
                     symmetric_cheby_diff_yz[i][ny] = symmetric_air_p_yz[i]*tmp3[i][ny];
                     symmetric_cheby_diff_zx[i][ny] = symmetric_air_p_zx[i]*tmp3[i][ny];
                     symmetric_cheby_diff_zy[i][ny] = symmetric_air_p_zy[i]*tmp_diff[i][ny];
                     symmetric_cheby_diff_zz[i][ny] = symmetric_air_p_zz[i]*tmp3[i][ny];
                         }

                     
                     reverse_quantity (i, force, symmetric_cheby_diff_xx,symmetric_air_p_grid_xx);
                     reverse_quantity (i, force, symmetric_cheby_diff_xy,symmetric_air_p_grid_xy);
                     reverse_quantity (i, force, symmetric_cheby_diff_xz,symmetric_air_p_grid_xz);
                     reverse_quantity (i, force, symmetric_cheby_diff_yx,symmetric_air_p_grid_yx);
                     reverse_quantity (i, force, symmetric_cheby_diff_yy,symmetric_air_p_grid_yy);
                     reverse_quantity (i, force, symmetric_cheby_diff_yz,symmetric_air_p_grid_yz);
                     reverse_quantity (i, force, symmetric_cheby_diff_zx,symmetric_air_p_grid_zx);
                     reverse_quantity (i, force, symmetric_cheby_diff_zy,symmetric_air_p_grid_zy);
                     reverse_quantity (i, force, symmetric_cheby_diff_zz,symmetric_air_p_grid_zz);

                    

                          
                             double preverse=0.0;
                             preverse = -20.0*pi*radius*radius*radius*b_ly*b_ly/(3.0*Reynolds);
                             
                             for(nx=0;nx<Nx;++nx)
                                for(ny=0;ny<Ny;++ny)
                                   for(nz=0;nz<Nz;++nz){
                                          symmetric_air_p_grid_xx[nx][ny][nz] = preverse*symmetric_air_p_grid_xx[nx][ny][nz];
                                          symmetric_air_p_grid_xy[nx][ny][nz] = preverse*symmetric_air_p_grid_xy[nx][ny][nz];
                                          symmetric_air_p_grid_xz[nx][ny][nz] = preverse*symmetric_air_p_grid_xz[nx][ny][nz];
                                          symmetric_air_p_grid_yx[nx][ny][nz] = preverse*symmetric_air_p_grid_yx[nx][ny][nz];
                                          symmetric_air_p_grid_yy[nx][ny][nz] = preverse*symmetric_air_p_grid_yy[nx][ny][nz];
                                          symmetric_air_p_grid_yz[nx][ny][nz] = preverse*symmetric_air_p_grid_yz[nx][ny][nz];
                                          symmetric_air_p_grid_zx[nx][ny][nz] = preverse*symmetric_air_p_grid_zx[nx][ny][nz];
                                          symmetric_air_p_grid_zy[nx][ny][nz] = preverse*symmetric_air_p_grid_zy[nx][ny][nz];
                                          symmetric_air_p_grid_zz[nx][ny][nz] = preverse*symmetric_air_p_grid_zz[nx][ny][nz];               
                                                 }

                            sx.setToZero();
                            sx.setState(Physical,Spectral);  

                for(int nx=0;nx<Nx;++nx)
                  for(int ny=0;ny<Ny;++ny) 
                      for(int nz=0;nz<Nz;++nz){
                            sx(nx,ny,nz,i3j(0,0))=symmetric_air_p_grid_xx[nx][ny][nz]; 
                            sx(nx,ny,nz,i3j(0,1))=symmetric_air_p_grid_xy[nx][ny][nz]; 
                            sx(nx,ny,nz,i3j(0,2))=symmetric_air_p_grid_xz[nx][ny][nz]; 
                            sx(nx,ny,nz,i3j(1,0))=symmetric_air_p_grid_yx[nx][ny][nz]; 
                            sx(nx,ny,nz,i3j(1,1))=symmetric_air_p_grid_yy[nx][ny][nz]; 
                            sx(nx,ny,nz,i3j(1,2))=symmetric_air_p_grid_yz[nx][ny][nz]; 
                            sx(nx,ny,nz,i3j(2,0))=symmetric_air_p_grid_zx[nx][ny][nz]; 
                            sx(nx,ny,nz,i3j(2,1))=symmetric_air_p_grid_zy[nx][ny][nz]; 
                            sx(nx,ny,nz,i3j(2,2))=symmetric_air_p_grid_zz[nx][ny][nz];
                                 }
                                  
                        sx.makeSpectral();

                     for(int j=0;j<force.Nd();++j)
                        for(int k=0;k<force.Nd();++k)
                           for(int mx=0;mx<force.Mx();++mx)
                               for(int mz=0;mz<force.Mz();++mz)
                                   for(int ny=0;ny<Ny;++ny)
                                       sx.cmplx(mx,ny,mz,i3j(j,k)) = factor*sx.cmplx(mx,ny,mz,i3j(j,k));

                        sx1.setState(Spectral,Spectral);  
                      
            for(int mx=0;mx<u.Mx();++mx)
               for(int mz=0;mz<u.Mz();++mz)
                  for(int ny=0;ny<Ny;++ny)
                     {
                       const channelflow::Complex Dx = sx.Dx(mx);
                       const channelflow::Complex Dz = sx.Dz(mz);
                        
              sx1.cmplx(mx,ny,mz,0) = -(Dx*sx.cmplx(mx,ny,mz,i3j(0,0))+ Dz*sx.cmplx(mx,ny,mz,i3j(0,2)) + sx.cmplx(mx,ny,mz,i3j(0,1)));
              sx1.cmplx(mx,ny,mz,1) = -(Dx*sx.cmplx(mx,ny,mz,i3j(1,0))+ Dz*sx.cmplx(mx,ny,mz,i3j(1,2)) + sx.cmplx(mx,ny,mz,i3j(1,1)));
              sx1.cmplx(mx,ny,mz,2) = -(Dx*sx.cmplx(mx,ny,mz,i3j(2,0))+ Dz*sx.cmplx(mx,ny,mz,i3j(2,2)) + sx.cmplx(mx,ny,mz,i3j(2,1)));

                
                     }
                      
           for(i=0;i<n_atom;++i){
          delete [] symmetric_cheby_diff_xx[i]; delete [] symmetric_cheby_diff_xy[i]; delete [] symmetric_cheby_diff_xz[i];
          delete [] symmetric_cheby_diff_yx[i]; delete [] symmetric_cheby_diff_yy[i]; delete [] symmetric_cheby_diff_yz[i];
          delete [] symmetric_cheby_diff_zx[i]; delete [] symmetric_cheby_diff_zy[i]; delete [] symmetric_cheby_diff_zz[i];    
                   }
          delete [] symmetric_cheby_diff_xx; delete [] symmetric_cheby_diff_xy; delete [] symmetric_cheby_diff_xz;
          delete [] symmetric_cheby_diff_yx; delete [] symmetric_cheby_diff_yy; delete [] symmetric_cheby_diff_yz; 
          delete [] symmetric_cheby_diff_zx; delete [] symmetric_cheby_diff_zy; delete [] symmetric_cheby_diff_zz;                

          for(int nx=0;nx<Nx;++nx){
           for(int ny=0;ny<Ny;++ny){
        delete [] symmetric_air_p_grid_xx[nx][ny];delete [] symmetric_air_p_grid_xy[nx][ny];delete [] symmetric_air_p_grid_xz[nx][ny];
        delete [] symmetric_air_p_grid_yx[nx][ny];delete [] symmetric_air_p_grid_yy[nx][ny];delete [] symmetric_air_p_grid_yz[nx][ny];
        delete [] symmetric_air_p_grid_zx[nx][ny];delete [] symmetric_air_p_grid_zy[nx][ny];delete [] symmetric_air_p_grid_zz[nx][ny];
                  }
        delete [] symmetric_air_p_grid_xx[nx]; delete [] symmetric_air_p_grid_xy[nx]; delete  [] symmetric_air_p_grid_xz[nx];
        delete [] symmetric_air_p_grid_yx[nx]; delete [] symmetric_air_p_grid_yy[nx]; delete  [] symmetric_air_p_grid_yz[nx];
        delete [] symmetric_air_p_grid_zx[nx]; delete [] symmetric_air_p_grid_zy[nx]; delete  [] symmetric_air_p_grid_zz[nx];
             }
        delete [] symmetric_air_p_grid_xx; delete [] symmetric_air_p_grid_xy; delete  [] symmetric_air_p_grid_xz;  
        delete [] symmetric_air_p_grid_yx; delete [] symmetric_air_p_grid_yy; delete  [] symmetric_air_p_grid_yz; 
        delete [] symmetric_air_p_grid_zx; delete [] symmetric_air_p_grid_zy; delete  [] symmetric_air_p_grid_zz;  
 
return;

}

//********************************************************************************************************************************************************//
//***************************************End of symmetric stress tensor calculation by interpolation******************************************************//
//********************************************************************************************************************************************************//

//********************************************************************************************************************************************************//
//****************************************Symmetric stress tensor by spectral method**********************************************************************//
//********************************************************************************************************************************************************//
/*void symmetric_stress_tensor_spectral_method(int i,FlowField& u,FlowField& sx,FlowField& sx1,const Vector& x_grid, const Vector&  z_grid)
 {
              int Mx = u.Mx();
              int Ny = u.Ny();
              int Mz = u.Mz();
              int Nx = u.Nx();
              int Nz = u.Nz();
              int ny = 0;
              
              Real a = u.a();   
              Real b = u.b();   
              Real Lx = u.Lx();
              Real Lz = u.Lz(); 
              Vector x = periodicpoints(Nx, Lx);
              Vector y = chebypoints(Ny,a,b);
              Vector z = periodicpoints(Nz, Lz);

             
           ChebyCoeff Ubase(Ny,a,b,Physical);
           double interpolated_sx[3][3];
           double symmetric_air_p_xx[n_atom],symmetric_air_p_xy[n_atom],symmetric_air_p_xz[n_atom],symmetric_air_p_yx[n_atom],symmetric_air_p_yy[n_atom],symmetric_air_p_yz[n_atom],
                  symmetric_air_p_zx[n_atom],symmetric_air_p_zy[n_atom],symmetric_air_p_zz[n_atom];

            int n_order=5;

            sigma=sqrt(sigma_sq);
		 double xp=p_x[i];
		 double yp=p_y[i];
		 double zp=p_z[i];
    
         for(ny=0;ny<Ny;++ny)
          Ubase[ny] =y(ny)-(a+b)/2.0;        
          
     sx.setState(Spectral,Spectral);
       sx.setToZero();  
          ComplexChebyCoeff g(Ny,a,b,Spectral);
          ComplexChebyCoeff gy(Ny,a,b,Spectral);  
             
           for(int mx=0;mx<Mx;++mx){
             for(int mz=0;mz<Mz;++mz){
               for(int j=0;j<3;++j){ 
                 for(ny=0;ny<Ny;++ny){
                       g.set(ny,u.cmplx(mx,ny,mz,j));                    
                        }
                          
                           diff(g,gy);
                          for(ny=0;ny<Ny;++ny){
                          sx.cmplx(mx,ny,mz,i3j(j,1))=gy[ny];  
                       } 
                      }
                      
                     for(ny=0;ny<Ny;++ny){ 
                       const channelflow::Complex Dx = u.Dx(mx);
                       const channelflow::Complex Dz = u.Dz(mz);
               sx.cmplx(mx,ny,mz,i3j(0,0)) = Dx*u.cmplx(mx,ny,mz,0);      
               sx.cmplx(mx,ny,mz,i3j(0,2)) = Dz*u.cmplx(mx,ny,mz,0);
               sx.cmplx(mx,ny,mz,i3j(1,0)) = Dx*u.cmplx(mx,ny,mz,1);                   //(i,j)--> i is the dir of velocity and j is the dir of derivatives              
               sx.cmplx(mx,ny,mz,i3j(1,2)) = Dz*u.cmplx(mx,ny,mz,1);
               sx.cmplx(mx,ny,mz,i3j(2,0)) = Dx*u.cmplx(mx,ny,mz,2);             
               sx.cmplx(mx,ny,mz,i3j(2,2)) = Dz*u.cmplx(mx,ny,mz,2);          
                     }
                  }
                
                }
       
         sx.makePhysical();
             
       interpolation_tensor(sx, x_grid, z_grid,xp,yp,zp, n_order, interpolated_sx);

       symmetric_air_p_xx[i] = interpolated_sx[0][ny][0];     
       symmetric_air_p_xy[i] = 0.5*(interpolated_sx[0][1]+interpolated_sx[1][0]);   //strain = E_ij = du_i/dx_j
       symmetric_air_p_xz[i] = 0.5*(interpolated_sx[0][2]+interpolated_sx[2][0]);     
       symmetric_air_p_yx[i] = 0.5*(interpolated_sx[1][0]+interpolated_sx[0][1]);
       symmetric_air_p_yy[i] = interpolated_sx[1][1];                                 //symmmetric part S_ij = 0.5(du_i/dx_j + du_j/dx_i)
       symmetric_air_p_yz[i] = 0.5*(interpolated_sx[1][2]+interpolated_sx[2][1]); 
       symmetric_air_p_zx[i] = 0.5*(interpolated_sx[2][0]+interpolated_sx[0][2]);     //antisymmetric part A_ij = 0.5(du_i/dx_j - du_j/dx_i)
       symmetric_air_p_zy[i] = 0.5*(interpolated_sx[2][1]+interpolated_sx[1][2]);    
       symmetric_air_p_zz[i] = interpolated_sx[2][2]; 

       
       sx.makeSpectral();
           
             double p=0.0;
             p = 20.0*pi*radius*radius*radius*vis_gas/3.0;

             for(ny=0;ny<Ny;++ny)
               for(int mx=0;mx<Mx;++mx)
                  for(int mz=0;mz<Mz;++mz)
                     tmp4[mx][ny][mz] *= (-p);
 
            for(i=0;i<n_atom;++i)
              for(ny=0;ny<Ny;++ny)
               for(int mx=0;mx<Mx;++mx)
                  for(int mz=0;mz<Mz;++mz){
                     sx.cmplx(mx,ny,mz,i3j(0,0))=symmetric_air_p_xx[i]*tmp4[mx][ny][mz];
                     sx.cmplx(mx,ny,mz,i3j(0,1))=symmetric_air_p_xy[i]*tmp4[mx][ny][mz];
                     sx.cmplx(mx,ny,mz,i3j(0,2))=symmetric_air_p_xz[i]*tmp4[mx][ny][mz];
                     sx.cmplx(mx,ny,mz,i3j(1,0))=sx.cmplx(mx,ny,mz,i3j(0,1));
                     sx.cmplx(mx,ny,mz,i3j(1,1))=symmetric_air_p_yy[i]*tmp4[mx][ny][mz];
                     sx.cmplx(mx,ny,mz,i3j(1,2))=symmetric_air_p_yz[i]*tmp4[mx][ny][mz];
                     sx.cmplx(mx,ny,mz,i3j(2,0))=sx.cmplx(mx,ny,mz,i3j(0,2));
                     sx.cmplx(mx,ny,mz,i3j(2,1))=sx.cmplx(mx,ny,mz,i3j(1,2));
                     sx.cmplx(mx,ny,mz,i3j(2,2))=symmetric_air_p_zz[i]*tmp4[mx][ny][mz];
                         }
                
             for(int k=0;k<sx.Nd()/3;++k)
               for(int j=0;j<sx.Nd()/3;++j)  
                for(ny=0;ny<Ny;++ny)
                  for(int mx=0;mx<Mx;++mx)
                    for(int mz=0;mz<Mz;++mz)             
                     sx.cmplx(mx,ny,mz,i3j(j,k)) /= (Nx*Nz*(Ny-1));


       sx1.setState(Spectral,Spectral);
           ComplexChebyCoeff h(Ny,a,b,Spectral);
           ComplexChebyCoeff hy(Ny,a,b,Spectral);
            for(int mx=0;mx<Mx;++mx){
               for(int mz=0;mz<Mz;++mz){
                  for(int ny=0;ny<Ny;++ny){
                       const channelflow::Complex Dx = sx.Dx(mx);
                       const channelflow::Complex Dz = sx.Dz(mz);
                        
              sx1.cmplx(mx,ny,mz,0) = -(Dx*sx.cmplx(mx,ny,mz,i3j(0,0))+Dz*sx.cmplx(mx,ny,mz,i3j(0,2)));
              sx1.cmplx(mx,ny,mz,1) = -(Dx*sx.cmplx(mx,ny,mz,i3j(1,0))+Dz*sx.cmplx(mx,ny,mz,i3j(1,2)));
              sx1.cmplx(mx,ny,mz,2) = -(Dx*sx.cmplx(mx,ny,mz,i3j(2,0))+Dz*sx.cmplx(mx,ny,mz,i3j(2,2)));
                           }
                       }
                  }
           for(int j=0;j<3;++j){   
            for(int mx=0;mx<Mx;++mx){
               for(int mz=0;mz<Mz;++mz){
                  for(int ny=0;ny<Ny;++ny){
                        h.set(ny,sx.cmplx(mx,ny,mz,i3j(j,1)));    
                      }
                      diff(h,hy);    
                   for(int ny=0;ny<Ny;++ny)
                        sx1.cmplx(mx,ny,mz,i) -= hy[ny];
                 } 
               }
             }
          
       sx1.makeSpectral();
     return;
       }
 
*/

//---------------------------------------------------------------------
  //  SUBROUTINE FOR COLLISION time CALCULATION
  //---------------------------------------------------------------
      double collision_time(int i,int j)
          {             
	    //FILE *fp56;   
	    // int     jx,jz;
	     
	     double  tmin= 1.e8;
	     double  vx_ij,vy_ij,vz_ij;
	     double  accx_ij,accy_ij, accz_ij, vel_dot_acc, accij_sq, vel_sq;
	     // double  rrx, rry, rrz;
	     double  rx, bx, bx1, rx_sq; 
	     double  ry, by, by1, ry_sq ;
	     double  rz,bz, bz1, rz_sq, b, b1, a_b1;
	     double  c,  ap1, bp1,cp1,dp1,ep1;
	     // double  projected_dist_x, projected_dist_y, projected_dist_z;
	     // double  projected_dist_sq;
	     

    double cubic_root(int i,int j,double bp1,double cp1,double dp1,double ep1);
    double quadratic_root(int i,int j,double cp1,double dp1,double ep1);
    double quartic_root(int i,int j,double ap1,double bp1,double cp1,double dp1,double ep1);

		 //  Difference in velocities in i,j pair
		 vx_ij= p_vxi-p_vx[j];
		 vy_ij= p_vyi-p_vy[j];
		 vz_ij= p_vzi-p_vz[j];
		 accx_ij=acc_x[i]-acc_x[j];
		 accy_ij=acc_y[i]-acc_y[j];
		 accz_ij=acc_z[i]-acc_z[j];
      
        vel_dot_acc=(vx_ij*accx_ij+ vy_ij* accy_ij+ vz_ij* accz_ij);  //  vel_dot_acc= VIJ_DOT_AIJ

		 accij_sq=(accx_ij*accx_ij+accy_ij*accy_ij+accz_ij*accz_ij);
      
		 vel_sq=vx_ij*vx_ij +vy_ij*vy_ij+vz_ij*vz_ij ;         //vel_sq=V_DOT_V
      

		 rx=p_xi-p_x[j];
		 ry=p_yi-p_y[j];
		 rz=p_zi-p_z[j];

		 if(rx > 0.50*b_lx) 
		   rx=rx-1.0*b_lx ; // !Closest image distance is that for collision
    
		 if (rz > 0.50*b_lz) 
		   rz=rz-1.0*b_lz;
      
		 if(rx < -0.50*b_lx)
		   rx=rx+1.0*b_lx;
   
		 if (rz < -0.50*b_lz)
		   rz=rz+1.0*b_lz;

		 bx=rx*vx_ij;
		 bx1=rx*accx_ij;
		 rx_sq=rx*rx;

		 by = ry*vy_ij;
		 by1=ry*accy_ij;
		 ry_sq= ry*ry;

		 bz = rz*vz_ij;
		 bz1= rz*accz_ij;
		 rz_sq= rz*rz;
		 
		 b= bx+ by +bz ;      //   b= R_DOT_U
		 b1=bx1+ by1 +bz1 ;   //   b1=R_DOT_A
		 a_b1=vel_sq+b1;


		 	 		 

		 /*
		 if(time_step>= 7374 && (i==869 && j==86))
	   
	      {
	    cout<<"::::::::::::in coll_time:::::::::::::::::::"<<endl;
	    c= rx_sq + ry_sq +rz_sq - sigma_sq;
 
	    ap1=accij_sq/4.0;
	    bp1=vel_dot_acc;
	    cp1=(vel_sq+b1);
	    dp1=(2.0*b); 
	    ep1=c;
	
   	      cout<<i<<" "<<j<<endl;
	      cout<<endl;
	printf("%40.30lf %40.30lf %40.30lf %40.30lf %40.30lf \n",ap1,bp1,cp1,dp1,ep1);
      cout<<"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<<endl;
      //cin.get();
	      }	

		 */
		 // if(time_step== 14248 && (i==1209 && j==242))cin.get();
		 
		 //  if((i==atom_a && j==atom_b)||(i==atom_b && j==atom_a))
		 //cout<<"##################b for atom_a&b=   "<<atom_a<<" "<<atom_b<<" "<<b<<endl;
		 

   //    if(b < 0.0 ||(projected_dist_sq -sigma_sq)<0.0)         //Eliminate i,j pair moving away from one another
 

		 //if(time_step>=33224 && i==1694 && j==464) cout<<"b="<<" "<<b<<endl;

	if(backward_index==1 && (i== overlaping_a||j==overlaping_a))
	cout<<"i and j are "<< i<<" "<<j <<" "<< cell_of_part[i]<<" "<<cell_of_part[j]<<endl;
		
 if(b < 0.0 ||(backward_index==1 && ((i==overlaping_a && j== overlaping_b)||(i==overlaping_b && j== overlaping_a)))) 
	    { 

            c= rx_sq + ry_sq +rz_sq - sigma_sq;

 if(c<0 && backward_index==1 && ((i==overlaping_a && j==overlaping_b)||(i==overlaping_b && j==overlaping_a)))
   //if(c<0 && backward_index==1 && ((i==previous_atom_a && j==previous_atom_b)||(i==previous_atom_b && j==previous_atom_a)))
{c=-c;       //--ibm
 cout<<" CHANGING SIGN OF C "<<endl;
 if(fabs(c)>1.e-16)
   {
   cout<<" colliding pair penetrated more in the last collision, so pausing"<<endl;
   cin.get();
     }
}


  //    if(c<0 && backward_index==1)c=-c; ///---alpha 

 // if(time_step>= 10854 &&  backward_index==1 && j==242)cout<<"------------   i and j-------- at 4 == "<<i<<" "<<j<<" "<<c<<endl;
	    ap1=accij_sq/4.0;
	    bp1=vel_dot_acc;
	    cp1=(vel_sq+b1);
	    dp1=(2.0*b); 
	    ep1=c;
	   


 	     if(c<0)
	      {
   cout<<"overlaping"<<c<<" "<<i<<" "<<j<<" "<< b <<" "<<overlaping_a <<" "<<overlaping_b <<" "<< backward_index <<endl;

   if(backward_index==1 && ((i==overlaping_a && j==overlaping_b)||(i==overlaping_b && j== overlaping_a)))
     cin.get();

   //above line means after executing the backward step also it is finding same overlaping pair till in overlap
   overlaping_a=i;
   overlaping_b=j;
   cout<<"overlaping_a and b= "<<overlaping_a<<" "<< overlaping_b<<" "<<"b=  "<<b<<endl;
   cout<<"sigma_sq=   "<<sigma_sq<<endl;

   backward_index=1;

   
   if( backward_index==1)
     {
cout<<"::::::::::::in coll_time:::::::::::::::::::"<<endl;
cout<<i<<" "<<j<<endl;
  cout<<endl;
	printf("%30.20lf %30.20lf %30.20lf %30.20lf %35.28lf \n",ap1,bp1,cp1,dp1,ep1);
cout<<"::::::::::::::::::::::::::::::::::::"<<endl;

     }
   
   backward_step(); // calling the function backward_step to sovle overlaping problem , that is including overlapping pair
   // in collision time calculation, which was missed in last calculation
   
   overlap=true;
   

   tmin=9999999.0; // just setting an arbitrary large number
   
   return(tmin);
   
	      }
 if( backward_index==1 && ((i==previous_atom_a && j==previous_atom_b)||(i==previous_atom_b && j==previous_atom_a)))	
   coll_time= quadratic_root(i,j,cp1,dp1,ep1);   //since for very closed particles cubic terms are negligble
 else
   {
    if(fabs(ap1/bp1)>1.e5) coll_time=quartic_root(i,j,ap1,bp1,cp1,dp1,ep1);	//accurecy reason
    else if(fabs(cp1/bp1)<1.e5)coll_time=cubic_root(i,j,bp1,cp1,dp1,ep1);
    else coll_time= quadratic_root(i,j,cp1,dp1,ep1);
   }


 if(coll_time<1.0e-9 && backward_index==1 && ((i==previous_atom_a && j==previous_atom_b)||(i==previous_atom_b && j==previous_atom_a)))
{
   cout<<coll_time<<endl;
   cout<<i<<" "<<j<<" "<< "verry low collision time"<< endl;
	      cout<<endl;
	printf("%40.30lf %40.30lf %40.30lf %40.30lf %40.30lf \n",ap1,bp1,cp1,dp1,ep1);
	cout<< "----------pausing-------------"<<endl;
	//cin.get();
 }    


/*if(fabs(cp1/bp1)>1.e6)
 coll_time= quadratic_root(i,j,cp1,dp1,ep1);
  else
 coll_time=cubic_root(i,j,bp1,cp1,dp1,ep1);	*/
		//calling cubic root

      
      
   //----------------added for debugging------
   
 //        if(time_step>= 19265 && i==3200 && j==234)

   if(backward_index==1 && ((i==overlaping_a && j== overlaping_b)||(i==overlaping_b && j== overlaping_a)))
	      {
	      cout<<"::::::::::::in coll_time:::::::::::::::::::"<<endl;

	
   	  
    cout<<i<<" "<<j<<endl;
	      cout<<endl;
	printf("%30.20lf %30.20lf %30.20lf %30.20lf %30.20lf %30.20lf\n",ap1,bp1,cp1,dp1,ep1,coll_time);
cout<<"::::::::::::::::::::::::::::::::::::"<<endl;
 cin.get();	
	      }


/*
				
	 if(time_step>=15971 && i==338 && j==3348)
     {
cout<<"::::::::::::in coll_time:::::::::::::::::::"<<endl;
cout<<i<<" "<<j<<endl;
  cout<<endl;
	printf("%30.20lf %30.20lf %30.20lf %30.20lf %35.28lf %30.20lf\n",ap1,bp1,cp1,dp1,ep1,coll_time);
cout<<"::::::::::::::::::::::::::::::::::::"<<endl;
if(time_step>= 15971)cin.get();	
     }

*/		
			 //----------------------------- 
      
   if(coll_time < tmin)
   tmin = coll_time; 

	    }
     
   return(tmin);
}

double wall_coll_time(int i,int j)
{
  // int sign;
  double  tmin_w=10.e6;   // initial setting of a large value
  double   sigma_by2, r_pw, rpw_acc, bw, discri,root;   
  
 sigma_by2=sqrt(sigma_sq)/2.0; ///sigma/2.0
 
 /*if(j==0)
 sign=-1;
 else
 sign=1;*/

 if(j==0)
  {
 r_pw=p_y[i]-p_y[n_atom+j];  //n_atom indicates left/bottom wall
 bw=r_pw*p_vy[i];
rpw_acc=r_pw*acc_y[i];
 discri= (p_vy[i]*p_vy[i])-2.0*acc_y[i]*(r_pw-sigma_by2);
 //if(bw<0.0 && discri >=0.0 )
 if((bw<0.0|| (rpw_acc<0.0 && abs(acc_y[i])>1.0e-14)) && discri >=0.0 )
    {
      if(abs(acc_y[i])>1.0e-14)
      root= (-p_vy[i]-sqrt(discri))/acc_y[i]; // for wall-particle collision quardratic equation is sufficient
      else
	root=-(r_pw-sigma_by2)/p_vy[i];

    if(root<=tmin_w)
      tmin_w=root;
    }
  }

 else
   {        // i.e for the other wall
  r_pw=p_y[i]-p_y[n_atom+j];  //n_atom indicates left wall
  bw=r_pw*p_vy[i];
  rpw_acc=r_pw*acc_y[i];
  discri= (p_vy[i]*p_vy[i])-2.0*acc_y[i]*(r_pw+sigma_by2);
  //if(bw<0.0 && discri >=0.0 )
  if((bw<0.0|| (rpw_acc<0.0 && abs(acc_y[i])>1.0e-14)) && discri >=0.0 )
    {

 if(abs(acc_y[i])>1.0e-14)
    root= (-p_vy[i]+sqrt(discri))/acc_y[i];
 else
	root=-(r_pw+sigma_by2)/p_vy[i];


    if(root<=tmin_w)
      tmin_w=root;
    }
   }


    return(tmin_w);



}




/*--------------------------------------------------------------------------------
 //			MOVE ALL THE SPHERE IN THE DURATION OF NEXT INSIDENCE
                  
 //-----------------------------------------------------------------------------*/
   void  advance_particle()
     {
       int i;

    for(i=0;i<n_atom;i++)
	{

	  /* before advancing storing the particle position and velocities to access in the next step */
	  
	   p_x_previous[i]= p_x[i];
	   p_y_previous[i]= p_y[i];
	   p_z_previous[i]= p_z[i];

	   p_vx_previous[i]=p_vx[i];
	   p_vy_previous[i]=p_vy[i];
	   p_vz_previous[i]=p_vz[i];
 
           p_angvelx_previous[i]=p_angvelx[i];
           p_angvely_previous[i]=p_angvely[i];
           p_angvelz_previous[i]=p_angvelz[i];
	   
	  p_x[i]=p_x[i]+p_vx[i]*t_next+0.50*acc_x[i]*(t_next*t_next); // Advance all spheres
          p_y[i]=p_y[i]+p_vy[i]*t_next+0.50*acc_y[i]*(t_next*t_next);
          p_z[i]=p_z[i]+p_vz[i]*t_next+0.50*acc_z[i]*(t_next*t_next);

              
//!!$ if(x(i)>2.*dlx .or.p_z(i)>2.*dlz)then
//!!$             print*,""
///!!$             print*,"########################"
//!!$             print*,""
//!!$              print*,i,"particle is crossing the (+ve)image box"
//!!$              print*,""
//!!$         !pause
//!!$              end if
//!!$              
//!!$           if(x(i)<-1.*dlx .or.p_z(i)<-1.*dlz)then
//!!$              print*,""
//!!$             print*,"########################"
//!!$             print*,""
//!!$              print*,i,"particle is crossing the (-ve)image box"
//!!$              print*,""
//!!$              print*,"########################"
//!!$           ! pause
//!!$              end if
//!!$ print*,"x,z,--before",i,x(i),z(i),dlx,dlz
//!!$         IF (X[i] > b_lx) X[i]=X[i]-b_lx           !apply periodic boundary condition
//!!$          !IF (p_y[i] > 1.D0) p_y[i]=p_y[i]-1.D0
//!!$          IF (p_z[i] > b_lz) p_z[i]=p_z[i]-b_lz
//!!$
//!!$          IF (X[i] < 0.D0) X[i]=X[i]+b_lx
//!!$         ! IF (p_y[i] < 0.D0) p_y[i]=p_y[i]+1.D0
//!!$          IF (p_z[i] < 0.D0) p_z[i]=p_z[i]+b_lz          
//!$print*,"x,z,-after",i,x(i),z(i),dlx,dlz
      //       ! print*,"x,z,--before",i,x(i),z(i),dlx,dlz

	  if(p_x[i] > b_lx) 
         p_x[i]=p_x[i]-(int(p_x[i]/b_lx)*b_lx) ;//apply periodic boundary condition
         
          if(p_z[i] > b_lz) 
	    p_z[i]=p_z[i]-(int(p_z[i]/b_lz)*b_lz);
	  
          if(p_x[i] < 0.0) 
	    p_x[i]=p_x[i]+(int(p_x[i]/b_lx)*b_lx+b_lx);
         
          if(p_z[i] < 0.0) 
	    p_z[i]=p_z[i]+(int(p_z[i]/b_lz)*b_lz+b_lz);
          
       
//!$         if( run_time>=1.2398E-002 .and. (i==97.or.i==98))then
//!!$print*,"***%%%%%%%%%%%%%%",i, X[i], p_y[i], p_z[i],tnext
	     
//!!$end if
          if((p_y[i]-(p_y[n_atom]+sqrt(sigma_sq)/2.0))<-10.e-9)
	    { 
	    cout<<" OVERLAP WITH BOTTOM WALL "<<i<<" "<<p_y[i]<<" "<<"allowed= "<<p_y[n_atom]+sqrt(sigma_sq)/2.0<<"  error= "<<(p_y[i]-(p_y[n_atom]+sqrt(sigma_sq)/2.0))<<endl;
	    cout<<"#################"<<"a[i]   "<<acc_y[i];
	  cin.get();
	    }
          if((-p_y[i]+(p_y[n_atom+1]-sqrt(sigma_sq)/2.0))<-10.0e-9 )
	    {
	    cout<<" OVERLAP WITH TOP WALL "<<i<<" "<<p_y[i]<<" "<<" allowed= "<<p_y[n_atom+1]-sqrt(sigma_sq)/2.0<<" error= " <<(-p_y[i]+(p_y[n_atom+1]-sqrt(sigma_sq)/2.0))<< endl;
	    cout<<"#################"<<"a[i]   "<<acc_y[i]<<" "<<(p_y[n_atom+1]-sqrt(sigma_sq)/2.0)<<endl;;
	  cin.get();
	    }

     if(t_next==del_tc)
	{		   
   if(i==atom_a && atom_b==n_atom && fabs(p_y[atom_a]-(p_y[n_atom]+ sqrt(sigma_sq)/2.0))>10.e-9)

     cout<<" COLIDING PARTICLE "<<atom_a<<" IS NOT REACHING TO THE BOTTOM WALL "<<"e rror= "<<fabs(p_y[atom_a]-(p_y[n_atom]+sqrt(sigma_sq)/2.0))<<endl;
               
           if(i==atom_a && atom_b==n_atom+1 && fabs(-p_y[atom_a]+(p_y[n_atom+1]-sqrt(sigma_sq)/2.0))>10.0e-9 )
			   
	     cout<<" COLIDING PARTICLE "<<atom_a<<" IS NOT REACHING TO THE TOP WALL " <<" error= "<<fabs(-p_y[atom_a]+(p_y[n_atom+1]-sqrt(sigma_sq)/2.0))<< endl;
		   }

          p_vx[i]=p_vx[i]+acc_x[i]*t_next;
          p_vy[i]=p_vy[i]+acc_y[i]*t_next;
          p_vz[i]=p_vz[i]+acc_z[i]*t_next;

          p_angvelx[i]=p_angvelx[i]+angaccx[i]*t_next;
          p_angvely[i]=p_angvely[i]+angaccy[i]*t_next;
          p_angvelz[i]=p_angvelz[i]+angaccz[i]*t_next;
            
	}




    return;
     }

/*-------------------------------------------------------------------------------------
   Bcakward movement of the particle before recalculating the collision time
   ----------------------------------------------------------------------------------*/

void  backward_step()
     {
       int i;
       Nof_backward_step=Nof_backward_step+1;
       cout<<" Executing backward time step "<<endl;
       
       run_time=run_time-t_next;
       time_step=time_step-1;
       
       if(event_noise_injection)t_noise=t_noise-del_tb;
 
    for(i=0;i<n_atom;i++)
	{
	  p_vx[i]=p_vx_previous[i];
          p_vy[i]=p_vy_previous[i];
          p_vz[i]=p_vz_previous[i];

	  p_x[i]=p_x_previous[i];
          p_y[i]=p_y_previous[i];
          p_z[i]=p_z_previous[i];


	}


    return;
     }

/*     For Linear Velocity collision         */
//----------------------------------------------------------------------------

//   DETERMINE CHANGE IN VELOCITY OF THE COLLIDING PAIR
//---------------------------------------------------------------------------    
   
 /*void post_coll_vel()
     
     {
       int i;
       double  rx, ry, rz,RR;
       double vx_ij,vy_ij,vz_ij,b, change_vx, change_vy, change_vz;// change in velocities after collision
       // double excess_value;
  
     if(atom_a<n_atom && atom_b<n_atom)  //For particle particle collision
       {
	 rx= p_x[atom_a]-p_x[atom_b];
	 ry= p_y[atom_a]-p_y[atom_b];
	 rz= p_z[atom_a]-p_z[atom_b];

     if(rx > 0.50*b_lx) 
       rx=rx-1.0*b_lx ; // !Closest image distance is that for collision
    
     if (rz > 0.50*b_lz) 
       rz=rz-1.0*b_lz;
      
     if(rx < -0.50*b_lx)
       rx=rx+1.0*b_lx;
   
     if (rz < -0.50*b_lz)
       rz=rz+1.0*b_lz;

     RR=rx*rx+ry*ry+rz*rz;
     
     //  cout<<"before change ---    "<<rx<<" "<<ry<<" "<<rz<<" "<<RR<<" "<<RR-sigma_sq<<endl;

        
     if( (RR-sigma_sq)<0.0) // if colliding-particle penetrates more for error in collision time it will adjust
       {
	 excess_value=fabs(RR-sigma_sq);
	 
	 if(rx>0.0)
	   p_x[atom_a]= p_x[atom_b]+sqrt(rx*rx+excess_value/3.0);
	   	  
	 else
	   p_x[atom_a]= p_x[atom_b]-sqrt(rx*rx+excess_value/3.0);
	  
	      
	 if(ry>0.0)
	   p_y[atom_a]= p_y[atom_b]+sqrt(ry*ry+excess_value/3.0); 
	 else
	   p_y[atom_a]= p_y[atom_b]-sqrt(ry*ry+excess_value/3.0); 
	
	      
	 
	 if(rz>0.0)
	   p_z[atom_a]= p_z[atom_b]+sqrt(rz*rz+excess_value/3.0);
	    
	 else
	   p_z[atom_a]= p_z[atom_b]-sqrt(rz*rz+excess_value/3.0);
	      

	 rx= p_x[atom_a]-p_x[atom_b];
	 ry= p_y[atom_a]-p_y[atom_b];
	 rz= p_z[atom_a]-p_z[atom_b];

     if(rx > 0.50*b_lx) 
       rx=rx-1.0*b_lx ; // !Closest image distance is that for collision
    
     if (rz > 0.50*b_lz) 
       rz=rz-1.0*b_lz;
      
     if(rx < -0.50*b_lx)
       rx=rx+1.0*b_lx;
   
     if (rz < -0.50*b_lz)
       rz=rz+1.0*b_lz;

     RR=rx*rx+ry*ry+rz*rz;
     //cout<<"after change ---    "<<rx<<" "<<ry<<" "<<rz<<" "<<RR<<" "<<RR-sigma_sq<<endl;
       }
     

     //cout<<"after change ---    "<<rx<<" "<<ry<<" "<<rz<<" "<<RR<<" "<<RR-sigma_sq<<endl;

   vx_ij = p_vx[atom_a]-p_vx[atom_b]; //! Differeces in sphere velocities
     vy_ij = p_vy[atom_a]-p_vy[atom_b];
     vz_ij = p_vz[atom_a]-p_vz[atom_b];
    
     b=(rx*vx_ij +ry*vy_ij +rz*vz_ij)/sigma_sq; // Project velocity differences in to line of centres
          
     

     change_vx = rx*b;  //       !Velocity change due to collision
     change_vy = ry*b;
     change_vz = rz*b;

    

     p_vx[atom_a]= p_vx[atom_a]-change_vx;
     p_vy[atom_a]= p_vy[atom_a]-change_vy;
     p_vz[atom_a]= p_vz[atom_a]-change_vz;

     p_vx[atom_b]= p_vx[atom_b]+change_vx;
     p_vy[atom_b]= p_vy[atom_b]+change_vy;
     p_vz[atom_b]= p_vz[atom_b]+change_vz;

     
       cout<<" ----------in change_vel------------"<<endl;
       cout<<endl;
       cout<<atom_a<<" "<<atom_b<<endl;
       cout<<RR-sigma_sq<<" "<<b<<" "<<change_vx<<endl;
       cout<<endl;
       cout<<" ----------in change_vel------------"<<endl;
       
       
     
        if(Tot_Nof_coll>=5528)
	  {
	    cout<< p_vx[atom_a] <<" "<<p_vy[atom_a] <<" "<<p_vz[atom_a] <<endl;
	    cout<< p_vx[atom_b] <<" "<< p_vy[atom_b]<<" "<< p_vy[atom_b]<<endl;

	    cout<<endl;
	    cout<<"-----------------------------"<<endl;
	    cout<< acc_x[atom_a]<<" "<< acc_y[atom_a]<<" "<<acc_z[atom_a]<<endl;

	  }
     
     


      
       }
	    
     else{   // for particle-wall collision
     
       p_vx[atom_a]=et* p_vx[atom_a]+(1-et)*p_vx[atom_b];
       p_vy[atom_a]=-en*p_vy[atom_a];
       p_vz[atom_a]=en*p_vz[atom_a];
     }
//! Update kinetic energy--------
   tot_KE=0.0;
        
   for(i=0;i<n_atom;i++)
          tot_KE= tot_KE+(p_vx[i]*p_vx[i]+p_vy[i]*p_vy[i]+p_vz[i]*p_vz[i]);
		      
		   tot_KE= 0.5*part_mass*tot_KE;


		   //cout<<"#################"<<run_time<<"   "<<tot_KE<<endl;

		   return;
		     		   
		   }*/

      /*          For Angular Velocity Collisions           */  

//----------------------------------------------------------------------------

//   DETERMINE CHANGE IN VELOCITY OF THE COLLIDING PAIR
//---------------------------------------------------------------------------    
   
 void post_coll_vel()
     
     {
       int i;
       double  rx, ry, rz,RR,b;
       double vx_ij,vy_ij,vz_ij, change_vx, change_vy, change_vz;// change in velocities after collision
       double angvelx_ij,angvely_ij,angvelz_ij, w,wx,wy,wz,change_angvelx, change_angvely, change_angvelz;// change in angular velocities after collision
       // double excess_value;
       double kappa=0.4;
          double rone = 1.0; double rtwo = 2.0;
    
    double etan = 0.0; double etat= 0.0;double etat1,kappa1;

   

     if(atom_a<n_atom && atom_b<n_atom)  //For particle particle collision
       {
	 rx= p_x[atom_a]-p_x[atom_b];
	 ry= p_y[atom_a]-p_y[atom_b];
	 rz= p_z[atom_a]-p_z[atom_b];

     if(rx > 0.50*b_lx) 
       rx=rx-1.0*b_lx ; // !Closest image distance is that for collision
    
     if (rz > 0.50*b_lz) 
       rz=rz-1.0*b_lz;
      
     if(rx < -0.50*b_lx)
       rx=rx+1.0*b_lx;
   
     if (rz < -0.50*b_lz)
       rz=rz+1.0*b_lz;

     RR=rx*rx+ry*ry+rz*rz;
     
     //  cout<<"before change ---    "<<rx<<" "<<ry<<" "<<rz<<" "<<RR<<" "<<RR-sigma_sq<<endl;

     /*   
     if( (RR-sigma_sq)<0.0) // if colliding-particle penetrates more for error in collision time it will adjust
       {
	 excess_value=fabs(RR-sigma_sq);
	 
	 if(rx>0.0)
	   p_x[atom_a]= p_x[atom_b]+sqrt(rx*rx+excess_value/3.0);
	   	  
	 else
	   p_x[atom_a]= p_x[atom_b]-sqrt(rx*rx+excess_value/3.0);
	  
	      
	 if(ry>0.0)
	   p_y[atom_a]= p_y[atom_b]+sqrt(ry*ry+excess_value/3.0); 
	 else
	   p_y[atom_a]= p_y[atom_b]-sqrt(ry*ry+excess_value/3.0); 
	
	      
	 
	 if(rz>0.0)
	   p_z[atom_a]= p_z[atom_b]+sqrt(rz*rz+excess_value/3.0);
	    
	 else
	   p_z[atom_a]= p_z[atom_b]-sqrt(rz*rz+excess_value/3.0);
	      

	 rx= p_x[atom_a]-p_x[atom_b];
	 ry= p_y[atom_a]-p_y[atom_b];
	 rz= p_z[atom_a]-p_z[atom_b];

     if(rx > 0.50*b_lx) 
       rx=rx-1.0*b_lx ; // !Closest image distance is that for collision
    
     if (rz > 0.50*b_lz) 
       rz=rz-1.0*b_lz;
      
     if(rx < -0.50*b_lx)
       rx=rx+1.0*b_lx;
   
     if (rz < -0.50*b_lz)
       rz=rz+1.0*b_lz;

     RR=rx*rx+ry*ry+rz*rz;
    //cout<<"after change ---    "<<rx<<" "<<ry<<" "<<rz<<" "<<RR<<" "<<RR-sigma_sq<<endl;
       }
     */

     //cout<<"after change ---    "<<rx<<" "<<ry<<" "<<rz<<" "<<RR<<" "<<RR-sigma_sq<<endl;

/*
   vx_ij = p_vx[atom_a]-p_vx[atom_b]; //! Differeces in sphere velocities
     vy_ij = p_vy[atom_a]-p_vy[atom_b];
     vz_ij = p_vz[atom_a]-p_vz[atom_b];
    
     b=(rx*vx_ij +ry*vy_ij +rz*vz_ij)/sigma_sq; // Project velocity differences in to line of centres
          
     

     change_vx = rx*b;  //       !Velocity change due to collision
     change_vy = ry*b;
     change_vz = rz*b;

    

     p_vx[atom_a]= p_vx[atom_a]-change_vx;
     p_vy[atom_a]= p_vy[atom_a]-change_vy;
     p_vz[atom_a]= p_vz[atom_a]-change_vz;

     p_vx[atom_b]= p_vx[atom_b]+change_vx;
     p_vy[atom_b]= p_vy[atom_b]+change_vy;
     p_vz[atom_b]= p_vz[atom_b]+change_vz;

     
       cout<<" ----------in change_vel------------"<<endl;
       cout<<endl;
       cout<<atom_a<<" "<<atom_b<<endl;
       cout<<RR-sigma_sq<<" "<<b<<" "<<change_vx<<endl;
       cout<<endl;
       cout<<" ----------in change_vel------------"<<endl;*/




   etan=(rone+en)/rtwo;
    etat1=(rone+et)/rtwo;
    kappa1=kappa/(rone+kappa);
    etat = etat1*kappa1;
    
     vx_ij = p_vx[atom_a]-p_vx[atom_b]; //! Differeces in sphere velocities
    vy_ij = p_vy[atom_a]-p_vy[atom_b];
     vz_ij = p_vz[atom_a]-p_vz[atom_b];
     
     angvelx_ij = (p_angvelx[atom_a]+p_angvelx[atom_b])*sigma/2.0;
     angvely_ij = (p_angvely[atom_a]+p_angvely[atom_b])*sigma/2.0;
     angvelz_ij = (p_angvelz[atom_a]+p_angvelz[atom_b])*sigma/2.0;
     
       rx=rx/sigma; ry=ry/sigma; rz=rz/sigma;  
    
      b=(rx*vx_ij +ry*vy_ij +rz*vz_ij); // Project velocity differences in to line of centres

      change_vx = etan*rx*b;
      change_vy = etan*ry*b;
      change_vz = etan*rz*b;

      change_vx = change_vx + etat*vx_ij- etat*rx*b;
      change_vy = change_vy + etat*vy_ij- etat*ry*b;
      change_vz = change_vz + etat*vz_ij- etat*rz*b;

      change_vx=change_vx-etat*(ry*angvelz_ij-rz*angvely_ij);
      change_vy=change_vy-etat*(rz*angvelx_ij-rx*angvelz_ij);
      change_vz=change_vz-etat*(rx*angvely_ij-ry*angvelx_ij);

     
     
     p_vx[atom_a]= p_vx[atom_a]-change_vx;
     p_vy[atom_a]= p_vy[atom_a]-change_vy;
     p_vz[atom_a]= p_vz[atom_a]-change_vz;

     p_vx[atom_b]= p_vx[atom_b]+change_vx;
     p_vy[atom_b]= p_vy[atom_b]+change_vy;
     p_vz[atom_b]= p_vz[atom_b]+change_vz;

     
     p_angvelx[atom_a]= p_angvelx[atom_a]-rtwo*(ry*change_vz-rz*change_vy)/(sigma*kappa);
     p_angvely[atom_a]= p_angvely[atom_a]-rtwo*(rz*change_vx-rx*change_vz)/(sigma*kappa);
     p_angvelz[atom_a]= p_angvelz[atom_a]-rtwo*(rx*change_vy-ry*change_vx)/(sigma*kappa);

     p_angvelx[atom_b]= p_angvelx[atom_b]-rtwo*(ry*change_vz-rz*change_vy)/(sigma*kappa);
     p_angvely[atom_b]= p_angvely[atom_b]-rtwo*(rz*change_vx-rx*change_vz)/(sigma*kappa);
     p_angvelz[atom_b]= p_angvelz[atom_b]-rtwo*(rx*change_vy-ry*change_vx)/(sigma*kappa);

/*
     wx=vx_ij-(ry*angvelz_ij-rz*angvely_ij);
     wy=vy_ij-(rz*angvelx_ij-rx*angvelz_ij);     
     wz=vz_ij-(rx*angvely_ij-ry*angvelx_ij);

      w=rx*wx+ry*wy+rz*wz; // Project velocity differences in to line of centres

     change_vx=et1*kappa*(wx-rx*w)/(1+kappa) + rx*w;

     change_vy=et1*kappa*(wy-ry*w)/(1+kappa)+ ry*w;

     change_vz=et1*kappa*(wz-rz*w)/(1+kappa) + rz*w;

     
change_angvelx=-rtwo*et1*(ry*wz-rz*wy)/(sigma*(kappa+1));
change_angvely=-rtwo*et1*(rz*wx-rx*wz)/(sigma*(kappa+1));
change_angvelz=-rtwo*et1*(rx*wy-ry*wx)/(sigma*(kappa+1));  

cout<<change_vx<<"     "<<change_vy<<"     "<<change_vz<<"    "<<change_angvelx<<"     "<<change_angvely<<endl;

     p_vx[atom_a]= p_vx[atom_a]-change_vx;
     p_vy[atom_a]= p_vy[atom_a]-change_vy;
     p_vz[atom_a]= p_vz[atom_a]-change_vz;

     p_vx[atom_b]= p_vx[atom_b]+change_vx;
     p_vy[atom_b]= p_vy[atom_b]+change_vy;
     p_vz[atom_b]= p_vz[atom_b]+change_vz;

     p_angvelx[atom_a]= p_angvelx[atom_a]+change_angvelx;
     p_angvely[atom_a]= p_angvely[atom_a]+change_angvely;
     p_angvelz[atom_a]= p_angvelz[atom_a]+change_angvelz;

     p_angvelx[atom_b]= p_angvelx[atom_b]+change_angvelx;
     p_angvely[atom_b]= p_angvely[atom_b]+change_angvely;
     p_angvelz[atom_b]= p_angvelz[atom_b]+change_angvelz;*/

  //   cout<<"change_vx="<<change_vx<<"    "<<"change_angvelx="<<change_angvelx<<endl;
    
          
/*
     
       cout<<" ----------in change_vel and change_angvel------------"<<endl;
       cout<<endl;
       cout<<atom_a<<" "<<atom_b<<endl;
       cout<<RR-sigma_sq<<" "<<change_vx<<"  "<<change_angvelx<<endl;
       cout<<endl;
       cout<<" ----------in change_vel and change_angvel------------"<<endl;
      
       */
     /* 
        if(Tot_Nof_coll>=5528)
	  {
	    cout<< p_vx[atom_a] <<" "<<p_vy[atom_a] <<" "<<p_vz[atom_a] <<endl;
	    cout<< p_vx[atom_b] <<" "<< p_vy[atom_b]<<" "<< p_vy[atom_b]<<endl;

	    cout<<endl;
	    cout<<"-----------------------------"<<endl;
	    cout<< acc_x[atom_a]<<" "<< acc_y[atom_a]<<" "<<acc_z[atom_a]<<endl;

	  }
     */
     


      
       }
	    
     else{   // for particle-wall collision
    
       p_vx[atom_a]=et* p_vx[atom_a]+(1-et)*p_vx[atom_b];
       p_vy[atom_a]=-en*p_vy[atom_a];
       p_vz[atom_a]=en*p_vz[atom_a];

    /*   p_vx[atom_a] = 5*(p_vx[atom_a] - 2*radius*p_angvelz[atom_a]/5)/7;
       p_vy[atom_a] = -en*p_vy[atom_a];
       p_vz[atom_a] = 5*(p_vz[atom_a] + 2*radius*p_angvelx[atom_a]/5)/7;

       p_angvelx[atom_a] = p_vz[atom_a]/radius;
       p_angvely[atom_a] = p_angvely[atom_a];
       p_angvelz[atom_a] = -p_vx[atom_a]/radius;*/

  /* double etw = 1.0;   double enw = 1.0;  
  
   etan=(rone+enw);
    etat1=(rone+etw);
    kappa1=kappa/(rone+kappa);
    etat = etat1*kappa1;          

     vx_ij = p_vx[atom_a]-p_vx[atom_b]; //! Differeces in sphere velocities
     vy_ij = p_vy[atom_a];
     vz_ij = p_vz[atom_a];
     
     angvelx_ij = (p_angvelx[atom_a])*sigma/2.0;
     angvely_ij = (p_angvely[atom_a])*sigma/2.0;
     angvelz_ij = (p_angvelz[atom_a])*sigma/2.0;
     
         rx= 0.0;
	 ry= (p_y[atom_b]-p_y[atom_a]);         
	 rz= 0.0;
     
      ry=ry/radius;  
   /*
     wx=vx_ij-(ry*angvelz_ij-rz*angvely_ij);
     wy=vy_ij-(rz*angvelx_ij-rx*angvelz_ij);     
     wz=vz_ij-(rx*angvely_ij-ry*angvelx_ij);
     
     w=rx*wx+ry*wy+rz*wz; // Project velocity differences in to line of centres


     b=(rx*vx_ij +ry*vy_ij +rz*vz_ij); // Project velocity differences in to line of centres
   
     change_vx = etan*rx*b;
      change_vy = etan*ry*b;
      change_vz = etan*rz*b;

      change_vx = change_vx + etat*vx_ij- etat*rx*b;
      change_vy = change_vy + etat*vy_ij- etat*ry*b;
      change_vz = change_vz + etat*vz_ij- etat*rz*b;

      change_vx=change_vx-etat*(ry*angvelz_ij-rz*angvely_ij);
      change_vy=change_vy-etat*(rz*angvelx_ij-rx*angvelz_ij);
      change_vz=change_vz-etat*(rx*angvely_ij-ry*angvelx_ij);




     p_vx[atom_a]= p_vx[atom_a]-change_vx;
     p_vy[atom_a]= p_vy[atom_a]-change_vy;
     p_vz[atom_a]= p_vz[atom_a]-change_vz;

   

  //   cout<<"change_vx="<<change_vx<<"    "<<"change_angvelx="<<change_angvelx<<endl;
     
     p_angvelx[atom_a]= p_angvelx[atom_a]-rtwo*(ry*change_vz-rz*change_vy)/(sigma*kappa);
     p_angvely[atom_a]= p_angvely[atom_a]-rtwo*(rz*change_vx-rx*change_vz)/(sigma*kappa);
     p_angvelz[atom_a]= p_angvelz[atom_a]-rtwo*(rx*change_vy-ry*change_vx)/(sigma*kappa);


  /* 
    change_vx=et1*kappa*(wx-rx*w)/(1+kappa) + rx*w;

     change_vy=et1*kappa*(wy-ry*w)/(1+kappa) + ry*w;

     change_vz=et1*kappa*(wz-rz*w)/(1+kappa) + rz*w;


change_angvelx=-2*et1*(ry*wz-rz*wy)/(sigma*(kappa+1));
change_angvely=-2*et1*(rz*wx-rx*wz)/(sigma*(kappa+1));
change_angvelz=-2*et1*(rx*wy-ry*wx)/(sigma*(kappa+1));  


     p_vx[atom_a]= p_vx[atom_a]+2*change_vx;
     p_vy[atom_a]= p_vy[atom_a]+2*change_vy;
     p_vz[atom_a]= p_vz[atom_a]+2*change_vz;
     
     p_angvelx[atom_a]= p_angvelx[atom_a]+2*change_angvelx;
     p_angvely[atom_a]= p_angvely[atom_a]+2*change_angvely;
     p_angvelz[atom_a]= p_angvelz[atom_a]+2*change_angvelz;*/


     }
     
//! Update kinetic energy--------
   tot_KE=0.0;
        
   for(i=0;i<n_atom;i++)
          tot_KE= tot_KE+(p_vx[i]*p_vx[i]+p_vy[i]*p_vy[i]+p_vz[i]*p_vz[i]);
		      
		   tot_KE= 0.5*part_mass*tot_KE;


		   //cout<<"#################"<<run_time<<"   "<<tot_KE<<endl;
//!  Update angular kinetic energy----------
   /*  tot_angularKE=0.0;
     
    for(i=0;i<n_atom;i++)
   tot_angularKE= tot_angularKE + (p_angvelx[i]*p_angvelx[i] + p_angvely[i]*p_angvely[i] + p_angvelz[i]*p_angvelz[i]);
              
          tot_angularKE= 0.5*moment_inertia*tot_angularKE;*/
                 
		   return;
		     		   
		   }

   //--------------------------------------------------------------------
   //    UPDATE/MAKE  COLLISION TABLE AFTER THE COLLISION
   //-----------------------------------------------------------------
    
void update_coll_table(const FlowField& u,const FlowField& omega,const Vector& x_grid, const Vector&  z_grid)

       {   
	 int i,j;
	 int i_cell,icell_zero, jcell_zero, jcell, nabour;
	 double  ctime_wall;

	restart_update:

	 overlap=false;

       for(i=0;i<n_atom;i++)
         tc[i] = 1.e10 ;   // Initially set elements of-
	                      // -table to arbitrarily heigh values
     
	   for(i= 0;i<n_atom;i++)
	     {
	     
	       drag(i,u,omega,x_grid,z_grid);                       //calling the function drag

	       /*
	       acc_x[i]=fdrag_x[i]/part_mass + g_x; // here part_mass== part_mass*, g_x==g_x*
	       acc_y[i]=fdrag_y[i]/part_mass + g_y;
 
	       acc_z[i]=fdrag_z[i]/part_mass + g_z;
	       */
	       acc_x[i]=fdrag_x[i] + g_x; // here fdrag is already acceleation
	       acc_y[i]=fdrag_y[i] + g_y;
 
	       acc_z[i]=fdrag_z[i] + g_z;

                 

	      
	     }
	
	   
	   links() ; // calling links to make the head of chain link list

	//   cout<<"link-- is completed"<<"time_step= "<<time_step<<" "<<endl;

	   icell_zero=Z*Z*(Z-1);

 for(i_cell=0;i_cell< Nof_lattice;i_cell++)

  {
    
      i= head[i_cell];

      
      
       repeat1:
       if ( i >= 0 )
            {
	       p_xi = p_x[i];
	       p_yi = p_y[i];
	       p_zi = p_z[i];
          
	       p_vxi= p_vx[i];
	       p_vyi= p_vy[i];
	       p_vzi= p_vz[i];

	       j = list[i];
	       
	         repeat2:  
	       if ( j >= 0 )

		 {
		   coll_time= collision_time(i,j);

		   if(overlap){
		     cout<<"is there overlap "<<overlap<<endl;
		     goto restart_update;
		   }
		   
		   if(coll_time < tc[i])
		  {
		tc[i]=coll_time;
		npair[i]=j;
		  }
                   if(coll_time < tc[j])
		     {
	           tc[j]=coll_time;
		   npair[j]=i;
		     }
		   
                      
		 
	       j = list[j];

                 goto repeat2;
		 }
 //         **********loop over neighbouring cell********

		 if(i_cell<Z*Z*(Z-1))  //------if the cell is not near the wall
		   {
		     jcell_zero = 13*(i_cell );

		     for(nabour=0;nabour<13;nabour++)

		  {

		  jcell = map [jcell_zero + nabour ];
		   j = head[jcell];
		  repeat3:   
		  if( j >= 0 ) 
		    {
		      coll_time= collision_time(i,j);
                if(overlap)goto restart_update;
	                
		      if(coll_time < tc[i])
			{
			  tc[i]=coll_time;
			  npair[i]=j;
			}
		            if(coll_time < tc[j])
			      {
				tc[j]=coll_time;
				npair[j]=i;
			      }
		   
                      
		    
			    
			    j = list[j];

			    goto repeat3;
		    }
		  }

		  if(i_cell<Z*Z)
		   {
		     j=0;  // wall at y=0
	   ctime_wall=wall_coll_time(i,j);  //!calling function to calculate time for particle--wall collision
		                       
	       if(ctime_wall < tc[i])        
                  {
		    tc[i]=ctime_wall;     
		    npair[i]=n_atom+j;
		  }
		   }
		  	  


	     	 i=list[i];
		 goto repeat1;
		   }
		 
		 else   //-------if the cell is near the wall
		   {
		     jcell_zero=Z*Z*(Z-1)*13;
		     jcell_zero=jcell_zero+(i_cell-icell_zero)*4;
		     
		     for(nabour=0;nabour<4;nabour++)

		       {
			 jcell = map [jcell_zero + nabour ];
			 j = head[jcell];
		  repeat4:   
		  if( j >= 0 ) 
		    {
		      
		      coll_time= collision_time(i,j);
		      if(overlap)goto restart_update;
		                     
                if(coll_time < tc[i])
		  {
		tc[i]=coll_time;
		npair[i]=j;
		  }
                   if(coll_time < tc[j])
		     {
	           tc[j]=coll_time;
		   npair[j]=i;
		     }
		   
                      
		    
	       j = list[j];
                 goto repeat4; 
		    }
		       }

      /* ------particle waqll collision---------*/

		     j=1;  // wall at y=0
 

	   ctime_wall=wall_coll_time(i,j);  //!calling function to calculate
		                        // time for particl--wall collision
	      if(ctime_wall < tc[i])        
                  {
		    tc[i]=ctime_wall;     
		    npair[i]=n_atom+j;
		  }


		   }


		 i=list[i];
		 goto repeat1;

	    }
  }


      backward_index=0;  
       return;
       }
           

/*-------------------------------------------------------------
 function links to make the head of chain link list
 ----------------------------------------------------------*/

 void links()              
{
 FILE *fp55;
  double celli, CELL;
  int   icell, i;
  // int   icell_242, icell_1209;   


      //C    ** ZERO head OF CHAIN ARRAY **

      for( icell = 0;icell< Nof_lattice;icell++)

	head[icell] = -1;


      celli = (double)Z ;
      CELL  = 1.0 / celli;

/*if( CELL< RCUT ) {

          cout<< "CELL SIZE TOO SMALL FOR CUTOFF "<<endl;
	exit(0);
      }
*/

      //   ** SORT ALL ATOMS **

      for(i=0;i<n_atom;i++)
	{

           icell = (int)( ( p_x[i]/b_lx ) * celli )
                       + (int)( ( p_y[i]/b_ly ) * celli ) * Z*Z
	                + (int) ( ( p_z[i]/b_lz)  * celli ) * Z ;

           list[i]     = head[icell];
           head[icell] = i;
	     }

      /*  if(Tot_Nof_coll>=5528){
cout<<"===========in links======================="<<endl;

  cout<<endl;

 icell_242 = (int)( ( p_x[242]/b_lx ) * celli )
                       + (int)( ( p_y[242]/b_ly ) * celli ) * Z*Z
	                + (int) ( ( p_z[242]/b_lz) * celli ) * Z ;

 icell_1209 = (int)( ( p_x[1209]/b_lx ) * celli )
                       + (int)( ( p_y[1209]/b_ly ) * celli ) * Z*Z
	                + (int) ( ( p_z[1209]/b_lz)  * celli ) * Z ;

 cout<< p_x[242]<<" "<<p_y[242]<<" "<< p_z[242]<<endl;
 cout<< p_vx[242]<<" "<<p_vy[242]<<" "<< p_vz[242]<<endl;
 cout<< acc_x[242]<<" "<<acc_y[242]<<" "<< acc_z[242]<<endl;
 cout<<endl;
 cout<< p_x[1209]<<" "<<p_y[1209]<<" "<< p_z[1209]<<endl;
 cout<< p_vx[1209]<<" "<<p_vy[1209]<<" "<< p_vz[1209]<<endl;
 cout<< acc_x[1209]<<" "<<acc_y[1209]<<" "<< acc_z[1209]<<endl;
 cout<<endl;

 cout<< icell_242<<" "<< icell_1209<<" "<<head[icell_242]<<endl;
 cout<<"================================================="<<endl;
 }*/

      //  if(time_step==5824)

      cell_of_part[i]=(int)( ( p_x[list[i]]/b_lx ) * celli )
                       + (int)( ( p_y[list[i]]/b_ly ) * celli ) * Z*Z
                          + (int) ( ( p_z[list[i]]/b_lz)  * celli ) * Z;



 if(time_step>=15970)
	{
 fp55=fopen("link1.txt","w");
 //cout<<" ------writing in the file"<<endl;
 for(i=0;i<n_atom;i++)
	{

  cell_of_part[i]=(int)( ( p_x[list[i]]/b_lx ) * celli )
                       + (int)( ( p_y[list[i]]/b_ly ) * celli ) * Z*Z
                          + (int) ( ( p_z[list[i]]/b_lz)  * celli ) * Z;
 
   fprintf(fp55,"%10d %10d %10d %10d\n",i,list[i], cell_of_part[i] , head[cell_of_part[i]]);
	  
   //cout<<cell_of_part[i]<<" "<<head[cell_of_part[i]]<<endl;
	}
 fclose(fp55);

 
	}

      return;
}



//-------------------------------------------------------------
  //   !GIVING BROWNIAN NOISE TO THE PARTICLES
//----------------------------------------------------------
     void adding_noise()
	 {
	   int i;
	   double yp;
	   //double diff_p[5]={0.0,0.0,0.0,0.0,0.0};
	   double diff_px, diff_py ,diff_pz,diff_pxy;
	   double sigma_u,sigma_v, sigma_w,corre_coeff_uv;
	  
	   double gita1,gita2,gita3;
	   double fatom;
	   /*double fnoise1_x, fnoise1_y, fnoise1_z;
	   double  xx, yy, zz,fatom;
	   double fnoise_x,fnoise_y,fnoise_z;
	   */
	   double gaussrand(void);
	   void   spline(double yp,double diff_p[]);

 for(i=0;i<n_atom;i++)
 {                 
                if(p_y[i]>= 0.50*b_ly)
                  yp=p_y[i];
                else
                   yp=b_ly-p_y[i];
		
		/*  
		spline(yp,diff_p);    //calling spline curve fitting
            		

             diff_px=diff_p[1];
             diff_py=diff_p[2];
             diff_pz=diff_p[3];
             diff_pxy=diff_p[4];
		*/
	     diff_px=const_diff_x;
             diff_py=const_diff_y;
             diff_pz=const_diff_z;
             diff_pxy=const_diff_xy;
	     

      sigma_u=sqrt(2*fabs(diff_px)*del_tb); // using fabs is not authentic but used here because some diffusivity is unexpectedly -ve due to the erros in model.
      sigma_v=sqrt(2*fabs(diff_py)*del_tb);
      sigma_w=sqrt(2*fabs(diff_pz)*del_tb);

      corre_coeff_uv=(diff_pxy/sqrt(fabs(diff_px*diff_py)));

      gita1=gaussrand();
      u_noise=sigma_u*gita1;
                        //gita1,2,3 are gaussian rand of zero mean and unit varience// 
      gita2=gaussrand();
      v_noise=sigma_v*(corre_coeff_uv*gita1+(sqrt(1.0-corre_coeff_uv*corre_coeff_uv))*gita2);

      gita3=gaussrand();
      w_noise=sigma_w*gita3;
      
     
      // cout<<"*************** in noise injection**********"<< p_vx[i]<<endl;
         p_vx[i]=p_vx[i]+u_noise;
         p_vy[i]=p_vy[i]+v_noise;
         p_vz[i]=p_vz[i]+w_noise;

	 //cout<<"*************** in noise injection**********"<< p_vx[i]<<endl;
      
	 // cout<<"in ADD NOISE()"<< "  "<<diff_px<< "  "<<diff_py<<" "<<diff_pz<< " "<<yp;

	 //	 cin.get();

 }
      
       

         tot_KE =0.0;
         fatom=(double)(n_atom);
      
			 for(i=0;i<n_atom;i++){

//!!$        p_vx[i]= p_vx[i]-(SUM_NOISE_X)/fatom  !Adjust linear velocities so that
//!!$        p_vy[i]= p_vy[i]-(SUM_NOISE_Y)/fatom   !total linear momentum is zero.
//!!$        p_vz[i]= p_vz[i]-(SUM_NOISE_z)/fatom
//!!$

            tot_KE=tot_KE + p_vx[i]* p_vx[i]+p_vy[i]* p_vy[i]+p_vz[i]* p_vz[i];
         
			 }

            tot_KE=0.50*part_mass*tot_KE;
            return;
	 }


double gaussrand(void)
{
  static int iset=0;
  static double gset;
  double gasdev,fac,rsq,v1,v2;
  //double ran1;
 
  if (iset==0) {

  redo:
     v1 = 2.0*((double)rand()/RAND_MAX)-1.0;
     v2 = 2.0*((double)rand()/RAND_MAX)-1.0;
     //cout<< v1<<" "<<v2;
    // v1 = 2*drand48()-1.0;
    //v2 = 2*drand48()-1.0;

    rsq = v1*v1 + v2*v2;
    if(rsq >= 1 || rsq ==0) goto redo;
        
    fac = sqrt(-2*log(rsq)/rsq);
    gset=v1*fac;
    gasdev=v2*fac;
    iset=1;
  }
  else {
    gasdev = gset;
    iset = 0;
  }
  
    
  return gasdev;
}


 //-------------------------------------------------------

    void spline(double yp,double diff_p[] )

	{
	  int npt_spline= N_diff+1; //total no of points
 
  int l,i,j,m;
  double diff[diff_nmax],u[diff_nmax], a[diff_nmax],d[diff_nmax],c[diff_nmax][diff_nmax],
	h[diff_nmax],D_diff[diff_nmax];
  double q1,q2,q3;

  void  gauss( int m, double c[diff_nmax][diff_nmax],double d[diff_nmax],double a[diff_nmax]);	
  //void  gauss(int , int, double ,double ,double);					       
  
       
// Compute distances between data points and function differences
  
// for(i=0;i<npt_spline;i++)
  // cout<<diff_x[i]<<"  "<<diff_y[i]<<"  "<<diff_z[i]<<" "<<y_ch[i]<<endl;


             for(l=1;l<5;l++)
	        {
		 if (l == 1)
		   {
                   for(i=0;i<npt_spline;i++)
                      diff[i]=diff_x[i];
		   }
                else if(l == 2)
		  {
                      for(i=0;i<npt_spline;i++)
                      diff[i]=diff_y[i];
		  }
		 else if(l==3)
		   { 
                      for(i=0;i<npt_spline;i++)
                      diff[i]=diff_z[i];
		   }
		 else
		   {
		     for(i=0;i<npt_spline;i++)
                      diff[i]=diff_xy[i];
		   }
		 for(i=1;i<npt_spline;i++)  // here i will be 1 to n-1
		 {
		   h[i]=y_ch[i]-y_ch[i-1];
		   D_diff[i]=diff[i]-diff[i-1];
		   
		   }
		
                
//  Initialise c matrix
                  
                  for(i=1;i<npt_spline-1;i++)
		    {
                    for(j=1;j<npt_spline-1;j++)
                        c[i][j]=0.0;
		    }
		 
//     Compute diagonal element of c

                       for(i=1;i<npt_spline-1;i++)
			 {
                           c[i][i]=2.0*(h[i]+h[i+1]);
                           //!write(*,*)" c(i,i)=", c(i,i)
			 }
//     Compute off diagonal elements of c

                          for(i=2;i<npt_spline-1;i++)
			    {
                              c[i-1][i]=h[i];
                              c[i][i-1]=h[i];
			    }   
 
//     Compute elements of d array
                              
			  for(i=1;i<npt_spline-1;i++)
			    {
			      d[i]=(D_diff[i+1]/h[i+1]-D_diff[i]/h[i])*6.0;
			      
			    }
//!     Compute elements of a using Gaussian elemination
//!     Change array subscripts from 2 to n-1 to 1 to n-1
//!before calling gauss

	  m=npt_spline-2;
                for(i=0;i<m;i++)
		  {
		    d[i]=d[i+1];
		    for(j=0;j<m;j++)
		      {
			c[i][j]=c[i+1][j+1];
			
			//	cout<<c[i][j]<<" "<<d[i]<<endl;   //---------------------
		      }
		  }
                   
		//    CALLING gauss--------------------------

                      gauss(m,c,d,a);
// compute the coefficients of natural cubic spline
                    
                 for(i=npt_spline-2;i>=1;i--)
		  {
		  a[i]=a[i-1];
		  }
						  
		 a[0]=0.0;
                 a[npt_spline-1]=0.0;

// Locate the domain of  Y_CHp

		 for(i=1;i<npt_spline;i++)
		   {if(yp<=y_ch[i])
		     break;
		   }
                         
		// 	i=1;
// 		       while( y_ch[i]>=yp)
// 				 i=i+1;
                        
	   //     Compute interpolation value at yp

	      u[i-1]=yp -y_ch[i-1];
	      u[i]=yp-y_ch[i];
	      q1=(h[i]*h[i]) * u[i]-(u[i]*u[i]*u[i]);
              q2=u[i-1]*u[i-1]*u[i-1]-(h[i]*h[i])*u[i-1];
              q3=diff[i]*u[i-1]-diff[i-1]*u[i];
              
				 //write(*,*)a(i-1),q1,a[i],q2,h[i],q3
             
               diff_p[l]=(a[i-1]*q1+a[i]*q2)/(6.0*h[i])+q3/h[i];

		}
	     return;
	}

//!-----------------------------------------------------------

//___________________________________________________________            

void gauss(int n,double a[diff_nmax][diff_nmax],double b[diff_nmax],double a1[diff_nmax])
  //int NMAX, n;
  //double a[][NMAX];
  //double b[NMAX];
  //double y_ch[NMAX];

{ 
    int k,i,j;
    double PIVOT,FACTOR,SUM;
    
//      Elimination begions
      for( k=0;k<n-1;k++)
      {
              
	PIVOT=a[k][k];

               
	for(i=k+1;i<n;i++)
	  {
                 
	    FACTOR = a[i][k]/PIVOT;
		    for(j=k+1;j<n;j++)
		      {
		      a[i][j] =a[i][j] - FACTOR * a[k][j];
		      }
                      
		    b[i]=b[i]-FACTOR*b[k];
			//write(*,*)"********",a(i,j), b[i]
 
			}
      }
   
 
    //      Back substitution begins
                       
    a1[n-1]=b[n-1]/a[n-1][n-1];
    //cout<<"in //gauss()"<< "  "<< a1[n]<< "  "<< a[n-1][n-1]<<"  "<< b[n-1]<< endl;
 			 for(k=n-2;k>=0;k--)
			   {
			    
			     SUM=0;
			     for(j=k+1;j<n;j++)
			     {
                             SUM=SUM +a[k][j]*a1[j];

			     }

                             a1[k]=(b[k]-SUM)/a[k][k];
                             //cout<<"in //gauss()"<< "  "<< a1[k]<<endl;
			   }

		       return;
  }


// -----------function for spatial interpolation of air field------ 



void interpolation(const FlowField& u, const Vector& x_grid, const Vector&  z_grid,double xp,double yp, double zp, int n_order, double u_interpolated[] )
{
//int n_order is the order of interpolation
int i,j,k;
int nx_part, nz_part;
int min_limit=((n_order+1)/2)-1;    // n_order is the order of poly int.
int max_limit=(n_order+1)/2;
int n_point=n_order+1;
int Nx_ = u.Nx();
int Ny_ = u.Ny();
int Nz_ = u.Nz();
int Nd_=u.Nd();
double a_=u.a();
double b_=u.b();


double *scratch_part_,*x,*z;
double **rdata_for_interpolation;
//double u_interpolated[3];              //3 is for 3 comp  of velocities
double dy[3];
double delta_x=x_grid[1]-x_grid[0];
double delta_z=z_grid[1]-z_grid[0];


int locate(const Vector&, int, double);
double eval_at_yp(double *scratch_part_,double yp, double a_, double b_,int Ny_);
void poli_interpolation_2d(double [], double [], double **, int , int , double , double , double *,double *);

scratch_part_= new double [Ny_];
x= new double [n_point];
z= new double [n_point];



rdata_for_interpolation=new double *[n_point];
for(i=0;i<n_point;i++) rdata_for_interpolation[i]=new double [n_point];

nx_part=locate(x_grid,Nx_,xp);
nz_part=locate(z_grid,Nz_,zp);

//cout<<" nx_part= "<<nx_part<<" "<<nz_part<<endl;

  for (int i=0; i<Nd_; ++i){ 
    for (int nx=nx_part-min_limit; nx<=nx_part+max_limit; ++nx){
   j=nx-(nx_part-min_limit);
   x[j]=nx*delta_x;
   int nx1=nx;
   if(nx1<0)nx1=nx+(Nx_-1);if(nx1>(Nx_-1))nx1=nx-(Nx_-1);        //--check -1 wont be there.

   	for (int nz=nz_part-min_limit; nz<=nz_part+max_limit; ++nz){
k=nz-(nz_part-min_limit);
z[k]=nz*delta_z;
int nz1=nz;
if(nz1<0)nz1=nz+(Nz_-1);if(nz1>(Nz_-1))nz1=nz-(Nz_-1);


            for (int ny=0; ny<Ny_; ++ny)
			{ 
			
	scratch_part_[ny]= u(nx1,ny,nz1,i);
			}    	// end loop for ny

   rdata_for_interpolation[j][k]=eval_at_yp(scratch_part_,yp,a_,b_,Ny_);

						   } 
						}

poli_interpolation_2d( x, z, rdata_for_interpolation,n_point,n_point,xp,zp,&u_interpolated[i],&dy[i] );

			   }

delete [] x;
delete [] z;
delete [] scratch_part_;

for(i=0;i<n_point;i++) delete [] rdata_for_interpolation[i];
delete [] rdata_for_interpolation;

}


void interpolation_tensor(FlowField& sx, const Vector& x_grid, const Vector&  z_grid,double xp,double yp, double zp, int n_order, double sx_interpolated[3][3] )
{
//int n_order is the order of interpolation
int i,j,k;
int nx_part, nz_part;
int min_limit=((n_order+1)/2)-1;    // n_order is the order of poly int.
int max_limit=(n_order+1)/2;
int n_point=n_order+1;
int Nx_ = sx.Nx();
int Ny_ = sx.Ny();
int Nz_ = sx.Nz();
int Nd_=sx.Nd();
double a_=sx.a();
double b_=sx.b();


double *scratch_part_,*x,*z;
double **rdata_for_interpolation;
//double u_interpolated[3];              //3 is for 3 comp  of velocities
double dy[3][3];
double delta_x=x_grid[1]-x_grid[0];
double delta_z=z_grid[1]-z_grid[0];


int locate(const Vector&, int, double);
double eval_at_yp(double *scratch_part_,double yp, double a_, double b_,int Ny_);
void poli_interpolation_2d(double [], double [], double **, int , int , double , double , double *,double *);

scratch_part_= new double [Ny_];
x= new double [n_point];
z= new double [n_point];



rdata_for_interpolation=new double *[n_point];
for(i=0;i<n_point;i++) rdata_for_interpolation[i]=new double [n_point];

nx_part=locate(x_grid,Nx_,xp);
nz_part=locate(z_grid,Nz_,zp);

//cout<<" nx_part= "<<nx_part<<" "<<nz_part<<endl;

  for (int i=0; i<Nd_/3; ++i){ 
    for (int f=0; f<Nd_/3; ++f){
          
       for (int nx=nx_part-min_limit; nx<=nx_part+max_limit; ++nx){
   j=nx-(nx_part-min_limit);
   x[j]=nx*delta_x;
   int nx1=nx;
   if(nx1<0)nx1=nx+(Nx_-1);if(nx1>(Nx_-1))nx1=nx-(Nx_-1);        //--check -1 wont be there.

   	for (int nz=nz_part-min_limit; nz<=nz_part+max_limit; ++nz){
k=nz-(nz_part-min_limit);
z[k]=nz*delta_z;
int nz1=nz;
if(nz1<0)nz1=nz+(Nz_-1);if(nz1>(Nz_-1))nz1=nz-(Nz_-1);


            for (int ny=0; ny<Ny_; ++ny)
			{
	scratch_part_[ny]= sx(nx1,ny,nz1,i3j(i,f));
			}    	// end loop for ny

   rdata_for_interpolation[j][k]=eval_at_yp(scratch_part_,yp,a_,b_,Ny_);

						   } 

						}
poli_interpolation_2d( x, z, rdata_for_interpolation,n_point,n_point,xp,zp,&sx_interpolated[i][f],&dy[i][f] );

			   }
                         }

delete [] x;
delete [] z;
delete [] scratch_part_;

for(i=0;i<n_point;i++) delete [] rdata_for_interpolation[i];
delete [] rdata_for_interpolation;

}


int locate(const Vector& xx, int n, double x)
{ 
int ju,jm,j1,j;
int ascnd;
j1=-1;
ju=n;
ascnd=(xx[n-1]>xx[0]);
while(ju-j1>1)
	{
	jm=(ju+j1) >>1;

	if(x>xx[jm] ==ascnd)
		j1=jm;
	else
		ju=jm;	


	}
j=j1;
return j;
}	



double eval_at_yp(double *data_part_, double x,double a_,double b_,int N_)
{
Real y = (2*x-a_-b_)/(b_-a_);
  Real y2 = 2*y;
  Real d=0.0;
  Real dd=0.0;
  for (int j=N_-1; j>0; --j) {
    Real sv=d;
    d = y2*d - dd + data_part_[j];
    dd=sv;
  }
  return (y*d - dd + data_part_[0]); // NR has 0.5*c[0], but that gives wrong results!
}


void poli_interpolation_2d(double x1a[], double x2a[], double **ya, int m, int n, double x1, double x2, double *y,double *dy)
{
void polint(double xa[], double ya[], int n, double x, double *y,double *dy);
int j;
double *ymtmp,y1,dy1;
ymtmp=new double [m];

for(j=0;j<m;j++)
{
polint(x2a,ya[j],n,x2,&ymtmp[j], &dy1);
}
polint(x1a,ymtmp,m,x1,&y1, &dy1);

*y=y1;
*dy=dy1;
delete [] ymtmp;
}

void polint(double xa[], double ya[], int n, double x, double *yy,double *dy)
{
  //***** take care---- here int n is not the poly_ order it is the no of input data point

int i,m,ns=0;
double den, dif, dift,ho,hp,w;
double *c,*d;
c= new double [n];
d= new double [n];

//for(i=0;i<n;i++)
//cout<<xa[i]<<endl;
//
// cin.get();
dif=fabs(x-xa[0]);

for(i=0;i<n;i++)
   {
if((dift=fabs(x-xa[i]))<dif)
	{
	ns=i;
	dif=dift;
	}	
	c[i]=ya[i];
	d[i]=ya[i];
     }
 *yy=ya[ns--];

  for(m=1;m<n;m++)
 {
  for(i=0;i<n-m;i++)
	{
	ho=xa[i]-x;
 hp=xa[i+m]-x;
 w=c[i+1]-d[i];
 if((den=ho-hp)==0.0)
 {cout<<"error in routin polint()"<< endl;
  cout<<i<<"   "<<m  <<"    "<<xa[i]<<"  "<<xa[i+m]<<endl;
 for(i=0;i<n;i++)
 cout<<xa[i]<<endl;
  exit(0);
  cin.get();
 }
 den=w/den;
 d[i]=hp*den;
 c[i]=ho*den;
 }
 *yy +=(*dy=(2*ns<(n-m)? c[ns+1] : d[ns--]));
 }
 delete [] d;
 delete [] c;
 }



void chebyshev_transform(const FlowField& force)
{


 for(int i=0;i<n_atom;++i){
               for(int ny=0;ny<force.Ny();++ny){
                    tmp3[i][ny] =2*cos(ny*acos((2*p_y[i]-(force.b()+force.a()))/(force.b()-force.a())));
                     }  
              tmp3[i][0] /= 2.0;
              tmp3[i][force.Ny()-1] /= 2.0;

        }

               

 return;              
}


void chebyshev_differentiation(const FlowField& force)
            
             {

               double a = force.a();
               double b = force.b();

            
              for(int i=0;i<n_atom;++i)
              {
               for(int ny=0;ny<force.Ny();++ny)
                     tmp_diff[i][ny] = 0.0;


                     tmp_diff[i][0] = 0.0;
                     tmp_diff[i][1] = 1.0;
                      
                     tmp_diff[i][2] = 2*2*cos(1*acos((2*p_y[i]-(b+a))/(b-a)));

                        for(int ny=3;ny<force.Ny();++ny)
                      tmp_diff[i][ny] = 2*ny*cos((ny-1)*acos((2*p_y[i]-(b+a))/(b-a))) + ny*tmp_diff[i][ny-2]/(ny-2);
                        
                  
     
                  for(int ny=0;ny<force.Ny();++ny)
                     tmp_diff2[i][ny] = 0.0;


                     tmp_diff2[i][0] = 0.0;
                     tmp_diff2[i][1] = 0.0;
                      
                     tmp_diff2[i][2] = 2*2* tmp_diff[i][1];

                        for(int ny=3;ny<force.Ny();++ny)
                      tmp_diff2[i][ny] = 2*ny*tmp_diff[i][ny-1] + ny*tmp_diff2[i][ny-2]/(ny-2);
                     }
                   return;
                 }
 


//************************************************************************************************//
//--------------------------Incorporating ReverseForce on fluid-----------------------------------//
//************************************************************************************************//

void reverse_quantity (int i, const FlowField& force, double **in, double ***out)

          {
              int nx=0;
              int ny=0;
              int nz=0;
              int s1=0;  int s3=0;
              int Nx = force.Nx();
              int Ny = force.Ny();
              int Nz = force.Nz();
              double delta_x=0;double delta_z=0;
              
              double x1=0;double x2=0;double z1=0;double z2=0;
              double a1=0;
              double totalreverseforce_p_x,totalreverseforce_grid_x,totalreverseforce_p_y,totalreverseforce_grid_y,
                    totalreverseforce_p_z,totalreverseforce_grid_z;
                         Real a = force.a();   
              Real b = force.b();   
              Real Lx = force.Lx();
              Real Lz = force.Lz(); 
              Vector x = periodicpoints(Nx, Lx);
              Vector y = chebypoints(Ny,a,b);
              Vector z = periodicpoints(Nz, Lz);
             
              
                   totalreverseforce_p_x=0;totalreverseforce_grid_x=0;totalreverseforce_p_y=0;totalreverseforce_grid_y=0;
                    totalreverseforce_p_z=0;totalreverseforce_grid_z=0;
                                             

                     for(nx=0;nx<Nx;++nx){
                        for(ny=0;ny<Ny;++ny){
                          for(nz=0;nz<Nz;++nz){
                               out[nx][ny][nz]=0;
                               out[nx][ny][nz]=0;
                               out[nx][ny][nz]=0;
                                        }
                                      }
                                   }

                     for(i=0;i<n_atom;++i)
                        {
                       for(nx=0;nx<Nx;++nx){
                         if(p_x[i]>x(nx) && p_x[i]<x(nx+1)){
                               
                                   x1 = p_x[i] - x(nx);
                                   x2 = x(nx+1) - p_x[i];
                                   delta_x = x(nx+1)-x(nx);
                                   s1=nx;
                                   
                                    }
              
                                  }

                    
                             for(nz=0;nz<Nz;++nz){
                              if(p_z[i]>z(nz) && p_z[i]<z(nz+1)){
                                   
                                     z1 = p_z[i] - z(nz);
                                     z2 = z(nz+1) - p_z[i];   
                                     delta_z = z(nz+1)-z(nz);
                                      s3=nz;
                                
                                        }
                                    } 
                                
                           for(ny=0;ny<Ny;++ny)
                             {
                             if (s1 == Nx-1  && s3 == Nz-1){
                                      nx=s1;nz=s3;     
                                              a1  = out[0][ny][0];
                             out[0][ny][0] = (x1*z1*in[i][ny])/(delta_x*delta_z);
                             out[0][ny][0] = out[0][ny][0]+a1;     
                                                
                                                    a1  = out[nx][ny][0];
                             out[nx][ny][0] = (x2*z1*in[i][ny])/(delta_x*delta_z);
                             out[nx][ny][0] = out[nx][ny][0]+a1;

                                                    a1  = out[0][ny][nz];
                             out[0][ny][nz] = (x1*z2*in[i][ny])/(delta_x*delta_z);
                             out[0][ny][nz] = out[0][ny][nz]+a1;

                                                  a1  = out[nx][ny][nz];
                             out[nx][ny][nz] = (x2*z2*in[i][ny])/(delta_x*delta_z);
                             out[nx][ny][nz] = out[nx][ny][nz]+a1;

                              
                                       }

                             else if (s1 == Nx-1){
                                           nz = s3;
                                      nx=s1;   
                                              a1  = out[0][ny][nz+1];
                             out[0][ny][nz+1] = (x1*z1*in[i][ny])/(delta_x*delta_z);
                             out[0][ny][nz+1] = out[0][ny][nz+1]+a1;     
                               

                                                    a1  = out[nx][ny][nz+1];
                             out[nx][ny][nz+1] = (x2*z1*in[i][ny])/(delta_x*delta_z);
                             out[nx][ny][nz+1] = out[nx][ny][nz+1]+a1;

                                                    a1  = out[0][ny][nz];
                             out[0][ny][nz] = (x1*z2*in[i][ny])/(delta_x*delta_z);
                             out[0][ny][nz] = out[0][ny][nz]+a1;


                                                  a1  = out[nx][ny][nz];
                             out[nx][ny][nz] = (x2*z2*in[i][ny])/(delta_x*delta_z);
                             out[nx][ny][nz] = out[nx][ny][nz]+a1;

                             

				 }
                                       else if (s3 == Nz-1){   
                                              nx=s1;nz=0;  
                                                
                                                    a1  = out[nx][ny][0];
                             out[nx][ny][0] = (x2*z1*in[i][ny])/(delta_x*delta_z);
                             out[nx][ny][0] = out[nx][ny][0]+a1;


                                                    a1  = out[nx+1][ny][0];
                             out[nx+1][ny][0] = (x1*z1*in[i][ny])/(delta_x*delta_z);
                             out[nx+1][ny][0] = out[nx+1][ny][0]+a1;

                                                
                                                    a1  = out[nx+1][ny][nz];
                             out[nx+1][ny][nz] = (x1*z2*in[i][ny])/(delta_x*delta_z);
                             out[nx+1][ny][nz] = out[nx+1][ny][nz]+a1;

                                                  a1  = out[nx][ny][nz];
                             out[nx][ny][nz] = (x2*z2*in[i][ny])/(delta_x*delta_z);
                             out[nx][ny][nz] = out[nx][ny][nz]+a1;

                                               
				}
                                        
                                         else {
                                                     nx=s1;nz=s3;

                                                     a1  = out[nx+1][ny][nz+1];
                             out[nx+1][ny][nz+1] = (x1*z1*in[i][ny])/(delta_x*delta_z);
                             out[nx+1][ny][nz+1] = out[nx+1][ny][nz+1]+a1;     
                                  
                                                    a1  = out[nx][ny][nz+1];
                             out[nx][ny][nz+1] = (x2*z1*in[i][ny])/(delta_x*delta_z);
                             out[nx][ny][nz+1] = out[nx][ny][nz+1]+a1;

                                                    a1  = out[nx+1][ny][nz];
                             out[nx+1][ny][nz] = (x1*z2*in[i][ny])/(delta_x*delta_z);
                             out[nx+1][ny][nz] = out[nx+1][ny][nz]+a1;

                                                
                                                  a1  = out[nx][ny][nz];
                             out[nx][ny][nz] = (x2*z2*in[i][ny])/(delta_x*delta_z);
                             out[nx][ny][nz] = out[nx][ny][nz]+a1;

                                              
                                     } 
                                 }
				} 
                              
                return;         
                      };


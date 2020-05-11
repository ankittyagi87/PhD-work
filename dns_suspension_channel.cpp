#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <complex>
#include <time.h>


#include "channelflow/flowfield.h"
#include "channelflow/vector.h"
#include "channelflow/chebyshev.h"
#include "channelflow/tausolver.h"
#include "channelflow/dns.h"
#include "channelflow/turbstats.h"
#include "channelflow/periodicfunc.h"
#include "channelflow/def"

#define Max_dim_position 4
using namespace std;
using namespace channelflow;


/*
#define diff_nmax 100
#define n_atom 4000
#define Z  10
#define Nof_lattice  Z * Z * Z 
#define map_size 13 * Nof_lattice 
#define dimension 3
#define Max_frame 2100
#define Max_y_grid 500
#define SEED 12345
typedef complex<double> Complex;
*/

//const int Max_dim_position=4;


double p_x[n_atom+2],p_y[n_atom+2],p_z[n_atom+2],//position
  p_vx[n_atom+2],p_vy[n_atom+2],p_vz[n_atom+2],p_v[n_atom+2],//velocity
  p_angvelx[n_atom+2],p_angvely[n_atom+2],p_angvelz[n_atom+2],//angular velocity
  p_x_previous[n_atom+2],p_y_previous[n_atom+2],p_z_previous[n_atom+2],//position at previous time step
  p_vx_previous[n_atom+2],p_vy_previous[n_atom+2],p_vz_previous[n_atom+2],//velocity at previous time step
  length_scale,time_scale,vel_scale,acc_scale,           //length scale , time scale, velocity scale         

  acc_x[n_atom+2],acc_y[n_atom+2],acc_z[n_atom+2],//acceleration
  angaccx[n_atom+2],angaccy[n_atom+2],angaccz[n_atom+2], //angular acceleration
  tc[n_atom],//collision related  
  et,en,ep,// restritution coefficient

  g_x,g_y,g_z,//gravitational const
  vx_max_air,vy_air,vz_air,// air velocity related
  den_gas,vis_gas,den_part,part_mass,sigma, sigma_sq, moment_inertia, radius, //air and particle properties

  b_lx,b_ly,b_lz, // Box length related
  fdrag_x[n_atom],fdrag_y[n_atom],fdrag_z[n_atom],//drag force
 
  diff_x[diff_nmax],diff_y[diff_nmax],diff_z[diff_nmax],diff_xy[diff_nmax],y_ch[diff_nmax],
  const_diff_x,const_diff_y,const_diff_z,const_diff_xy, //for giving const diffusivity
  tau_vp,tau_x,tau_y,tau_z, // diffusivity and integral time scale
  time_max,t_coll,t_noise,del_tc,del_tb, del_tb_actual,sampling_time,// run time related

   del_ts,t_sample,t_equb, // sampling time and equilibrium time
  del_x,del_y,del_z; //length of the grid
 
  // --------------variable used in function prperty()------------------------

 double  vx_inst[Max_y_grid], vy_inst[Max_y_grid], vz_inst[Max_y_grid],  //500 is for max y_grid
         vx_mean[Max_y_grid],vy_mean[Max_y_grid],vz_mean[Max_y_grid], Nofp_mean[Max_y_grid];
double  pv_dash[n_atom+2], pvx_dash[n_atom+2],pvy_dash[n_atom+2],pvz_dash[n_atom+2],pvx_dash_sq[Max_y_grid],
pvy_dash_sq[Max_y_grid],pvz_dash_sq[Max_y_grid],pvxy_dash[Max_y_grid];
 double sigma_xx[Max_y_grid],sigma_yy[Max_y_grid],sigma_zz[Max_y_grid], tau_xy[Max_y_grid];

double vx_vx_inst[Max_y_grid],vy_vy_inst[Max_y_grid],vz_vz_inst[Max_y_grid],vx_vy_inst[Max_y_grid],
  sigma_xx_1[Max_y_grid],sigma_yy_1[Max_y_grid],sigma_zz_1[Max_y_grid],sigma_xy_1[Max_y_grid];

 int   Nofp_inst[Max_y_grid];
//---------------------------------------------------------------------
  
 double del_vx,del_vy,del_vz,del_v,p_vx_max,p_vy_max,p_vz_max,p_v_max, //velocity distribution related
  del_r;// radial distribution function related


 double sum_B,tot_KE,run_time ,t_next, sum_KE, t_zero,B_zero;
 double eta;
 double  p_xi, p_yi, p_zi, p_vxi, p_vyi, p_vzi,coll_time;

 double  Total_Rep_x, Total_Rep_y, Total_Rep_z;
  int    count_Rep;

 // --------------variable used in angular function prperty()------------------------

 double  angvelx_inst[Max_y_grid], angvely_inst[Max_y_grid], angvelz_inst[Max_y_grid],  //500 is for max y_grid
         angvelx_mean[Max_y_grid],angvely_mean[Max_y_grid],angvelz_mean[Max_y_grid];
double  pangvel_dash[n_atom+2], pangvelx_dash[n_atom+2],pangvely_dash[n_atom+2],pangvelz_dash[n_atom+2],pangvelx_dash_sq[Max_y_grid],
pangvely_dash_sq[Max_y_grid],pangvelz_dash_sq[Max_y_grid],pangvelxy_dash[Max_y_grid];
 double angsigma_xx[Max_y_grid],angsigma_yy[Max_y_grid],angsigma_zz[Max_y_grid], angtau_xy[Max_y_grid];

double angvelx_angvelx_inst[Max_y_grid],angvely_angvely_inst[Max_y_grid],angvelz_angvelz_inst[Max_y_grid],angvelx_angvely_inst[Max_y_grid],
  angsigma_xx_1[Max_y_grid],angsigma_yy_1[Max_y_grid],angsigma_zz_1[Max_y_grid],angsigma_xy_1[Max_y_grid];

//---------------------------------------------------------------------
  
 double del_angvelx,del_angvely,del_angvelz,del_angvel,p_angvelx_max,p_angvely_max,p_angvelz_max,p_angvel_max; //angular velocity distribution related
  

 double  p_angvelxi, p_angvelyi, p_angvelzi;

//----------------------- parameters used in function final_property()----------------

 double  sum_vx[Max_y_grid],sum_vy[Max_y_grid],sum_vz[Max_y_grid],
         mean_vx[Max_y_grid], mean_vy[Max_y_grid],
         mean_vz[Max_y_grid];

//double  vx_dash[Max_frame][Max_y_grid], vy_dash[Max_frame][Max_y_grid],
//   vz_dash[Max_frame][Max_y_grid], Np_dash[Max_frame][Max_y_grid];

  double sum_sigma_xx[Max_y_grid],sum_sigma_yy[Max_y_grid],sum_sigma_zz[Max_y_grid],
          sum_tau_xy[Max_y_grid],sum_np_sq[Max_y_grid];

  double  mean_sigma_xx[Max_y_grid],mean_sigma_yy[Max_y_grid],mean_sigma_zz[Max_y_grid],
          mean_tau_xy[Max_y_grid],mean_np_sq[Max_y_grid];

  int sum_Nofp[Max_y_grid]; 

//----------------------- parameters used in function angular final_property()----------------

 double  sum_angvelx[Max_y_grid],sum_angvely[Max_y_grid],sum_angvelz[Max_y_grid],
         mean_angvelx[Max_y_grid], mean_angvely[Max_y_grid],
         mean_angvelz[Max_y_grid];

//double  angvelx_dash[Max_frame][Max_y_grid], angvely_dash[Max_frame][Max_y_grid],
//   angvelz_dash[Max_frame][Max_y_grid], Np_dash[Max_frame][Max_y_grid];

  double sum_angsigma_xx[Max_y_grid],sum_angsigma_yy[Max_y_grid],sum_angsigma_zz[Max_y_grid],
          sum_angtau_xy[Max_y_grid];

  double  mean_angsigma_xx[Max_y_grid],mean_angsigma_yy[Max_y_grid],mean_angsigma_zz[Max_y_grid],
          mean_angtau_xy[Max_y_grid];
 

//---------------------------parameter to calculate distribution function-----

  int npair[n_atom],//pair of colliding particle
      Nof_vx[Max_dim_position][2000],Nof_vy[Max_dim_position][2000],Nof_vz[Max_dim_position][2000],
    Nof_v[Max_dim_position][2000], Ngof_r[1000][1000], N_atom_count[1000],//velocity dist. related
      Ndel_x,Ndel_y,Ndel_z;// number of grid
  double fof_vx_sq[Max_dim_position][2000],fof_vy_sq[Max_dim_position][2000],fof_vz_sq[Max_dim_position][2000];
  int tot_atom_count[Max_dim_position]; //used in property()


//---------------------------parameters to calculate velocity distribution function for free flight------

int Nof_vx_free_flight[Max_dim_position][2000],Nof_vy_free_flight[Max_dim_position][2000],
  Nof_vz_free_flight[Max_dim_position][2000];
double fof_vx_free_flight_sq[Max_dim_position][2000],fof_vy_free_flight_sq[Max_dim_position][2000],
  fof_vz_free_flight_sq[Max_dim_position][2000];
int tot_atom_vel_dist_free_flight[Max_dim_position];

 void vel_dist_free_flight();
 void vel_distri_func_free_flight();


//---------------------------parameter to calculate angular distribution function-----

  int Nof_angvelx[Max_dim_position][2000],Nof_angvely[Max_dim_position][2000],Nof_angvelz[Max_dim_position][2000],
    Nof_angvel[Max_dim_position][2000]; //velocity dist. related
     
  double fof_angvelx_sq[Max_dim_position][2000],fof_angvely_sq[Max_dim_position][2000],fof_angvelz_sq[Max_dim_position][2000];
  


//---------------------------parameters to calculate angular velocity distribution function for free flight------

int Nof_angvelx_free_flight[Max_dim_position][2000],Nof_angvely_free_flight[Max_dim_position][2000],
  Nof_angvelz_free_flight[Max_dim_position][2000];
double fof_angvelx_free_flight_sq[Max_dim_position][2000],fof_angvely_free_flight_sq[Max_dim_position][2000],
  fof_angvelz_free_flight_sq[Max_dim_position][2000];
int tot_atom_angvel_dist_free_flight[Max_dim_position];

 void angular_vel_dist_free_flight();
 void angular_vel_distri_func_free_flight();

//---------------------------parameters to calculate acceleration  distribution function-----------------

  int  Nof_accx[Max_dim_position][2000], Nof_accy[Max_dim_position][2000], Nof_accz[Max_dim_position][2000],
    Nof_accxy[Max_dim_position][2000], Ndel_accx, Ndel_accy,Ndel_accz,Ndel_accxy;

  double fof_accx_sq[Max_dim_position][2000],fof_accy_sq[Max_dim_position][2000],fof_accz_sq[Max_dim_position][2000],
    fof_accxy_sq[Max_dim_position][2000];

  int tot_atom_acc_dist[Max_dim_position];

double del_accx,del_accy,del_accz,del_accxy,accx_max,accy_max,accz_max,accxy_max;

//---------------------------parameters to calculate acceleration distribution function from air------
int index_acc_dist_from_air;
void acc_dist_from_air(const FlowField& u,const Vector& x_grid, const Vector&  z_grid);
void acc_distri_func_from_air() ;

int Nof_accx_from_air[Max_dim_position][2000],Nof_accy_from_air[Max_dim_position][2000],
  Nof_accz_from_air[Max_dim_position][2000],Nof_accxy_from_air[Max_dim_position][2000];

double fof_accx_from_air_sq[Max_dim_position][2000],fof_accy_from_air_sq[Max_dim_position][2000];
double fof_accz_from_air_sq[Max_dim_position][2000],fof_accxy_from_air_sq[Max_dim_position][2000];
int Ndel_accx_from_air,Ndel_accy_from_air, Ndel_accz_from_air, Ndel_accxy_from_air;

 double accx_max_from_air,accy_max_from_air,accz_max_from_air,accxy_max_from_air;
 double del_accx_from_air,del_accy_from_air,del_accz_from_air,del_accxy_from_air;
//---------------------------------------------------------------------------------------------------

//---------------------------parameters to calculate angular acceleration  distribution function-----------------

  int  Nof_angaccx[Max_dim_position][2000], Nof_angaccy[Max_dim_position][2000], Nof_angaccz[Max_dim_position][2000],
    Nof_angaccxy[Max_dim_position][2000], Ndel_angaccx, Ndel_angaccy,Ndel_angaccz,Ndel_angaccxy;

  double fof_angaccx_sq[Max_dim_position][2000],fof_angaccy_sq[Max_dim_position][2000],fof_angaccz_sq[Max_dim_position][2000],
    fof_angaccxy_sq[Max_dim_position][2000];

  int tot_atom_angacc_dist[Max_dim_position];

double del_angaccx,del_angaccy,del_angaccz,del_angaccxy,angaccx_max,angaccy_max,angaccz_max,angaccxy_max;

//---------------------------parameters to calculate angular acceleration distribution function from air------
int index_angacc_dist_from_air;
void angular_acc_dist_from_air(const FlowField& omega,const Vector& x_grid, const Vector&  z_grid);
void angular_acc_distri_func_from_air();

int Nof_angaccx_from_air[Max_dim_position][2000],Nof_angaccy_from_air[Max_dim_position][2000],
  Nof_angaccz_from_air[Max_dim_position][2000],Nof_angaccxy_from_air[Max_dim_position][2000];

double fof_angaccx_from_air_sq[Max_dim_position][2000],fof_angaccy_from_air_sq[Max_dim_position][2000];
double fof_angaccz_from_air_sq[Max_dim_position][2000],fof_angaccxy_from_air_sq[Max_dim_position][2000];
int Ndel_angaccx_from_air,Ndel_angaccy_from_air, Ndel_angaccz_from_air, Ndel_angaccxy_from_air;

 double angaccx_max_from_air,angaccy_max_from_air,angaccz_max_from_air,angaccxy_max_from_air;
 double del_angaccx_from_air,del_angaccy_from_air,del_angaccz_from_air,del_angaccxy_from_air;
//----------------------------------------


//-------------------------parameters to calculate air angular velocity distribution function-----------------
void air_vel_fluc_dist(const FlowField& u,const Vector& x_grid,const Vector& z_grid, const Vector& y_grid);
  void air_vel_distri_func() ;

 int Nof_sample_airvel_dist[Max_dim_position];
 int Nof_air_velx[Max_dim_position][2000], Nof_air_vely[Max_dim_position][2000], Nof_air_velz[Max_dim_position][2000];
 int Ndel_air_velx,Ndel_air_vely,Ndel_air_velz;


 double max_air_velx,max_air_vely,max_air_velz;
 double del_air_velx_dist,del_air_vely_dist,del_air_velz_dist;
//------------------------------------------------------------------------------------------
  int Ndel_v,Ndel_vx,Ndel_vy,Ndel_vz, // numver of grid in velocity space
      n_unit,n_sample,sampls, k_write,atom_a,atom_b,
      sum_Nofp_inst[Max_frame][2000],
      N_rdels,
      Tot_Nof_coll,eqb_coll,max_coll,Tot_Nof_p_p_coll,Nof_p_w_coll,
      Nof_p_p_coll[Max_y_grid]={0}; // particle-particle and particle-wall collision



//-------------------------parameters to calculate air angular velocity distribution function-----------------
//void air_angular_vel_fluc_dist(const FlowField& omega,const Vector& x_grid,const Vector& z_grid, const Vector& y_grid);
  void air_angular_vel_distri_func() ;

 int Nof_sample_airangvel_dist[Max_dim_position];
 int Nof_air_angvelx[Max_dim_position][2000], Nof_air_angvely[Max_dim_position][2000], Nof_air_angvelz[Max_dim_position][2000];
 int Ndel_air_angvelx,Ndel_air_angvely,Ndel_air_angvelz;


 double max_air_angvelx,max_air_angvely,max_air_angvelz;
 double del_air_angvelx_dist,del_air_angvely_dist,del_air_angvelz_dist;

  int Ndel_angvel, Ndel_angvelx, Ndel_angvely,Ndel_angvelz;   //number of grid in angular velocity space

//--------------------------------------------------------------------------------------------------
  int time_step,time_step_write;
  int  N_diff;

 
  int  backward_index;
  int  Nof_backward_step;
  int  overlaping_a, overlaping_b;
  int  previous_atom_a, previous_atom_b;
  int  index_free_flight;
  int  data_saving_index;
  int  restart_for_part_index;

extern int index_air_vel_corre;


  double del_data_saving_time,starting_time,data_saving_time;

/*************Cheby diff***********************/
double **tmp3, **tmp_diff, **tmp_diff2;
double **frev_cheby_x, **frev_cheby_y, **frev_cheby_z; 
double **frev_cheby_diff2_x, **frev_cheby_diff2_y, **frev_cheby_diff2_z;

double symmetric_air_p_xx[n_atom],symmetric_air_p_xy[n_atom],symmetric_air_p_xz[n_atom],
                  symmetric_air_p_yx[n_atom],symmetric_air_p_yy[n_atom],symmetric_air_p_yz[n_atom],
                  symmetric_air_p_zx[n_atom],symmetric_air_p_zy[n_atom],symmetric_air_p_zz[n_atom];
//double ***freverse_x2_cheby_diff2, ***freverse_y2_cheby_diff2, ***freverse_z2_cheby_diff2;
/* ------parameters used for lattice cell---*/

  int  list[n_atom], head[Nof_lattice], map[map_size];
  int cell_of_part[n_atom+2];
  void    maping_neighbour();
  void    links ( ); 
  int     index_cell(int, int, int);
  double u_noise,v_noise,w_noise; // parameter used for perturbation of particle velocity

  const bool pause_here = true;

  bool overlap;
  bool event_noise_injection;
  bool event_collision;

  void   read_parameters();
  void   initial_position();
  void   initial_velocity(const FlowField& u,const Vector& x_grid, const Vector&  z_grid);
  void   initial_position_velocity_angularvelocity_restart(const string& );
   void   initial_angular_velocity(const FlowField& omega,const Vector& x_grid, const Vector& z_grid);
 
  void   drag(int , const FlowField& u,const FlowField& omega, const Vector& x_grid, const Vector&  z_grid);
  void reverseforce (int i, const FlowField& force);

 void chebyshev_transform(const FlowField& force);
 void chebyshev_differentiation(const FlowField& force);
 void reverseforce_spectral(int i, FlowField& force);

 

  double collision_time(int i,int j);

  void  backward_step();

  int n_order=5; // --5 th order lagrangian interpolation

  void interpolation(const FlowField& , const Vector& ,const Vector& ,double ,double , double , int, double[]);

  double Reynolds=2500.0;   //Reynolds number based on .5*(difference of the wall velocity)*channel half width
  //channelflow::Complex ****tmp4;
  //void Fourier_chebhyshev_transform(int i, const FlowField& u);
  

int main()
{

  time_t start,end;
  time(& start);

  FILE *fp11;
  FILE *fp12;
  FILE *fp110;
  FILE *fp111;

  int i, ii,ny;
  int ij;

 double wall_coll_time(int i,int j);
 void	advance_particle();
 void   post_coll_vel();
 void	update_coll_table(const FlowField& u,const FlowField& omega, const Vector& x_grid, const Vector&  z_grid);
 void   set_counter();
 void   distri_func();
 void   property();
 void   final_property();
 void   adding_noise();
 void   acc_property();   
 void   acc_distri_func();

void   angular_distri_func();
 void   angular_property();
 void   angular_final_property();                                                                         // angular
 void   angular_acc_property();   
 void   angular_acc_distri_func();

 void call_part_velo_corre();
 void call_part_acc_corre();
void call_fluid_vel_corre(const FlowField& ,const Vector& ,const Vector& );

void call_part_angular_velo_corre();
 void call_part_angular_acc_corre();                                                                      //angular
void call_fluid_angular_vel_corre(const FlowField& ,const Vector& ,const Vector& );
//--------------------parameter for dns------------------

 
   read_parameters(); // calling initial paremeter
  
   sigma=sqrt(sigma_sq);
   part_mass=3.1415926540*(sigma*sigma*sigma)*den_part/6.0;     //here sigma is nondimensional but not den_part
    moment_inertia = 0.40*(part_mass*radius*radius);
   cout<<"part_mass="<<part_mass<<endl;

   // initialization of the  parameters.
   run_time=starting_time;
   time_step=0;
   Tot_Nof_coll=0;      
   Tot_Nof_p_p_coll=0;
   Nof_p_w_coll=0;	
   t_noise=starting_time;
   Nof_backward_step=0;
   
  const int Nx=128;
  const int Ny=65;
  const int Nz=64;
 


  const Real Lx = b_lx; // 4.1448
  const Real Lz =b_lz; // 2.0724


  const Real a=0.0;
  const Real b=b_ly;
  //const Real bsub=-0.8;

  cout<<" L"<<Lx<<" "<<Lz<<" "<<b<<endl;
  const Real nu = 1.0/Reynolds;

  // const Real CFLmin = 0.70;
  // const Real CFLmax = 0.90;
  // const Real dtmin  = 0.0025;
  // const Real dtmax  = 0.05;
  //const Real dt =.0025;     //----p.s.g  

  const Real dt= del_tb;
  //const Real T0 = 0.0;   // integration start time
  // const Real T1 = 100;   // grow turbulence from perturbations: T0 < t <T1
  // const Real T2 = 300;   // take statistics: T1 <t <T2
  //const Real dT = 1.0;   // plotting and statistics interval
  const int n = 1;       //--p.s.g


  ofstream outfile ("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/time_1.bf");

  DNSFlags flags;
  flags.timestepping = CNAB2;//SBDF3;
  flags.nonlinearity = Rotational;
  flags.dealiasing   = DealiasXZ; 
  flags.constraint   = BulkVelocity;
  flags.Ubulk = -2*nu;
 
  const int kxmax=3;     // maximum Fourier mode for perturbations
  const int kzmax=3;
  //const Real perturbMag = 0.01;
  //const Real decay = 0.5; // spectral decay of perturbation profiles

  const char sp= ' ';
  const char nl= '\n';

  cout<<"delta_t= "<<dt<<endl;
  //cin.get();

//void ascii_save(int nx_part, int nz_part,const FlowField& u,  const Vector& x,const Vector& y,const Vector& z, const string& filebase) ; 

  cout << setprecision(6);

  Vector x = periodicpoints(Nx, Lx);
  Vector y = chebypoints(Ny,a,b);
  Vector z = periodicpoints(Nz, Lz);
  x.save("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/x");
  y.save("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/y");
  z.save("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/z");

 ChebyTransform trans(Ny);
 ChebyCoeff U(Ny,a,b,Physical);

  for (int ny=0; ny<Ny; ++ny) 
    U[ny] = (1.0 - square(abs(y[ny]-(b+a)/2.0)/((b-a)/2.0)));   
    U.save("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/Ubase");
    U.makeSpectral(trans);

  //FlowField u(Nx,Ny,Nz,3,Lx,Lz,a,b);
  // FlowField q(Nx,Ny,Nz,1,Lx,Lz,a,b);
  //u.addPerturbations(kxmax,kzmax,perturbMag,decay);

  //FlowField u("/global/scratch/chepart/couette_dns/run_v110/u_start199");
 // FlowField q("/global/scratch/chepart/couette_dns/run_v110/q_start199");

    /*if(restart_for_part_index==1)
    {
      FlowField u("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/u_start"+i2s(int(starting_time)));
      FlowField q("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/u_start"+i2s(int(starting_time)));
    }
    else{
      FlowField u("/global/scratch/chepart/couette_dns/run1/u_start159");
      FlowField q("/global/scratch/chepart/couette_dns/run1/q_start159");
      }*/

    string input_u_file,input_q_file;

    if(restart_for_part_index==1)
  { 
    input_u_file="/localscratch/cheantya/multipole/density_2000/volfrac_0.001/u_start_"+i2s(int(starting_time));
    input_q_file="/localscratch/cheantya/multipole/density_2000/volfrac_0.001/q_start_"+i2s(int(starting_time));
       }

    else
      {  input_u_file="/localscratch/cheantya/multipole/density_2000/volfrac_0.001/u_start160";
         input_q_file="/localscratch/cheantya/multipole/density_2000/volfrac_0.001/q_start160";
      }

    FlowField u(input_u_file);
    FlowField q(input_q_file);
    FlowField omega(Nx,Ny,Nz,3,Lx,Lz,a,b);
    FlowField e(Nx,Ny,Nz,1,Lx,Lz,a,b);
    FlowField force(Nx,Ny,Nz,3,Lx,Lz,a,b);
    FlowField tmp(Nx,Ny,Nz,12,Lx,Lz,a,b);
   // reverseforce(i,force);

 cout << "div(u)/L2Norm(u)  == " << divNorm(u)/L2Norm(u) << endl;  
  
  cout << "optimizing FFTW..." << flush;
  fftw_loadwisdom();
  u.optimizeFFTW();
  fftw_savewisdom();
  cout << "done" << endl;
  



  DNS dns(u, U, nu, dt, flags);
  dns.reset_Ubulk(2.0/3.0);               //--- very imp parameter... obtained by integration laminar profile

   

  ofstream modestream("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/fmodes.asc");
  ofstream dragstream("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/drags.asc");

  modestream << "% ";
  for (int kx=-kxmax; kx<=kxmax; ++kx) 
    for (int kz=0; kz<=kzmax; ++kz) 
      modestream << kx << ',' << kz << sp;
  modestream << nl;
  
  TurbStats stats(U, nu);

 
  initial_position();  //calling the function to get initial particle position

 fp110=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/part_initial_position.txt","w");
 for(i=0;i<n_atom;i++)
   fprintf(fp110,"%15d %30.20lf %30.20lf %30.20lf \n",i,p_x[i],p_y[i],p_z[i]);//--- writing the initial position in the file

  fclose(fp110);
    


       
      /* 
       if (pause_here) {
     cout << "hit return to continue..." << flush;
      char buff[10];
      cin.getline(buff,10);
      cin.get();
      	}
      */

  if(restart_for_part_index==1)
    {
initial_position_velocity_angularvelocity_restart("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/part_restart_"+i2s(int(starting_time)));
    }

  else
    {

       FlowField ucurl = u;
 
       ucurl += U;

      curl(ucurl,omega);
      
      u.makePhysical_xz();
      omega.makePhysical_xz();

      initial_velocity(u,x,z);    // calling the function to get initial velocity --- whent no restart
      initial_angular_velocity(omega,x,z);   
   
      u.makeSpectral_xz();      
      omega.makeSpectral_xz();
    
    }
  
      maping_neighbour();  //maping the neighbour of each lattice cell
     
      //------------BUILD THE TABLE OF COLLISION--------------------

       cout<<"     START to  BUILD THE TABLE OF COLLISION         "<<endl;

      
      FlowField ucurl = u;

       ucurl += U;

      curl(ucurl,omega); 

      u.makePhysical_xz();
      omega.makePhysical_xz();

      update_coll_table(u,omega,x,z); //here it is building collision table
      u.makeSpectral_xz();
      omega.makeSpectral_xz();

      cout<<"  END   "<<endl;

      fp11=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/time_ke_1.txt","w");
      fp12=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/repeated_coll_1.txt","w");

     
 if (eqb_coll >0){
   
   
   while(run_time< time_max){

     //----------------------------------------------------------
     //  IDENTIFY THE COLLIDING PAIR and  SHORTEST COLLISION time
     //----------------------------------------------------------

     //if(run_time>t_equb-2.0)del_tb=del_tb_actual;//  modifying del_tb ,2 is fixed arbitrarily

     //if(time_step>4000)del_tb=del_tb_actual;


  
  
     if((data_saving_time-run_time)<1e-6 && data_saving_index==1)
       {

       cout << "saving flowfields..." << endl;
      u.binarySave("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/u_start_"+i2s(int(run_time+0.5)));
      FlowField ucurl = u;
      curl(ucurl,omega);
      omega.binarySave("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/omega_start_"+i2s(int(run_time+0.5)));
      q.binarySave("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/q_start_"+i2s(int(run_time+0.5)));
      cout << "done" << endl;
      

   string filename("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/part_restart_"+i2s(int(run_time+0.5)));
  filename += ".asc";
  ofstream os1(filename.c_str(),ios::out);
  os1 << scientific << setprecision(10);
  for(ij=0;ij<n_atom;ij++)
    os1<<ij<<" "<<p_x[ij]<<"  "<<p_y[ij]<<"  "<<p_z[ij]<<"  "<<p_vx[ij]<<"  "<<p_vy[ij]<<"  "<<p_vz[ij]<<"   "<<p_angvelx[ij]<<"  "<<p_angvely[ij]<<"  "<<p_angvelz[ij]<<endl;

  os1.close();

 data_saving_time=data_saving_time+ del_data_saving_time;
    }



     event_noise_injection= false;
     event_collision=false;

 del_tc=1.e8;
 for( ii=0;ii<n_atom;ii++)
  {
   if(tc[ii]< del_tc)
                  {
	       
		del_tc=tc[ii];
	        atom_a=ii;
		  }
   }
               
 atom_b=npair[atom_a];

 
  
  // Check for negative collision 

   for( ii=0;ii<n_atom;ii++)
	       {
	   
	       if(tc[ii]< 0.0)
		 {
		cout<<" COLLITION time NEGATIVE  "<<endl;
		cout<<"    "<<endl;
		cout<<ii<<" "<<npair[ii]<<" "<<tc[ii]<<endl;
	        exit(1);
	                     }
			     }

	     t_coll=run_time + del_tc;
	     
	     

	   if( t_sample <t_noise && t_sample < t_coll)
	     
	     {   
	                           // looop for sampling the data
	      t_next= t_sample-run_time;
	      cout<<"tnext="<<t_next;
	      advance_particle(); //calling function to move particle 
            
	     FlowField u_dummy=u;
    
            FlowField ucurl = u_dummy;   

            ucurl += U;            
  
           curl(ucurl,omega);	 

           u_dummy.makePhysical_xz();                 
          omega.makePhysical_xz();
                    
	  update_coll_table(u_dummy,omega,x,z);              // updating callision table
	  u_dummy.makeSpectral_xz();
          omega.makeSpectral_xz();

	  if(n_sample==0){
	   set_counter();
	  	  }

	  cout<<" %%%%%%%%% Now sampling the data %%%%%%%%%%%%%%%%%%"<<endl;

	  

        property();   // calling property() for calculation of diff. property
        angular_property();

	run_time=run_time+t_next;
	time_step= time_step+1;
	t_sample=t_sample+del_ts;   //t_sample is the  next sampling time


       stats.addData(u,tmp);

      stats.msave("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/uu");
      stats.msave("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/uustar", true);
      
      Real ustar = stats.ustar();
      Vector yp = stats.yplus();
      yp.save("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/yp");

      ChebyCoeff Umean = stats.U();
      Umean.makeSpectral(trans);
      ChebyCoeff Umeany = diff(Umean);
      ChebyCoeff Umeanyy = diff2(Umean);

      Umean.makePhysical(trans);
      Umeany.makePhysical(trans);
      Umeanyy.makePhysical(trans);

      Umean.save("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/Umean");
      Umeany.save("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/Umeany");
      Umeanyy.save("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/Umeanyy");      

      Umean /= ustar;
      Umean.save("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/Uplus");

      ChebyCoeff ubase = stats.ubase();
      ubase.save("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/ubase");

      ChebyCoeff uv = stats.uv();
      uv.save("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/uv");
      
      uv.makeSpectral(trans);
      ChebyCoeff uvy = diff(uv);
      uv.makePhysical(trans);
      uvy.makePhysical(trans);
      uvy.save("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/uvy");

      save(ustar, "/localscratch/cheantya/multipole/density_2000/volfrac_0.001/ustar");
      save(nu, "/localscratch/cheantya/multipole/density_2000/volfrac_0.001/nu");
      FlowField e = energy(u,U);
      
      e.makeSpectral();
      e.asciiSave("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/energy");

      
      omega.setState(Spectral,Spectral);
      FlowField rot_energy = energy(omega);
      rot_energy.asciiSave("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/rot_energy");
      
     
	}
    

      else
	{ 
  	  if(t_coll > t_noise)
	    {                              // loop for integrating the air field
                        
	      t_next= t_noise-run_time;  //time difference between present time and next event of noise injection

	      advance_particle();

	      run_time=run_time+t_next;      
	      time_step= time_step+1;
	      
	      event_noise_injection=true;
    
	 cout<<"********integrating air field ******   "<< t_next<<" "<<time_step<<" "<<endl;
	          
            // Fourier_chebhyshev_transform(i,u); 
                frev_cheby_x = new double *[n_atom];frev_cheby_y = new double *[n_atom];frev_cheby_z = new double *[n_atom];
            frev_cheby_diff2_x = new double *[n_atom];frev_cheby_diff2_y = new double *[n_atom];frev_cheby_diff2_z = new double *[n_atom];
            for(i=0;i<n_atom;++i)
             {
            frev_cheby_x[i] = new double [Ny];frev_cheby_y[i] = new double [Ny];frev_cheby_z[i] = new double [Ny];
            frev_cheby_diff2_x[i] = new double [Ny];frev_cheby_diff2_y[i] = new double [Ny];frev_cheby_diff2_z[i] = new double [Ny];
             }
                
             tmp3 = new double *[n_atom];
                for(i=0;i<n_atom;++i)
                  tmp3[i] = new double [force.Ny()];

               tmp_diff2 = new double *[n_atom];
                  for(i=0;i<n_atom;++i)
                     tmp_diff2[i] = new double [force.Ny()];

 
               tmp_diff = new double *[n_atom];
                  for(i=0;i<n_atom;++i)
                     tmp_diff[i] = new double [force.Ny()];
            
              chebyshev_transform(force);
              chebyshev_differentiation(force);
              reverseforce_spectral(i,force);
              
               dns.advance(u, q, n);   //---calling the function to advance flowfield
             
               for(i=0;i<n_atom;++i){
	          delete [] tmp3[i];delete [] tmp_diff[i];delete [] tmp_diff2[i];
                    }
                  delete [] tmp3;delete [] tmp_diff;delete [] tmp_diff2;

                for(i=0;i<n_atom;++i){  
                  delete [] frev_cheby_x[i]; delete [] frev_cheby_y[i];  delete [] frev_cheby_z[i];
                  delete [] frev_cheby_diff2_x[i]; delete [] frev_cheby_diff2_y[i];  delete [] frev_cheby_diff2_z[i];
                   }
	          delete [] frev_cheby_x; delete [] frev_cheby_y;  delete [] frev_cheby_z;
                  delete [] frev_cheby_diff2_x; delete [] frev_cheby_diff2_y;  delete [] frev_cheby_diff2_z;

                  

	     // cout << "   CFL ==  " << dns_cfl << " " <<dns_l2norm2<<" "<<dns_divnorm2<<endl;
    
	      outfile << run_time << "  " << dns.CFL() << "  " << L2Norm2(u)  <<" " << divNorm2(u)<<endl; 

           t_noise=t_noise+del_tb;  //---- t_noise is the next airfield integration time
          
	   FlowField u_dummy=u;
	   FlowField u1_dummy=u;

           FlowField ucurl = u_dummy;

           ucurl += U;           

           curl(ucurl,omega);

           u_dummy.makePhysical_xz();
	   u1_dummy.makePhysical();

     	   
           omega.makePhysical_xz();

	   update_coll_table(u_dummy,omega,x,z);   // ----calling function to update colision table 
             
	   index_free_flight=1;
	   if( index_free_flight==1 && run_time>t_equb+2.0*del_ts)
	     {

	   acc_property();
           angular_acc_property();

     if(index_acc_dist_from_air==1){


   acc_dist_from_air(u_dummy,x,z);// to calculate distribution function  of u'/tau_v and v_p at mean free time
   angular_acc_dist_from_air(omega,x,z); // to calculate angular distribution function 
  
      air_vel_fluc_dist(u1_dummy,x, z, y);}
      //air_angular_vel_fluc_dist(omega1_dummy,x,z,y);

	    vel_dist_free_flight();
            angular_vel_dist_free_flight();
 
	    call_part_acc_corre();
            call_part_angular_acc_corre();
if(index_air_vel_corre==1){call_fluid_vel_corre(u_dummy,x,z); call_fluid_angular_vel_corre(omega,x,z);}
	    call_part_velo_corre();
            call_part_angular_velo_corre();
            
	     }
	    
      u_dummy.makeSpectral_xz();
      u1_dummy.makeSpectral();
      omega.makeSpectral_xz();
	   index_free_flight=0;

	    }
 
	  else
	    {               	      // loop for allowing the particle to collide
          
	     t_next=del_tc ;
     
	  //    cout<< " in main "<<atom_b<<"  "<<p_vx[atom_b]<<endl;

           advance_particle();        // Move all the particles during the next occurence
           run_time=run_time+t_next;
	   time_step= time_step+1;

	   event_collision=true;

	//   cout<<"-------------------------"<<endl;
   //    cout<<endl;

	cout<<time_step<<" "<<run_time<<" "<<tot_KE<<" "<< atom_a<<" "<<atom_b<<" "<<Tot_Nof_coll<<" "<<Tot_Nof_p_p_coll<<" "<< Nof_p_w_coll<<" "<<t_next<<endl;
/*
	   cout<<endl;
           cout<<"------------Linear Velocities--------------------"<<endl;
	   cout<<p_vx[atom_a]<<" "<<p_vy[atom_a]<<" " <<p_vz[atom_a]<<endl;
	   cout<<p_vx[atom_b]<<" "<<p_vy[atom_b]<<" " <<p_vz[atom_b]<<endl;
           cout<<"------------Angular Velocities--------------------"<<endl;
           cout<<p_angvelx[atom_a]<<"   "<<p_angvely[atom_a]<<"   "<<p_angvelz[atom_a]<<endl;
           cout<<p_angvelx[atom_b]<<"   "<<p_angvely[atom_b]<<"   "<<p_angvelz[atom_b]<<endl;

	   cout<<endl;
	   cout<< "-----------------------"<<endl;
*/
       if((atom_a == previous_atom_a && atom_b==previous_atom_b)||(atom_a== previous_atom_b && atom_b==previous_atom_a))
	 fprintf(fp12,"%15d %8d %8d %8d %8d\n",time_step,atom_a,atom_b,previous_atom_a ,previous_atom_b);


           post_coll_vel();     //Determine change in velocity of the colliding pair
	   
	   FlowField u_dummy=u; 
           FlowField ucurl = u_dummy;
          
           ucurl += U;

           curl(ucurl,omega);
         
           u_dummy.makePhysical_xz();          
           omega.makePhysical_xz();

	   update_coll_table(u_dummy,omega,x,z);  // Update the collision table

	   u_dummy.makeSpectral_xz();
           omega.makeSpectral_xz(); 
 

	   previous_atom_a=atom_a; // storing the pair to recall after next updating
	   previous_atom_b=atom_b;
         
         if(atom_a>=n_atom ||atom_b>=n_atom)
	  Nof_p_w_coll=Nof_p_w_coll+1;
	 else
	   {
	  Tot_Nof_p_p_coll= Tot_Nof_p_p_coll+1; 

	   ny=(int)((((p_y[atom_a]+p_y[atom_b])/2.0-(sigma/2.0))/del_y));
	   Nof_p_p_coll[ny]= Nof_p_p_coll[ny]+1;

	   
	   }
 		

	    Tot_Nof_coll=Tot_Nof_coll+1;
        
	    }
	}

 
     
	 //    cout<<run_time<<"    "<<tot_KE<<"  "<<Nof_backward_step<<endl;

	     if(fmod((double)time_step,(double)time_step_write)==0){
      
  fprintf(fp11,"%15d %30.20lf %30.20lf %8d %8d %8d %8d\n",time_step,run_time,tot_KE,Tot_Nof_coll,Tot_Nof_p_p_coll,Nof_p_w_coll,Nof_backward_step);
   }
      
                 
   }
 }

 if(run_time>t_equb)
   {

 final_property(); // average final property

 angular_final_property();         //average angular final property
 
 distri_func();   //velosity and radial distribution function  

 angular_distri_func();      //angular velocity distribution

 acc_distri_func();

 angular_acc_distri_func();     //angular acceleration distribution
 

 if(index_acc_dist_from_air==1){
acc_distri_func_from_air() ;
air_vel_distri_func();

angular_acc_distri_func_from_air();
air_angular_vel_distri_func();
 }

 vel_distri_func_free_flight();
 angular_vel_distri_func_free_flight(); 

 fp111=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/part_final_position.txt","w");
 for(i=0;i<n_atom;i++)
   fprintf(fp111,"%15d %30.20lf %30.20lf %30.20lf \n",i,p_x[i],p_y[i],p_z[i]);
  fclose(fp111);

   }

  /*
   u.makePhysical();
    u.saveSlice(2,0,45,"/localscratch/cheantya/multipole/density_2000/volfrac_0.001/u_at_xy");
    u.saveSlice(2,1,45,"/localscratch/cheantya/multipole/density_2000/volfrac_0.001/v_at_xy");
    u.saveSlice(2,2,45,"/localscratch/cheantya/multipole/density_2000/volfrac_0.001/w_at_xy");

    u.saveSlice(1,0,27,"/localscratch/cheantya/multipole/density_2000/volfrac_0.001/u_at_xz_centre");
    u.saveSlice(1,1,27,"/localscratch/cheantya/multipole/density_2000/volfrac_0.001/v_at_xz_centre");
    u.saveSlice(1,2,27,"/localscratch/cheantya/multipole/density_2000/volfrac_0.001/w_at_xz_centre");
    
    u.saveSlice(1,0,1,"/localscratch/cheantya/multipole/density_2000/volfrac_0.001/u_at_xz_nearwall");
    u.saveSlice(1,1,1,"/localscratch/cheantya/multipole/density_2000/volfrac_0.001/v_at_xz_nearwall");
    u.saveSlice(1,2,1,"/localscratch/cheantya/multipole/density_2000/volfrac_0.001/w_at_xz_nearwall");

    u.makeSpectral();

    curl(u, omega);
    omega.saveSlice(1,0,27,"/localscratch/cheantya/multipole/density_2000/volfrac_0.001/omegax_xz_centre");
    omega.saveSlice(1,1,27,"/localscratch/cheantya/multipole/density_2000/volfrac_0.001/omegay_xz_centre");
    omega.saveSlice(1,2,27,"/localscratch/cheantya/multipole/density_2000/volfrac_0.001/omegaz_xz_centre");

    omega.saveSlice(1,0,1,"/localscratch/cheantya/multipole/density_2000/volfrac_0.001/omegax_xz_nearwall");
    omega.saveSlice(1,1,1,"/localscratch/cheantya/multipole/density_2000/volfrac_0.001/omegay_xz_nearwall");
    omega.saveSlice(1,2,1,"/localscratch/cheantya/multipole/density_2000/volfrac_0.001/omegaz_xz_nearwall");

    omega.makePhysical();
          omega.asciiSave("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/omega_ascii");
  */


 cout<<"******* Programme has ended successfully **************"<<endl; 
 time( &end);
 double diff_=difftime(end,start);
 cout<< " "<<endl;
 cout<<" average Rep_X,Rep_y and Rep_z are   "<< Total_Rep_x/count_Rep<<"   "<< Total_Rep_y/count_Rep<<"  "<< Total_Rep_z/count_Rep;
 cout<<  "  "<<endl;
 cout<<" time taken for the simulation = "<<diff_<<endl;  
 

 return(0);
}                             //end of the main program

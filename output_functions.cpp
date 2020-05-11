#include "channelflow/flowfield.h"
#include "channelflow/vector.h"
#include "channelflow/def"

#define Max_dim_position 4
using namespace channelflow;

//extern  int Max_dim_position;

extern int Ndel_v,Ndel_vx,Ndel_vy,Ndel_vz;
extern int Ndel_angvel,Ndel_angvelx,Ndel_angvely,Ndel_angvelz;                                //angular
extern int Ndel_x,Ndel_y,Ndel_z;
extern int Tot_Nof_coll,Tot_Nof_p_p_coll,Nof_p_w_coll,sampls, n_sample;

extern int   Nof_vx[][2000],Nof_vy[][2000],Nof_vz[][2000],Nof_v[][2000], Ngof_r[][1000], 
N_atom_count[], Nof_p_p_coll[];

extern int Nof_angvelx[][2000], Nof_angvely[][2000],Nof_angvelz[][2000],Nof_angvel[][2000];           //angular

double Nof_v_inst[Max_dim_position][2000],Nof_vx_inst[Max_dim_position][2000], 
  Nof_vy_inst[Max_dim_position][2000], Nof_vz_inst[Max_dim_position][2000];

double Nof_angvel_inst[Max_dim_position][2000],Nof_angvelx_inst[Max_dim_position][2000],               //angular
  Nof_angvely_inst[Max_dim_position][2000], Nof_angvelz_inst[Max_dim_position][2000];

extern double   fof_vx_sq[][2000],fof_vy_sq[][2000],fof_vz_sq[][2000];
extern double   fof_angvelx_sq[][2000],fof_angvely_sq[][2000],fof_angvelz_sq[][2000];                 //angular

extern double sum_KE,sigma,b_lx,b_ly,b_lz;
extern double p_xi, p_yi, p_zi,p_x[n_atom+2],p_y[n_atom+2],p_z[n_atom+2];
extern double p_vx[n_atom+2],p_vy[n_atom+2],p_vz[n_atom+2],p_v[n_atom+2];          //particle velocities
extern double  p_angvelx[n_atom+2],p_angvely[n_atom+2],p_angvelz[n_atom+2],p_angvel[n_atom+2];// particle angular velocities
extern double acc_x[n_atom+2],acc_y[n_atom+2],acc_z[n_atom+2];                       //acceleration
extern double angaccx[n_atom+2],angaccy[n_atom+2],angaccz[n_atom+2];                //particle angular velocities

extern double del_r;
extern double del_vx,del_vy,del_vz,del_v,p_vx_max,p_vy_max,p_vz_max,p_v_max;
extern double del_angvelx,del_angvely,del_angvelz,del_angvel,p_angvelx_max,p_angvely_max,p_angvelz_max,p_angvel_max;                  //angular
extern double del_x,del_y,del_z;
extern double del_accx,del_accy,del_accz,del_accxy,accx_max,accy_max,accz_max,accxy_max;
extern double del_angaccx,del_angaccy,del_angaccz,del_angaccxy,angaccx_max,angaccy_max,angaccz_max,angaccxy_max;                 //angular
 
extern double  y_min_for_distribution_1,y_max_for_distribution_1,
  y_min_for_distribution_2,y_max_for_distribution_2,
  y_min_for_distribution_3,y_max_for_distribution_3;

// --------------variable used in function prperty()------------------------

extern double  vx_inst[Max_y_grid], vy_inst[Max_y_grid], vz_inst[Max_y_grid],  //500 is for max y_grid
         vx_mean[Max_y_grid],vy_mean[Max_y_grid],vz_mean[Max_y_grid], Nofp_mean[Max_y_grid];
double	Nofp_inst_sq[Max_y_grid];

extern double  angvelx_inst[Max_y_grid], angvely_inst[Max_y_grid], angvelz_inst[Max_y_grid],  //500 is for max y_grid
         angvelx_mean[Max_y_grid],angvely_mean[Max_y_grid],angvelz_mean[Max_y_grid];// Nofp_mean[Max_y_grid];                //angular
//double	Nofp_inst_sq[Max_y_grid];

 double symmetric_p_xx_mean[Max_y_grid],
    symmetric_p_xy_mean[Max_y_grid],
    symmetric_p_xz_mean[Max_y_grid],
    symmetric_p_yy_mean[Max_y_grid],
    symmetric_p_yz_mean[Max_y_grid],
    symmetric_p_zz_mean[Max_y_grid],

    symmetric_particle_corr_xx_mean[Max_y_grid],
    symmetric_particle_corr_xy_mean[Max_y_grid],
    symmetric_particle_corr_xz_mean[Max_y_grid],
    symmetric_particle_corr_yy_mean[Max_y_grid],
    symmetric_particle_corr_yz_mean[Max_y_grid],
    symmetric_particle_corr_zz_mean[Max_y_grid],

    symmetric_p_xx_inst[Max_y_grid],
    symmetric_p_xy_inst[Max_y_grid],
    symmetric_p_xz_inst[Max_y_grid],
    symmetric_p_yy_inst[Max_y_grid],
    symmetric_p_yz_inst[Max_y_grid],
    symmetric_p_zz_inst[Max_y_grid];

extern double symmetric_air_p_xx[n_atom],symmetric_air_p_xy[n_atom],symmetric_air_p_xz[n_atom],
                  symmetric_air_p_yx[n_atom],symmetric_air_p_yy[n_atom],symmetric_air_p_yz[n_atom],
                  symmetric_air_p_zx[n_atom],symmetric_air_p_zy[n_atom],symmetric_air_p_zz[n_atom];

   double angvel_corr_x_mean[Max_y_grid];
   double angvel_corr_y_mean[Max_y_grid];
   double angvel_corr_z_mean[Max_y_grid];

extern double pv_dash[n_atom+2], pvx_dash[n_atom+2],pvy_dash[n_atom+2],pvz_dash[n_atom+2],pvx_dash_sq[Max_y_grid],
pvy_dash_sq[Max_y_grid],pvz_dash_sq[Max_y_grid],pvxy_dash[Max_y_grid];
extern double sigma_xx[Max_y_grid],sigma_yy[Max_y_grid],sigma_zz[Max_y_grid], tau_xy[Max_y_grid];

extern double pangvel_dash[n_atom+2], pangvelx_dash[n_atom+2],pangvely_dash[n_atom+2],pangvelz_dash[n_atom+2],pangvelx_dash_sq[Max_y_grid],
pangvely_dash_sq[Max_y_grid],pangvelz_dash_sq[Max_y_grid],pangvelxy_dash[Max_y_grid];                                                           //angular
extern double angsigma_xx[Max_y_grid],angsigma_yy[Max_y_grid],angsigma_zz[Max_y_grid], angtau_xy[Max_y_grid];

extern double vx_vx_inst[Max_y_grid],vy_vy_inst[Max_y_grid],vz_vz_inst[Max_y_grid],vx_vy_inst[Max_y_grid],
sigma_xx_1[Max_y_grid],sigma_yy_1[Max_y_grid],sigma_zz_1[Max_y_grid],sigma_xy_1[Max_y_grid];

extern double angvelx_angvelx_inst[Max_y_grid],angvely_angvely_inst[Max_y_grid],angvelz_angvelz_inst[Max_y_grid],angvelx_angvely_inst[Max_y_grid],   //angular
angsigma_xx_1[Max_y_grid],angsigma_yy_1[Max_y_grid],angsigma_zz_1[Max_y_grid],angsigma_xy_1[Max_y_grid];

double avg_stress_x[Max_y_grid],avg_stress_y[Max_y_grid],avg_stress_z[Max_y_grid],avg_stress_xy[Max_y_grid];
double var_sigma_xx[Max_y_grid], var_sigma_yy[Max_y_grid],var_sigma_zz[Max_y_grid],var_tau_xy[Max_y_grid];

double avg_angstress_x[Max_y_grid],avg_angstress_y[Max_y_grid],avg_angstress_z[Max_y_grid],avg_angstress_xy[Max_y_grid];
double var_angsigma_xx[Max_y_grid], var_angsigma_yy[Max_y_grid],var_angsigma_zz[Max_y_grid],var_angtau_xy[Max_y_grid];                          //angular

double  pvx_dash_4[Max_y_grid], pvy_dash_4[Max_y_grid], pvz_dash_4[Max_y_grid], pvxy_dash_4[Max_y_grid];
double  avg_pvx_dash_4[Max_y_grid], avg_pvy_dash_4[Max_y_grid], avg_pvz_dash_4[Max_y_grid], avg_pvxy_dash_4[Max_y_grid];

double  pangvelx_dash_4[Max_y_grid], pangvely_dash_4[Max_y_grid], pangvelz_dash_4[Max_y_grid], pangvelxy_dash_4[Max_y_grid];
double  avg_pangvelx_dash_4[Max_y_grid], avg_pangvely_dash_4[Max_y_grid], avg_pangvelz_dash_4[Max_y_grid], avg_pangvelxy_dash_4[Max_y_grid];  //angular

double  pvx_dash_3[Max_y_grid], pvy_dash_3[Max_y_grid], pvz_dash_3[Max_y_grid];
double  avg_pvx_dash_3[Max_y_grid], avg_pvy_dash_3[Max_y_grid], avg_pvz_dash_3[Max_y_grid];

double  pangvelx_dash_3[Max_y_grid], pangvely_dash_3[Max_y_grid], pangvelz_dash_3[Max_y_grid];
double  avg_pangvelx_dash_3[Max_y_grid], avg_pangvely_dash_3[Max_y_grid], avg_pangvelz_dash_3[Max_y_grid];                        //angular

double accx_dash_2[Max_y_grid],accy_dash_2[Max_y_grid],accz_dash_2[Max_y_grid];
double accx_dash_3[Max_y_grid],accy_dash_3[Max_y_grid],accz_dash_3[Max_y_grid];
double accx_dash_4[Max_y_grid],accy_dash_4[Max_y_grid],accz_dash_4[Max_y_grid];

double angaccx_dash_2[Max_y_grid],angaccy_dash_2[Max_y_grid],angaccz_dash_2[Max_y_grid];
double angaccx_dash_3[Max_y_grid],angaccy_dash_3[Max_y_grid],angaccz_dash_3[Max_y_grid];                                           //angular
double angaccx_dash_4[Max_y_grid],angaccy_dash_4[Max_y_grid],angaccz_dash_4[Max_y_grid];

double avg_accx_dash_2[Max_y_grid],avg_accy_dash_2[Max_y_grid],avg_accz_dash_2[Max_y_grid];
double avg_accx_dash_3[Max_y_grid],avg_accy_dash_3[Max_y_grid],avg_accz_dash_3[Max_y_grid];
double avg_accx_dash_4[Max_y_grid],avg_accy_dash_4[Max_y_grid],avg_accz_dash_4[Max_y_grid];

double avg_angaccx_dash_2[Max_y_grid],avg_angaccy_dash_2[Max_y_grid],avg_angaccz_dash_2[Max_y_grid];
double avg_angaccx_dash_3[Max_y_grid],avg_angaccy_dash_3[Max_y_grid],avg_angaccz_dash_3[Max_y_grid];                          //angular
double avg_angaccx_dash_4[Max_y_grid],avg_angaccy_dash_4[Max_y_grid],avg_angaccz_dash_4[Max_y_grid];

double  mean_accx_sq[Max_y_grid], mean_accy_sq[Max_y_grid], mean_accz_sq[Max_y_grid];
double accx_mean[Max_y_grid],accy_mean[Max_y_grid],accz_mean[Max_y_grid];
double mean_accx_dash_2[Max_y_grid],mean_accy_dash_2[Max_y_grid],mean_accz_dash_2[Max_y_grid];

double  mean_angaccx_sq[Max_y_grid], mean_angaccy_sq[Max_y_grid], mean_angaccz_sq[Max_y_grid];
double angaccx_mean[Max_y_grid],angaccy_mean[Max_y_grid],angaccz_mean[Max_y_grid];                                           //angular
double mean_angaccx_dash_2[Max_y_grid],mean_angaccy_dash_2[Max_y_grid],mean_angaccz_dash_2[Max_y_grid];


double accx_dash_from_air_2[Max_y_grid],accy_dash_from_air_2[Max_y_grid],accz_dash_from_air_2[Max_y_grid];
double accx_dash_from_air_3[Max_y_grid],accy_dash_from_air_3[Max_y_grid],accz_dash_from_air_3[Max_y_grid];
double accx_dash_from_air_4[Max_y_grid],accy_dash_from_air_4[Max_y_grid],accz_dash_from_air_4[Max_y_grid];

double angaccx_dash_from_air_2[Max_y_grid],angaccy_dash_from_air_2[Max_y_grid],angaccz_dash_from_air_2[Max_y_grid];
double angaccx_dash_from_air_3[Max_y_grid],angaccy_dash_from_air_3[Max_y_grid],angaccz_dash_from_air_3[Max_y_grid];             //angular
double angaccx_dash_from_air_4[Max_y_grid],angaccy_dash_from_air_4[Max_y_grid],angaccz_dash_from_air_4[Max_y_grid];

double avg_accx_dash_from_air_2[Max_y_grid],avg_accy_dash_from_air_2[Max_y_grid],avg_accz_dash_from_air_2[Max_y_grid];
double avg_accx_dash_from_air_3[Max_y_grid],avg_accy_dash_from_air_3[Max_y_grid],avg_accz_dash_from_air_3[Max_y_grid];
double avg_accx_dash_from_air_4[Max_y_grid],avg_accy_dash_from_air_4[Max_y_grid],avg_accz_dash_from_air_4[Max_y_grid];

double avg_angaccx_dash_from_air_2[Max_y_grid],avg_angaccy_dash_from_air_2[Max_y_grid],avg_angaccz_dash_from_air_2[Max_y_grid];
double avg_angaccx_dash_from_air_3[Max_y_grid],avg_angaccy_dash_from_air_3[Max_y_grid],avg_angaccz_dash_from_air_3[Max_y_grid];         //angular
double avg_angaccx_dash_from_air_4[Max_y_grid],avg_angaccy_dash_from_air_4[Max_y_grid],avg_angaccz_dash_from_air_4[Max_y_grid];

double avg_vx_dash_air_2[Max_y_grid],avg_vy_dash_air_2[Max_y_grid], avg_vz_dash_air_2[Max_y_grid];
double avg_vx_dash_air_3[Max_y_grid],avg_vy_dash_air_3[Max_y_grid], avg_vz_dash_air_3[Max_y_grid];
double avg_vx_dash_air_4[Max_y_grid],avg_vy_dash_air_4[Max_y_grid], avg_vz_dash_air_4[Max_y_grid];

double avg_angvelx_dash_air_2[Max_y_grid],avg_angvely_dash_air_2[Max_y_grid], avg_angvelz_dash_air_2[Max_y_grid];
double avg_angvelx_dash_air_3[Max_y_grid],avg_angvely_dash_air_3[Max_y_grid], avg_angvelz_dash_air_3[Max_y_grid];                     //angular
double avg_angvelx_dash_air_4[Max_y_grid],avg_angvely_dash_air_4[Max_y_grid], avg_angvelz_dash_air_4[Max_y_grid];

extern int   Nofp_inst[Max_y_grid];
int atom_count_inst[Max_dim_position];
int   n_sample_acc_dist;
int   n_sample_angacc_dist;
int  Nofp_from_air_[Max_y_grid];
int no_of_sample_for_air_vel_fluc[Max_y_grid];
int no_of_sample_for_air_angvel_fluc[Max_y_grid];

//---------------------------------------------------------------------

extern double  sampling_time,t_zero,run_time,tot_KE;
extern double sum_B,eta;

extern int sum_Nofp[],N_rdels;
extern int  tot_atom_count[];
//int  ang_tot_atom_count[];

extern double  Total_Rep_x, Total_Rep_y, Total_Rep_z;
extern int    count_Rep;

//-----------------------------------------------

extern void set_counter();
extern void  property();
extern void final_property();
extern void distri_func();

extern void angular_property();
extern void angular_final_property();                                    //angular
extern void angular_distri_func();                          

extern void  acc_property();
extern void  call_part_velo_corre();
extern void  call_part_acc_corre();

extern void  angular_acc_property();
extern void  call_part_angular_velo_corre();                            //angular
extern void  call_part_angular_acc_corre();
  
extern void vel_dist_free_flight();
extern void vel_distri_func_free_flight();
extern void acc_dist_from_air(const FlowField& u,const Vector& x_grid, const Vector&  z_grid);
extern void acc_distri_func_from_air() ;

extern void angular_vel_dist_free_flight();
extern void angular_vel_distri_func_free_flight();
extern void angular_acc_dist_from_air(const FlowField& omega,const Vector& x_grid, const Vector&  z_grid);                  //angular
extern void angular_acc_distri_func_from_air() ;

extern void air_vel_fluc_dist(const FlowField& u,const Vector& x_grid,const Vector& z_grid, const Vector& y_grid);
extern void air_vel_distri_func() ;

void air_angular_vel_fluc_dist(const FlowField& omega,const Vector& x_grid,const Vector& z_grid, const Vector& y_grid);
extern void air_angular_vel_distri_func() ;                                                                                          //angular

void set_acc_corre();
void func_part_acc_corre();
void set_vel_corre();
void func_part_vel_corre();

void set_angular_acc_corre();
void func_part_angular_acc_corre();
void set_angular_vel_corre();
void func_part_angular_vel_corre();
extern void interpolation(const FlowField& , const Vector& ,const Vector& ,double ,double , double , int, double[]);

//---------------------------parameters to calculate velocity distribution function for free flight------
void air_vel_fluc_dist(const FlowField& u,const Vector& x_grid,const Vector& z_grid, const Vector& y_grid);
  void air_vel_distri_func() ;

extern int Nof_vx_free_flight[][2000],Nof_vy_free_flight[][2000],Nof_vz_free_flight[][2000];
extern double fof_vx_free_flight_sq[][2000],fof_vy_free_flight_sq[][2000],fof_vz_free_flight_sq[][2000];

extern int tot_atom_vel_dist_free_flight[];
int tot_atom_count_vel_free_flight[Max_dim_position];
int n_sample_vel_dist_free_flight;

//---------------------------parameters to calculate angular velocity distribution function for free flight------------
void air_angular_vel_fluc_dist(const FlowField& omega,const Vector& x_grid,const Vector& z_grid, const Vector& y_grid);
  void air_angular_vel_distri_func() ;

extern int Nof_angvelx_free_flight[][2000],Nof_angvely_free_flight[][2000],Nof_angvelz_free_flight[][2000];
extern double fof_angvelx_free_flight_sq[][2000],fof_angvely_free_flight_sq[][2000],fof_angvelz_free_flight_sq[][2000];

extern int tot_atom_angvel_dist_free_flight[];
int tot_atom_count_angvel_free_flight[Max_dim_position];
int n_sample_angvel_dist_free_flight;

//--------------------parameters for - acceleration distribution function-----
extern int tot_atom_acc_dist[Max_dim_position];
extern int   Nof_accx[][2000], Nof_accy[][2000], Nof_accz[][2000],Nof_accxy[][2000],
             Ndel_accx, Ndel_accy,Ndel_accz,Ndel_accxy;

extern double fof_accx_sq[][2000],fof_accy_sq[][2000],fof_accz_sq[][2000],fof_accxy_sq[][2000];

extern double g_x,g_y,g_z;

extern double del_accx,del_accy,del_accz,del_accxy,accx_max,accy_max,accz_max,accxy_max;
double *accx_dash,*accy_dash,*accz_dash;
double  Nof_accx_inst[Max_dim_position][2000],Nof_accy_inst[Max_dim_position][2000], 
Nof_accz_inst[Max_dim_position][2000], Nof_accxy_inst[Max_dim_position][2000];

//--------------------parameters for - angular acceleration distribution function-------
extern int tot_atom_angacc_dist[Max_dim_position];
extern int   Nof_angaccx[][2000], Nof_angaccy[][2000], Nof_angaccz[][2000],Nof_angaccxy[][2000],
             Ndel_angaccx, Ndel_angaccy,Ndel_angaccz,Ndel_angaccxy;

extern double fof_angaccx_sq[][2000],fof_angaccy_sq[][2000],fof_angaccz_sq[][2000],fof_angaccxy_sq[][2000];

extern double g_x,g_y,g_z;

extern double del_angaccx,del_angaccy,del_angaccz,del_angaccxy,angaccx_max,angaccy_max,angaccz_max,angaccxy_max;
double *angaccx_dash,*angaccy_dash,*angaccz_dash;
double  Nof_angaccx_inst[Max_dim_position][2000],Nof_angaccy_inst[Max_dim_position][2000], 
Nof_angaccz_inst[Max_dim_position][2000], Nof_angaccxy_inst[Max_dim_position][2000];

//--------------------parameters for acceleration distribution function from air velocity fluctuation---


 int  n_sample_acc_dist_from_air;
int tot_atom_acc_dist_from_air[Max_dim_position];
extern int index_acc_dist_from_air;


extern  int Nof_accx_from_air[][2000],Nof_accy_from_air[][2000],Nof_accz_from_air[][2000],Nof_accxy_from_air[][2000];
extern  int Ndel_accx_from_air,Ndel_accy_from_air, Ndel_accz_from_air, Ndel_accxy_from_air;
extern  double accx_max_from_air,accy_max_from_air,accz_max_from_air,accxy_max_from_air;
extern  double del_accx_from_air,del_accy_from_air,del_accz_from_air,del_accxy_from_air; 
extern  double fof_accx_from_air_sq[][2000],fof_accy_from_air_sq[][2000];
extern  double fof_accz_from_air_sq[][2000],fof_accxy_from_air_sq[][2000];

//--------------------parameters for angular acceleration distribution function from air angular velocity fluctuation---


 int  n_sample_angacc_dist_from_air;
int tot_atom_angacc_dist_from_air[Max_dim_position];
extern int index_ang_acc_dist_from_air;


extern  int Nof_angaccx_from_air[][2000],Nof_angaccy_from_air[][2000],Nof_angaccz_from_air[][2000],Nof_angaccxy_from_air[][2000];
extern  int Ndel_angaccx_from_air,Ndel_angaccy_from_air, Ndel_angaccz_from_air, Ndel_angaccxy_from_air;
extern  double angaccx_max_from_air,angaccy_max_from_air,angaccz_max_from_air,angaccxy_max_from_air;
extern  double del_angaccx_from_air,del_angaccy_from_air,del_angaccz_from_air,del_angaccxy_from_air; 
extern  double fof_angaccx_from_air_sq[][2000],fof_angaccy_from_air_sq[][2000];
extern  double fof_angaccz_from_air_sq[][2000],fof_angaccxy_from_air_sq[][2000];

//---------------------parameters for distribution function of air velocity fluctuation-----------------
extern int Nof_sample_airvel_dist[Max_dim_position];
extern int Nof_air_velx[][2000], Nof_air_vely[][2000], Nof_air_velz[][2000];
extern int Ndel_air_velx,Ndel_air_vely,Ndel_air_velz;


extern double max_air_velx,max_air_vely,max_air_velz;
extern double del_air_velx_dist,del_air_vely_dist,del_air_velz_dist;
 int index_air_vel_fluc_dist;

//---------------------parameters for distribution function of air angular velocity fluctuation-----------------
extern int Nof_sample_airangvel_dist[Max_dim_position];
extern int Nof_air_angvelx[][2000], Nof_air_angvely[][2000], Nof_air_angvelz[][2000];
extern int Ndel_air_angvelx,Ndel_air_angvely,Ndel_air_angvelz;


extern double max_air_angvelx,max_air_angvely,max_air_angvelz;
extern double del_air_angvelx_dist,del_air_angvely_dist,del_air_angvelz_dist;
 int index_air_angvel_fluc_dist;

//---------------------parameter for particle acceleration correlation------------

int Nmax_part_acc_corre,Nof_part_corre;
int N_count_acc_corre;

int Nmax_part_vel_corre,Nof_part_vel_corre;
int N_count_vel_corre;
int n_del_y_acc_corre;
int index_part_y_post[n_atom+2];

double acc_initial_time;
double  vel_corre_initial_time;
 
int *part_acc_corre;
int *part_vel_corre;

double *accx_dash_initial,*accy_dash_initial,*accz_dash_initial;
//double *norm_x,*norm_y,*norm_z,*norm_xy,*norm_yx;
double *pvx_dash_initial, *pvy_dash_initial, *pvz_dash_initial;

double *avg_norm_x,*avg_norm_y,*avg_norm_z,*avg_norm_xy,*avg_norm_yx;
double avg_vel_norm_x,avg_vel_norm_y,avg_vel_norm_z,avg_vel_norm_xy,avg_vel_norm_yx;
double del_y_acc_corre;

void spline(double *x,double *y, int n, double *y2);
int locate(const double xx[], int n, double x);
double splint(double *xa, double *ya, double *y2a, int n, double x);
void intake_average();

int index_average_from_input;
int Max_dns_grid;

int index_stress_varience,index_fourth_moment;

double *y2_part_velo,*y2_air_velo;
double *y_part_grid,*y_air_grid;
double *avg_vx_air,avg_vx_part[Max_y_grid],avg_vy_part[Max_y_grid],avg_vz_part[Max_y_grid];
extern double tau_vp;

int min_Nop_for_corre;
int index_return_from_part_acc_corre;
int index_return_from_corre;
double y_min_for_corre,y_max_for_corre;    

//   extern int N_count_acc_corre;
//   extern void set_acc_corre();
//   extern void part_acc_corre();

//---------------------parameter for particle angular acceleration correlation------------

int Nmax_part_angacc_corre,Nof_part_ang_corre;
int N_count_angacc_corre;

int Nmax_part_angvel_corre,Nof_part_angvel_corre;
int N_count_angvel_corre;
int n_del_y_angacc_corre;
int index_part_ang_y_post[n_atom+2];

double angacc_initial_time;
double  angvel_corre_initial_time;
 
int *part_angacc_corre;
int *part_angvel_corre;

double *angaccx_dash_initial,*angaccy_dash_initial,*angaccz_dash_initial;
//double *norm_x,*norm_y,*norm_z,*norm_xy,*norm_yx;
double *pangvelx_dash_initial, *pangvely_dash_initial, *pangvelz_dash_initial;

double *angavg_norm_x,*angavg_norm_y,*angavg_norm_z,*angavg_norm_xy,*angavg_norm_yx;
double avg_angvel_norm_x,avg_angvel_norm_y,avg_angvel_norm_z,avg_angvel_norm_xy,avg_angvel_norm_yx;
double del_y_angacc_corre;


//--------------

double *avg_angvelx_air,avg_angvelx_part[Max_y_grid],avg_angvely_part[Max_y_grid],avg_angvelz_part[Max_y_grid];

int index_return_from_part_angacc_corre;
int angular_index_return_from_corre;
//double ang_y_min_for_corre,y_max_for_corre; 

//---------------------parameters for air velocity  correlation at particle position  ------------

int N_count_fluid_vel_corre,n_del_y_air_vel_corre;

int part_y_pos_for_air_vel_corre[n_atom+2];
int consider_part_air_vel_corre[n_atom+2];
int index_return_from_air_vel_corre;

extern int index_air_vel_corre;
int Nmax_part_air_vel_corre;

double del_y_air_vel_corre,air_vel_initial_time;

double vx_dash_air_vel_corre[n_atom+2],vy_dash_air_vel_corre[n_atom+2],vz_dash_air_vel_corre[n_atom+2];
double *vx_dash_air_initial,*vy_dash_air_initial,*vz_dash_air_initial;
double *avg_norm_air_x,*avg_norm_air_y,*avg_norm_air_z,*avg_norm_air_xy,*avg_norm_air_yx;

extern void call_fluid_vel_corre(const FlowField& ,const Vector& ,const Vector& );

void set_fluid_vel_corre();
void func_fluid_vel_corre();
//----------------------------------------------------------------------------------------------

//---------------------parameters for air angular velocity  correlation at particle position  ------------

int N_count_fluid_angvel_corre,n_del_y_air_angvel_corre;

int part_y_pos_for_air_angvel_corre[n_atom+2];
int consider_part_air_angvel_corre[n_atom+2];
int index_return_from_air_angvel_corre;

extern int index_air_angvel_corre;
int Nmax_part_air_angvel_corre;

double del_y_air_angvel_corre,air_angvel_initial_time;

double angvelx_dash_air_angvel_corre[n_atom+2],angvely_dash_air_angvel_corre[n_atom+2],angvelz_dash_air_angvel_corre[n_atom+2];
double *angvelx_dash_air_initial,*angvely_dash_air_initial,*angvelz_dash_air_initial;
double *angavg_norm_air_x,*angavg_norm_air_y,*angavg_norm_air_z,*angavg_norm_air_xy,*angavg_norm_air_yx;

extern void angular_call_fluid_vel_corre(const FlowField& ,const Vector& ,const Vector& );

void set_fluid_angular_vel_corre();
void func_fluid_angular_vel_corre();
//----------------------------------------------------------------------------------------------


int position_index;
int max_position_index;



FILE *fp18;
FILE *fp19;
FILE *fp20;
FILE *fp21;
FILE *fp211;
FILE *fp212;
FILE *fp213;
FILE *fp214;
FILE *fp315;
FILE *fp316;

FILE *fp1018;
FILE *fp1019;
FILE *fp1020;
FILE *fp1021;
FILE *fp10211;
FILE *fp10212;                     //angular file save
FILE *fp10213;
FILE *fp10214;
FILE *fp10315;
FILE *fp10316;
void  set_counter()
    {	
      int i,j,ny;
      //int nx,nz, N_shell;
   string rough_line;
      double rough;

	 sum_B= 0.0;
	 sampls=0;
	 sum_KE=0.0;
        t_zero=run_time;
	Tot_Nof_coll=0;
	Tot_Nof_p_p_coll=0;
	Nof_p_w_coll=0;

	
	n_sample_acc_dist=0;
        n_sample_angacc_dist=0;                                           //angular

  n_sample_acc_dist_from_air=0; 
  n_sample_vel_dist_free_flight=0;

  n_sample_angacc_dist_from_air=0;
  n_sample_angvel_dist_free_flight=0;                                       //angular

  for(j=0;j<Max_dim_position;j++) { Nof_sample_airvel_dist[j]=0;   Nof_sample_airangvel_dist[j]=0;}


 for(i=0;i<Max_dim_position;i++)
   {
     tot_atom_count[i]=0;
     tot_atom_acc_dist[i]=0;
     tot_atom_count_vel_free_flight[i]=0;
     tot_atom_acc_dist_from_air[i]=0;

     //ang_tot_atom_count[i]=0;
     tot_atom_angacc_dist[i]=0;
     tot_atom_count_angvel_free_flight[i]=0;                                    //angular
     tot_atom_angacc_dist_from_air[i]=0;
   }
 



 for(j=0;j<Max_dim_position;j++)
   for(i=0;i<=Ndel_v;i++) Nof_v[j][i]=0;

  for(j=0;j<Max_dim_position;j++)
   for(i=0;i<=Ndel_angvel;i++) Nof_angvel[j][i]=0;                              //angular

 for(j=0;j<Max_dim_position;j++)
 for(i=0;i<=Ndel_vx;i++){
   Nof_vx[j][i]=0;
   Nof_vx_free_flight[j][i]=0;
   fof_vx_sq[j][i]=0;
   fof_vx_free_flight_sq[j][i]=0;
    }

 for(j=0;j<Max_dim_position;j++)
 for(i=0;i<=Ndel_angvelx;i++){
   Nof_angvelx[j][i]=0;
   Nof_angvelx_free_flight[j][i]=0;                                           //angular
   fof_angvelx_sq[j][i]=0;
   fof_angvelx_free_flight_sq[j][i]=0;
    }

 for(j=0;j<Max_dim_position;j++)
  for(i=0;i<=Ndel_vy;i++){
  Nof_vy[j][i]=0;
  Nof_vy_free_flight[j][i]=0;
  fof_vy_sq[j][i]=0;
  fof_vy_free_flight_sq[j][i]=0;
 }

 for(j=0;j<Max_dim_position;j++)
  for(i=0;i<=Ndel_angvely;i++){
  Nof_angvely[j][i]=0;
  Nof_angvely_free_flight[j][i]=0;                                         //angular
  fof_angvely_sq[j][i]=0;
  fof_angvely_free_flight_sq[j][i]=0;
 }

 for(j=0;j<Max_dim_position;j++)
 for(i=0;i<=Ndel_vz;i++){
   Nof_vz[j][i]=0;
   Nof_vz_free_flight[j][i]=0;
   fof_vz_sq[j][i]=0;
   fof_vz_free_flight_sq[j][i]=0;
    }

for(j=0;j<Max_dim_position;j++)
 for(i=0;i<=Ndel_angvelz;i++){
   Nof_angvelz[j][i]=0;
   Nof_angvelz_free_flight[j][i]=0;                                        //angular
   fof_angvelz_sq[j][i]=0;
   fof_angvelz_free_flight_sq[j][i]=0;
    }

  for(j=0;j<Max_dim_position;j++)
    for(i=0;i<=Ndel_accx;i++){ Nof_accx[j][i]=0;  fof_accx_sq[j][i]=0;}

for(j=0;j<Max_dim_position;j++)
    for(i=0;i<=Ndel_angaccx;i++){ Nof_angaccx[j][i]=0;  fof_angaccx_sq[j][i]=0;}                //angular

 for(j=0;j<Max_dim_position;j++)
 for(i=0;i<=Ndel_accy;i++) 
{Nof_accy[j][i]=0; fof_accy_sq[j][i]=0;}

for(j=0;j<Max_dim_position;j++)
 for(i=0;i<=Ndel_angaccy;i++) 
{Nof_angaccy[j][i]=0; fof_angaccy_sq[j][i]=0;}                                    //angular

  for(j=0;j<Max_dim_position;j++)  
 for(i=0;i<=Ndel_accz;i++) 
{Nof_accz[j][i]=0; fof_accz_sq[j][i]=0;}

for(j=0;j<Max_dim_position;j++)  
 for(i=0;i<=Ndel_angaccz;i++)                                                   //angular 
{Nof_angaccz[j][i]=0; fof_angaccz_sq[j][i]=0;}

 for(j=0;j<Max_dim_position;j++) 
 for(i=0;i<=Ndel_accxy;i++) 
{Nof_accxy[j][i]=0; fof_accxy_sq[j][i]=0;}

for(j=0;j<Max_dim_position;j++) 
 for(i=0;i<=Ndel_angaccxy;i++)                                                //angular 
{Nof_angaccxy[j][i]=0; fof_angaccxy_sq[j][i]=0;}

 for(j=0;j<Max_dim_position;j++) 
 for(i=0;i<=Ndel_accx_from_air;i++)
 {Nof_accx_from_air[j][i]=0;fof_accx_from_air_sq[j][i]=0;}

for(j=0;j<Max_dim_position;j++) 
 for(i=0;i<=Ndel_angaccx_from_air;i++)                                          //angular
 {Nof_angaccx_from_air[j][i]=0;fof_angaccx_from_air_sq[j][i]=0;}

 for(j=0;j<Max_dim_position;j++) 
 for(i=0;i<=Ndel_accy_from_air;i++)
{Nof_accy_from_air[j][i]=0;fof_accy_from_air_sq[j][i]=0;}

for(j=0;j<Max_dim_position;j++) 
 for(i=0;i<=Ndel_angaccy_from_air;i++)                                            //angular
{Nof_angaccy_from_air[j][i]=0;fof_angaccy_from_air_sq[j][i]=0;}

 for(j=0;j<Max_dim_position;j++) 
 for(i=0;i<=Ndel_accz_from_air;i++)
{Nof_accz_from_air[j][i]=0;fof_accz_from_air_sq[j][i]=0;}

for(j=0;j<Max_dim_position;j++) 
 for(i=0;i<=Ndel_angaccz_from_air;i++)                                              //angular
{Nof_angaccz_from_air[j][i]=0;fof_angaccz_from_air_sq[j][i]=0;}

 for(j=0;j<Max_dim_position;j++) 
 for(i=0;i<=Ndel_accxy_from_air;i++)
{Nof_accxy_from_air[j][i]=0;fof_accxy_from_air_sq[j][i]=0;}

for(j=0;j<Max_dim_position;j++) 
 for(i=0;i<=Ndel_angaccxy_from_air;i++)                                            //angular
{Nof_angaccxy_from_air[j][i]=0;fof_angaccxy_from_air_sq[j][i]=0;}

  for(j=0;j<Max_dim_position;j++) for(i=0;i<=Ndel_air_velx;i++)Nof_air_velx[j][i]=0;
  for(j=0;j<Max_dim_position;j++) for(i=0;i<=Ndel_air_vely;i++)Nof_air_vely[j][i]=0;
  for(j=0;j<Max_dim_position;j++)  for(i=0;i<=Ndel_air_velz;i++)Nof_air_velz[j][i]=0;

  for(j=0;j<Max_dim_position;j++) for(i=0;i<=Ndel_air_angvelx;i++)Nof_air_angvelx[j][i]=0;
  for(j=0;j<Max_dim_position;j++) for(i=0;i<=Ndel_air_angvely;i++)Nof_air_angvely[j][i]=0;               //angular
  for(j=0;j<Max_dim_position;j++)  for(i=0;i<=Ndel_air_angvelz;i++)Nof_air_angvelz[j][i]=0;

 for(ny=0;ny<Ndel_y;ny++)
      for(j=1;j<=N_rdels;j++)
	Ngof_r[ny][j]=0;

 for(ny=0;ny<Ndel_y;ny++)
   N_atom_count[ny]=0;
 

 for(ny=0;ny<Ndel_y;ny++)
  Nof_p_p_coll[ny]=0;         
                     
 for(ny=0;ny<Ndel_y;ny++)
  {
    vx_mean[ny]=0;
    vy_mean[ny]=0;
    vz_mean[ny]=0;

    symmetric_p_xx_mean[ny]=0.0;
    symmetric_p_xy_mean[ny]=0.0;
    symmetric_p_xz_mean[ny]=0.0;
    symmetric_p_yy_mean[ny]=0.0;
    symmetric_p_yz_mean[ny]=0.0;
    symmetric_p_zz_mean[ny]=0.0;

    symmetric_particle_corr_xx_mean[ny] = 0.0;
    symmetric_particle_corr_xy_mean[ny] = 0.0;
    symmetric_particle_corr_xz_mean[ny] = 0.0;
    symmetric_particle_corr_yy_mean[ny] = 0.0;
    symmetric_particle_corr_yz_mean[ny] = 0.0;
    symmetric_particle_corr_zz_mean[ny] = 0.0;

    angvelx_mean[ny]=0;
    angvely_mean[ny]=0;                              //angular
    angvelz_mean[ny]=0;

    angvel_corr_x_mean[ny]=0;
    angvel_corr_y_mean[ny]=0;                              //angular
    angvel_corr_z_mean[ny]=0;

    Nofp_mean[ny]=0;
 Nofp_inst_sq[ny]=0.0;

 sigma_xx[ny]=0;
 sigma_yy[ny]=0;
 sigma_zz[ny]=0;   	    
 tau_xy[ny]=0;

 sigma_xx_1[ny]=0.0;
 sigma_yy_1[ny]=0.0;
 sigma_zz_1[ny]=0.0;
 sigma_xy_1[ny]=0.0;

 var_sigma_xx[ny]=0;
 var_sigma_yy[ny]=0;
 var_sigma_zz[ny]=0;
 var_tau_xy[ny]=0;


 angsigma_xx[ny]=0;
 angsigma_yy[ny]=0;                            //angular
 angsigma_zz[ny]=0;   	    
 angtau_xy[ny]=0;

 angsigma_xx_1[ny]=0.0;
 angsigma_yy_1[ny]=0.0;
 angsigma_zz_1[ny]=0.0;                       //angular
 angsigma_xy_1[ny]=0.0;

 var_angsigma_xx[ny]=0;
 var_angsigma_yy[ny]=0;                      //angular
 var_angsigma_zz[ny]=0;
 var_angtau_xy[ny]=0;


 avg_pvx_dash_4[ny]=0;
 avg_pvy_dash_4[ny]=0;
 avg_pvz_dash_4[ny]=0;
 avg_pvxy_dash_4[ny]=0;

 avg_pvx_dash_3[ny]=0;
 avg_pvy_dash_3[ny]=0;
 avg_pvz_dash_3[ny]=0;

 avg_accx_dash_2[ny]=0;
 avg_accy_dash_2[ny]=0;
 avg_accz_dash_2[ny]=0;

mean_accx_sq[ny]=0.0;
 mean_accy_sq[ny]=0.0;
 mean_accz_sq[ny]=0.0;

 accx_mean[ny]=0.0;
accy_mean[ny]=0.0;
accz_mean[ny]=0.0;

 avg_accx_dash_3[ny]=0;
 avg_accy_dash_3[ny]=0;
 avg_accz_dash_3[ny]=0;

 avg_accx_dash_4[ny]=0;
 avg_accy_dash_4[ny]=0;
 avg_accz_dash_4[ny]=0;

 avg_accx_dash_from_air_2[ny]=0;
 avg_accy_dash_from_air_2[ny]=0;
 avg_accz_dash_from_air_2[ny]=0;

 avg_accx_dash_from_air_3[ny]=0;
 avg_accy_dash_from_air_3[ny]=0;
 avg_accz_dash_from_air_3[ny]=0;

 avg_accx_dash_from_air_4[ny]=0;
 avg_accy_dash_from_air_4[ny]=0;
 avg_accz_dash_from_air_4[ny]=0;


//------------for angular velocities and acceleration
 avg_pangvelx_dash_4[ny]=0;
 avg_pangvely_dash_4[ny]=0;
 avg_pangvelz_dash_4[ny]=0;
 avg_pangvelxy_dash_4[ny]=0;

 avg_pangvelx_dash_3[ny]=0;
 avg_pangvely_dash_3[ny]=0;
 avg_pangvelz_dash_3[ny]=0;

 avg_angaccx_dash_2[ny]=0;
 avg_angaccy_dash_2[ny]=0;
 avg_angaccz_dash_2[ny]=0;

mean_angaccx_sq[ny]=0.0;
 mean_angaccy_sq[ny]=0.0;
 mean_angaccz_sq[ny]=0.0;

 angaccx_mean[ny]=0.0;
angaccy_mean[ny]=0.0;
angaccz_mean[ny]=0.0;


 avg_angaccx_dash_3[ny]=0;
 avg_angaccy_dash_3[ny]=0;
 avg_angaccz_dash_3[ny]=0;

 avg_angaccx_dash_4[ny]=0;
 avg_angaccy_dash_4[ny]=0;
 avg_angaccz_dash_4[ny]=0;

 avg_angaccx_dash_from_air_2[ny]=0;
 avg_angaccy_dash_from_air_2[ny]=0;
 avg_angaccz_dash_from_air_2[ny]=0;

 avg_angaccx_dash_from_air_3[ny]=0;
 avg_angaccy_dash_from_air_3[ny]=0;
 avg_angaccz_dash_from_air_3[ny]=0;

 avg_angaccx_dash_from_air_4[ny]=0;
 avg_angaccy_dash_from_air_4[ny]=0;
 avg_angaccz_dash_from_air_4[ny]=0;
//------------------------------------------
  }

 for(ny=0;ny<Max_dns_grid;ny++)
{
  avg_vx_dash_air_2[ny]=0;
 avg_vy_dash_air_2[ny]=0;
 avg_vz_dash_air_2[ny]=0;

 avg_vx_dash_air_3[ny]=0; 
 avg_vy_dash_air_3[ny]=0; 
 avg_vz_dash_air_3[ny]=0; 

 avg_vx_dash_air_4[ny]=0; 
 avg_vy_dash_air_4[ny]=0; 
 avg_vz_dash_air_4[ny]=0; 

 avg_vx_dash_air_2[ny]=0;
 avg_vy_dash_air_2[ny]=0;
 avg_vz_dash_air_2[ny]=0;

 avg_vx_dash_air_3[ny]=0; 
 avg_vy_dash_air_3[ny]=0; 
 avg_vz_dash_air_3[ny]=0; 

 avg_vx_dash_air_4[ny]=0; 
 avg_vy_dash_air_4[ny]=0; 
 avg_vz_dash_air_4[ny]=0; 


//--------------angular--------------
  avg_angvelx_dash_air_2[ny]=0;
 avg_angvely_dash_air_2[ny]=0;
 avg_angvelz_dash_air_2[ny]=0;

 avg_angvelx_dash_air_3[ny]=0; 
 avg_angvely_dash_air_3[ny]=0; 
 avg_angvelz_dash_air_3[ny]=0; 

 avg_angvelx_dash_air_4[ny]=0; 
 avg_angvely_dash_air_4[ny]=0; 
 avg_angvelz_dash_air_4[ny]=0; 

 avg_angvelx_dash_air_2[ny]=0;
 avg_angvely_dash_air_2[ny]=0;
 avg_angvelz_dash_air_2[ny]=0;

 avg_angvelx_dash_air_3[ny]=0; 
 avg_angvely_dash_air_3[ny]=0; 
 avg_angvelz_dash_air_3[ny]=0; 

 avg_angvelx_dash_air_4[ny]=0; 
 avg_angvely_dash_air_4[ny]=0; 
 avg_angvelz_dash_air_4[ny]=0;
//------------------------------------------

   no_of_sample_for_air_vel_fluc[ny]=0;
   no_of_sample_for_air_angvel_fluc[ny]=0;

  }

 Total_Rep_x=0.0;
 Total_Rep_y=0.0;
 Total_Rep_z=0.0;
 count_Rep=0;

 if(index_stress_varience==1)
   {
ifstream infile_avg_part_stress("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/stress_0.txt");
getline(infile_avg_part_stress,rough_line);

 for(ny=0;ny<Ndel_y;ny++)
   {
     infile_avg_part_stress>>rough>>rough>>rough>>avg_stress_x[ny]>>rough>>avg_stress_y[ny]>> rough>>avg_stress_z[ny]>>rough>>avg_stress_xy[ny];
     
     //cout<<"*******  "<<avg_stress_x[ny]<<"  "<<avg_stress_y[ny]<<"  "<<avg_stress_z[ny]<<" "<<avg_stress_xy[ny]<<endl;
     // cin.get();
     }
   infile_avg_part_stress.close();

  ifstream infile_avg_part_angstress("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angstress_0.txt");
getline(infile_avg_part_angstress,rough_line);

 for(ny=0;ny<Ndel_y;ny++)
   {
     infile_avg_part_angstress>>rough>>rough>>rough>>avg_angstress_x[ny]>>rough>>avg_angstress_y[ny]>> rough>>avg_angstress_z[ny]>>rough>>avg_angstress_xy[ny];
     
     //cout<<"*******  "<<avg_angstress_x[ny]<<"  "<<avg_angstress_y[ny]<<"  "<<avg_angstress_z[ny]<<" "<<avg_angstress_xy[ny]<<endl;
     // cin.get();
     }
 infile_avg_part_angstress.close();
   }

 if(index_fourth_moment==1 || index_stress_varience==1)
   {
 
ifstream infile_avg_part_velo("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_0.txt");
  getline(infile_avg_part_velo,rough_line);
  getline(infile_avg_part_velo,rough_line);
 
 for(ny=0;ny<Ndel_y;ny++)
   {
     infile_avg_part_velo>>rough>>rough>>avg_vx_part[ny]>>avg_vy_part[ny]>>avg_vz_part[ny]>>rough>>rough;
     
     // cout<<"*******  "<<avg_vx_part[ny]<<" "<<avg_vy_part[ny]<<" "<<avg_vz_part[ny]<<endl;
     //cin.get();
     }

 infile_avg_part_velo.close();  



  ifstream infile_avg_part_angvelo("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angvelocity_0.txt");
  getline(infile_avg_part_angvelo,rough_line);
  getline(infile_avg_part_angvelo,rough_line);
 
 for(ny=0;ny<Ndel_y;ny++)
   {
     infile_avg_part_angvelo>>rough>>rough>>avg_angvelx_part[ny]>>avg_angvely_part[ny]>>avg_angvelz_part[ny]>>rough>>rough;
     
     // cout<<"*******  "<<avg_angvelx_part[ny]<<" "<<avg_angvely_part[ny]<<" "<<avg_angvelz_part[ny]<<endl;
     //cin.get();
     }

 infile_avg_part_angvelo.close();

   }

 if(index_acc_dist_from_air==1|| index_air_vel_corre==1||index_air_vel_fluc_dist==1)intake_average();



    	     return;
    }





//-----------------------------------------------------------------------
   // CALCULATION OF THE PROPERTIES
   //--------------------------------------------------------------------
  void  property()
     {   
       int i,j, ipi, nm1,ny,kk,jj;
       int ishell ;
      
       int N_shell,N_shellx, N_shelly,N_shellz;
       
       // int k,N_rdels,

       //double  TELASP,BELASP, VIRIAL;
       
      
       double ave_KE,rx, ry, rz,  r_ij;
       double  offset,offset_x, offset_y , offset_z ;
       
       //double fof_v, fatom,HH, HINST;


       // double pi = 3.141592650;
       //FILE *fp2;

       //FILE *fp3;
 
                        
	 sampls=sampls+1;
	 n_sample=n_sample+1;
	 // Accumulate running average kinetic energy
                        
  sum_KE = sum_KE + tot_KE;
  ave_KE = sum_KE/(double)(n_sample);
                       
                    
 //-----------------------------------------------
 // Storing couner for Radial distribution function g(r)
 //------------------------------------------------                        
            
   nm1=n_atom-1;
		  for( i=0;i<nm1;i++)
		    {
		      p_xi=p_x[i];
		      p_yi=p_y[i];
		      p_zi=p_z[i];
		      ipi=i+1;
                    
		      for(j=ipi;j<n_atom;j++)
		      {
                      rx=p_xi-p_x[j];
                      ry=p_yi-p_y[j];
                      rz=p_zi-p_z[j];

                     
                      
              if( rx > 0.50*b_lx)
		rx= rx - 1.0*b_lx ; // Apply minimum image criteria
             
             if( rz > 0.50*b_lz)
	       rz= rz - 1.0*b_lz;
                      
              if(rx < -0.50*b_lx)
		rx=rx + 1.0*b_lx;
              
              if (rz < -0.50*b_lz) 
		rz=rz + 1.0*b_lz;
	      
	       ny=(int)((((p_y[i]+p_y[j])/2.0)-(sigma/2.0))/del_y);
            
	       r_ij = sqrt(rx*rx +rz*rz);   //g(r) is calculated at different y-plane

                      if(r_ij<= 0.50*b_lx) 
			{
			  ishell =(int)(r_ij/del_r + 0.50);
                         
		    Ngof_r[ny][ishell]=Ngof_r[ny][ishell]+1;   //ny is used to get at different grid point
	       N_atom_count[ny]=N_atom_count[ny]+1;        //  N_atom_count is used to store no of atom used as origin   

					   }
                         
		      }
		    }
                        
	 
//-------------------------------------------------------------                  
 //  Calculation for numberdensity , particle velocities and particle stress
//-------------------------------------------------------------- 

 j=n_sample;

    for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
    Nofp_inst[ny]=0;
    vx_inst[ny]=0.0;
    vy_inst[ny]=0.0;
    vz_inst[ny]=0.0;

    symmetric_p_xx_inst[ny]=0.0;
    symmetric_p_xy_inst[ny]=0.0;
    symmetric_p_xz_inst[ny]=0.0;
    symmetric_p_yy_inst[ny]=0.0;
    symmetric_p_yz_inst[ny]=0.0;
    symmetric_p_zz_inst[ny]=0.0;
    
    pvx_dash_sq[ny]=0.0;
    pvy_dash_sq[ny]=0.0;
    pvz_dash_sq[ny]=0.0;
    pvxy_dash[ny]=0.0;

    vx_vx_inst[ny]=0.0;
    vy_vy_inst[ny]=0.0;
    vz_vz_inst[ny]=0.0;
    vx_vy_inst[ny]=0.0;

    
           }

    for(i=0;i<n_atom;i++)// summing up all the velocities in each the y-grid for all the particle 
    {
      ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

      // cout<<"in property"<<p_y[i]<<" "<<ny<<endl;


              Nofp_inst[ny]= Nofp_inst[ny]+ 1;

              vx_inst[ny]= vx_inst[ny]+p_vx[i];
              vy_inst[ny]= vy_inst[ny]+p_vy[i];  
              vz_inst[ny]= vz_inst[ny]+p_vz[i];
             
              symmetric_p_xx_inst[ny] = symmetric_p_xx_inst[ny] + symmetric_air_p_xx[i]; 
              symmetric_p_xy_inst[ny] = symmetric_p_xy_inst[ny] + symmetric_air_p_xy[i];  
              symmetric_p_xz_inst[ny] = symmetric_p_xz_inst[ny] + symmetric_air_p_xz[i];   
              symmetric_p_yy_inst[ny] = symmetric_p_yy_inst[ny] + symmetric_air_p_yy[i];   
              symmetric_p_yz_inst[ny] = symmetric_p_yz_inst[ny] + symmetric_air_p_yz[i];        
              symmetric_p_zz_inst[ny] = symmetric_p_zz_inst[ny] + symmetric_air_p_zz[i]; 
      }

    for(ny=0;ny<Ndel_y;ny++)                // for calculating density fluctuation(error bar)
      {
	Nofp_inst_sq[ny]= Nofp_inst_sq[ny]+Nofp_inst[ny]*Nofp_inst[ny];
      }


// --------------newly added for stress callculation by <u^2>-Uavg^2 folmulea

 for(i=0;i<n_atom;i++)// summing up all the velocities in each the y-grid for all the particle 
    {
      ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

      // cout<<"in property"<<p_y[i]<<" "<<ny<<endl;


      //       Nofp_inst[ny]= Nofp_inst[ny]+ 1;

              vx_vx_inst[ny]= vx_vx_inst[ny]+p_vx[i]*p_vx[i];
              vy_vy_inst[ny]= vy_vy_inst[ny]+p_vy[i]*p_vy[i];  
              vz_vz_inst[ny]= vz_vz_inst[ny]+p_vz[i]*p_vz[i];
	      vx_vy_inst[ny]= vx_vy_inst[ny]+p_vx[i]*p_vy[i];  


      }

 for(ny=0;ny<Ndel_y;ny++) // doing the particle average in each y-grid-in each sample frame
       {
    if(Nofp_inst[ny]>0.0)
      {
 vx_vx_inst[ny]= vx_vx_inst[ny]/Nofp_inst[ny];
 vy_vy_inst[ny]= vy_vy_inst[ny]/Nofp_inst[ny];
 vz_vz_inst[ny]= vz_vz_inst[ny]/Nofp_inst[ny];
 vx_vy_inst[ny]= vx_vy_inst[ny]/Nofp_inst[ny];
      }

    sigma_xx_1[ny]=sigma_xx_1[ny]+ vx_vx_inst[ny];
    sigma_yy_1[ny]=sigma_yy_1[ny]+ vy_vy_inst[ny];
    sigma_zz_1[ny]=sigma_zz_1[ny]+ vz_vz_inst[ny];
    sigma_xy_1[ny]=sigma_xy_1[ny]+ vx_vy_inst[ny];

       }

 //------------------back to the old stress calculation by calculation fluc over particle avg----------

  for(ny=0;ny<Ndel_y;ny++) // doing the particle average in each y-grid-in each sample frame
       {
    if(Nofp_inst[ny]>0.0)
      {
     vx_inst[ny]= vx_inst[ny]/Nofp_inst[ny];  // vx_inst[ny] is basically particle averaged inst in time.
     vy_inst[ny]= vy_inst[ny]/Nofp_inst[ny];
     vz_inst[ny]= vz_inst[ny]/Nofp_inst[ny];

     symmetric_p_xx_inst[ny] = symmetric_p_xx_inst[ny]/Nofp_inst[ny];
     symmetric_p_xy_inst[ny] = symmetric_p_xy_inst[ny]/Nofp_inst[ny];
     symmetric_p_xz_inst[ny] = symmetric_p_xz_inst[ny]/Nofp_inst[ny];
     symmetric_p_yy_inst[ny] = symmetric_p_yy_inst[ny]/Nofp_inst[ny];
     symmetric_p_yz_inst[ny] = symmetric_p_yz_inst[ny]/Nofp_inst[ny];
     symmetric_p_zz_inst[ny] = symmetric_p_zz_inst[ny]/Nofp_inst[ny];

      }
     vx_mean[ny]=vx_mean[ny]+vx_inst[ny];    // data storing for the time averaging but take care 
     vy_mean[ny]=vy_mean[ny]+vy_inst[ny];   //---- that this quantities are not average
     vz_mean[ny]=vz_mean[ny]+vz_inst[ny];

     symmetric_p_xx_mean[ny] = symmetric_p_xx_mean[ny] + symmetric_p_xx_inst[ny];
     symmetric_p_xy_mean[ny] = symmetric_p_xy_mean[ny] + symmetric_p_xy_inst[ny];
     symmetric_p_xz_mean[ny] = symmetric_p_xz_mean[ny] + symmetric_p_xz_inst[ny];
     symmetric_p_yy_mean[ny] = symmetric_p_yy_mean[ny] + symmetric_p_yy_inst[ny];
     symmetric_p_yz_mean[ny] = symmetric_p_yz_mean[ny] + symmetric_p_yz_inst[ny];
     symmetric_p_zz_mean[ny] = symmetric_p_zz_mean[ny] + symmetric_p_zz_inst[ny];

     Nofp_mean[ny]=Nofp_mean[ny]+Nofp_inst[ny];
    
       }


 for(i=0;i<n_atom;i++)
   {
 ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
 pvx_dash[i]=p_vx[i]-vx_inst[ny];  // calculating fluctuation for each particle
 pvy_dash[i]=p_vy[i]-vy_inst[ny];
 pvz_dash[i]=p_vz[i]-vz_inst[ny];
 pv_dash[i]=sqrt(p_vx[i]*p_vx[i]+p_vy[i]*p_vy[i]+p_vz[i]*p_vz[i])
             -sqrt(vx_inst[ny]*vx_inst[ny]+vy_inst[ny]*vy_inst[ny]+vz_inst[ny]*vz_inst[ny]);

 pvx_dash_sq[ny]=pvx_dash_sq[ny]+pvx_dash[i]*pvx_dash[i];
 pvy_dash_sq[ny]=pvy_dash_sq[ny]+pvy_dash[i]*pvy_dash[i];
 pvz_dash_sq[ny]=pvz_dash_sq[ny]+pvz_dash[i]*pvz_dash[i];
 pvxy_dash[ny]=pvxy_dash[ny]+pvx_dash[i]*pvy_dash[i];
   }

 for(ny=0;ny<Ndel_y;ny++) 
  {
    if(Nofp_inst[ny]>0.0)
      {
 pvx_dash_sq[ny]= pvx_dash_sq[ny]/Nofp_inst[ny];// doing the particle average i
 pvy_dash_sq[ny]= pvy_dash_sq[ny]/Nofp_inst[ny];
 pvz_dash_sq[ny]= pvz_dash_sq[ny]/Nofp_inst[ny];
 pvxy_dash[ny]  = pvxy_dash[ny]/Nofp_inst[ny];
      }
    sigma_xx[ny]=sigma_xx[ny]+ pvx_dash_sq[ny]; // data storing for the time averaging
    sigma_yy[ny]=sigma_yy[ny]+ pvy_dash_sq[ny];
    sigma_zz[ny]=sigma_zz[ny]+ pvz_dash_sq[ny];
    tau_xy[ny]=tau_xy[ny]+ pvxy_dash[ny];

    }


 /*-----------------------------------------------------------------------
Calculation of the third and the fourth moment, considering fluctuations over particle
average, because it is found to give the same result if we calculate the fluctuation over 
particle + time average.... that routine is also implemented the programme
------------------------------------------------------------------------------------*/

 for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
   
    pvx_dash_4[ny]=0.0;
    pvy_dash_4[ny]=0.0;
    pvz_dash_4[ny]=0.0;
    pvxy_dash_4[ny]=0.0;

    pvx_dash_3[ny]=0.0;
    pvy_dash_3[ny]=0.0;
    pvz_dash_3[ny]=0.0;
   
        }
      
      for(i=0;i<n_atom;i++)
	{
	  ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
	  
	  pvx_dash_4[ny]=pvx_dash_4[ny]+pvx_dash[i]*pvx_dash[i]*pvx_dash[i]*pvx_dash[i];
	  pvy_dash_4[ny]=pvy_dash_4[ny]+pvy_dash[i]*pvy_dash[i]*pvy_dash[i]*pvy_dash[i];
	  pvz_dash_4[ny]=pvz_dash_4[ny]+pvz_dash[i]*pvz_dash[i]*pvz_dash[i]*pvz_dash[i];
	
	  pvx_dash_3[ny]=pvx_dash_3[ny]+pvx_dash[i]*pvx_dash[i]*pvx_dash[i];
	  pvy_dash_3[ny]=pvy_dash_3[ny]+pvy_dash[i]*pvy_dash[i]*pvy_dash[i];
	  pvz_dash_3[ny]=pvz_dash_3[ny]+pvz_dash[i]*pvz_dash[i]*pvz_dash[i];

	     }

      for(ny=0;ny<Ndel_y;ny++) 
	{
        if(Nofp_inst[ny]>0.0)
	  {
	    pvx_dash_4[ny]= pvx_dash_4[ny]/Nofp_inst[ny];// doing the particle average i
	    pvy_dash_4[ny]= pvy_dash_4[ny]/Nofp_inst[ny];
	    pvz_dash_4[ny]= pvz_dash_4[ny]/Nofp_inst[ny];
	    
	    pvx_dash_3[ny]= pvx_dash_3[ny]/Nofp_inst[ny];// doing the particle average i
	    pvy_dash_3[ny]= pvy_dash_3[ny]/Nofp_inst[ny];
	    pvz_dash_3[ny]= pvz_dash_3[ny]/Nofp_inst[ny];

	  }

	avg_pvx_dash_4[ny]=avg_pvx_dash_4[ny]+pvx_dash_4[ny];
	avg_pvy_dash_4[ny]=avg_pvy_dash_4[ny]+pvy_dash_4[ny];
	avg_pvz_dash_4[ny]=avg_pvz_dash_4[ny]+pvz_dash_4[ny];

	avg_pvx_dash_3[ny]=avg_pvx_dash_3[ny]+pvx_dash_3[ny];
	avg_pvy_dash_3[ny]=avg_pvy_dash_3[ny]+pvy_dash_3[ny];
	avg_pvz_dash_3[ny]=avg_pvz_dash_3[ny]+pvz_dash_3[ny];
	

	}

      

      //----------------------------------------------------------


 //------------------------------------------------------------------------
      //   Storing counter for  Velocity distribution  calculation
      //--------------------------------------------------------------------------
               for(jj=0;jj<Max_dim_position;jj++)
		  for(i=0;i<Ndel_v;i++)
		    Nof_v_inst[jj][i]=0;

	       for(jj=0;jj<Max_dim_position;jj++)
		  for(i=0;i<Ndel_vx;i++)
		    Nof_vx_inst[jj][i]=0;
      
	       for(jj=0;jj<Max_dim_position;jj++)
		  for(i=0;i<Ndel_vy;i++)
		    Nof_vy_inst[jj][i]=0;

	       for(jj=0;jj<Max_dim_position;jj++)
		  for(i=0;i<Ndel_vz;i++)
		    Nof_vz_inst[jj][i]=0;

	       for(jj=0;jj<Max_dim_position;jj++)
		  atom_count_inst[jj]=0;
           
            offset =1.5+p_v_max/del_v;
            offset_x =1.5+p_vx_max/del_vx;
            offset_y =1.5+p_vy_max/del_vy;
            offset_z =1.5+p_vz_max/del_vz;


	    
          for(i=0;i<n_atom;i++)
           {
	     position_index=0;    //initialization of position_index

	     if(p_y[i]>y_min_for_distribution_1 && p_y[i]<y_max_for_distribution_1)position_index=1;
	     if(p_y[i]>y_min_for_distribution_2 && p_y[i]<y_max_for_distribution_2)position_index=2;
	     if(p_y[i]>y_min_for_distribution_3 && p_y[i]<y_max_for_distribution_3)position_index=3;
	     

	     if(position_index>0)
	       {                                                
	       
	if((pv_dash[i]>= -p_v_max && pv_dash[i]<=p_v_max) &&(pvx_dash[i]>=-p_vx_max && pvx_dash[i]<=p_vx_max)&&
           (pvy_dash[i]>=-p_vy_max && pvy_dash[i]<=p_vy_max)&& (pvz_dash[i]>=-p_vz_max && pvz_dash[i]<=p_vz_max))  
				                // to avoid exceeding arrey subscript accidentally
				 {
     // p_v[i]= sqrt(p_vx[i]*p_vx[i]+p_vy[i]*p_vy[i]+p_vz[i]*p_vz[i]);

	        N_shell=(int)(pv_dash[i]/del_v+offset);
		Nof_v[position_index][N_shell]= Nof_v[position_index][N_shell]+1;
		Nof_v_inst[position_index][N_shell]= Nof_v_inst[position_index][N_shell]+1;

             N_shellx=(int)(pvx_dash[i]/del_vx+offset_x);
             Nof_vx_inst[position_index][N_shellx]= Nof_vx_inst[position_index][N_shellx]+1;
	     Nof_vx[position_index][N_shellx]= Nof_vx[position_index][N_shellx] +1 ; //!Nofvx_s SEPERATE X COMPONENT OF Nof_vx

            	                                                     //--COMPONENT OF Nof_vx

	  N_shelly=(int)(pvy_dash[i]/del_vy+offset_y);
	  Nof_vy_inst[position_index][N_shelly]=Nof_vy_inst[position_index][N_shelly]+1;
	  Nof_vy[position_index][N_shelly]= Nof_vy[position_index][N_shelly]+1;
             
	  N_shellz=(int)(pvz_dash[i]/del_vz+offset_z);
	  Nof_vz_inst[position_index][N_shellz]=Nof_vz_inst[position_index][N_shellz]+1;
	  Nof_vz[position_index][N_shellz]= Nof_vz[position_index][N_shellz]+1;
	  
	  atom_count_inst[position_index]=atom_count_inst[position_index]+1;
	  tot_atom_count[position_index]=tot_atom_count[position_index]+1;
   }
             }
	   }

	  for(position_index=1;position_index<=max_position_index;position_index++) 
	    	  for( kk=1;kk<=Ndel_vx;kk++)           // k should start from 1 not from 0  
                 {
		   if(atom_count_inst[position_index]>0)	  
    fof_vx_sq[position_index][kk]=fof_vx_sq[position_index][kk]+(double)(Nof_vx_inst[position_index][kk]/(atom_count_inst[position_index]*del_vx))*
(double)(Nof_vx_inst[position_index][kk]/(atom_count_inst[position_index]*del_vx));
		 }

 for(position_index=1;position_index<=max_position_index;position_index++) 
	   for( kk=1;kk<=Ndel_vy;kk++)           // k should start from 1 not from 0  
                 {
	if(atom_count_inst[position_index]>0)	   
   fof_vy_sq[position_index][kk]=fof_vy_sq[position_index][kk]+(double)( Nof_vy_inst[position_index][kk]/(atom_count_inst[position_index]*del_vy))*
(double)( Nof_vy_inst[position_index][kk]/(atom_count_inst[position_index]*del_vy));
		 }

 for(position_index=1;position_index<=max_position_index;position_index++) 
	   for( kk=1;kk<=Ndel_vz;kk++)           // k should start from 1 not from 0  
                 {
		if(atom_count_inst[position_index]>0)   
   fof_vz_sq[position_index][kk]=fof_vz_sq[position_index][kk]+ (double)( Nof_vz_inst[position_index][kk]/(atom_count_inst[position_index]*del_vz))*
(double)( Nof_vz_inst[position_index][kk]/(atom_count_inst[position_index]*del_vz));
		 }



 /*------ Calculation of the varience of the stress-------------------*/


 if(index_stress_varience==1)  // this will give the idea of the difference of particle average to the particle+time average.
                               // but the error bar is calculated in different way, saving 20 ensemble of 1000 sample each
  {

 for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
   
    pvx_dash_sq[ny]=0.0;
    pvy_dash_sq[ny]=0.0;
    pvz_dash_sq[ny]=0.0;
    pvxy_dash[ny]=0.0;

        }

 for(i=0;i<n_atom;i++)
   {
 ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
 pvx_dash[i]=p_vx[i]-avg_vx_part[ny];  // calculating fluctuation for each particle
 pvy_dash[i]=p_vy[i]-avg_vy_part[ny];
 pvz_dash[i]=p_vz[i]-avg_vz_part[ny];
 
 pvx_dash_sq[ny]=pvx_dash_sq[ny]+pvx_dash[i]*pvx_dash[i];
 pvy_dash_sq[ny]=pvy_dash_sq[ny]+pvy_dash[i]*pvy_dash[i];
 pvz_dash_sq[ny]=pvz_dash_sq[ny]+pvz_dash[i]*pvz_dash[i];
 pvxy_dash[ny]=pvxy_dash[ny]+pvx_dash[i]*pvy_dash[i];
   }

 for(ny=0;ny<Ndel_y;ny++) 
  {

    if(Nofp_inst[ny]>0.0)
      {
 pvx_dash_sq[ny]= pvx_dash_sq[ny]/Nofp_inst[ny];// doing the particle average i
 pvy_dash_sq[ny]= pvy_dash_sq[ny]/Nofp_inst[ny];
 pvz_dash_sq[ny]= pvz_dash_sq[ny]/Nofp_inst[ny];
 pvxy_dash[ny]  = pvxy_dash[ny]/Nofp_inst[ny];
      }
    
  }

 for(ny=0;ny<Ndel_y;ny++) 
  {
    var_sigma_xx[ny]= var_sigma_xx[ny]+ (pvx_dash_sq[ny]-avg_stress_x[ny])*(pvx_dash_sq[ny]-avg_stress_x[ny]);
    var_sigma_yy[ny]= var_sigma_yy[ny]+ (pvy_dash_sq[ny]-avg_stress_y[ny])*(pvy_dash_sq[ny]-avg_stress_y[ny]);
    var_sigma_zz[ny]= var_sigma_zz[ny]+ (pvz_dash_sq[ny]-avg_stress_z[ny])*(pvz_dash_sq[ny]-avg_stress_z[ny]);
    var_tau_xy[ny]= var_tau_xy[ny]    + (pvxy_dash[ny]-avg_stress_xy[ny])* ( pvxy_dash[ny]-avg_stress_xy[ny]);

  }
									    

  }

 /*-------------------------------------------------------------------------------------    */
 /*---- Calculation of the fourth moment of the particle velocity fluctuation-by loading the avg. particle vel.
//------*/ 
 /*---------------------------------------------------------------------------------------*/

 if(index_fourth_moment==1)
   {
      for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
   
    pvx_dash_4[ny]=0.0;
    pvy_dash_4[ny]=0.0;
    pvz_dash_4[ny]=0.0;
    pvxy_dash_4[ny]=0.0;

        }
      
      for(i=0;i<n_atom;i++)
	{
	  ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
	  
	  pvx_dash[i]=p_vx[i]-avg_vx_part[ny];  // calculating fluctuation for each particle
	  pvy_dash[i]=p_vy[i]-avg_vy_part[ny];
	  pvz_dash[i]=p_vz[i]-avg_vz_part[ny];

	  pvx_dash_4[ny]=pvx_dash_4[ny]+pvx_dash[i]*pvx_dash[i]*pvx_dash[i]*pvx_dash[i];
	  pvy_dash_4[ny]=pvy_dash_4[ny]+pvy_dash[i]*pvy_dash[i]*pvy_dash[i]*pvy_dash[i];
	  pvz_dash_4[ny]=pvz_dash_4[ny]+pvz_dash[i]*pvz_dash[i]*pvz_dash[i]*pvz_dash[i];
	  pvxy_dash_4[ny]= pvxy_dash_4[ny]+pvx_dash[i]*pvx_dash[i]*pvy_dash[i]*pvy_dash[i];
	     }

      for(ny=0;ny<Ndel_y;ny++) 
	{
        if(Nofp_inst[ny]>0.0)
	  {
	    pvx_dash_4[ny]= pvx_dash_4[ny]/Nofp_inst[ny];// doing the particle average i
	    pvy_dash_4[ny]= pvy_dash_4[ny]/Nofp_inst[ny];
	    pvz_dash_4[ny]= pvz_dash_4[ny]/Nofp_inst[ny];
	    pvxy_dash_4[ny]  = pvxy_dash_4[ny]/Nofp_inst[ny];
	  }

	avg_pvx_dash_4[ny]=avg_pvx_dash_4[ny]+pvx_dash_4[ny];
	avg_pvy_dash_4[ny]=avg_pvy_dash_4[ny]+pvy_dash_4[ny];
	avg_pvz_dash_4[ny]=avg_pvz_dash_4[ny]+pvz_dash_4[ny];
	avg_pvxy_dash_4[ny]=avg_pvxy_dash_4[ny]+pvxy_dash_4[ny];
	}

   }




     

	  /*		
			     
 fp3=fopen("/home/psg/dns_out/channel/with_pp_col/volfrac_0.0004/run_density_rough_2000/all_data_1.txt","a");
 fprintf(fp3," sample_No= %d\n",n_sample);

        for(ny=0;ny<Ndel_y;ny++)             
 
	fprintf(fp3,"%8d %20.14lf %20.14lf %20.14lf  %20d  \n",
 ny,vx_inst[ny] , vy_inst[ny],vz_inst[ny],Nofp_inst[ny]);
	 fclose(fp3);
	  */	    
             return;

     }



//--------------------------------------------------------------
//     Functional properties of angular velocity distribution
//--------------------------------------------------------------
void angular_property(){

      int i,ny,kk,jj;
       int N_shell,N_shellx, N_shelly,N_shellz;
       
       // int k,N_rdels,

       //double  TELASP,BELASP, VIRIAL;
       
      
       double ave_KE;
       double  offset,offset_x, offset_y , offset_z ;
       
       //double fof_v, fatom,HH, HINST;


       // double pi = 3.141592650;
       //FILE *fp2;

       //FILE *fp3;
 
                        
	 sampls=sampls+1;
	 n_sample=n_sample+1;
	 // Accumulate running average kinetic energy
                        
  sum_KE = sum_KE + tot_KE;
  ave_KE = sum_KE/(double)(n_sample);
int j;
j=n_sample;

    for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous  angular velocities
       {
    Nofp_inst[ny]=0;
    angvelx_inst[ny]=0.0;
    angvely_inst[ny]=0.0;
    angvelz_inst[ny]=0.0;
    
    pangvelx_dash_sq[ny]=0.0;
    pangvely_dash_sq[ny]=0.0;
    pangvelz_dash_sq[ny]=0.0;
    pangvelxy_dash[ny]=0.0;

    angvelx_angvelx_inst[ny]=0.0;
    angvely_angvely_inst[ny]=0.0;
    angvelz_angvelz_inst[ny]=0.0;
    angvelx_angvely_inst[ny]=0.0;

    
           }

    for(i=0;i<n_atom;i++)// summing up all the angvelocities in each the y-grid for all the particle 
    {
      ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

      // cout<<"in property"<<p_y[i]<<" "<<ny<<endl;


              Nofp_inst[ny]= Nofp_inst[ny]+ 1;

              angvelx_inst[ny]= angvelx_inst[ny]+p_angvelx[i];
              angvely_inst[ny]= angvely_inst[ny]+p_angvely[i];  
              angvelz_inst[ny]= angvelz_inst[ny]+p_angvelz[i];
      }

    for(ny=0;ny<Ndel_y;ny++)                // for calculating density fluctuation(error bar)
      {
	Nofp_inst_sq[ny]= Nofp_inst_sq[ny]+Nofp_inst[ny]*Nofp_inst[ny];
      }


// --------------newly added for stress callculation by <u^2>-Uavg^2 folmulea

 for(i=0;i<n_atom;i++)// summing up all the velocities in each the y-grid for all the particle 
    {
      ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

      // cout<<"in property"<<p_y[i]<<" "<<ny<<endl;


      //       Nofp_inst[ny]= Nofp_inst[ny]+ 1;

              angvelx_angvelx_inst[ny]= angvelx_angvelx_inst[ny]+p_angvelx[i]*p_angvelx[i];
              angvely_angvely_inst[ny]= angvely_angvely_inst[ny]+p_angvely[i]*p_angvely[i];  
              angvelz_angvelz_inst[ny]= angvelz_angvelz_inst[ny]+p_angvelz[i]*p_angvelz[i];
	      angvelx_angvely_inst[ny]= angvelx_angvely_inst[ny]+p_angvelx[i]*p_angvely[i];  


      }

 for(ny=0;ny<Ndel_y;ny++) // doing the particle average in each y-grid-in each sample frame
       {
    if(Nofp_inst[ny]>0.0)
      {
 angvelx_angvelx_inst[ny]= angvelx_angvelx_inst[ny]/Nofp_inst[ny];
 angvely_angvely_inst[ny]= angvely_angvely_inst[ny]/Nofp_inst[ny];
 angvelz_angvelz_inst[ny]= angvelz_angvelz_inst[ny]/Nofp_inst[ny];
 angvelx_angvely_inst[ny]= angvelx_angvely_inst[ny]/Nofp_inst[ny];
      }

    angsigma_xx_1[ny]=angsigma_xx_1[ny]+ angvelx_angvelx_inst[ny];
    angsigma_yy_1[ny]=angsigma_yy_1[ny]+ angvely_angvely_inst[ny];
    angsigma_zz_1[ny]=angsigma_zz_1[ny]+ angvelz_angvelz_inst[ny];
    angsigma_xy_1[ny]=angsigma_xy_1[ny]+ angvelx_angvely_inst[ny];

       }

 //------------------back to the old stress calculation by calculation fluc over particle avg----------

  for(ny=0;ny<Ndel_y;ny++) // doing the particle average in each y-grid-in each sample frame
       {
    if(Nofp_inst[ny]>0.0)
      {
     angvelx_inst[ny]= angvelx_inst[ny]/Nofp_inst[ny];  // vx_inst[ny] is basically particle averaged inst in time.
     angvely_inst[ny]= angvely_inst[ny]/Nofp_inst[ny];
     angvelz_inst[ny]= angvelz_inst[ny]/Nofp_inst[ny];
      }
     angvelx_mean[ny]=angvelx_mean[ny]+angvelx_inst[ny];    // data storing for the time aangveleraging but take care 
     angvely_mean[ny]=angvely_mean[ny]+angvely_inst[ny];   //---- that this quantities are not average
     angvelz_mean[ny]=angvelz_mean[ny]+angvelz_inst[ny];
     Nofp_mean[ny]=Nofp_mean[ny]+Nofp_inst[ny];
    
       }


 for(i=0;i<n_atom;i++)
   {
 ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
 pangvelx_dash[i]=p_angvelx[i]-angvelx_inst[ny];  // calculating fluctuation for each particle
 pangvely_dash[i]=p_angvely[i]-angvely_inst[ny];
 pangvelz_dash[i]=p_angvelz[i]-angvelz_inst[ny];
 pangvel_dash[i]=sqrt(p_angvelx[i]*p_angvelx[i]+p_angvely[i]*p_angvely[i]+p_angvelz[i]*p_angvelz[i])
             -sqrt(angvelx_inst[ny]*angvelx_inst[ny]+angvely_inst[ny]*angvely_inst[ny]+angvelz_inst[ny]*angvelz_inst[ny]);

 pangvelx_dash_sq[ny]=pangvelx_dash_sq[ny]+pangvelx_dash[i]*pangvelx_dash[i];
 pangvely_dash_sq[ny]=pangvely_dash_sq[ny]+pangvely_dash[i]*pangvely_dash[i];
 pangvelz_dash_sq[ny]=pangvelz_dash_sq[ny]+pangvelz_dash[i]*pangvelz_dash[i];
 pangvelxy_dash[ny]=pangvelxy_dash[ny]+pangvelx_dash[i]*pangvely_dash[i];
   }

 for(ny=0;ny<Ndel_y;ny++) 
  {
    if(Nofp_inst[ny]>0.0)
      {
 pangvelx_dash_sq[ny]= pangvelx_dash_sq[ny]/Nofp_inst[ny];// doing the particle average i
 pangvely_dash_sq[ny]= pangvely_dash_sq[ny]/Nofp_inst[ny];
 pangvelz_dash_sq[ny]= pangvelz_dash_sq[ny]/Nofp_inst[ny];
 pangvelxy_dash[ny]  = pangvelxy_dash[ny]/Nofp_inst[ny];
      }
    sigma_xx[ny]=sigma_xx[ny]+ pangvelx_dash_sq[ny]; // data storing for the time aangveleraging
    sigma_yy[ny]=sigma_yy[ny]+ pangvely_dash_sq[ny];
    sigma_zz[ny]=sigma_zz[ny]+ pangvelz_dash_sq[ny];
    tau_xy[ny]=tau_xy[ny]+ pangvelxy_dash[ny];

    }


 /*-----------------------------------------------------------------------
Calculation of the third and the fourth moment, considering fluctuations over particle
average, because it is found to give the same result if we calculate the fluctuation over 
particle + time average.... that routine is also implemented the programme
------------------------------------------------------------------------------------*/

 for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
   
    pangvelx_dash_4[ny]=0.0;
    pangvely_dash_4[ny]=0.0;
    pangvelz_dash_4[ny]=0.0;
    pangvelxy_dash_4[ny]=0.0;

    pangvelx_dash_3[ny]=0.0;
    pangvely_dash_3[ny]=0.0;
    pangvelz_dash_3[ny]=0.0;
   
        }
      
      for(i=0;i<n_atom;i++)
	{
	  ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
	  
	  pangvelx_dash_4[ny]=pangvelx_dash_4[ny]+pangvelx_dash[i]*pangvelx_dash[i]*pangvelx_dash[i]*pangvelx_dash[i];
	  pangvely_dash_4[ny]=pangvely_dash_4[ny]+pangvely_dash[i]*pangvely_dash[i]*pangvely_dash[i]*pangvely_dash[i];
	  pangvelz_dash_4[ny]=pangvelz_dash_4[ny]+pangvelz_dash[i]*pangvelz_dash[i]*pangvelz_dash[i]*pangvelz_dash[i];
	
	  pangvelx_dash_3[ny]=pangvelx_dash_3[ny]+pangvelx_dash[i]*pangvelx_dash[i]*pangvelx_dash[i];
	  pangvely_dash_3[ny]=pangvely_dash_3[ny]+pangvely_dash[i]*pangvely_dash[i]*pangvely_dash[i];
	  pangvelz_dash_3[ny]=pangvelz_dash_3[ny]+pangvelz_dash[i]*pangvelz_dash[i]*pangvelz_dash[i];

	     }

      for(ny=0;ny<Ndel_y;ny++) 
	{
        if(Nofp_inst[ny]>0.0)
	  {
	    pangvelx_dash_4[ny]= pangvelx_dash_4[ny]/Nofp_inst[ny];// doing the particle average i
	    pangvely_dash_4[ny]= pangvely_dash_4[ny]/Nofp_inst[ny];
	    pangvelz_dash_4[ny]= pangvelz_dash_4[ny]/Nofp_inst[ny];
	    
	    pangvelx_dash_3[ny]= pangvelx_dash_3[ny]/Nofp_inst[ny];// doing the particle average i
	    pangvely_dash_3[ny]= pangvely_dash_3[ny]/Nofp_inst[ny];
	    pangvelz_dash_3[ny]= pangvelz_dash_3[ny]/Nofp_inst[ny];

	  }

	avg_pangvelx_dash_4[ny]=avg_pangvelx_dash_4[ny]+pangvelx_dash_4[ny];
	avg_pangvely_dash_4[ny]=avg_pangvely_dash_4[ny]+pangvely_dash_4[ny];
	avg_pangvelz_dash_4[ny]=avg_pangvelz_dash_4[ny]+pangvelz_dash_4[ny];

	avg_pangvelx_dash_3[ny]=avg_pangvelx_dash_3[ny]+pangvelx_dash_3[ny];
	avg_pangvely_dash_3[ny]=avg_pangvely_dash_3[ny]+pangvely_dash_3[ny];
	avg_pangvelz_dash_3[ny]=avg_pangvelz_dash_3[ny]+pangvelz_dash_3[ny];
	

	}

//------------------------------------------------------------------------
      //   Storing counter for  Angular Velocity distribution  calculation
      //--------------------------------------------------------------------------
               for(jj=0;jj<Max_dim_position;jj++)
		  for(i=0;i<Ndel_angvel;i++)
		    Nof_angvel_inst[jj][i]=0;

	       for(jj=0;jj<Max_dim_position;jj++)
		  for(i=0;i<Ndel_angvelx;i++)
		    Nof_angvelx_inst[jj][i]=0;
      
	       for(jj=0;jj<Max_dim_position;jj++)
		  for(i=0;i<Ndel_angvely;i++)
		    Nof_angvely_inst[jj][i]=0;

	       for(jj=0;jj<Max_dim_position;jj++)
		  for(i=0;i<Ndel_angvelz;i++)
		    Nof_angvelz_inst[jj][i]=0;

	       for(jj=0;jj<Max_dim_position;jj++)
		  atom_count_inst[jj]=0;
           
            offset =1.5+p_angvel_max/del_angvel;
            offset_x =1.5+p_angvelx_max/del_angvelx;
            offset_y =1.5+p_angvely_max/del_angvely;
            offset_z =1.5+p_angvelz_max/del_angvelz;


	    
          for(i=0;i<n_atom;i++)
           {
	     position_index=0;    //initialization of position_index

	     if(p_y[i]>y_min_for_distribution_1 && p_y[i]<y_max_for_distribution_1)position_index=1;
	     if(p_y[i]>y_min_for_distribution_2 && p_y[i]<y_max_for_distribution_2)position_index=2;
	     if(p_y[i]>y_min_for_distribution_3 && p_y[i]<y_max_for_distribution_3)position_index=3;
	     

	     if(position_index>0)
	       {                                                
	       
	if((pangvel_dash[i]>= -p_angvel_max && pangvel_dash[i]<=p_angvel_max) &&(pangvelx_dash[i]>=-p_angvelx_max && pangvelx_dash[i]<=p_angvelx_max)&&
           (pangvely_dash[i]>=-p_angvely_max && pangvely_dash[i]<=p_angvely_max)&& (pangvelz_dash[i]>=-p_angvelz_max && pangvelz_dash[i]<=p_angvelz_max))  
				                // to avoid exceeding arrey subscript accidentally
				 {
     // p_angvel[i]= sqrt(p_angvelx[i]*p_angvelx[i]+p_angvely[i]*p_angvely[i]+p_angvelz[i]*p_angvelz[i]);

	        N_shell=(int)(pangvel_dash[i]/del_angvel+offset);
		Nof_angvel[position_index][N_shell]= Nof_angvel[position_index][N_shell]+1;
		Nof_angvel_inst[position_index][N_shell]= Nof_angvel_inst[position_index][N_shell]+1;

             N_shellx=(int)(pangvelx_dash[i]/del_angvelx+offset_x);
             Nof_angvelx_inst[position_index][N_shellx]= Nof_angvelx_inst[position_index][N_shellx]+1;
	     Nof_angvelx[position_index][N_shellx]= Nof_angvelx[position_index][N_shellx] +1 ; //!Nofangvelx_s SEPERATE X COMPONENT OF Nof_angvelx

            	                                                     //--COMPONENT OF Nof_angvelx

	  N_shelly=(int)(pangvely_dash[i]/del_angvely+offset_y);
	  Nof_angvely_inst[position_index][N_shelly]=Nof_angvely_inst[position_index][N_shelly]+1;
	  Nof_angvely[position_index][N_shelly]= Nof_angvely[position_index][N_shelly]+1;
             
	  N_shellz=(int)(pangvelz_dash[i]/del_angvelz+offset_z);
	  Nof_angvelz_inst[position_index][N_shellz]=Nof_angvelz_inst[position_index][N_shellz]+1;
	  Nof_angvelz[position_index][N_shellz]= Nof_angvelz[position_index][N_shellz]+1;
	  
	  atom_count_inst[position_index]=atom_count_inst[position_index]+1;
	  tot_atom_count[position_index]=tot_atom_count[position_index]+1;
   }
             }
	   }

	  for(position_index=1;position_index<=max_position_index;position_index++) 
	    	  for( kk=1;kk<=Ndel_angvelx;kk++)           // k should start from 1 not from 0  
                 {
		   if(atom_count_inst[position_index]>0)	  
    fof_angvelx_sq[position_index][kk]=fof_angvelx_sq[position_index][kk]+(double)(Nof_angvelx_inst[position_index][kk]/(atom_count_inst[position_index]*del_angvelx))*
(double)(Nof_angvelx_inst[position_index][kk]/(atom_count_inst[position_index]*del_angvelx));
		 }

 for(position_index=1;position_index<=max_position_index;position_index++) 
	   for( kk=1;kk<=Ndel_angvely;kk++)           // k should start from 1 not from 0  
                 {
	if(atom_count_inst[position_index]>0)	   
   fof_angvely_sq[position_index][kk]=fof_angvely_sq[position_index][kk]+(double)( Nof_angvely_inst[position_index][kk]/(atom_count_inst[position_index]*del_angvely))*
(double)( Nof_angvely_inst[position_index][kk]/(atom_count_inst[position_index]*del_angvely));
		 }

 for(position_index=1;position_index<=max_position_index;position_index++) 
	   for( kk=1;kk<=Ndel_angvelz;kk++)           // k should start from 1 not from 0  
                 {
		if(atom_count_inst[position_index]>0)   
   fof_angvelz_sq[position_index][kk]=fof_angvelz_sq[position_index][kk]+ (double)( Nof_angvelz_inst[position_index][kk]/(atom_count_inst[position_index]*del_angvelz))*
(double)( Nof_angvelz_inst[position_index][kk]/(atom_count_inst[position_index]*del_angvelz));
		 }



 /*------ Calculation of the varience of the stress-------------------*/


 if(index_stress_varience==1)  // this will give the idea of the difference of particle average to the particle+time average.
                               // but the error bar is calculated in different way, saving 20 ensemble of 1000 sample each
  {

 for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
   
    pangvelx_dash_sq[ny]=0.0;
    pangvely_dash_sq[ny]=0.0;
    pangvelz_dash_sq[ny]=0.0;
    pangvelxy_dash[ny]=0.0;

        }

 for(i=0;i<n_atom;i++)
   {
 ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
 pangvelx_dash[i]=p_angvelx[i]-avg_angvelx_part[ny];  // calculating fluctuation for each particle
 pangvely_dash[i]=p_angvely[i]-avg_angvely_part[ny];
 pangvelz_dash[i]=p_angvelz[i]-avg_angvelz_part[ny];
 
 pangvelx_dash_sq[ny]=pangvelx_dash_sq[ny]+pangvelx_dash[i]*pangvelx_dash[i];
 pangvely_dash_sq[ny]=pangvely_dash_sq[ny]+pangvely_dash[i]*pangvely_dash[i];
 pangvelz_dash_sq[ny]=pangvelz_dash_sq[ny]+pangvelz_dash[i]*pangvelz_dash[i];
 pangvelxy_dash[ny]=pangvelxy_dash[ny]+pangvelx_dash[i]*pangvely_dash[i];
   }

 for(ny=0;ny<Ndel_y;ny++) 
  {

    if(Nofp_inst[ny]>0.0)
      {
 pangvelx_dash_sq[ny]= pangvelx_dash_sq[ny]/Nofp_inst[ny];// doing the particle average i
 pangvely_dash_sq[ny]= pangvely_dash_sq[ny]/Nofp_inst[ny];
 pangvelz_dash_sq[ny]= pangvelz_dash_sq[ny]/Nofp_inst[ny];
 pangvelxy_dash[ny]  = pangvelxy_dash[ny]/Nofp_inst[ny];
      }
    
  }

 for(ny=0;ny<Ndel_y;ny++) 
  {
    var_angsigma_xx[ny]= var_angsigma_xx[ny]+ (pangvelx_dash_sq[ny]-avg_angstress_x[ny])*(pangvelx_dash_sq[ny]-avg_angstress_x[ny]);
    var_angsigma_yy[ny]= var_angsigma_yy[ny]+ (pangvely_dash_sq[ny]-avg_angstress_y[ny])*(pangvely_dash_sq[ny]-avg_angstress_y[ny]);
    var_angsigma_zz[ny]= var_angsigma_zz[ny]+ (pangvelz_dash_sq[ny]-avg_angstress_z[ny])*(pangvelz_dash_sq[ny]-avg_angstress_z[ny]);
    var_angtau_xy[ny]= var_angtau_xy[ny]    + (pangvelxy_dash[ny]-avg_angstress_xy[ny])* ( pangvelxy_dash[ny]-avg_angstress_xy[ny]);

  }
									    

  }

 /*-------------------------------------------------------------------------------------    */
 /*---- Calculation of the fourth moment of the particle velocity fluctuation-by loading the avg. particle vel.
//------*/ 
 /*---------------------------------------------------------------------------------------*/

 if(index_fourth_moment==1)
   {
      for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
   
    pangvelx_dash_4[ny]=0.0;
    pangvely_dash_4[ny]=0.0;
    pangvelz_dash_4[ny]=0.0;
    pangvelxy_dash_4[ny]=0.0;

        }
      
      for(i=0;i<n_atom;i++)
	{
	  ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
	  
	  pangvelx_dash[i]=p_angvelx[i]-avg_angvelx_part[ny];  // calculating fluctuation for each particle
	  pangvely_dash[i]=p_angvely[i]-avg_angvely_part[ny];
	  pangvelz_dash[i]=p_angvelz[i]-avg_angvelz_part[ny];

	  pangvelx_dash_4[ny]=pangvelx_dash_4[ny]+pangvelx_dash[i]*pangvelx_dash[i]*pangvelx_dash[i]*pangvelx_dash[i];
	  pangvely_dash_4[ny]=pangvely_dash_4[ny]+pangvely_dash[i]*pangvely_dash[i]*pangvely_dash[i]*pangvely_dash[i];
	  pangvelz_dash_4[ny]=pangvelz_dash_4[ny]+pangvelz_dash[i]*pangvelz_dash[i]*pangvelz_dash[i]*pangvelz_dash[i];
	  pangvelxy_dash_4[ny]= pangvelxy_dash_4[ny]+pangvelx_dash[i]*pangvelx_dash[i]*pangvely_dash[i]*pangvely_dash[i];
	     }

      for(ny=0;ny<Ndel_y;ny++) 
	{
        if(Nofp_inst[ny]>0.0)
	  {
	    pangvelx_dash_4[ny]= pangvelx_dash_4[ny]/Nofp_inst[ny];// doing the particle average i
	    pangvely_dash_4[ny]= pangvely_dash_4[ny]/Nofp_inst[ny];
	    pangvelz_dash_4[ny]= pangvelz_dash_4[ny]/Nofp_inst[ny];
	    pangvelxy_dash_4[ny]  = pangvelxy_dash_4[ny]/Nofp_inst[ny];
	  }

	avg_pangvelx_dash_4[ny]=avg_pangvelx_dash_4[ny]+pangvelx_dash_4[ny];
	avg_pangvely_dash_4[ny]=avg_pangvely_dash_4[ny]+pangvely_dash_4[ny];
	avg_pangvelz_dash_4[ny]=avg_pangvelz_dash_4[ny]+pangvelz_dash_4[ny];
	avg_pangvelxy_dash_4[ny]=avg_pangvelxy_dash_4[ny]+pangvelxy_dash_4[ny];
	}

   }




     

	  /*		
			     
 fp3=fopen("/home/psg/dns_out/channel/with_pp_col/volfrac_0.0004/run_density_rough_2000/all_data_1.txt","a");
 fprintf(fp3," sample_No= %d\n",n_sample);

        for(ny=0;ny<Ndel_y;ny++)             
 
	fprintf(fp3,"%8d %20.14lf %20.14lf %20.14lf  %20d  \n",
 ny,angvelx_inst[ny] , angvely_inst[ny],angvelz_inst[ny],Nofp_inst[ny]);
	 fclose(fp3);
	  */	    
             return;

     }

//-------------------------------------------------------------
// Average and fluctuation calculation of different properties
//-------------------------------------------------------------
void final_property()
{

 int ny;
 double  y_position;
 double  Var_Nofp_[Max_y_grid];
 double symmetric_norm[Max_y_grid], symmetric_norm_corr[Max_y_grid];

   
  FILE *fp4;
  FILE *fp5;
  FILE *fp51;
  FILE *fp52;

sampling_time=run_time-t_zero;

  
//------ for grid dependent property calculation----------


 for(ny=0;ny<Ndel_y;ny++) // loop for doing time averaging
   { 
     Nofp_mean[ny]=Nofp_mean[ny]/(double)n_sample; 

     vx_mean[ny]=vx_mean[ny]/(double)n_sample;    
     vy_mean[ny]=vy_mean[ny]/(double)n_sample;  
     vz_mean[ny]=vz_mean[ny]/(double)n_sample;
     
     symmetric_p_xx_mean[ny] = symmetric_p_xx_mean[ny]/(double)n_sample;
     symmetric_p_xy_mean[ny] = symmetric_p_xy_mean[ny]/(double)n_sample;
     symmetric_p_xz_mean[ny] = symmetric_p_xz_mean[ny]/(double)n_sample;
     symmetric_p_yy_mean[ny] = symmetric_p_yy_mean[ny]/(double)n_sample;
     symmetric_p_yz_mean[ny] = symmetric_p_yz_mean[ny]/(double)n_sample;
     symmetric_p_zz_mean[ny] = symmetric_p_zz_mean[ny]/(double)n_sample;

     symmetric_norm[ny] = sqrt(symmetric_p_xx_mean[ny]*symmetric_p_xx_mean[ny] + 2*symmetric_p_xy_mean[ny]*symmetric_p_xy_mean[ny] + 2*symmetric_p_xz_mean[ny]*symmetric_p_xz_mean[ny] + symmetric_p_yy_mean[ny]*symmetric_p_yy_mean[ny] + 2*symmetric_p_yz_mean[ny]*symmetric_p_yz_mean[ny] + symmetric_p_yy_mean[ny]*symmetric_p_yy_mean[ny]);

     symmetric_particle_corr_xx_mean[ny] = Nofp_mean[ny]*symmetric_p_xx_mean[ny];
     symmetric_particle_corr_xy_mean[ny] = Nofp_mean[ny]*symmetric_p_xy_mean[ny];
     symmetric_particle_corr_xz_mean[ny] = Nofp_mean[ny]*symmetric_p_xz_mean[ny];
     symmetric_particle_corr_yy_mean[ny] = Nofp_mean[ny]*symmetric_p_yy_mean[ny];
     symmetric_particle_corr_yz_mean[ny] = Nofp_mean[ny]*symmetric_p_yz_mean[ny];
     symmetric_particle_corr_zz_mean[ny] = Nofp_mean[ny]*symmetric_p_zz_mean[ny];

     symmetric_norm_corr[ny] = Nofp_mean[ny]*symmetric_norm[ny];

			 sigma_xx[ny]=sigma_xx[ny]/(double)n_sample;
			 sigma_yy[ny]=sigma_yy[ny]/(double)n_sample;
			 sigma_zz[ny]=sigma_zz[ny]/(double)n_sample;

			 tau_xy[ny]=tau_xy[ny]/(double)n_sample;

			 sigma_xx_1[ny]=sigma_xx_1[ny]/(double)n_sample- vx_mean[ny]* vx_mean[ny];
			 sigma_yy_1[ny]=sigma_yy_1[ny]/(double)n_sample- vy_mean[ny]* vy_mean[ny];
			 sigma_zz_1[ny]=sigma_zz_1[ny]/(double)n_sample- vz_mean[ny]* vz_mean[ny];
			 sigma_xy_1[ny]=sigma_xy_1[ny]/(double)n_sample- vx_mean[ny]* vy_mean[ny];

			 Var_Nofp_[ny]=Nofp_inst_sq[ny]/(double)n_sample- Nofp_mean[ny]* Nofp_mean[ny];
			 
			 }

 if(index_stress_varience==1)
   {
      for(ny=0;ny<Ndel_y;ny++) // loop for doing time averaging
	{ 
	   var_sigma_xx[ny]=var_sigma_xx[ny]/(double)n_sample;
	   var_sigma_yy[ny]=var_sigma_yy[ny]/(double)n_sample;
	   var_sigma_zz[ny]=var_sigma_zz[ny]/(double)n_sample;
 	   var_tau_xy[ny]=var_tau_xy[ny]/(double)n_sample;
	}

   }

  if(index_fourth_moment==1)
    {
      for(ny=0;ny<Ndel_y;ny++)          // loop for doing time averaging
       {
	 avg_pvx_dash_4[ny]=avg_pvx_dash_4[ny]/(double)n_sample;
	 avg_pvy_dash_4[ny]=avg_pvy_dash_4[ny]/(double)n_sample;
	 avg_pvz_dash_4[ny]=avg_pvz_dash_4[ny]/(double)n_sample;
	 


	 avg_pvx_dash_3[ny]=avg_pvx_dash_3[ny]/(double)n_sample;
	 avg_pvy_dash_3[ny]=avg_pvy_dash_3[ny]/(double)n_sample;
	 avg_pvz_dash_3[ny]=avg_pvz_dash_3[ny]/(double)n_sample;
      }
   
    }

for(ny=0;ny<Ndel_y;ny++) 
    {

    accx_mean[ny]=accx_mean[ny]/(double)n_sample;
    accy_mean[ny]=accy_mean[ny]/(double)n_sample;
    accz_mean[ny]=accz_mean[ny]/(double)n_sample;

    mean_accx_dash_2[ny]= mean_accx_sq[ny]/(double)n_sample-accx_mean[ny]*accx_mean[ny];
    mean_accy_dash_2[ny]= mean_accy_sq[ny]/(double)n_sample-accy_mean[ny]*accy_mean[ny];
    mean_accz_dash_2[ny]= mean_accz_sq[ny]/(double)n_sample-accz_mean[ny]*accz_mean[ny];
    }



  for(ny=0;ny<Ndel_y;ny++) 
    {
      avg_accx_dash_2[ny]=avg_accx_dash_2[ny]/(double)n_sample;
      avg_accy_dash_2[ny]=avg_accy_dash_2[ny]/(double)n_sample;
      avg_accz_dash_2[ny]=avg_accz_dash_2[ny]/(double)n_sample;
      

      avg_accx_dash_3[ny]=avg_accx_dash_3[ny]/(double)n_sample;
      avg_accy_dash_3[ny]=avg_accy_dash_3[ny]/(double)n_sample;
      avg_accz_dash_3[ny]=avg_accz_dash_3[ny]/(double)n_sample;
      
      avg_accx_dash_4[ny]=avg_accx_dash_4[ny]/(double)n_sample;
      avg_accy_dash_4[ny]=avg_accy_dash_4[ny]/(double)n_sample;
      avg_accz_dash_4[ny]=avg_accz_dash_4[ny]/(double)n_sample;
      
    }



fp4=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_2.txt","w");
fprintf(fp4," %s Sampling duration =%20.14lf\n","%", sampling_time);
fprintf(fp4," %s Grid_index  Grid point Vx    Vy    Vz   No_of_Part  var_no_of_P   No of Collision \n ","%"); 
  

 for(ny=0;ny<Ndel_y;ny++)
   {
   y_position=(ny+0.5)*del_y+(0.5*sigma);

fprintf(fp4,"%8d %20.14lf %20.14lf %20.14lf %20.14lf  %20.14lf %20.14lf %15d\n",
ny, y_position,vx_mean[ny] , vy_mean[ny],vz_mean[ny],Nofp_mean[ny], Var_Nofp_[ny], Nof_p_p_coll[ny]);
 
   }

fclose(fp4);

fp51=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/symmetric_particle_2.txt","w");
fprintf(fp51,"%s Grid_index     Grid_point       Symm_xx    Symm_xy    Symm_xz   Symm_yy   Symm_yz   Symm_zz  Symm_norm  No_of_Part  \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

fprintf(fp51,"%8d %20.14lf %20.14lf %20.14lf %20.14lf  %20.14lf %20.14lf %20.14lf %20.14lf %20.14lf\n",
ny, y_position,symmetric_p_xx_mean[ny] ,symmetric_p_xy_mean[ny], symmetric_p_xz_mean[ny],symmetric_p_yy_mean[ny],symmetric_p_yz_mean[ny],symmetric_p_zz_mean[ny],symmetric_norm[ny],
Nofp_mean[ny]);
   }
 
fclose(fp51);

fp52=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/symmetric_particle_corr_2.txt","w");
fprintf(fp52,"%s Grid_index     Grid_point       Symm_corr_xx    Symm_corr_xy    Symm_corr_xz   Symm_corr_yy   Symm_corr_yz   Symm_corr_zz   Symm_norm_corr  \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

fprintf(fp52,"%8d %20.14lf %20.14lf %20.14lf %20.14lf  %20.14lf %20.14lf %20.14lf %20.14lf \n",
ny, y_position,symmetric_particle_corr_xx_mean[ny] ,symmetric_particle_corr_xy_mean[ny], symmetric_particle_corr_xz_mean[ny],symmetric_particle_corr_yy_mean[ny],symmetric_particle_corr_yz_mean[ny],symmetric_particle_corr_zz_mean[ny],symmetric_norm_corr[ny]);
   }
 
fclose(fp52);

fp5=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/stress_2.txt","w");
fprintf(fp5,"%s Grid_index     Grid_point       sigma_x  sigmax_1      sigma_y sigmay_1       sigma_z sigma_z_1  tau_xy tauxy_1  \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

fprintf(fp5,"%8d %20.14lf %20.14lf %20.14lf %20.14lf  %20.14lf %20.14lf %20.14lf %20.14lf %20.14lf \n",
ny, y_position,sigma_xx[ny] ,sigma_xx_1[ny], sigma_yy[ny],sigma_yy_1[ny],sigma_zz[ny],sigma_zz_1[ny],tau_xy[ny],sigma_xy_1[ny]);
   }
 
fclose(fp5);

 fp18=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/corre_part_history.txt","a");
 fprintf(fp18 ,"%20.14lf %20.14lf %20.14lf %20.14lf  \n",run_time,p_x[n_atom/2],p_y[n_atom/2],
	   p_z[n_atom/2]);
   fclose(fp18);

 if(index_stress_varience==1)
   {
 fp20=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/varience_stress.txt","w");
 fprintf(fp20,"%s Grid_index     Grid_point       var_sigma_x1      var_sigma_y1        var_sigma_z_1       var_tauxy_1  \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

 fprintf(fp20,"%8d %20.14lf %20.14lf %20.14lf %20.14lf  %20.14lf  \n",
	ny,    y_position,    var_sigma_xx[ny] ,     var_sigma_yy[ny],       var_sigma_zz[ny],     var_tau_xy[ny]);
   }

 fclose(fp20);
   }


 fp21=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/vel_moment_4.txt","w");
 fprintf(fp21,"%s Grid_index     Grid_point       moment_4_x     moment_4_y       moment_4_z       \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp21,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,     avg_pvx_dash_4[ny] ,  avg_pvy_dash_4[ny] ,  avg_pvz_dash_4[ny]  );
   }
 fclose(fp21);
  

fp211=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/vel_moment_3.txt","w");
 fprintf(fp211,"%s Grid_index     Grid_point       moment_3_x     moment_3_y       moment_3_z       \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp211,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,     avg_pvx_dash_3[ny] ,  avg_pvy_dash_3[ny] ,  avg_pvz_dash_3[ny]  );
   }
 fclose(fp211);


fp212=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_moment_2.txt","w");
 fprintf(fp212,"%s Grid_index     Grid_point       moment_2_x     moment_2_y       moment_2_z       \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp212,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,     avg_accx_dash_2[ny] ,  avg_accy_dash_2[ny] ,  avg_accz_dash_2[ny]  );
   }
 fclose(fp212);


fp213=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_moment_3.txt","w");
 fprintf(fp213,"%s Grid_index     Grid_point       moment_3_x     moment_3_y       moment_3_z       \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp213,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,     avg_accx_dash_3[ny] ,  avg_accy_dash_3[ny] ,  avg_accz_dash_3[ny]  );
   }
 fclose(fp213);

fp214=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_moment_4.txt","w");
 fprintf(fp214,"%s Grid_index     Grid_point       moment_4_x     moment_4_y       moment_4_z       \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp214,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,     avg_accx_dash_4[ny] ,  avg_accy_dash_4[ny] ,  avg_accz_dash_4[ny]  );
   }
 fclose(fp214);


fp315=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_fluctuation.txt","w");
 fprintf(fp315,"%s Grid_index     Grid_point       acc_fluc_x_sq     acc_fluc_y_sq     acc_fluc_z_sq    \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp315,"%8d %20.14lf %.14e   %.14e   %.14e  \n",
	  ny,    y_position,     mean_accx_dash_2[ny] ,  mean_accy_dash_2[ny] ,  mean_accz_dash_2[ny]  );
   }
 fclose(fp315);

fp316=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/mean_acc.txt","w");
 fprintf(fp316,"%s Grid_index     Grid_point         mean_acc_x         mean_acc_y       mean_acc_z    \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp316,"%8d %20.14lf %.14e   %.14e   %.14e  \n",
	  ny,    y_position,    accx_mean[ny]-g_x ,  accy_mean[ny]-g_y ,  accz_mean[ny]-g_z );
   }
 fclose(fp316);






 //cout<<"----it is at ----00"<<endl;
 return;
}




//--------------------------------------------------------------------------------------
// Average and fluctuation calculation of different properties of angular velocities
//--------------------------------------------------------------------------------------
void angular_final_property()
{

 int ny;
 double  y_position;
 double  Var_Nofp_[Max_y_grid];
 double angvel_norm[Max_y_grid]; 
 double angvel_norm_corr[Max_y_grid];
   
  FILE *fp104;
  FILE *fp1041;
  FILE *fp105;
 
sampling_time=run_time-t_zero;

  
//------ for grid dependent property calculation----------


 for(ny=0;ny<Ndel_y;ny++) // loop for doing time averaging
   { 
     Nofp_mean[ny]=Nofp_mean[ny]/(double)n_sample; 

     angvelx_mean[ny]=angvelx_mean[ny]/(double)n_sample;    
     angvely_mean[ny]=angvely_mean[ny]/(double)n_sample;  
     angvelz_mean[ny]=angvelz_mean[ny]/(double)n_sample;

    
     angvel_corr_x_mean[ny] = Nofp_mean[ny]*angvelx_mean[ny];
     angvel_corr_y_mean[ny] = Nofp_mean[ny]*angvely_mean[ny];
     angvel_corr_z_mean[ny] = Nofp_mean[ny]*angvelz_mean[ny];

     angvel_norm[ny] = sqrt(angvelx_mean[ny]*angvelx_mean[ny] + angvely_mean[ny]*angvely_mean[ny] + angvelz_mean[ny]*angvelz_mean[ny]);

     angvel_norm_corr[ny] = Nofp_mean[ny]*angvel_norm[ny];
     
			 angsigma_xx[ny]=angsigma_xx[ny]/(double)n_sample;
			 angsigma_yy[ny]=angsigma_yy[ny]/(double)n_sample;
			 angsigma_zz[ny]=angsigma_zz[ny]/(double)n_sample;

			 tau_xy[ny]=tau_xy[ny]/(double)n_sample;

			 angsigma_xx_1[ny]=angsigma_xx_1[ny]/(double)n_sample- angvelx_mean[ny]* angvelx_mean[ny];
			 angsigma_yy_1[ny]=angsigma_yy_1[ny]/(double)n_sample- angvely_mean[ny]* angvely_mean[ny];
			 angsigma_zz_1[ny]=angsigma_zz_1[ny]/(double)n_sample- angvelz_mean[ny]* angvelz_mean[ny];
			 angsigma_xy_1[ny]=angsigma_xy_1[ny]/(double)n_sample- angvelx_mean[ny]* angvely_mean[ny];

			 Var_Nofp_[ny]=Nofp_inst_sq[ny]/(double)n_sample- Nofp_mean[ny]* Nofp_mean[ny];
			 
			 }

 if(index_stress_varience==1)
   {
      for(ny=0;ny<Ndel_y;ny++) // loop for doing time averaging
	{ 
	   var_angsigma_xx[ny]=var_angsigma_xx[ny]/(double)n_sample;
	   var_angsigma_yy[ny]=var_angsigma_yy[ny]/(double)n_sample;
	   var_angsigma_zz[ny]=var_angsigma_zz[ny]/(double)n_sample;
 	   var_angtau_xy[ny]=var_angtau_xy[ny]/(double)n_sample;
	}

   }

  if(index_fourth_moment==1)
    {
      for(ny=0;ny<Ndel_y;ny++)          // loop for doing time averaging
       {
	 avg_pangvelx_dash_4[ny]=avg_pangvelx_dash_4[ny]/(double)n_sample;
	 avg_pangvely_dash_4[ny]=avg_pangvely_dash_4[ny]/(double)n_sample;
	 avg_pangvelz_dash_4[ny]=avg_pangvelz_dash_4[ny]/(double)n_sample;
	 


	 avg_pangvelx_dash_3[ny]=avg_pangvelx_dash_3[ny]/(double)n_sample;
	 avg_pangvely_dash_3[ny]=avg_pangvely_dash_3[ny]/(double)n_sample;
	 avg_pangvelz_dash_3[ny]=avg_pangvelz_dash_3[ny]/(double)n_sample;
      }
   
    }

for(ny=0;ny<Ndel_y;ny++) 
    {

    angaccx_mean[ny]=angaccx_mean[ny]/(double)n_sample;
    angaccy_mean[ny]=angaccy_mean[ny]/(double)n_sample;
    angaccz_mean[ny]=angaccz_mean[ny]/(double)n_sample;

    mean_angaccx_dash_2[ny]= mean_angaccx_sq[ny]/(double)n_sample-angaccx_mean[ny]*angaccx_mean[ny];
    mean_angaccy_dash_2[ny]= mean_angaccy_sq[ny]/(double)n_sample-angaccy_mean[ny]*angaccy_mean[ny];
    mean_angaccz_dash_2[ny]= mean_angaccz_sq[ny]/(double)n_sample-angaccz_mean[ny]*angaccz_mean[ny];
    }



  for(ny=0;ny<Ndel_y;ny++) 
    {
      avg_angaccx_dash_2[ny]=avg_angaccx_dash_2[ny]/(double)n_sample;
      avg_angaccy_dash_2[ny]=avg_angaccy_dash_2[ny]/(double)n_sample;
      avg_angaccz_dash_2[ny]=avg_angaccz_dash_2[ny]/(double)n_sample;
      

      avg_angaccx_dash_3[ny]=avg_angaccx_dash_3[ny]/(double)n_sample;
      avg_angaccy_dash_3[ny]=avg_angaccy_dash_3[ny]/(double)n_sample;
      avg_angaccz_dash_3[ny]=avg_angaccz_dash_3[ny]/(double)n_sample;
      
      avg_angaccx_dash_4[ny]=avg_angaccx_dash_4[ny]/(double)n_sample;
      avg_angaccy_dash_4[ny]=avg_angaccy_dash_4[ny]/(double)n_sample;
      avg_angaccz_dash_4[ny]=avg_angaccz_dash_4[ny]/(double)n_sample;
      
    }



fp104=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angvelocity_2.txt","w");
fprintf(fp104," %s Sampling duration =%20.14lf\n","%", sampling_time);
fprintf(fp104," %s Grid_index  Grid point angVx    angVy    angVz   angvel_norm No_of_Part  var_no_of_P   No of Collision \n ","%"); 
  

 for(ny=0;ny<Ndel_y;ny++)
   {
   y_position=(ny+0.5)*del_y+(0.5*sigma);

fprintf(fp104,"%8d %20.14lf %20.14lf %20.14lf %20.14lf %20.14lf  %20.14lf %20.14lf %15d\n",
ny, y_position,angvelx_mean[ny] , angvely_mean[ny],angvelz_mean[ny],angvel_norm[ny], Nofp_mean[ny], Var_Nofp_[ny], Nof_p_p_coll[ny]);
 
   }

fclose(fp104);

fp1041=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angvelocity_particle_corr_2.txt","w");
fprintf(fp1041," %s Grid_index  Grid point angvel_corr_x_mean    angvel_corr_y_mean    angvel_corr_z_mean   angvel_norm_corr \n ","%"); 
  

 for(ny=0;ny<Ndel_y;ny++)
   {
   y_position=(ny+0.5)*del_y+(0.5*sigma);

fprintf(fp1041,"%8d %20.14lf %20.14lf %20.14lf %20.14lf %20.14lf  \n",
ny, y_position,angvel_corr_x_mean[ny] , angvel_corr_y_mean[ny],angvel_corr_z_mean[ny], angvel_norm_corr[ny]);
 
   }

fclose(fp1041);


fp105=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angstress_2.txt","w");
fprintf(fp105,"%s Grid_index     Grid_point       angsigma_x  angsigmax_1      angsigma_y angsigmay_1       angsigma_z angsigma_z_1  angtau_xy angtauxy_1  \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

fprintf(fp105,"%8d %20.14lf %20.14lf %20.14lf %20.14lf  %20.14lf %20.14lf %20.14lf %20.14lf %20.14lf \n",
ny, y_position,angsigma_xx[ny] ,angsigma_xx_1[ny], angsigma_yy[ny],angsigma_yy_1[ny],angsigma_zz[ny],angsigma_zz_1[ny],angtau_xy[ny],angsigma_xy_1[ny]);
   }
 
fclose(fp105);

 if(index_stress_varience==1)
   {
 fp1020=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/varience_angstress.txt","w");
 fprintf(fp1020,"%s Grid_index     Grid_point       angvar_sigma_x1      angvar_sigma_y1        angvar_sigma_z_1       angvar_tauxy_1  \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

 fprintf(fp1020,"%8d %20.14lf %20.14lf %20.14lf %20.14lf  %20.14lf  \n",
	ny,    y_position,    var_angsigma_xx[ny] ,     var_angsigma_yy[ny],       var_angsigma_zz[ny],     var_angtau_xy[ny]);
   }

 fclose(fp1020);
   }


 fp1021=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angvel_moment_4.txt","w");
 fprintf(fp1021,"%s Grid_index     Grid_point       moment_4_x     moment_4_y       moment_4_z       \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp1021,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,     avg_pangvelx_dash_4[ny] ,  avg_pangvely_dash_4[ny] ,  avg_pangvelz_dash_4[ny]  );
   }
 fclose(fp1021);
  

fp10211=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angvel_moment_3.txt","w");
 fprintf(fp10211,"%s Grid_index     Grid_point       moment_3_x     moment_3_y       moment_3_z       \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp10211,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,     avg_pangvelx_dash_3[ny] ,  avg_pangvely_dash_3[ny] ,  avg_pangvelz_dash_3[ny]  );
   }
 fclose(fp10211);


fp10212=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_moment_2.txt","w");
 fprintf(fp10212,"%s Grid_index     Grid_point       moment_2_x     moment_2_y       moment_2_z       \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp10212,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,     avg_angaccx_dash_2[ny] ,  avg_angaccy_dash_2[ny] ,  avg_angaccz_dash_2[ny]  );
   }
 fclose(fp10212);


fp10213=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_moment_3.txt","w");
 fprintf(fp10213,"%s Grid_index     Grid_point       moment_3_x     moment_3_y       moment_3_z       \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp10213,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,     avg_angaccx_dash_3[ny] ,  avg_angaccy_dash_3[ny] ,  avg_angaccz_dash_3[ny]  );
   }
 fclose(fp10213);

fp10214=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_moment_4.txt","w");
 fprintf(fp10214,"%s Grid_index     Grid_point       moment_4_x     moment_4_y       moment_4_z       \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp10214,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,     avg_angaccx_dash_4[ny] ,  avg_angaccy_dash_4[ny] ,  avg_angaccz_dash_4[ny]  );
   }
 fclose(fp10214);


fp10315=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_fluctuation.txt","w");
 fprintf(fp10315,"%s Grid_index     Grid_point       acc_fluc_x_sq     acc_fluc_y_sq     acc_fluc_z_sq    \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp10315,"%8d %20.14lf %.14e   %.14e   %.14e  \n",
	  ny,    y_position,     mean_angaccx_dash_2[ny] ,  mean_angaccy_dash_2[ny] ,  mean_angaccz_dash_2[ny]  );
   }
 fclose(fp10315);

fp10316=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/mean_angacc.txt","w");
 fprintf(fp10316,"%s Grid_index     Grid_point         mean_acc_x         mean_acc_y       mean_acc_z    \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp10316,"%8d %20.14lf %.14e   %.14e   %.14e  \n",
	  ny,    y_position,    angaccx_mean[ny] ,  angaccy_mean[ny] ,  angaccz_mean[ny] );
   }
 fclose(fp10316);






 //cout<<"----it is at ----00"<<endl;
 return;
}

//-------------------------------------------------------------------
    
// CALCULATION OF VELOCITY DISTRIBUTION FUNCTION
//---------------------------------------------------------------
     
           void distri_func() 
          {
              int k, j,ny ; 
              double fnorm,p_v;
	      double vx, fof_v,fof_vx,vy,fof_vy,vz,fof_vz;
	      double  origins,RDELSI;
	      double gof_r[Max_y_grid],radius;
	      double  area_shell;
	      double density_per_area;
	      double radius_by_sigma;
              double pi =3.14159265;
	      double var_fof_vx,var_fof_vy,var_fof_vz;

	     

	      FILE *fp6;
	      FILE *fp7;
	      FILE *fp71;
	      FILE *fp72;




	      FILE *fp8;
	      FILE *fp81;
	      FILE *fp82;


	      FILE *fp9;
	      FILE *fp91;
	      FILE *fp92;

	      FILE *fp10;
  fp6=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_center.txt","w");
  fp7=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_center_x.txt","w");
  fp71=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_middle_x.txt","w");
  fp72=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_wall_x.txt","w");


  fp8=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_center_y.txt","w");
  fp81=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_middle_y.txt","w");
  fp82=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_wall_y.txt","w");

  fp9=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_center_z.txt","w");
  fp91=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_middle_z.txt","w");
  fp92=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_wall_z.txt","w");

	  
	                                           //n_sample = (Tot_Nof_coll- eqb_coll)/n_sample
	    // fnorm =(double)(n_atom*n_sample);
                                                    //tmp_ave = 2.0*sum_KE/(3.0*fnorm);
  fnorm =(double)(tot_atom_count[1]);                  //   fact = 1.0/sqrt(2.0*pi*tmp_ave);
                                                       //HMB= 0.0;
    
//! Loop over veolocity shells each of thickness del_v

      for(k=1;k<=Ndel_v;k++)   // k should start frm 1 not from 0  
            {
      p_v= -p_v_max+ del_v*(k-1);
        
      //Sampled time average distribution

      fof_v = (double)(Nof_v[1][k])/(del_v* fnorm);
      fprintf(fp6,"%20.14lf %20.14lf\n",p_v,fof_v);
	   }

	       cout<<"****************"<<n_sample<<endl;
	      

	    for(k=1;k<=Ndel_vx;k++)           // k should start from 1 not from 0  
                 {
		   vx= -p_vx_max+ del_vx*(k-1);

		   fof_vx= (double)( Nof_vx[1][k])/((double)tot_atom_count[1]*del_vx);
		   var_fof_vx=(fof_vx_sq[1][k]/(double)n_sample)-(fof_vx*fof_vx);
		   fprintf(fp7,"%20.14lf %20.14lf %20.14lf\n",vx,fof_vx, var_fof_vx);

		   fof_vx= (double)( Nof_vx[2][k])/((double)tot_atom_count[2]*del_vx);
		   var_fof_vx=(fof_vx_sq[2][k]/(double)n_sample)-(fof_vx*fof_vx);
		   fprintf(fp71,"%20.14lf %20.14lf %20.14lf\n",vx,fof_vx, var_fof_vx);

		   fof_vx= (double)( Nof_vx[3][k])/((double)tot_atom_count[3]*del_vx);
		   var_fof_vx=(fof_vx_sq[3][k]/(double)n_sample)-(fof_vx*fof_vx);
		   fprintf(fp72,"%20.14lf %20.14lf %20.14lf\n",vx,fof_vx, var_fof_vx);

		 }
        

	   for(k=1;k<=Ndel_vy;k++)                   // k should start from 1 not from 0  
		 {
	    vy= -p_vy_max+ del_vy*(k-1);

	    fof_vy= (double)(Nof_vy[1][k])/((double)tot_atom_count[1]*del_vy);
	    var_fof_vy=(fof_vy_sq[1][k]/(double)n_sample)-(fof_vy*fof_vy);
            fprintf(fp8,"%20.14lf %20.14lf %20.14lf\n",vy,fof_vy,var_fof_vy);

	    fof_vy= (double)(Nof_vy[2][k])/((double)tot_atom_count[2]*del_vy);
	    var_fof_vy=(fof_vy_sq[2][k]/(double)n_sample)-(fof_vy*fof_vy);
            fprintf(fp81,"%20.14lf %20.14lf %20.14lf\n",vy,fof_vy,var_fof_vy);

	    fof_vy= (double)(Nof_vy[3][k])/((double)tot_atom_count[3]*del_vy);
	    var_fof_vy=(fof_vy_sq[3][k]/(double)n_sample)-(fof_vy*fof_vy);
            fprintf(fp82,"%20.14lf %20.14lf %20.14lf\n",vy,fof_vy,var_fof_vy);



	   }

      for(k=1;k<=Ndel_vz;k++)
         {
	   vz= -p_vz_max+ del_vz*(k-1);

	   fof_vz= (double)(Nof_vz[1][k])/((double)tot_atom_count[1]*del_vz);
	   var_fof_vz=(fof_vz_sq[1][k]/(double)n_sample)-(fof_vz*fof_vz);
	   fprintf(fp9 ,"%20.14lf %20.14lf %20.14lf\n",vz,fof_vz,var_fof_vz);

	   fof_vz= (double)(Nof_vz[2][k])/((double)tot_atom_count[2]*del_vz);
	   var_fof_vz=(fof_vz_sq[2][k]/(double)n_sample)-(fof_vz*fof_vz);
	   fprintf(fp91 ,"%20.14lf %20.14lf %20.14lf\n",vz,fof_vz,var_fof_vz);

	   fof_vz= (double)(Nof_vz[3][k])/((double)tot_atom_count[3]*del_vz);
	   var_fof_vz=(fof_vz_sq[3][k]/(double)n_sample)-(fof_vz*fof_vz);
	   fprintf(fp92 ,"%20.14lf %20.14lf %20.14lf\n",vz,fof_vz,var_fof_vz);
	  }
       

//----------------------------------------------------
//     Normalization of  counters for radial distribution function
//----------------------------------------------------

				  

      //double density = 6.0*eta/pi;
 //  IF(MOD(Tot_Nof_coll,k_write)==0)THEN

   fp10=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/radial_dist_2.txt","w");
      
	 
   // ! origins = (Tot_Nof_coll / KSAMPL)*(n_atom/2)
   
   // !Loop over radial shells
   RDELSI = del_r/sigma;

     
  for(ny=0;ny<Ndel_y;ny=ny+20)
    {
      
      fprintf(fp10," at ny=%d\n",ny);
    origins = (int)( N_atom_count[ny]);

    for(j=1;j<=N_rdels;j++)
     {
      
         if(Ngof_r[ny][j] > 0)
	   {
	     //RADIUS= RDELSI*(double)j;
	     radius=del_r*(double)j;
	     //vol_shell=4.0* pi * RDELSI * RADIUS*RADIUS+(pi*RDELSI*RDELSI*RDELSI)/3.0;
	     area_shell= 2.0*pi*radius*del_r;
	     density_per_area=(Nofp_mean[ny]/del_y)*sigma,
	     //gof_r[ny] = (double)(Ngof_r[ny][j])/(density*origins*vol_shell);
	     gof_r[ny] = (double)(Ngof_r[ny][j])/(density_per_area*origins*area_shell);


	   // Correct first value of g(r); the first sampling interval is hahf of rdel
	   //            because the spheres are inpenitrable

           // if(RADIUS < (1.+ RDELSI))
	      if(radius<(sigma+del_r))
	      {
		//RADIUS = 1.0 + 0.5 * RDELSI;
		radius=sigma+0.5*del_r;
		gof_r[ny] = 2. *gof_r[ny];
              }
	      radius_by_sigma=radius/sigma;

	     
	fprintf(fp10,"%d %d %20.14lf %20.14lf \n", j, Ngof_r[ny][j],radius_by_sigma,gof_r[ny] ) ;  
	   }
	
     }

    }     

 	  return;   
	  }                   
	   

    //-------------------------------------------------------------------
    
// CALCULATION OF ANGULAR VELOCITY DISTRIBUTION FUNCTION
//---------------------------------------------------------------
     
           void angular_distri_func() 
          {
              int k ; 
              double fnorm,p_angvel;
	      double angvelx, fof_angvel,fof_angvelx,angvely,fof_angvely,angvelz,fof_angvelz;
	     
	      double var_fof_angvelx,var_fof_angvely,var_fof_angvelz;

	     

	      FILE *fp106;
	      FILE *fp107;
	      FILE *fp1071;
	      FILE *fp1072;




	      FILE *fp108;
	      FILE *fp1081;
	      FILE *fp1082;


	      FILE *fp109;
	      FILE *fp1091;
	      FILE *fp1092;

	      
  fp106=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_center.txt","w");
  fp107=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_center_x.txt","w");
  fp1071=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_middle_x.txt","w");
  fp1072=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_wall_x.txt","w");


  fp108=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_center_y.txt","w");
  fp1081=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_middle_y.txt","w");
  fp1082=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_wall_y.txt","w");

  fp109=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_center_z.txt","w");
  fp1091=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_middle_z.txt","w");
  fp1092=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_wall_z.txt","w");

	  
	                                           //n_sample = (Tot_Nof_coll- eqb_coll)/n_sample
	    // fnorm =(double)(n_atom*n_sample);
                                                    //tmp_ave = 2.0*sum_KE/(3.0*fnorm);
  fnorm =(double)(tot_atom_count[1]);                  //   fact = 1.0/sqrt(2.0*pi*tmp_ave);
                                                       //HMB= 0.0;
    
//! Loop over angular veolocity shells each of thickness del_angvel

      for(k=1;k<=Ndel_angvel;k++)   // k should start frm 1 not from 0  
            {
      p_angvel= -p_angvel_max+ del_angvel*(k-1);
        
      //Sampled time average distribution

      fof_angvel = (double)(Nof_angvel[1][k])/(del_angvel* fnorm);
      fprintf(fp106,"%20.14lf %20.14lf\n",p_angvel,fof_angvel);
	   }

	       cout<<"****************"<<n_sample<<endl;
	      

	    for(k=1;k<=Ndel_angvelx;k++)           // k should start from 1 not from 0  
                 {
		   angvelx= -p_angvelx_max+ del_angvelx*(k-1);

		   fof_angvelx= (double)( Nof_angvelx[1][k])/((double)tot_atom_count[1]*del_angvelx);
		   var_fof_angvelx=(fof_angvelx_sq[1][k]/(double)n_sample)-(fof_angvelx*fof_angvelx);
		   fprintf(fp107,"%20.14lf %20.14lf %20.14lf\n",angvelx,fof_angvelx, var_fof_angvelx);

		   fof_angvelx= (double)( Nof_angvelx[2][k])/((double)tot_atom_count[2]*del_angvelx);
		   var_fof_angvelx=(fof_angvelx_sq[2][k]/(double)n_sample)-(fof_angvelx*fof_angvelx);
		   fprintf(fp1071,"%20.14lf %20.14lf %20.14lf\n",angvelx,fof_angvelx, var_fof_angvelx);

		   fof_angvelx= (double)( Nof_angvelx[3][k])/((double)tot_atom_count[3]*del_angvelx);
		   var_fof_angvelx=(fof_angvelx_sq[3][k]/(double)n_sample)-(fof_angvelx*fof_angvelx);
		   fprintf(fp1072,"%20.14lf %20.14lf %20.14lf\n",angvelx,fof_angvelx, var_fof_angvelx);

		 }
        

	   for(k=1;k<=Ndel_angvely;k++)                   // k should start from 1 not from 0  
		 {
	    angvely= -p_angvely_max+ del_angvely*(k-1);

	    fof_angvely= (double)(Nof_angvely[1][k])/((double)tot_atom_count[1]*del_angvely);
	    var_fof_angvely=(fof_angvely_sq[1][k]/(double)n_sample)-(fof_angvely*fof_angvely);
            fprintf(fp108,"%20.14lf %20.14lf %20.14lf\n",angvely,fof_angvely,var_fof_angvely);

	    fof_angvely= (double)(Nof_angvely[2][k])/((double)tot_atom_count[2]*del_angvely);
	    var_fof_angvely=(fof_angvely_sq[2][k]/(double)n_sample)-(fof_angvely*fof_angvely);
            fprintf(fp1081,"%20.14lf %20.14lf %20.14lf\n",angvely,fof_angvely,var_fof_angvely);

	    fof_angvely= (double)(Nof_angvely[3][k])/((double)tot_atom_count[3]*del_angvely);
	    var_fof_angvely=(fof_angvely_sq[3][k]/(double)n_sample)-(fof_angvely*fof_angvely);
            fprintf(fp1082,"%20.14lf %20.14lf %20.14lf\n",angvely,fof_angvely,var_fof_angvely);



	   }

      for(k=1;k<=Ndel_angvelz;k++)
         {
	   angvelz= -p_angvelz_max+ del_angvelz*(k-1);

	   fof_angvelz= (double)(Nof_angvelz[1][k])/((double)tot_atom_count[1]*del_angvelz);
	   var_fof_angvelz=(fof_angvelz_sq[1][k]/(double)n_sample)-(fof_angvelz*fof_angvelz);
	   fprintf(fp109 ,"%20.14lf %20.14lf %20.14lf\n",angvelz,fof_angvelz,var_fof_angvelz);

	   fof_angvelz= (double)(Nof_angvelz[2][k])/((double)tot_atom_count[2]*del_angvelz);
	   var_fof_angvelz=(fof_angvelz_sq[2][k]/(double)n_sample)-(fof_angvelz*fof_angvelz);
	   fprintf(fp1091 ,"%20.14lf %20.14lf %20.14lf\n",angvelz,fof_angvelz,var_fof_angvelz);

	   fof_angvelz= (double)(Nof_angvelz[3][k])/((double)tot_atom_count[3]*del_angvelz);
	   var_fof_angvelz=(fof_angvelz_sq[3][k]/(double)n_sample)-(fof_angvelz*fof_angvelz);
	   fprintf(fp1092 ,"%20.14lf %20.14lf %20.14lf\n",angvelz,fof_angvelz,var_fof_angvelz);
	  }
       
      
     }  



    void vel_dist_free_flight()
     {

       int i, ny,kk,jj;
       int N_shellx, N_shelly,N_shellz;
       double  offset,offset_x, offset_y , offset_z ;
     

       n_sample_vel_dist_free_flight=n_sample_vel_dist_free_flight+1;


    for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
    Nofp_inst[ny]=0;
    vx_inst[ny]=0.0;
    vy_inst[ny]=0.0;
    vz_inst[ny]=0.0;
       }  


for(i=0;i<n_atom;i++)// summing up all the velocities in each the y-grid for all the particle 
    {
      ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

      // cout<<"in property"<<p_y[i]<<" "<<ny<<endl;


              Nofp_inst[ny]= Nofp_inst[ny]+ 1;

              vx_inst[ny]= vx_inst[ny]+p_vx[i];
              vy_inst[ny]= vy_inst[ny]+p_vy[i];  
              vz_inst[ny]= vz_inst[ny]+p_vz[i];
      }

    
    
    for(ny=0;ny<Ndel_y;ny++) // doing the particle average in each y-grid-in each sample frame
       {
    if(Nofp_inst[ny]>0.0)
      {
     vx_inst[ny]= vx_inst[ny]/Nofp_inst[ny];  // vx_inst[ny] is basically particle averaged inst in time.
     vy_inst[ny]= vy_inst[ny]/Nofp_inst[ny];
     vz_inst[ny]= vz_inst[ny]/Nofp_inst[ny];
      }
       }

 for(i=0;i<n_atom;i++)
   {
 ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
 pvx_dash[i]=p_vx[i]-vx_inst[ny];  // calculating fluctuation for each particle
 pvy_dash[i]=p_vy[i]-vy_inst[ny];
 pvz_dash[i]=p_vz[i]-vz_inst[ny];

   }

            for(jj=0;jj<Max_dim_position;jj++)
                for(i=0;i<Ndel_v;i++)
		    Nof_v_inst[jj][i]=0;

             for(jj=0;jj<Max_dim_position;jj++)
		  for(i=0;i<Ndel_vx;i++)
		    Nof_vx_inst[jj][i]=0;

	     for(jj=0;jj<Max_dim_position;jj++)
		  for(i=0;i<Ndel_vy;i++)
		    Nof_vy_inst[jj][i]=0;

	     for(jj=0;jj<Max_dim_position;jj++)
		  for(i=0;i<Ndel_vz;i++)
		    Nof_vz_inst[jj][i]=0;

		 for(jj=0;jj<Max_dim_position;jj++) atom_count_inst[jj]=0;
           
            offset =1.5+p_v_max/del_v;
            offset_x =1.5+p_vx_max/del_vx;
            offset_y =1.5+p_vy_max/del_vy;
            offset_z =1.5+p_vz_max/del_vz;


  for(i=0;i<n_atom;i++)
           {
	     position_index=0;    //initialization of position_index

	     if(p_y[i]>y_min_for_distribution_1 && p_y[i]<y_max_for_distribution_1)position_index=1;
	     if(p_y[i]>y_min_for_distribution_2 && p_y[i]<y_max_for_distribution_2)position_index=2;
	     if(p_y[i]>y_min_for_distribution_3 && p_y[i]<y_max_for_distribution_3)position_index=3;
	     

	     if(position_index>0)
	     
	         {                                                
	       
 if((pvx_dash[i]>=-p_vx_max && pvx_dash[i]<=p_vx_max)&&(pvy_dash[i]>=-p_vy_max && pvy_dash[i]<=p_vy_max)&& 
(pvz_dash[i]>=-p_vz_max && pvz_dash[i]<=p_vz_max))  // to avoid exceeding arrey subscript accidentally
   {
     
	     N_shellx=(int)(pvx_dash[i]/del_vx+offset_x);
             Nof_vx_inst[position_index][N_shellx]= Nof_vx_inst[position_index][N_shellx]+1;
	     Nof_vx_free_flight[position_index][N_shellx]= Nof_vx_free_flight[position_index][N_shellx] +1 ; 

	  N_shelly=(int)(pvy_dash[i]/del_vy+offset_y);
	  Nof_vy_inst[position_index][N_shelly]=Nof_vy_inst[position_index][N_shelly]+1;
	  Nof_vy_free_flight[position_index][N_shelly]= Nof_vy_free_flight[position_index][N_shelly]+1;
             
	  N_shellz=(int)(pvz_dash[i]/del_vz+offset_z);
	  Nof_vz_inst[position_index][N_shellz]=Nof_vz_inst[position_index][N_shellz]+1;
	  Nof_vz_free_flight[position_index][N_shellz]= Nof_vz_free_flight[position_index][N_shellz]+1;
	  
	  atom_count_inst[position_index]=atom_count_inst[position_index]+1;
	  tot_atom_count_vel_free_flight[position_index]=tot_atom_count_vel_free_flight[position_index]+1;
   }
             }
	   }

 for(position_index=1;position_index<=max_position_index;position_index++) 
  for( kk=1;kk<=Ndel_vx;kk++)           // k should start from 1 not from 0  
                 {
		   if(atom_count_inst[position_index]>0)	  
      fof_vx_free_flight_sq[position_index][kk]=fof_vx_free_flight_sq[position_index][kk]+(double)(Nof_vx_inst[position_index][kk]/
(atom_count_inst[position_index]*del_vx))*(double)(Nof_vx_inst[position_index][kk]/(atom_count_inst[position_index]*del_vx));
		 }


 for(position_index=1;position_index<=max_position_index;position_index++) 
	   for( kk=1;kk<=Ndel_vy;kk++)           // k should start from 1 not from 0  
                 {
	if(atom_count_inst[position_index]>0)	   
   fof_vy_free_flight_sq[position_index][kk]=fof_vy_free_flight_sq[position_index][kk]+(double)( Nof_vy_inst[position_index][kk]/
(atom_count_inst[position_index]*del_vy))*(double)( Nof_vy_inst[position_index][kk]/(atom_count_inst[position_index]*del_vy));
		 }


 for(position_index=1;position_index<=max_position_index;position_index++) 
	   for( kk=1;kk<=Ndel_vz;kk++)           // k should start from 1 not from 0  
                 {
		if(atom_count_inst[position_index]>0)   
   fof_vz_free_flight_sq[position_index][kk]=fof_vz_free_flight_sq[position_index][kk]+ (double)( Nof_vz_inst[position_index][kk]/
(atom_count_inst[position_index]*del_vz))*(double)( Nof_vz_inst[position_index][kk]/(atom_count_inst[position_index]*del_vz));
		 }

	   return;
     }



 void vel_distri_func_free_flight()
{ 
              int k; 
              double fnorm;
	      double vx, fof_vx,vy,fof_vy,vz,fof_vz;
	      double var_fof_vx,var_fof_vy,var_fof_vz;

	     

	      FILE *fp40;
	      FILE *fp401;
	      FILE *fp402;

	      FILE *fp41;
	      FILE *fp411;
	      FILE *fp412;


	      FILE *fp42;
	      FILE *fp421;
	      FILE *fp422;
	    
 fp40=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_free_flight_center_x.txt","w");
 fp401=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_free_flight_middle_x.txt","w");
 fp402=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_free_flight_wall_x.txt","w");

  fp41=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_free_flight_center_y.txt","w");
 fp411=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_free_flight_middle_y.txt","w");
 fp412=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_free_flight_wall_y.txt","w");

  fp42=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_free_flight_center_z.txt","w");
  fp421=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_free_flight_middle_z.txt","w");
  fp422=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_dist_free_flight_wall_z.txt","w");

   fnorm =(double)(tot_atom_count_vel_free_flight[1]);                  
               

	    for(k=1;k<=Ndel_vx;k++)           // k should start from 1 not from 0-- this is right only 
                 {
		   vx= -p_vx_max+ del_vx*(k-1);

		   fof_vx= (double)( Nof_vx_free_flight[1][k])/((double)(tot_atom_count_vel_free_flight[1])*del_vx);
		   var_fof_vx=(fof_vx_free_flight_sq[1][k]/(double)n_sample_vel_dist_free_flight)-(fof_vx*fof_vx);
		   fprintf(fp40,"%20.14lf %20.14lf %20.14lf\n",vx,fof_vx, var_fof_vx);

		   fof_vx= (double)( Nof_vx_free_flight[2][k])/((double)(tot_atom_count_vel_free_flight[2])*del_vx);
		   var_fof_vx=(fof_vx_free_flight_sq[2][k]/(double)n_sample_vel_dist_free_flight)-(fof_vx*fof_vx);
		   fprintf(fp401,"%20.14lf %20.14lf %20.14lf\n",vx,fof_vx, var_fof_vx);

		   fof_vx= (double)( Nof_vx_free_flight[3][k])/((double)(tot_atom_count_vel_free_flight[3])*del_vx);
		   var_fof_vx=(fof_vx_free_flight_sq[3][k]/(double)n_sample_vel_dist_free_flight)-(fof_vx*fof_vx);
		   fprintf(fp402,"%20.14lf %20.14lf %20.14lf\n",vx,fof_vx, var_fof_vx);

		 }
        

	   for(k=1;k<=Ndel_vy;k++)                   // k should start from 1 not from 0  
		 {
		   vy= -p_vy_max+ del_vy*(k-1);

		   fof_vy= (double)(Nof_vy_free_flight[1][k])/((double)(tot_atom_count_vel_free_flight[1])*del_vy);
		   var_fof_vy=(fof_vy_free_flight_sq[1][k]/(double)n_sample_vel_dist_free_flight)-(fof_vy*fof_vy);
	    	   fprintf(fp41,"%20.14lf %20.14lf %20.14lf\n",vy,fof_vy,var_fof_vy);

		   fof_vy= (double)(Nof_vy_free_flight[2][k])/((double)(tot_atom_count_vel_free_flight[2])*del_vy);
		   var_fof_vy=(fof_vy_free_flight_sq[2][k]/(double)n_sample_vel_dist_free_flight)-(fof_vy*fof_vy);
	    	   fprintf(fp411,"%20.14lf %20.14lf %20.14lf\n",vy,fof_vy,var_fof_vy);

		   fof_vy= (double)(Nof_vy_free_flight[3][k])/((double)(tot_atom_count_vel_free_flight[3])*del_vy);
		   var_fof_vy=(fof_vy_free_flight_sq[3][k]/(double)n_sample_vel_dist_free_flight)-(fof_vy*fof_vy);
	    	   fprintf(fp412,"%20.14lf %20.14lf %20.14lf\n",vy,fof_vy,var_fof_vy);
		   
		 }



	   for(k=1;k<=Ndel_vz;k++)
	     {
	       vz= -p_vz_max+ del_vz*(k-1);

	       fof_vz= (double)(Nof_vz_free_flight[1][k])/((double)tot_atom_count_vel_free_flight[1]*del_vz);
	       var_fof_vz=(fof_vz_free_flight_sq[1][k]/ (double)n_sample_vel_dist_free_flight)-(fof_vz*fof_vz);
	       fprintf(fp42 ,"%20.14lf %20.14lf %20.14lf\n",vz,fof_vz,var_fof_vz);

	       fof_vz= (double)(Nof_vz_free_flight[2][k])/((double)tot_atom_count_vel_free_flight[2]*del_vz);
	       var_fof_vz=(fof_vz_free_flight_sq[2][k]/(double)n_sample_vel_dist_free_flight)-(fof_vz*fof_vz);
	       fprintf(fp421 ,"%20.14lf %20.14lf %20.14lf\n",vz,fof_vz,var_fof_vz);

	       fof_vz= (double)(Nof_vz_free_flight[3][k])/((double)tot_atom_count_vel_free_flight[3]*del_vz);
	       var_fof_vz=(fof_vz_free_flight_sq[3][k]/(double)n_sample_vel_dist_free_flight)-(fof_vz*fof_vz);
	       fprintf(fp422 ,"%20.14lf %20.14lf %20.14lf\n",vz,fof_vz,var_fof_vz);

	     }
       
     }

void angular_vel_dist_free_flight()
     {

       int i, ny,kk,jj;
       int N_shellx, N_shelly,N_shellz;
       double  offset,offset_x, offset_y , offset_z ;
     

       n_sample_angvel_dist_free_flight=n_sample_angvel_dist_free_flight+1;


    for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
    Nofp_inst[ny]=0;
    angvelx_inst[ny]=0.0;
    angvely_inst[ny]=0.0;
    angvelz_inst[ny]=0.0;
       }  


for(i=0;i<n_atom;i++)// summing up all the velocities in each the y-grid for all the particle 
    {
      ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

      // cout<<"in property"<<p_y[i]<<" "<<ny<<endl;


              Nofp_inst[ny]= Nofp_inst[ny]+ 1;

              angvelx_inst[ny]= angvelx_inst[ny]+p_angvelx[i];
              angvely_inst[ny]= angvely_inst[ny]+p_angvely[i];  
              angvelz_inst[ny]= angvelz_inst[ny]+p_angvelz[i];
      }

    
    
    for(ny=0;ny<Ndel_y;ny++) // doing the particle average in each y-grid-in each sample frame
       {
    if(Nofp_inst[ny]>0.0)
      {
     angvelx_inst[ny]= angvelx_inst[ny]/Nofp_inst[ny];  // vx_inst[ny] is basically particle averaged inst in time.
     angvely_inst[ny]= angvely_inst[ny]/Nofp_inst[ny];
     angvelz_inst[ny]= angvelz_inst[ny]/Nofp_inst[ny];
      }
       }

 for(i=0;i<n_atom;i++)
   {
 ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
 pangvelx_dash[i]=p_angvelx[i]-angvelx_inst[ny];  // calculating fluctuation for each particle
 pangvely_dash[i]=p_angvely[i]-angvely_inst[ny];
 pangvelz_dash[i]=p_angvelz[i]-angvelz_inst[ny];

   }

            for(jj=0;jj<Max_dim_position;jj++)
                for(i=0;i<Ndel_angvel;i++)
		    Nof_angvel_inst[jj][i]=0;

             for(jj=0;jj<Max_dim_position;jj++)
		  for(i=0;i<Ndel_angvelx;i++)
		    Nof_angvelx_inst[jj][i]=0;

	     for(jj=0;jj<Max_dim_position;jj++)
		  for(i=0;i<Ndel_angvely;i++)
		    Nof_angvely_inst[jj][i]=0;

	     for(jj=0;jj<Max_dim_position;jj++)
		  for(i=0;i<Ndel_angvelz;i++)
		    Nof_angvelz_inst[jj][i]=0;

		 for(jj=0;jj<Max_dim_position;jj++) atom_count_inst[jj]=0;
           
            offset =1.5+p_angvel_max/del_angvel;
            offset_x =1.5+p_angvelx_max/del_angvelx;
            offset_y =1.5+p_angvely_max/del_angvely;
            offset_z =1.5+p_angvelz_max/del_angvelz;


  for(i=0;i<n_atom;i++)
           {
	     position_index=0;    //initialization of position_index

	     if(p_y[i]>y_min_for_distribution_1 && p_y[i]<y_max_for_distribution_1)position_index=1;
	     if(p_y[i]>y_min_for_distribution_2 && p_y[i]<y_max_for_distribution_2)position_index=2;
	     if(p_y[i]>y_min_for_distribution_3 && p_y[i]<y_max_for_distribution_3)position_index=3;
	     

	     if(position_index>0)
	     
	         {                                                
	       
 if((pangvelx_dash[i]>=-p_angvelx_max && pangvelx_dash[i]<=p_angvelx_max)&&(pangvely_dash[i]>=-p_angvely_max && pangvely_dash[i]<=p_angvely_max)&& 
(pangvelz_dash[i]>=-p_angvelz_max && pangvelz_dash[i]<=p_angvelz_max))  // to aangveloid exceeding arrey subscript accidentally
   {
     
	     N_shellx=(int)(pangvelx_dash[i]/del_angvelx+offset_x);
             Nof_angvelx_inst[position_index][N_shellx]= Nof_angvelx_inst[position_index][N_shellx]+1;
	     Nof_angvelx_free_flight[position_index][N_shellx]= Nof_angvelx_free_flight[position_index][N_shellx] +1 ; 

	  N_shelly=(int)(pangvely_dash[i]/del_angvely+offset_y);
	  Nof_angvely_inst[position_index][N_shelly]=Nof_angvely_inst[position_index][N_shelly]+1;
	  Nof_angvely_free_flight[position_index][N_shelly]= Nof_angvely_free_flight[position_index][N_shelly]+1;
             
	  N_shellz=(int)(pangvelz_dash[i]/del_angvelz+offset_z);
	  Nof_angvelz_inst[position_index][N_shellz]=Nof_angvelz_inst[position_index][N_shellz]+1;
	  Nof_angvelz_free_flight[position_index][N_shellz]= Nof_angvelz_free_flight[position_index][N_shellz]+1;
	  
	  atom_count_inst[position_index]=atom_count_inst[position_index]+1;
	  tot_atom_count_angvel_free_flight[position_index]=tot_atom_count_angvel_free_flight[position_index]+1;
   }
             }
	   }

 for(position_index=1;position_index<=max_position_index;position_index++) 
  for( kk=1;kk<=Ndel_angvelx;kk++)           // k should start from 1 not from 0  
                 {
		   if(atom_count_inst[position_index]>0)	  
      fof_angvelx_free_flight_sq[position_index][kk]=fof_angvelx_free_flight_sq[position_index][kk]+(double)(Nof_angvelx_inst[position_index][kk]/
(atom_count_inst[position_index]*del_angvelx))*(double)(Nof_angvelx_inst[position_index][kk]/(atom_count_inst[position_index]*del_angvelx));
		 }


 for(position_index=1;position_index<=max_position_index;position_index++) 
	   for( kk=1;kk<=Ndel_angvely;kk++)           // k should start from 1 not from 0  
                 {
	if(atom_count_inst[position_index]>0)	   
   fof_angvely_free_flight_sq[position_index][kk]=fof_angvely_free_flight_sq[position_index][kk]+(double)( Nof_angvely_inst[position_index][kk]/
(atom_count_inst[position_index]*del_angvely))*(double)( Nof_angvely_inst[position_index][kk]/(atom_count_inst[position_index]*del_angvely));
		 }


 for(position_index=1;position_index<=max_position_index;position_index++) 
	   for( kk=1;kk<=Ndel_angvelz;kk++)           // k should start from 1 not from 0  
                 {
		if(atom_count_inst[position_index]>0)   
   fof_angvelz_free_flight_sq[position_index][kk]=fof_angvelz_free_flight_sq[position_index][kk]+ (double)( Nof_angvelz_inst[position_index][kk]/
(atom_count_inst[position_index]*del_angvelz))*(double)( Nof_angvelz_inst[position_index][kk]/(atom_count_inst[position_index]*del_angvelz));
		 }

	   return;
     }



 void angular_vel_distri_func_free_flight()
{ 
              int k; 
              double fnorm;
	      double angvelx, fof_angvelx,angvely,fof_angvely,angvelz,fof_angvelz;
	      double var_fof_angvelx,var_fof_angvely,var_fof_angvelz;

	     

	      FILE *fp1040;
	      FILE *fp10401;
	      FILE *fp10402;

	      FILE *fp1041;
	      FILE *fp10411;
	      FILE *fp10412;


	      FILE *fp1042;
	      FILE *fp10421;
	      FILE *fp10422;
	    
 fp1040=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_free_flight_center_x.txt","w");
 fp10401=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_free_flight_middle_x.txt","w");
 fp10402=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_free_flight_wall_x.txt","w");

  fp1041=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_free_flight_center_y.txt","w");
 fp10411=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_free_flight_middle_y.txt","w");
 fp10412=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_free_flight_wall_y.txt","w");

  fp1042=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_free_flight_center_z.txt","w");
  fp10421=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_free_flight_middle_z.txt","w");
  fp10422=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angular_velocity_dist_free_flight_wall_z.txt","w");

   fnorm =(double)(tot_atom_count_vel_free_flight[1]);                  
               

	    for(k=1;k<=Ndel_angvelx;k++)           // k should start from 1 not from 0-- this is right only 
                 {
		   angvelx= -p_angvelx_max+ del_angvelx*(k-1);

		   fof_angvelx= (double)( Nof_angvelx_free_flight[1][k])/((double)(tot_atom_count_angvel_free_flight[1])*del_angvelx);
		   var_fof_angvelx=(fof_angvelx_free_flight_sq[1][k]/(double)n_sample_angvel_dist_free_flight)-(fof_angvelx*fof_angvelx);
		   fprintf(fp1040,"%20.14lf %20.14lf %20.14lf\n",angvelx,fof_angvelx, var_fof_angvelx);

		   fof_angvelx= (double)( Nof_angvelx_free_flight[2][k])/((double)(tot_atom_count_angvel_free_flight[2])*del_angvelx);
		   var_fof_angvelx=(fof_angvelx_free_flight_sq[2][k]/(double)n_sample_angvel_dist_free_flight)-(fof_angvelx*fof_angvelx);
		   fprintf(fp10401,"%20.14lf %20.14lf %20.14lf\n",angvelx,fof_angvelx, var_fof_angvelx);

		   fof_angvelx= (double)( Nof_angvelx_free_flight[3][k])/((double)(tot_atom_count_angvel_free_flight[3])*del_angvelx);
		   var_fof_angvelx=(fof_angvelx_free_flight_sq[3][k]/(double)n_sample_angvel_dist_free_flight)-(fof_angvelx*fof_angvelx);
		   fprintf(fp10402,"%20.14lf %20.14lf %20.14lf\n",angvelx,fof_angvelx, var_fof_angvelx);

		 }
        

	   for(k=1;k<=Ndel_angvely;k++)                   // k should start from 1 not from 0  
		 {
		   angvely= -p_angvely_max+ del_angvely*(k-1);

		   fof_angvely= (double)(Nof_angvely_free_flight[1][k])/((double)(tot_atom_count_angvel_free_flight[1])*del_angvely);
		   var_fof_angvely=(fof_angvely_free_flight_sq[1][k]/(double)n_sample_angvel_dist_free_flight)-(fof_angvely*fof_angvely);
	    	   fprintf(fp1041,"%20.14lf %20.14lf %20.14lf\n",angvely,fof_angvely,var_fof_angvely);

		   fof_angvely= (double)(Nof_angvely_free_flight[2][k])/((double)(tot_atom_count_angvel_free_flight[2])*del_angvely);
		   var_fof_angvely=(fof_angvely_free_flight_sq[2][k]/(double)n_sample_angvel_dist_free_flight)-(fof_angvely*fof_angvely);
	    	   fprintf(fp10411,"%20.14lf %20.14lf %20.14lf\n",angvely,fof_angvely,var_fof_angvely);

		   fof_angvely= (double)(Nof_angvely_free_flight[3][k])/((double)(tot_atom_count_angvel_free_flight[3])*del_angvely);
		   var_fof_angvely=(fof_angvely_free_flight_sq[3][k]/(double)n_sample_angvel_dist_free_flight)-(fof_angvely*fof_angvely);
	    	   fprintf(fp10412,"%20.14lf %20.14lf %20.14lf\n",angvely,fof_angvely,var_fof_angvely);
		   
		 }



	   for(k=1;k<=Ndel_angvelz;k++)
	     {
	       angvelz= -p_angvelz_max+ del_angvelz*(k-1);

	       fof_angvelz= (double)(Nof_angvelz_free_flight[1][k])/((double)tot_atom_count_angvel_free_flight[1]*del_angvelz);
	       var_fof_angvelz=(fof_angvelz_free_flight_sq[1][k]/ (double)n_sample_angvel_dist_free_flight)-(fof_angvelz*fof_angvelz);
	       fprintf(fp1042 ,"%20.14lf %20.14lf %20.14lf\n",angvelz,fof_angvelz,var_fof_angvelz);

	       fof_angvelz= (double)(Nof_angvelz_free_flight[2][k])/((double)tot_atom_count_angvel_free_flight[2]*del_angvelz);
	       var_fof_angvelz=(fof_angvelz_free_flight_sq[2][k]/(double)n_sample_angvel_dist_free_flight)-(fof_angvelz*fof_angvelz);
	       fprintf(fp10421 ,"%20.14lf %20.14lf %20.14lf\n",angvelz,fof_angvelz,var_fof_angvelz);

	       fof_angvelz= (double)(Nof_angvelz_free_flight[3][k])/((double)tot_atom_count_angvel_free_flight[3]*del_angvelz);
	       var_fof_angvelz=(fof_angvelz_free_flight_sq[3][k]/(double)n_sample_angvel_dist_free_flight)-(fof_angvelz*fof_angvelz);
	       fprintf(fp10422 ,"%20.14lf %20.14lf %20.14lf\n",angvelz,fof_angvelz,var_fof_angvelz);

	     }
       
     }


void  acc_property()
     {   
       int ny,kk,jj;      
       int i,N_shell_accx, N_shell_accy,N_shell_accz,N_shell_accxy;
       double  offset_accx, offset_accy , offset_accz ,offset_accxy;
       
       int *Nofp;
       double *accx_inst,*accy_inst,*accz_inst;
        double *accx_inst_sq, *accy_inst_sq, *accz_inst_sq;
      

       int atom_acc_inst_dist[Max_dim_position];

       n_sample_acc_dist= n_sample_acc_dist+1;


       
       Nofp=new int[Ndel_y];

       accx_inst=new double[Ndel_y];   
       accy_inst=new double[Ndel_y];
       accz_inst=new double[Ndel_y];

        accx_inst_sq=new double[Ndel_y];
	accy_inst_sq=new double[Ndel_y];
	accz_inst_sq=new double[Ndel_y];
       


       accx_dash=new double[n_atom];
       accy_dash=new double[n_atom];
       accz_dash=new double[n_atom];
       
       
   for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
    Nofp[ny]=0;
    accx_inst[ny]=0.0;
    accy_inst[ny]=0.0;
    accz_inst[ny]=0.0;

    accx_inst_sq[ny]=0.0;
    accy_inst_sq[ny]=0.0;
    accz_inst_sq[ny]=0.0;

       }

 for(i=0;i<n_atom;i++)// summing up all the acceleration in each the y-grid for all the particle 
    {
      ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

      // cout<<"in property"<<p_y[i]<<" "<<ny<<endl;


              Nofp[ny]= Nofp[ny]+ 1;

              accx_inst[ny]= accx_inst[ny]+acc_x[i];
              accy_inst[ny]= accy_inst[ny]+acc_y[i];  
              accz_inst[ny]= accz_inst[ny]+acc_z[i];
      }


for(i=0;i<n_atom;i++)// summing up all the acceleration in each the y-grid for all the particle 
    {
      ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

              accx_inst_sq[ny]= accx_inst_sq[ny]+acc_x[i]*acc_x[i];
              accy_inst_sq[ny]= accy_inst_sq[ny]+acc_y[i]*acc_y[i];  
              accz_inst_sq[ny]= accz_inst_sq[ny]+acc_z[i]*acc_z[i];
      }


for(ny=0;ny<Ndel_y;ny++) // doing the particle average in each y-grid-in each sample frame
       {
    if(Nofp[ny]>0.0)
      {
	accx_inst_sq[ny]= accx_inst_sq[ny]/Nofp[ny];
	accy_inst_sq[ny]= accy_inst_sq[ny]/Nofp[ny];
	accz_inst_sq[ny]= accz_inst_sq[ny]/Nofp[ny];
      }

    mean_accx_sq[ny]=mean_accx_sq[ny]+accx_inst_sq[ny];
    mean_accy_sq[ny]=mean_accy_sq[ny]+accy_inst_sq[ny];
    mean_accz_sq[ny]=mean_accz_sq[ny]+accz_inst_sq[ny];

       }
 
for(ny=0;ny<Ndel_y;ny++) // doing the particle average in each y-grid-in each sample frame
       {
    if(Nofp[ny]>0.0)
      {
     accx_inst[ny]= accx_inst[ny]/Nofp[ny];  // accx_inst[ny] is basically particle averaged inst in time.
     accy_inst[ny]= accy_inst[ny]/Nofp[ny];
     accz_inst[ny]= accz_inst[ny]/Nofp[ny];
      }
    
    accx_mean[ny]=accx_mean[ny]+accx_inst[ny];
    accy_mean[ny]=accy_mean[ny]+accy_inst[ny];
    accz_mean[ny]=accz_mean[ny]+accz_inst[ny];
       }



 for(i=0;i<n_atom;i++)
   {
 ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
 accx_dash[i]=acc_x[i]-accx_inst[ny];                // calculating fluctuation for each particle
 accy_dash[i]=acc_y[i]-accy_inst[ny];
 accz_dash[i]=acc_z[i]-accz_inst[ny];
 
    }




 /*----------------------------------------------------------------
Calculation of the 2 nd 3 rd and the 4 th moment of the acceleration
-------------------------------------------------------------------------*/



 for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
   
    accx_dash_2[ny]=0.0;
    accy_dash_2[ny]=0.0;
    accz_dash_2[ny]=0.0;

    accx_dash_3[ny]=0.0;
    accy_dash_3[ny]=0.0;
    accz_dash_3[ny]=0.0;

    accx_dash_4[ny]=0.0;
    accy_dash_4[ny]=0.0;
    accz_dash_4[ny]=0.0;
     
        }


 for(i=0;i<n_atom;i++)
	{
	  ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
	  
	   accx_dash_2[ny]=accx_dash_2[ny]+accx_dash[i]*accx_dash[i];
	   accy_dash_2[ny]=accy_dash_2[ny]+accy_dash[i]*accy_dash[i];
	   accz_dash_2[ny]=accz_dash_2[ny]+accz_dash[i]*accz_dash[i];

	   accx_dash_3[ny]=accx_dash_3[ny]+accx_dash[i]*accx_dash[i]*accx_dash[i];
	   accy_dash_3[ny]=accy_dash_3[ny]+accy_dash[i]*accy_dash[i]*accy_dash[i];
	   accz_dash_3[ny]=accz_dash_3[ny]+accz_dash[i]*accz_dash[i]*accz_dash[i];

	   accx_dash_4[ny]=accx_dash_4[ny]+accx_dash[i]*accx_dash[i]*accx_dash[i]*accx_dash[i];
	   accy_dash_4[ny]=accy_dash_4[ny]+accy_dash[i]*accy_dash[i]*accy_dash[i]*accy_dash[i];
	   accz_dash_4[ny]=accz_dash_4[ny]+accz_dash[i]*accz_dash[i]*accz_dash[i]*accz_dash[i];

	     }
     
 for(ny=0;ny<Ndel_y;ny++) 
	{
        if(Nofp[ny]>0.0)
	  {
	    accx_dash_2[ny]=accx_dash_2[ny]/Nofp[ny];
	    accy_dash_2[ny]=accy_dash_2[ny]/Nofp[ny];
	    accz_dash_2[ny]= accz_dash_2[ny]/Nofp[ny];

	    accx_dash_3[ny]=accx_dash_3[ny]/Nofp[ny];
	    accy_dash_3[ny]=accy_dash_3[ny]/Nofp[ny];
	    accz_dash_3[ny]=accz_dash_3[ny]/Nofp[ny];

	    accx_dash_4[ny]=accx_dash_4[ny]/Nofp[ny];
	    accy_dash_4[ny]=accy_dash_4[ny]/Nofp[ny];
	    accz_dash_4[ny]=accz_dash_4[ny]/Nofp[ny];

	  }

	avg_accx_dash_2[ny]=avg_accx_dash_2[ny]+accx_dash_2[ny];
	avg_accy_dash_2[ny]=avg_accy_dash_2[ny]+accy_dash_2[ny];
	avg_accz_dash_2[ny]=avg_accz_dash_2[ny]+accz_dash_2[ny];

	avg_accx_dash_3[ny]=avg_accx_dash_3[ny]+accx_dash_3[ny];
	avg_accy_dash_3[ny]=avg_accy_dash_3[ny]+accy_dash_3[ny];
	avg_accz_dash_3[ny]=avg_accz_dash_3[ny]+accz_dash_3[ny];

	avg_accx_dash_4[ny]=avg_accx_dash_4[ny]+accx_dash_4[ny];
	avg_accy_dash_4[ny]=avg_accy_dash_4[ny]+accy_dash_4[ny];
	avg_accz_dash_4[ny]=avg_accz_dash_4[ny]+accz_dash_4[ny];

	}

 /*-------------calculation starts for distribution function----------------*/
       
  
            offset_accx =1.5+ accx_max/del_accx;
            offset_accy =1.5+ accy_max/del_accy;
            offset_accz =1.5+ accz_max/del_accz;
            offset_accxy=1.5+ accxy_max/del_accxy;

         

 for(jj=0;jj<Max_dim_position;jj++) 
 for(i=0;i<=Ndel_accx;i++) Nof_accx_inst[jj][i]=0;

 for(jj=0;jj<Max_dim_position;jj++)
 for(i=0;i<=Ndel_accy;i++) Nof_accy_inst[jj][i]=0;

 for(jj=0;jj<Max_dim_position;jj++)      
 for(i=0;i<=Ndel_accz;i++) Nof_accz_inst[jj][i]=0;

 for(jj=0;jj<Max_dim_position;jj++)
 for(i=0;i<=Ndel_accxy;i++) Nof_accxy_inst[jj][i]=0;

 for(jj=0;jj<Max_dim_position;jj++)atom_acc_inst_dist[jj]=0;


 for(i=0;i<n_atom;i++)
           {
	     position_index=0;    //initialization of position_index

	     if(p_y[i]>y_min_for_distribution_1 && p_y[i]<y_max_for_distribution_1)position_index=1;
	     if(p_y[i]>y_min_for_distribution_2 && p_y[i]<y_max_for_distribution_2)position_index=2;
	     if(p_y[i]>y_min_for_distribution_3 && p_y[i]<y_max_for_distribution_3)position_index=3;
	     

	     if(position_index>0)
	     
	       {  
 if((accx_dash[i]>=-accx_max && accx_dash[i]<=accx_max) && (accy_dash[i]>=-accy_max && accy_dash[i]<=accy_max)
&& (accz_dash[i]>=-accz_max && accz_dash[i]<=accz_max)&& (accx_dash[i]*accy_dash[i]>=-accxy_max && accx_dash[i]*accy_dash[i]<=accxy_max))  // to avoid exceeding arrey subscript accidentally
   {
		 N_shell_accx=(int)(accx_dash[i]/del_accx+offset_accx);
                 Nof_accx[position_index][N_shell_accx]= Nof_accx[position_index][N_shell_accx] +1 ;
		 Nof_accx_inst[position_index][N_shell_accx]= Nof_accx_inst[position_index][N_shell_accx] +1 ;

		 N_shell_accy=(int)(accy_dash[i]/del_accy + offset_accy);
                 Nof_accy[position_index][N_shell_accy]= Nof_accy[position_index][N_shell_accy] +1 ; 
		 Nof_accy_inst[position_index][N_shell_accy]= Nof_accy_inst[position_index][N_shell_accy] +1 ; 

		 N_shell_accz=(int)(accz_dash[i]/del_accz + offset_accz);
                 Nof_accz[position_index][N_shell_accz]= Nof_accz[position_index][N_shell_accz] +1 ;
		 Nof_accz_inst[position_index][N_shell_accz]= Nof_accz_inst[position_index][N_shell_accz] +1 ;
		 
		 N_shell_accxy=(int)(accx_dash[i]*accy_dash[i]/del_accxy + offset_accxy);
                 Nof_accxy[position_index][N_shell_accxy]= Nof_accxy[position_index][N_shell_accxy] +1 ; 
		 Nof_accxy_inst[position_index][N_shell_accxy]= Nof_accxy_inst[position_index][N_shell_accxy] +1 ; 

		 atom_acc_inst_dist[position_index]= atom_acc_inst_dist[position_index]+1;
	   tot_atom_acc_dist[position_index]= tot_atom_acc_dist[position_index] +1;
   }
	       }

	   }

  for(position_index=1;position_index<=max_position_index;position_index++) 
 for( kk=1;kk<=Ndel_accx;kk++) 
   { 
     fof_accx_sq[position_index][kk]=fof_accx_sq[position_index][kk]+(double)(Nof_accx_inst[position_index][kk]/
(atom_acc_inst_dist[position_index]*del_accx))*(double)(Nof_accx_inst[position_index][kk]/(atom_acc_inst_dist[position_index]*del_accx));
     // if(kk==173){cout<<"++++++++++++++++++++++++   "<<fof_accx_sq[kk]<<"    "<<" "<<atom_acc_inst_dist<<" "<<" "<<del_accx<<" "<<Nof_accx[kk]<<endl;
     //  cin.get();}

       }

 for(position_index=1;position_index<=max_position_index;position_index++) 
 for( kk=1;kk<=Ndel_accy;kk++) 
   { 
     fof_accy_sq[position_index][kk]=fof_accy_sq[position_index][kk]+(double)(Nof_accy_inst[position_index][kk]/
(atom_acc_inst_dist[position_index]*del_accy))*(double)(Nof_accy_inst[position_index][kk]/(atom_acc_inst_dist[position_index]*del_accy));
       }

 for(position_index=1;position_index<=max_position_index;position_index++) 
 for( kk=1;kk<=Ndel_accz;kk++) 
   { 
     fof_accz_sq[position_index][kk]=fof_accz_sq[position_index][kk]+(double)(Nof_accz_inst[position_index][kk]/
(atom_acc_inst_dist[position_index]*del_accz))*(double)(Nof_accz_inst[position_index][kk]/(atom_acc_inst_dist[position_index]*del_accz));
       }

 for(position_index=1;position_index<=max_position_index;position_index++) 
 for( kk=1;kk<=Ndel_accxy;kk++) 
   { 
     fof_accxy_sq[position_index][kk]=fof_accxy_sq[position_index][kk]+(double)(Nof_accxy_inst[position_index][kk]/
(atom_acc_inst_dist[position_index]*del_accxy))*(double)(Nof_accxy_inst[position_index][kk]/(atom_acc_inst_dist[position_index]*del_accxy));
       }



	       delete [] Nofp;

	       delete []accx_inst ; 
	       delete []accy_inst ;
	       delete []accz_inst ;

	       delete [] accx_dash ;
	       delete [] accy_dash ;
	       delete [] accz_dash ;

	        		


     }



       void acc_distri_func() 
          {
              int k ; 
              double fnorm;
	      double accx, fof_accx,accy,fof_accy,accz,fof_accz,accxy,fof_accxy;
	      double  var_fof_accx,var_fof_accy,var_fof_accz,var_fof_accxy;

	      FILE *fp11;
	      FILE *fp111;
	       FILE *fp112;

	      FILE *fp12;
	      FILE *fp121;
	      FILE *fp122;

	      FILE *fp13;
	      FILE *fp131;
	      FILE *fp132;

	      FILE *fp14;
	      FILE *fp141;
	      FILE *fp142;



	      
    fp11=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_center_x.txt","w");
    fp111=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_middle_x.txt","w");
    fp112=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_wall_x.txt","w");

    fp12=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_center_y.txt","w");
    fp121=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_middle_y.txt","w");
    fp122=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_wall_y.txt","w");


    fp13=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_center_z.txt","w");
    fp131=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_middle_z.txt","w");
    fp132=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_wall_z.txt","w");

    fp14=fopen("/home/ankit/Desktop/dns/channel/examples/volfrac_0.0004/acc_dist_center_xy.txt","w");
    fp141=fopen("/home/ankit/Desktop/dns/channel/examples/volfrac_0.0004/acc_dist_middle_xy.txt","w");
    fp142=fopen("/home/ankit/Desktop/dns/channel/examples/volfrac_0.0004/acc_dist_wall_xy.txt","w");	   
                                                    
     fnorm =(double)(tot_atom_acc_dist[1]);               
                                                   
      
//! Loop over veolocity shells each of thickness del_v

	  

	       for(k=1;k<=Ndel_accx;k++)           // k should start from 1 not from 0  
                 {
		   accx= -accx_max+ del_accx*(k-1);

		   fof_accx= (double)( Nof_accx[1][k])/((double)(tot_atom_acc_dist[1])*del_accx);
		   var_fof_accx=(fof_accx_sq[1][k]/(double)n_sample_acc_dist)- (fof_accx* fof_accx);
		   fprintf(fp11,"%20.14lf %20.14lf %20.14lf\n",accx,fof_accx,var_fof_accx);

		   fof_accx= (double)( Nof_accx[2][k])/((double)(tot_atom_acc_dist[2])*del_accx);
		   var_fof_accx=(fof_accx_sq[2][k]/(double)n_sample_acc_dist)- (fof_accx* fof_accx);
		   fprintf(fp111,"%20.14lf %20.14lf %20.14lf\n",accx,fof_accx,var_fof_accx);

		   fof_accx= (double)( Nof_accx[3][k])/((double)(tot_atom_acc_dist[3])*del_accx);
		   var_fof_accx=(fof_accx_sq[3][k]/(double)n_sample_acc_dist)- (fof_accx* fof_accx);
		   fprintf(fp112,"%20.14lf %20.14lf %20.14lf\n",accx,fof_accx,var_fof_accx);



		 }
        
	       for(k=1;k<=Ndel_accy;k++)                   // k should start from 1 not from 0  
		 {
		   accy= -accy_max+ del_accy*(k-1);

		   fof_accy= (double)(Nof_accy[1][k])/((double)(tot_atom_acc_dist[1])*del_accy);
		   var_fof_accy= (fof_accy_sq[1][k]/(double)n_sample_acc_dist)- (fof_accy* fof_accy);
		   fprintf(fp12,"%20.14lf %20.14lf  %20.14lf\n",accy,fof_accy,var_fof_accy);

		   fof_accy= (double)(Nof_accy[2][k])/((double)(tot_atom_acc_dist[2])*del_accy);
		   var_fof_accy= (fof_accy_sq[2][k]/(double)n_sample_acc_dist)- (fof_accy* fof_accy);
		   fprintf(fp121,"%20.14lf %20.14lf  %20.14lf\n",accy,fof_accy,var_fof_accy);

		   fof_accy= (double)(Nof_accy[3][k])/((double)(tot_atom_acc_dist[3])*del_accy);
		   var_fof_accy= (fof_accy_sq[3][k]/(double)n_sample_acc_dist)- (fof_accy* fof_accy);
		   fprintf(fp122,"%20.14lf %20.14lf  %20.14lf\n",accy,fof_accy,var_fof_accy);

		 }

	       for(k=1;k<=Ndel_accz;k++)
		 {
		   accz= -accz_max+ del_accz*(k-1);

		   fof_accz= (double)(Nof_accz[1][k])/((double)(tot_atom_acc_dist[1])*del_accz);
		   var_fof_accz= (fof_accz_sq[1][k]/(double)n_sample_acc_dist)- (fof_accz* fof_accz);
		   fprintf(fp13 ,"%20.14lf %20.14lf  %20.14lf\n",accz,fof_accz, var_fof_accz);

		   fof_accz= (double)(Nof_accz[2][k])/((double)(tot_atom_acc_dist[2])*del_accz);
		   var_fof_accz= (fof_accz_sq[2][k]/(double)n_sample_acc_dist)- (fof_accz* fof_accz);
		   fprintf(fp131 ,"%20.14lf %20.14lf  %20.14lf\n",accz,fof_accz, var_fof_accz);

		   fof_accz= (double)(Nof_accz[3][k])/((double)(tot_atom_acc_dist[3])*del_accz);
		   var_fof_accz= (fof_accz_sq[3][k]/(double)n_sample_acc_dist)- (fof_accz* fof_accz);
		   fprintf(fp132 ,"%20.14lf %20.14lf  %20.14lf\n",accz,fof_accz, var_fof_accz);
		   
		 }
	       
	       for(k=1;k<=Ndel_accxy;k++)
		 {
		   accxy= -accxy_max+ del_accxy*(k-1);

		   fof_accxy= (double)(Nof_accxy[1][k])/((double)(tot_atom_acc_dist[1])*del_accxy);
		   var_fof_accxy= (fof_accxy_sq[1][k]/(double)n_sample_acc_dist)- (fof_accxy* fof_accxy);
		   fprintf(fp14 ,"%20.14lf %20.14lf  %20.14lf\n",accxy,fof_accxy,var_fof_accxy);
	 
		   fof_accxy= (double)(Nof_accxy[2][k])/((double)(tot_atom_acc_dist[2])*del_accxy);
		   var_fof_accxy= (fof_accxy_sq[2][k]/(double)n_sample_acc_dist)- (fof_accxy* fof_accxy);
		   fprintf(fp141 ,"%20.14lf %20.14lf  %20.14lf\n",accxy,fof_accxy,var_fof_accxy);
		   
		   fof_accxy= (double)(Nof_accxy[3][k])/((double)(tot_atom_acc_dist[3])*del_accxy);
		   var_fof_accxy= (fof_accxy_sq[3][k]/(double)n_sample_acc_dist)- (fof_accxy* fof_accxy);
		   fprintf(fp142 ,"%20.14lf %20.14lf %20.14lf",accxy,fof_accxy,var_fof_accxy);
		   
		 
		 }

	       return;
	  }


void  angular_acc_property()
     {   
       int ny,kk,jj;      
       int i,N_shell_angaccx, N_shell_angaccy,N_shell_angaccz,N_shell_angaccxy;
       double  offset_angaccx, offset_angaccy , offset_angaccz ,offset_angaccxy;
       
       int *Nofp;
       double *angaccx_inst,*angaccy_inst,*angaccz_inst;
        double *angaccx_inst_sq, *angaccy_inst_sq, *angaccz_inst_sq;
      

       int atom_angacc_inst_dist[Max_dim_position];

       n_sample_angacc_dist= n_sample_angacc_dist+1;


       
       Nofp=new int[Ndel_y];

       angaccx_inst=new double[Ndel_y];   
       angaccy_inst=new double[Ndel_y];
       angaccz_inst=new double[Ndel_y];

        angaccx_inst_sq=new double[Ndel_y];
	angaccy_inst_sq=new double[Ndel_y];
	angaccz_inst_sq=new double[Ndel_y];
       


       angaccx_dash=new double[n_atom];
       angaccy_dash=new double[n_atom];
       angaccz_dash=new double[n_atom];
       
       
   for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
    Nofp[ny]=0;
    angaccx_inst[ny]=0.0;
    angaccy_inst[ny]=0.0;
    angaccz_inst[ny]=0.0;

    angaccx_inst_sq[ny]=0.0;
    angaccy_inst_sq[ny]=0.0;
    angaccz_inst_sq[ny]=0.0;

       }

 for(i=0;i<n_atom;i++)// summing up all the acceleration in each the y-grid for all the particle 
    {
      ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

      // cout<<"in property"<<p_y[i]<<" "<<ny<<endl;


              Nofp[ny]= Nofp[ny]+ 1;

              angaccx_inst[ny]= angaccx_inst[ny]+angaccx[i];
              angaccy_inst[ny]= angaccy_inst[ny]+angaccy[i];  
              angaccz_inst[ny]= angaccz_inst[ny]+angaccz[i];
      }


for(i=0;i<n_atom;i++)// summing up all the acceleration in each the y-grid for all the particle 
    {
      ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

              angaccx_inst_sq[ny]= angaccx_inst_sq[ny]+angaccx[i]*angaccx[i];
              angaccy_inst_sq[ny]= angaccy_inst_sq[ny]+angaccy[i]*angaccy[i];  
              angaccz_inst_sq[ny]= angaccz_inst_sq[ny]+angaccz[i]*angaccz[i];
      }


for(ny=0;ny<Ndel_y;ny++) // doing the particle average in each y-grid-in each sample frame
       {
    if(Nofp[ny]>0.0)
      {
	angaccx_inst_sq[ny]= angaccx_inst_sq[ny]/Nofp[ny];
	angaccy_inst_sq[ny]= angaccy_inst_sq[ny]/Nofp[ny];
	angaccz_inst_sq[ny]= angaccz_inst_sq[ny]/Nofp[ny];
      }

    mean_angaccx_sq[ny]=mean_angaccx_sq[ny]+angaccx_inst_sq[ny];
    mean_angaccy_sq[ny]=mean_angaccy_sq[ny]+angaccy_inst_sq[ny];
    mean_angaccz_sq[ny]=mean_angaccz_sq[ny]+angaccz_inst_sq[ny];

       }
 
for(ny=0;ny<Ndel_y;ny++) // doing the particle average in each y-grid-in each sample frame
       {
    if(Nofp[ny]>0.0)
      {
     angaccx_inst[ny]= angaccx_inst[ny]/Nofp[ny];  // angaccx_inst[ny] is basically particle averaged inst in time.
     angaccy_inst[ny]= angaccy_inst[ny]/Nofp[ny];
     angaccz_inst[ny]= angaccz_inst[ny]/Nofp[ny];
      }
    
    angaccx_mean[ny]=angaccx_mean[ny]+angaccx_inst[ny];
    angaccy_mean[ny]=angaccy_mean[ny]+angaccy_inst[ny];
    angaccz_mean[ny]=angaccz_mean[ny]+angaccz_inst[ny];
       }



 for(i=0;i<n_atom;i++)
   {
 ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
 angaccx_dash[i]=angaccx[i]-angaccx_inst[ny];                // calculating fluctuation for each particle
 angaccy_dash[i]=angaccy[i]-angaccy_inst[ny];
 angaccz_dash[i]=angaccz[i]-angaccz_inst[ny];
 
    }




 /*----------------------------------------------------------------
Calculation of the 2 nd 3 rd and the 4 th moment of the angacceleration
-------------------------------------------------------------------------*/



 for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
   
    angaccx_dash_2[ny]=0.0;
    angaccy_dash_2[ny]=0.0;
    angaccz_dash_2[ny]=0.0;

    angaccx_dash_3[ny]=0.0;
    angaccy_dash_3[ny]=0.0;
    angaccz_dash_3[ny]=0.0;

    angaccx_dash_4[ny]=0.0;
    angaccy_dash_4[ny]=0.0;
    angaccz_dash_4[ny]=0.0;
     
        }


 for(i=0;i<n_atom;i++)
	{
	  ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
	  
	   angaccx_dash_2[ny]=angaccx_dash_2[ny]+angaccx_dash[i]*angaccx_dash[i];
	   angaccy_dash_2[ny]=angaccy_dash_2[ny]+angaccy_dash[i]*angaccy_dash[i];
	   angaccz_dash_2[ny]=angaccz_dash_2[ny]+angaccz_dash[i]*angaccz_dash[i];

	   angaccx_dash_3[ny]=angaccx_dash_3[ny]+angaccx_dash[i]*angaccx_dash[i]*angaccx_dash[i];
	   angaccy_dash_3[ny]=angaccy_dash_3[ny]+angaccy_dash[i]*angaccy_dash[i]*angaccy_dash[i];
	   angaccz_dash_3[ny]=angaccz_dash_3[ny]+angaccz_dash[i]*angaccz_dash[i]*angaccz_dash[i];

	   angaccx_dash_4[ny]=angaccx_dash_4[ny]+angaccx_dash[i]*angaccx_dash[i]*angaccx_dash[i]*angaccx_dash[i];
	   angaccy_dash_4[ny]=angaccy_dash_4[ny]+angaccy_dash[i]*angaccy_dash[i]*angaccy_dash[i]*angaccy_dash[i];
	   angaccz_dash_4[ny]=angaccz_dash_4[ny]+angaccz_dash[i]*angaccz_dash[i]*angaccz_dash[i]*angaccz_dash[i];

	     }
     
 for(ny=0;ny<Ndel_y;ny++) 
	{
        if(Nofp[ny]>0.0)
	  {
	    angaccx_dash_2[ny]=angaccx_dash_2[ny]/Nofp[ny];
	    angaccy_dash_2[ny]=angaccy_dash_2[ny]/Nofp[ny];
	    angaccz_dash_2[ny]= angaccz_dash_2[ny]/Nofp[ny];

	    angaccx_dash_3[ny]=angaccx_dash_3[ny]/Nofp[ny];
	    angaccy_dash_3[ny]=angaccy_dash_3[ny]/Nofp[ny];
	    angaccz_dash_3[ny]=angaccz_dash_3[ny]/Nofp[ny];

	    angaccx_dash_4[ny]=angaccx_dash_4[ny]/Nofp[ny];
	    angaccy_dash_4[ny]=angaccy_dash_4[ny]/Nofp[ny];
	    angaccz_dash_4[ny]=angaccz_dash_4[ny]/Nofp[ny];

	  }

	avg_angaccx_dash_2[ny]=avg_angaccx_dash_2[ny]+angaccx_dash_2[ny];
	avg_angaccy_dash_2[ny]=avg_angaccy_dash_2[ny]+angaccy_dash_2[ny];
	avg_angaccz_dash_2[ny]=avg_angaccz_dash_2[ny]+angaccz_dash_2[ny];

	avg_angaccx_dash_3[ny]=avg_angaccx_dash_3[ny]+angaccx_dash_3[ny];
	avg_angaccy_dash_3[ny]=avg_angaccy_dash_3[ny]+angaccy_dash_3[ny];
	avg_angaccz_dash_3[ny]=avg_angaccz_dash_3[ny]+angaccz_dash_3[ny];

	avg_angaccx_dash_4[ny]=avg_angaccx_dash_4[ny]+angaccx_dash_4[ny];
	avg_angaccy_dash_4[ny]=avg_angaccy_dash_4[ny]+angaccy_dash_4[ny];
	avg_angaccz_dash_4[ny]=avg_angaccz_dash_4[ny]+angaccz_dash_4[ny];

	}

 /*-------------calculation starts for distribution function----------------*/
       
  
            offset_angaccx =1.5+ angaccx_max/del_angaccx;
            offset_angaccy =1.5+ angaccy_max/del_angaccy;
            offset_angaccz =1.5+ angaccz_max/del_angaccz;
            offset_angaccxy=1.5+ angaccxy_max/del_angaccxy;

         

 for(jj=0;jj<Max_dim_position;jj++) 
 for(i=0;i<=Ndel_angaccx;i++) Nof_angaccx_inst[jj][i]=0;

 for(jj=0;jj<Max_dim_position;jj++)
 for(i=0;i<=Ndel_angaccy;i++) Nof_angaccy_inst[jj][i]=0;

 for(jj=0;jj<Max_dim_position;jj++)      
 for(i=0;i<=Ndel_angaccz;i++) Nof_angaccz_inst[jj][i]=0;

 for(jj=0;jj<Max_dim_position;jj++)
 for(i=0;i<=Ndel_angaccxy;i++) Nof_angaccxy_inst[jj][i]=0;

 for(jj=0;jj<Max_dim_position;jj++)atom_angacc_inst_dist[jj]=0;


 for(i=0;i<n_atom;i++)
           {
	     position_index=0;    //initialization of position_index

	     if(p_y[i]>y_min_for_distribution_1 && p_y[i]<y_max_for_distribution_1)position_index=1;
	     if(p_y[i]>y_min_for_distribution_2 && p_y[i]<y_max_for_distribution_2)position_index=2;
	     if(p_y[i]>y_min_for_distribution_3 && p_y[i]<y_max_for_distribution_3)position_index=3;
	     

	     if(position_index>0)
	     
	       {  
 if((angaccx_dash[i]>=-angaccx_max && angaccx_dash[i]<=angaccx_max) && (angaccy_dash[i]>=-angaccy_max && angaccy_dash[i]<=angaccy_max)
&& (angaccz_dash[i]>=-angaccz_max && angaccz_dash[i]<=angaccz_max)&& (angaccx_dash[i]*angaccy_dash[i]>=-angaccxy_max && angaccx_dash[i]*angaccy_dash[i]<=angaccxy_max))  // to avoid exceeding arrey subscript accidentally
   {
		 N_shell_angaccx=(int)(angaccx_dash[i]/del_angaccx+offset_angaccx);
                 Nof_angaccx[position_index][N_shell_angaccx]= Nof_angaccx[position_index][N_shell_angaccx] +1 ;
		 Nof_angaccx_inst[position_index][N_shell_angaccx]= Nof_angaccx_inst[position_index][N_shell_angaccx] +1 ;

		 N_shell_angaccy=(int)(angaccy_dash[i]/del_angaccy + offset_angaccy);
                 Nof_angaccy[position_index][N_shell_angaccy]= Nof_angaccy[position_index][N_shell_angaccy] +1 ; 
		 Nof_angaccy_inst[position_index][N_shell_angaccy]= Nof_angaccy_inst[position_index][N_shell_angaccy] +1 ; 

		 N_shell_angaccz=(int)(angaccz_dash[i]/del_angaccz + offset_angaccz);
                 Nof_angaccz[position_index][N_shell_angaccz]= Nof_angaccz[position_index][N_shell_angaccz] +1 ;
		 Nof_angaccz_inst[position_index][N_shell_angaccz]= Nof_angaccz_inst[position_index][N_shell_angaccz] +1 ;
		 
		 N_shell_angaccxy=(int)(angaccx_dash[i]*angaccy_dash[i]/del_angaccxy + offset_angaccxy);
                 Nof_angaccxy[position_index][N_shell_angaccxy]= Nof_angaccxy[position_index][N_shell_angaccxy] +1 ; 
		 Nof_angaccxy_inst[position_index][N_shell_angaccxy]= Nof_angaccxy_inst[position_index][N_shell_angaccxy] +1 ; 

		 atom_angacc_inst_dist[position_index]= atom_angacc_inst_dist[position_index]+1;
	   tot_atom_angacc_dist[position_index]= tot_atom_angacc_dist[position_index] +1;
   }
	       }

	   }

  for(position_index=1;position_index<=max_position_index;position_index++) 
 for( kk=1;kk<=Ndel_angaccx;kk++) 
   { 
     fof_angaccx_sq[position_index][kk]=fof_angaccx_sq[position_index][kk]+(double)(Nof_angaccx_inst[position_index][kk]/
(atom_angacc_inst_dist[position_index]*del_angaccx))*(double)(Nof_angaccx_inst[position_index][kk]/(atom_angacc_inst_dist[position_index]*del_angaccx));
     // if(kk==173){cout<<"++++++++++++++++++++++++   "<<fof_angaccx_sq[kk]<<"    "<<" "<<atom_angacc_inst_dist<<" "<<" "<<del_angaccx<<" "<<Nof_angaccx[kk]<<endl;
     //  cin.get();}

       }

 for(position_index=1;position_index<=max_position_index;position_index++) 
 for( kk=1;kk<=Ndel_angaccy;kk++) 
   { 
     fof_angaccy_sq[position_index][kk]=fof_angaccy_sq[position_index][kk]+(double)(Nof_angaccy_inst[position_index][kk]/
(atom_angacc_inst_dist[position_index]*del_angaccy))*(double)(Nof_angaccy_inst[position_index][kk]/(atom_angacc_inst_dist[position_index]*del_angaccy));
       }

 for(position_index=1;position_index<=max_position_index;position_index++) 
 for( kk=1;kk<=Ndel_angaccz;kk++) 
   { 
     fof_angaccz_sq[position_index][kk]=fof_angaccz_sq[position_index][kk]+(double)(Nof_angaccz_inst[position_index][kk]/
(atom_angacc_inst_dist[position_index]*del_angaccz))*(double)(Nof_angaccz_inst[position_index][kk]/(atom_angacc_inst_dist[position_index]*del_angaccz));
       }

 for(position_index=1;position_index<=max_position_index;position_index++) 
 for( kk=1;kk<=Ndel_angaccxy;kk++) 
   { 
     fof_angaccxy_sq[position_index][kk]=fof_angaccxy_sq[position_index][kk]+(double)(Nof_angaccxy_inst[position_index][kk]/
(atom_angacc_inst_dist[position_index]*del_angaccxy))*(double)(Nof_angaccxy_inst[position_index][kk]/(atom_angacc_inst_dist[position_index]*del_angaccxy));
       }



	       delete [] Nofp;

	       delete []angaccx_inst ; 
	       delete []angaccy_inst ;
	       delete []angaccz_inst ;

	       delete [] angaccx_dash ;
	       delete [] angaccy_dash ;
	       delete [] angaccz_dash ;

	        		


     }



       void angular_acc_distri_func() 
          {
              int k ; 
              double fnorm;
	      double angaccx, fof_angaccx,angaccy,fof_angaccy,angaccz,fof_angaccz,angaccxy,fof_angaccxy;
	      double  var_fof_angaccx,var_fof_angaccy,var_fof_angaccz,var_fof_angaccxy;

	      FILE *fp1011;
	      FILE *fp10111;
	       FILE *fp10112;

	      FILE *fp1012;
	      FILE *fp10121;
	      FILE *fp10122;

	      FILE *fp1013;
	      FILE *fp10131;
	      FILE *fp10132;

	      FILE *fp1014;
	      FILE *fp10141;
	      FILE *fp10142;



	      
    fp1011=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_center_x.txt","w");
    fp10111=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_middle_x.txt","w");
    fp10112=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_wall_x.txt","w");

    fp1012=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_center_y.txt","w");
    fp10121=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_middle_y.txt","w");
    fp10122=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_wall_y.txt","w");


    fp1013=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_center_z.txt","w");
    fp10131=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_middle_z.txt","w");
    fp10132=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_wall_z.txt","w");

    fp1014=fopen("/home/ankit/Desktop/dns/channel/examples/volfrac_0.0004/angacc_dist_center_xy.txt","w");
    fp10141=fopen("/home/ankit/Desktop/dns/channel/examples/volfrac_0.0004/angacc_dist_middle_xy.txt","w");
    fp10142=fopen("/home/ankit/Desktop/dns/channel/examples/volfrac_0.0004/angacc_dist_wall_xy.txt","w");	   
                                                    
     fnorm =(double)(tot_atom_angacc_dist[1]);               
                                                   
      
//! Loop over veolocity shells each of thickness del_v

	  

	       for(k=1;k<=Ndel_angaccx;k++)           // k should start from 1 not from 0  
                 {
		   angaccx= -angaccx_max+ del_angaccx*(k-1);

		   fof_angaccx= (double)( Nof_angaccx[1][k])/((double)(tot_atom_angacc_dist[1])*del_angaccx);
		   var_fof_angaccx=(fof_angaccx_sq[1][k]/(double)n_sample_angacc_dist)- (fof_angaccx* fof_angaccx);
		   fprintf(fp1011,"%20.14lf %20.14lf %20.14lf\n",angaccx,fof_angaccx,var_fof_angaccx);

		   fof_angaccx= (double)( Nof_angaccx[2][k])/((double)(tot_atom_angacc_dist[2])*del_angaccx);
		   var_fof_angaccx=(fof_angaccx_sq[2][k]/(double)n_sample_angacc_dist)- (fof_angaccx* fof_angaccx);
		   fprintf(fp10111,"%20.14lf %20.14lf %20.14lf\n",angaccx,fof_angaccx,var_fof_angaccx);

		   fof_angaccx= (double)( Nof_angaccx[3][k])/((double)(tot_atom_angacc_dist[3])*del_angaccx);
		   var_fof_angaccx=(fof_angaccx_sq[3][k]/(double)n_sample_angacc_dist)- (fof_angaccx* fof_angaccx);
		   fprintf(fp10112,"%20.14lf %20.14lf %20.14lf\n",angaccx,fof_angaccx,var_fof_angaccx);



		 }
        
	       for(k=1;k<=Ndel_angaccy;k++)                   // k should start from 1 not from 0  
		 {
		   angaccy= -angaccy_max+ del_angaccy*(k-1);

		   fof_angaccy= (double)(Nof_angaccy[1][k])/((double)(tot_atom_angacc_dist[1])*del_angaccy);
		   var_fof_angaccy= (fof_angaccy_sq[1][k]/(double)n_sample_angacc_dist)- (fof_angaccy* fof_angaccy);
		   fprintf(fp1012,"%20.14lf %20.14lf  %20.14lf\n",angaccy,fof_angaccy,var_fof_angaccy);

		   fof_angaccy= (double)(Nof_angaccy[2][k])/((double)(tot_atom_angacc_dist[2])*del_angaccy);
		   var_fof_angaccy= (fof_angaccy_sq[2][k]/(double)n_sample_angacc_dist)- (fof_angaccy* fof_angaccy);
		   fprintf(fp10121,"%20.14lf %20.14lf  %20.14lf\n",angaccy,fof_angaccy,var_fof_angaccy);

		   fof_angaccy= (double)(Nof_angaccy[3][k])/((double)(tot_atom_angacc_dist[3])*del_angaccy);
		   var_fof_angaccy= (fof_angaccy_sq[3][k]/(double)n_sample_angacc_dist)- (fof_angaccy* fof_angaccy);
		   fprintf(fp10122,"%20.14lf %20.14lf  %20.14lf\n",angaccy,fof_angaccy,var_fof_angaccy);

		 }

	       for(k=1;k<=Ndel_angaccz;k++)
		 {
		   angaccz= -angaccz_max+ del_angaccz*(k-1);

		   fof_angaccz= (double)(Nof_angaccz[1][k])/((double)(tot_atom_angacc_dist[1])*del_angaccz);
		   var_fof_angaccz= (fof_angaccz_sq[1][k]/(double)n_sample_angacc_dist)- (fof_angaccz* fof_angaccz);
		   fprintf(fp1013 ,"%20.14lf %20.14lf  %20.14lf\n",angaccz,fof_angaccz, var_fof_angaccz);

		   fof_angaccz= (double)(Nof_angaccz[2][k])/((double)(tot_atom_angacc_dist[2])*del_angaccz);
		   var_fof_angaccz= (fof_angaccz_sq[2][k]/(double)n_sample_angacc_dist)- (fof_angaccz* fof_angaccz);
		   fprintf(fp10131 ,"%20.14lf %20.14lf  %20.14lf\n",angaccz,fof_angaccz, var_fof_angaccz);

		   fof_angaccz= (double)(Nof_angaccz[3][k])/((double)(tot_atom_angacc_dist[3])*del_angaccz);
		   var_fof_angaccz= (fof_angaccz_sq[3][k]/(double)n_sample_angacc_dist)- (fof_angaccz* fof_angaccz);
		   fprintf(fp10132 ,"%20.14lf %20.14lf  %20.14lf\n",angaccz,fof_angaccz, var_fof_angaccz);
		   
		 }
	       
	       for(k=1;k<=Ndel_angaccxy;k++)
		 {
		   angaccxy= -angaccxy_max+ del_angaccxy*(k-1);

		   fof_angaccxy= (double)(Nof_angaccxy[1][k])/((double)(tot_atom_angacc_dist[1])*del_angaccxy);
		   var_fof_angaccxy= (fof_angaccxy_sq[1][k]/(double)n_sample_angacc_dist)- (fof_angaccxy* fof_angaccxy);
		   fprintf(fp1014 ,"%20.14lf %20.14lf  %20.14lf\n",angaccxy,fof_angaccxy,var_fof_angaccxy);
	 
		   fof_angaccxy= (double)(Nof_angaccxy[2][k])/((double)(tot_atom_angacc_dist[2])*del_angaccxy);
		   var_fof_angaccxy= (fof_angaccxy_sq[2][k]/(double)n_sample_angacc_dist)- (fof_angaccxy* fof_angaccxy);
		   fprintf(fp10141 ,"%20.14lf %20.14lf  %20.14lf\n",angaccxy,fof_angaccxy,var_fof_angaccxy);
		   
		   fof_angaccxy= (double)(Nof_angaccxy[3][k])/((double)(tot_atom_angacc_dist[3])*del_angaccxy);
		   var_fof_angaccxy= (fof_angaccxy_sq[3][k]/(double)n_sample_angacc_dist)- (fof_angaccxy* fof_angaccxy);
		   fprintf(fp10142 ,"%20.14lf %20.14lf %20.14lf",angaccxy,fof_angaccxy,var_fof_angaccxy);
		   
		 
		 }

	       return;
	  }


 void air_vel_fluc_dist(const FlowField& u1,const Vector& x_grid,const Vector& z_grid, const Vector& y_grid)
  {
    int nx,ny,nz; 
    int N_shell_air_velx, N_shell_air_vely,N_shell_air_velz;

   
       double a_=u1.a();
       double b_=u1.b();
       double ubase;
       double air_velx_dash,air_vely_dash,air_velz_dash;
             
       double  offset_air_velx, offset_air_vely , offset_air_velz;


     for(nx=0;nx<u1.Nx();nx++)
	  for(nz=0;nz<u1.Nz();nz++)
	    for(ny=0;ny<u1.Ny();ny++)
	     {	      	  
	 ubase = (1.0 - square(abs(y_grid[ny]-(b_+a_)/2.0)/((b_-a_)/2.0))); /// chennel flow
	 air_velx_dash=u1(nx,ny,nz,0)+ubase-avg_vx_air[ny];
	 air_vely_dash=u1(nx,ny,nz,1);
	 air_velz_dash=u1(nx,ny,nz,2);

	 /*if(ny==30)
{cout<<" in output function part vel_fluc_air "<<  air_velx_dash << " "<<avg_vx_air[ny] << endl;
											      
cin.get();}
*/

	 avg_vx_dash_air_2[ny]= avg_vx_dash_air_2[ny]+air_velx_dash*air_velx_dash;
	 avg_vy_dash_air_2[ny]= avg_vy_dash_air_2[ny]+air_vely_dash*air_vely_dash;
	 avg_vz_dash_air_2[ny]= avg_vz_dash_air_2[ny]+air_velz_dash*air_velz_dash;

	 avg_vx_dash_air_3[ny]= avg_vx_dash_air_3[ny]+air_velx_dash*air_velx_dash*air_velx_dash;
	 avg_vy_dash_air_3[ny]= avg_vy_dash_air_3[ny]+air_vely_dash*air_vely_dash*air_vely_dash;
	 avg_vz_dash_air_3[ny]= avg_vz_dash_air_3[ny]+air_velz_dash*air_velz_dash*air_velz_dash;

	 avg_vx_dash_air_4[ny]= avg_vx_dash_air_4[ny]+air_velx_dash*air_velx_dash*air_velx_dash*air_velx_dash;
	 avg_vy_dash_air_4[ny]= avg_vy_dash_air_4[ny]+air_vely_dash*air_vely_dash*air_vely_dash*air_vely_dash;
	 avg_vz_dash_air_4[ny]= avg_vz_dash_air_4[ny]+air_velz_dash*air_velz_dash*air_velz_dash*air_velz_dash;

	 no_of_sample_for_air_vel_fluc[ny]= no_of_sample_for_air_vel_fluc[ny]+1;

	     }

      
       offset_air_velx = 1.5+ max_air_velx/del_air_velx_dist;
       offset_air_vely = 1.5+ max_air_vely/del_air_vely_dist;
       offset_air_velz = 1.5+ max_air_velz/del_air_velz_dist;
       

       for(nx=0;nx<u1.Nx();nx++)
	  for(nz=0;nz<u1.Nz();nz++)
	    for(ny=0;ny<u1.Ny();ny++)
	      {	      
	      position_index=0;    //initialization of position_index
	     if(y_grid[ny]>y_min_for_distribution_1 && y_grid[ny]<y_max_for_distribution_1)position_index=1;
	     if(y_grid[ny]>y_min_for_distribution_2 && y_grid[ny]<y_max_for_distribution_2)position_index=2;
	     if(y_grid[ny]>y_min_for_distribution_3 && y_grid[ny]<y_max_for_distribution_3)position_index=3;
	     if(position_index>0)
	       { 

      	 ubase = (1.0 - square(abs(y_grid[ny]-(b_+a_)/2.0)/((b_-a_)/2.0))); /// chennel flow
	 air_velx_dash=u1(nx,ny,nz,0)+ubase-avg_vx_air[ny];
	 air_vely_dash=u1(nx,ny,nz,1);
	 air_velz_dash=u1(nx,ny,nz,2);
	 
	 if((air_velx_dash>-max_air_velx && air_velx_dash<max_air_velx) && (air_vely_dash>-max_air_vely && air_vely_dash<max_air_vely)
	    && (air_velz_dash>-max_air_velx && air_velz_dash<max_air_velz))
	   {

	 N_shell_air_velx=(int)(air_velx_dash/del_air_velx_dist+ offset_air_velx);
	 Nof_air_velx[position_index][N_shell_air_velx]= Nof_air_velx[position_index][N_shell_air_velx]+1; 

	 N_shell_air_vely=(int)(air_vely_dash/del_air_vely_dist+ offset_air_vely);
	 Nof_air_vely[position_index][N_shell_air_vely]= Nof_air_vely[position_index][N_shell_air_vely]+1; 

	 N_shell_air_velz=(int)(air_velz_dash/del_air_velz_dist+ offset_air_velz);
	 Nof_air_velz[position_index][N_shell_air_velz]= Nof_air_velz[position_index][N_shell_air_velz]+1; 

	 Nof_sample_airvel_dist[position_index]=Nof_sample_airvel_dist[position_index]+1;


	 //cout<<"  Nof_sample_airvel_dist "<<  Nof_sample_airvel_dist<<endl;
	

	   }
	       }

	      }
  }

  void air_vel_distri_func() //rutine to calculate  dstribution (PDF) of air velocity fluctuation 
    {
  int k,ny;
  double fnorm;
  double air_vx, fof_air_vx, air_vy, fof_air_vy,air_vz,fof_air_vz;
  double y_position;
              FILE *fp301;
	      FILE *fp3011;
	      FILE *fp3012;

	      FILE *fp302;
	      FILE *fp3021;
	      FILE *fp3022;

	      FILE *fp303;
	      FILE *fp3031;
	      FILE *fp3032;

	       FILE *fp3033;
	       FILE *fp3034;
	       FILE *fp3035;
	      
 fp301=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_center_x.txt","w");
 fp3011=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_middle_x.txt","w");
 fp3012=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_wall_x.txt","w");

 fp302=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_center_y.txt","w");
 fp3021=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_middle_y.txt","w");
 fp3022=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_wall_y.txt","w");

 fp303=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_center_z.txt","w");
 fp3031=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_middle_z.txt","w");
 fp3032=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_wall_z.txt","w");	   

 fnorm =(double)(Nof_sample_airvel_dist[1]);       


              for(k=1;k<=Ndel_air_velx;k++)           // k should start from 1 not from 0  
                 {
		   air_vx= -max_air_velx + del_air_velx_dist*(k-1);
		  
		   fof_air_vx= (double)(Nof_air_velx[1][k])/((double)(Nof_sample_airvel_dist[1])*del_air_velx_dist);
		   fprintf(fp301,"%20.14lf %20.14lf  \n",air_vx,fof_air_vx );

		   fof_air_vx= (double)(Nof_air_velx[2][k])/((double)(Nof_sample_airvel_dist[2])*del_air_velx_dist);
		   fprintf(fp3011,"%20.14lf %20.14lf  \n",air_vx,fof_air_vx );

		  fof_air_vx= (double)(Nof_air_velx[3][k])/((double)(Nof_sample_airvel_dist[3])*del_air_velx_dist);
		   fprintf(fp3012,"%20.14lf %20.14lf  \n",air_vx,fof_air_vx ); 

		 }
        
	       for(k=1;k<=Ndel_air_vely;k++)                   // k should start from 1 not from 0  
		 {
		   air_vy= -max_air_vely + del_air_vely_dist*(k-1);

		   fof_air_vy= (double)(Nof_air_vely[1][k])/((double)(Nof_sample_airvel_dist[1])*del_air_vely_dist);
		   fprintf(fp302,"%20.14lf %20.14lf  \n",air_vy,fof_air_vy );

		   fof_air_vy= (double)(Nof_air_vely[2][k])/((double)(Nof_sample_airvel_dist[2])*del_air_vely_dist);
		   fprintf(fp3021,"%20.14lf %20.14lf  \n",air_vy,fof_air_vy );

		   fof_air_vy= (double)(Nof_air_vely[3][k])/((double)(Nof_sample_airvel_dist[3])*del_air_vely_dist);
		   fprintf(fp3022,"%20.14lf %20.14lf  \n",air_vy,fof_air_vy );


		 }

	       for(k=1;k<=Ndel_air_velz;k++)
		 {
		   air_vz= -max_air_velz + del_air_velz_dist*(k-1);

		   fof_air_vz= (double)(Nof_air_velz[1][k])/((double)(Nof_sample_airvel_dist[1])*del_air_velz_dist);
		   fprintf(fp303,"%20.14lf %20.14lf  \n",air_vz,fof_air_vz );

		   fof_air_vz= (double)(Nof_air_velz[2][k])/((double)(Nof_sample_airvel_dist[2])*del_air_velz_dist);
		   fprintf(fp3031,"%20.14lf %20.14lf  \n",air_vz,fof_air_vz );

		   fof_air_vz= (double)(Nof_air_velz[3][k])/((double)(Nof_sample_airvel_dist[3])*del_air_velz_dist);
		   fprintf(fp3032,"%20.14lf %20.14lf  \n",air_vz,fof_air_vz );

   
		   
		 }
	       for(ny=0;ny<Max_dns_grid;ny++)
		 {
	 avg_vx_dash_air_2[ny]= avg_vx_dash_air_2[ny]/no_of_sample_for_air_vel_fluc[ny];
	 avg_vy_dash_air_2[ny]= avg_vy_dash_air_2[ny]/no_of_sample_for_air_vel_fluc[ny];
	 avg_vz_dash_air_2[ny]= avg_vz_dash_air_2[ny]/no_of_sample_for_air_vel_fluc[ny];

	 avg_vx_dash_air_3[ny]= avg_vx_dash_air_3[ny]/no_of_sample_for_air_vel_fluc[ny];
	 avg_vy_dash_air_3[ny]= avg_vy_dash_air_3[ny]/no_of_sample_for_air_vel_fluc[ny];
	 avg_vz_dash_air_3[ny]= avg_vz_dash_air_3[ny]/no_of_sample_for_air_vel_fluc[ny];

	 avg_vx_dash_air_4[ny]= avg_vx_dash_air_4[ny]/no_of_sample_for_air_vel_fluc[ny];
	 avg_vy_dash_air_4[ny]= avg_vy_dash_air_4[ny]/no_of_sample_for_air_vel_fluc[ny];
	 avg_vz_dash_air_4[ny]= avg_vz_dash_air_4[ny]/no_of_sample_for_air_vel_fluc[ny];

		 }


 fp3033=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/vel_moment_air_2.txt","w");
 fprintf(fp3033,"%s Grid_index     Grid_point       moment_2_x     moment_2_y       moment_2_z       \n ","%"); 
 for(ny=0;ny<Max_dns_grid;ny++)
   {
     y_position=y_air_grid[Max_dns_grid-1-ny];  // be careful about the indexing of y_air_grid , check intake_average() for clarification
  fprintf(fp3033,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,    avg_vx_dash_air_2[ny], avg_vy_dash_air_2[ny] , avg_vz_dash_air_2[ny] );
   }
  fclose(fp3033);

 fp3034=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/vel_moment_air_3.txt","w");
 fprintf(fp3034,"%s Grid_index     Grid_point       moment_3_x     moment_3_y       moment_3_z       \n ","%"); 
 for(ny=0;ny<Max_dns_grid;ny++)
   {
     y_position=y_air_grid[Max_dns_grid-1-ny];  
  fprintf(fp3034,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,    avg_vx_dash_air_3[ny], avg_vy_dash_air_3[ny] , avg_vz_dash_air_3[ny] );
   }
  fclose(fp3034);

 fp3035=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/vel_moment_air_4.txt","w");
 fprintf(fp3035,"%s Grid_index     Grid_point       moment_4_x     moment_4_y       moment_4_z       \n ","%"); 
 for(ny=0;ny<Max_dns_grid;ny++)
   {
     y_position=y_air_grid[Max_dns_grid-1-ny];   
  fprintf(fp3035,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,    avg_vx_dash_air_4[ny], avg_vy_dash_air_4[ny] , avg_vz_dash_air_4[ny] );
   }
  fclose(fp3035);


	       return;
    }


void angular_air_vel_fluc_dist(const FlowField& omega1,const Vector& x_grid,const Vector& z_grid, const Vector& y_grid)
  {
    int nx,ny,nz; 
    int N_shell_air_angvelx, N_shell_air_angvely,N_shell_air_angvelz;

   
      // double a_=omega1.a();
      // double b_=omega1.b();
       
       double air_angvelx_dash,air_angvely_dash,air_angvelz_dash;
             
       double  offset_air_angvelx, offset_air_angvely , offset_air_angvelz;


     for(nx=0;nx<omega1.Nx();nx++)
	  for(nz=0;nz<omega1.Nz();nz++)
	    for(ny=0;ny<omega1.Ny();ny++)
	     {	      	  
	 //ubase = (1.0 - square(abs(y_grid[ny]-(b_+a_)/2.0)/((b_-a_)/2.0))); /// chennel flow
	 air_angvelx_dash=omega1(nx,ny,nz,0);//-avg_omegax_air[ny];
	 air_angvely_dash=omega1(nx,ny,nz,1);
	 air_angvelz_dash=omega1(nx,ny,nz,2);

	 /*if(ny==30)
{cout<<" in output function part vel_fluc_air "<<  air_velx_dash << " "<<avg_vx_air[ny] << endl;
											      
cin.get();}
*/

	 avg_angvelx_dash_air_2[ny]= avg_angvelx_dash_air_2[ny]+air_angvelx_dash*air_angvelx_dash;
	 avg_angvely_dash_air_2[ny]= avg_angvely_dash_air_2[ny]+air_angvely_dash*air_angvely_dash;
	 avg_angvelz_dash_air_2[ny]= avg_angvelz_dash_air_2[ny]+air_angvelz_dash*air_angvelz_dash;

	 avg_angvelx_dash_air_3[ny]= avg_angvelx_dash_air_3[ny]+air_angvelx_dash*air_angvelx_dash*air_angvelx_dash;
	 avg_angvely_dash_air_3[ny]= avg_angvely_dash_air_3[ny]+air_angvely_dash*air_angvely_dash*air_angvely_dash;
	 avg_angvelz_dash_air_3[ny]= avg_angvelz_dash_air_3[ny]+air_angvelz_dash*air_angvelz_dash*air_angvelz_dash;

	 avg_angvelx_dash_air_4[ny]= avg_angvelx_dash_air_4[ny]+air_angvelx_dash*air_angvelx_dash*air_angvelx_dash*air_angvelx_dash;
	 avg_angvely_dash_air_4[ny]= avg_angvely_dash_air_4[ny]+air_angvely_dash*air_angvely_dash*air_angvely_dash*air_angvely_dash;
	 avg_angvelz_dash_air_4[ny]= avg_angvelz_dash_air_4[ny]+air_angvelz_dash*air_angvelz_dash*air_angvelz_dash*air_angvelz_dash;

	 no_of_sample_for_air_angvel_fluc[ny]= no_of_sample_for_air_angvel_fluc[ny]+1;

	     }

      
       offset_air_angvelx = 1.5+ max_air_angvelx/del_air_angvelx_dist;
       offset_air_angvely = 1.5+ max_air_angvely/del_air_angvely_dist;
       offset_air_angvelz = 1.5+ max_air_angvelz/del_air_angvelz_dist;
       

       for(nx=0;nx<omega1.Nx();nx++)
	  for(nz=0;nz<omega1.Nz();nz++)
	    for(ny=0;ny<omega1.Ny();ny++)
	      {	      
	      position_index=0;    //initialization of position_index
	     if(y_grid[ny]>y_min_for_distribution_1 && y_grid[ny]<y_max_for_distribution_1)position_index=1;
	     if(y_grid[ny]>y_min_for_distribution_2 && y_grid[ny]<y_max_for_distribution_2)position_index=2;
	     if(y_grid[ny]>y_min_for_distribution_3 && y_grid[ny]<y_max_for_distribution_3)position_index=3;
	     if(position_index>0)
	       { 

      	 //ubase = (1.0 - square(abs(y_grid[ny]-(b_+a_)/2.0)/((b_-a_)/2.0))); /// chennel flow
	 air_angvelx_dash=omega1(nx,ny,nz,0)-avg_angvelx_air[ny];
	 air_angvely_dash=omega1(nx,ny,nz,1);
	 air_angvelz_dash=omega1(nx,ny,nz,2);
	 
	 if((air_angvelx_dash>-max_air_angvelx && air_angvelx_dash<max_air_angvelx) && (air_angvely_dash>-max_air_angvely && air_angvely_dash<max_air_angvely)
	    && (air_angvelz_dash>-max_air_angvelz && air_angvelz_dash<max_air_angvelz))
	   {

	 N_shell_air_angvelx=(int)(air_angvelx_dash/del_air_angvelx_dist+ offset_air_angvelx);
	 Nof_air_angvelx[position_index][N_shell_air_angvelx]= Nof_air_angvelx[position_index][N_shell_air_angvelx]+1; 

	 N_shell_air_angvely=(int)(air_angvely_dash/del_air_angvely_dist+ offset_air_angvely);
	 Nof_air_angvely[position_index][N_shell_air_angvely]= Nof_air_angvely[position_index][N_shell_air_angvely]+1; 

	 N_shell_air_angvelz=(int)(air_angvelz_dash/del_air_angvelz_dist+ offset_air_angvelz);
	 Nof_air_angvelz[position_index][N_shell_air_angvelz]= Nof_air_angvelz[position_index][N_shell_air_angvelz]+1; 

	 Nof_sample_airangvel_dist[position_index]=Nof_sample_airangvel_dist[position_index]+1;


	 //cout<<"  Nof_sample_airangvel_dist "<<  Nof_sample_airangvel_dist<<endl;
	

	   }
	       }

	      }
  }

  void air_angular_vel_distri_func() //rutine to calculate  dstribution (PDF) of air velocity fluctuation 
    {
  int k,ny;
  double fnorm;
  double air_angvelx, fof_air_angvelx, air_angvely, fof_air_angvely,air_angvelz,fof_air_angvelz;
  double y_position;
              FILE *fp10301;
	      FILE *fp103011;
	      FILE *fp103012;

	      FILE *fp10302;
	      FILE *fp103021;
	      FILE *fp103022;

	      FILE *fp10303;
	      FILE *fp103031;
	      FILE *fp103032;

	       FILE *fp103033;
	       FILE *fp103034;
	       FILE *fp103035;
	      
 fp10301=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_center_x.txt","w");
 fp103011=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_middle_x.txt","w");
 fp103012=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_wall_x.txt","w");

 fp10302=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_center_y.txt","w");
 fp103021=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_middle_y.txt","w");
 fp103022=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_wall_y.txt","w");

 fp10303=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_center_z.txt","w");
 fp103031=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_middle_z.txt","w");
 fp103032=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_dist_wall_z.txt","w");	   

 fnorm =(double)(Nof_sample_airangvel_dist[1]);       


              for(k=1;k<=Ndel_air_angvelx;k++)           // k should start from 1 not from 0  
                 {
		   air_angvelx= -max_air_angvelx + del_air_angvelx_dist*(k-1);
		  
		   fof_air_angvelx= (double)(Nof_air_angvelx[1][k])/((double)(Nof_sample_airangvel_dist[1])*del_air_angvelx_dist);
		   fprintf(fp10301,"%20.14lf %20.14lf  \n",air_angvelx,fof_air_angvelx );

		   fof_air_angvelx= (double)(Nof_air_angvelx[2][k])/((double)(Nof_sample_airangvel_dist[2])*del_air_angvelx_dist);
		   fprintf(fp103011,"%20.14lf %20.14lf  \n",air_angvelx,fof_air_angvelx );

		  fof_air_angvelx= (double)(Nof_air_angvelx[3][k])/((double)(Nof_sample_airangvel_dist[3])*del_air_angvelx_dist);
		   fprintf(fp103012,"%20.14lf %20.14lf  \n",air_angvelx,fof_air_angvelx ); 

		 }
        
	       for(k=1;k<=Ndel_air_angvely;k++)                   // k should start from 1 not from 0  
		 {
		   air_angvely= -max_air_angvely + del_air_angvely_dist*(k-1);

		   fof_air_angvely= (double)(Nof_air_angvely[1][k])/((double)(Nof_sample_airangvel_dist[1])*del_air_angvely_dist);
		   fprintf(fp10302,"%20.14lf %20.14lf  \n",air_angvely,fof_air_angvely );

		   fof_air_angvely= (double)(Nof_air_angvely[2][k])/((double)(Nof_sample_airangvel_dist[2])*del_air_angvely_dist);
		   fprintf(fp103021,"%20.14lf %20.14lf  \n",air_angvely,fof_air_angvely );

		   fof_air_angvely= (double)(Nof_air_angvely[3][k])/((double)(Nof_sample_airangvel_dist[3])*del_air_angvely_dist);
		   fprintf(fp103022,"%20.14lf %20.14lf  \n",air_angvely,fof_air_angvely );


		 }

	       for(k=1;k<=Ndel_air_angvelz;k++)
		 {
		   air_angvelz= -max_air_angvelz + del_air_angvelz_dist*(k-1);

		   fof_air_angvelz= (double)(Nof_air_angvelz[1][k])/((double)(Nof_sample_airangvel_dist[1])*del_air_angvelz_dist);
		   fprintf(fp10303,"%20.14lf %20.14lf  \n",air_angvelz,fof_air_angvelz );

		   fof_air_angvelz= (double)(Nof_air_angvelz[2][k])/((double)(Nof_sample_airangvel_dist[2])*del_air_angvelz_dist);
		   fprintf(fp103031,"%20.14lf %20.14lf  \n",air_angvelz,fof_air_angvelz );

		   fof_air_angvelz= (double)(Nof_air_angvelz[3][k])/((double)(Nof_sample_airangvel_dist[3])*del_air_angvelz_dist);
		   fprintf(fp103032,"%20.14lf %20.14lf  \n",air_angvelz,fof_air_angvelz );

   
		   
		 }
	       for(ny=0;ny<Max_dns_grid;ny++)
		 {
	 avg_angvelx_dash_air_2[ny]= avg_angvelx_dash_air_2[ny]/no_of_sample_for_air_vel_fluc[ny];
	 avg_angvely_dash_air_2[ny]= avg_angvely_dash_air_2[ny]/no_of_sample_for_air_angvel_fluc[ny];
	 avg_angvelz_dash_air_2[ny]= avg_angvelz_dash_air_2[ny]/no_of_sample_for_air_angvel_fluc[ny];

	 avg_angvelx_dash_air_3[ny]= avg_angvelx_dash_air_3[ny]/no_of_sample_for_air_angvel_fluc[ny];
	 avg_angvely_dash_air_3[ny]= avg_angvely_dash_air_3[ny]/no_of_sample_for_air_angvel_fluc[ny];
	 avg_angvelz_dash_air_3[ny]= avg_angvelz_dash_air_3[ny]/no_of_sample_for_air_angvel_fluc[ny];

	 avg_angvelx_dash_air_4[ny]= avg_angvelx_dash_air_4[ny]/no_of_sample_for_air_angvel_fluc[ny];
	 avg_angvely_dash_air_4[ny]= avg_angvely_dash_air_4[ny]/no_of_sample_for_air_angvel_fluc[ny];
	 avg_angvelz_dash_air_4[ny]= avg_angvelz_dash_air_4[ny]/no_of_sample_for_air_angvel_fluc[ny];

		 }


 fp103033=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angvel_moment_air_2.txt","w");
 fprintf(fp103033,"%s Grid_index     Grid_point       moment_2_x     moment_2_y       moment_2_z       \n ","%"); 
 for(ny=0;ny<Max_dns_grid;ny++)
   {
     y_position=y_air_grid[Max_dns_grid-1-ny];  // be careful about the indexing of y_air_grid , check intake_average() for clarification
  fprintf(fp103033,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,    avg_angvelx_dash_air_2[ny], avg_angvely_dash_air_2[ny] , avg_angvelz_dash_air_2[ny] );
   }
  fclose(fp103033);

 fp103034=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angvel_moment_air_3.txt","w");
 fprintf(fp103034,"%s Grid_index     Grid_point       moment_3_x     moment_3_y       moment_3_z       \n ","%"); 
 for(ny=0;ny<Max_dns_grid;ny++)
   {
     y_position=y_air_grid[Max_dns_grid-1-ny];  
  fprintf(fp103034,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,    avg_angvelx_dash_air_3[ny], avg_angvely_dash_air_3[ny] , avg_angvelz_dash_air_3[ny] );
   }
  fclose(fp103034);

 fp103035=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angvel_moment_air_4.txt","w");
 fprintf(fp103035,"%s Grid_index     Grid_point       moment_4_x     moment_4_y       moment_4_z       \n ","%"); 
 for(ny=0;ny<Max_dns_grid;ny++)
   {
     y_position=y_air_grid[Max_dns_grid-1-ny];   
  fprintf(fp103035,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,    avg_angvelx_dash_air_4[ny], avg_angvely_dash_air_4[ny] , avg_angvelz_dash_air_4[ny] );
   }
  fclose(fp103035);


	       return;
    }



    void acc_dist_from_air(const FlowField& u,const Vector& x_grid,const Vector& z_grid)
         {
        int ny,kk,jj; 
	int i,N_shell_accx, N_shell_accy,N_shell_accz,N_shell_accxy;

       int atom_acc_inst_dist[Max_dim_position];

       double interpolated_u[3];
       int n_order=5;
       double a_=u.a();
       double b_=u.b();
       double ubase;
       double vx_dash_air,vy_dash_air, vz_dash_air;
       double  avg_vx_air_at_py;
       double xp,yp,zp;
       double  offset_accx_from_air, offset_accy_from_air , offset_accz_from_air ,offset_accxy_from_air;


       n_sample_acc_dist_from_air=n_sample_acc_dist_from_air+1;
       
       accx_dash=new double[n_atom];
       accy_dash=new double[n_atom];
       accz_dash=new double[n_atom];

       
        for(i=0;i<n_atom;i++)
           {
	     //position_index=0;    //initialization of position_index

	     //if(p_y[i]>y_min_for_distribution_1 && p_y[i]<y_max_for_distribution_1)position_index=1;
	     //if(p_y[i]>y_min_for_distribution_2 && p_y[i]<y_max_for_distribution_2)position_index=2;
	     //if(p_y[i]>y_min_for_distribution_3 && p_y[i]<y_max_for_distribution_3)position_index=3;
	     

	     //if(position_index>0)
	     //{  
		  xp=p_x[i];
		  yp=p_y[i];
		  zp=p_z[i];
		 
    	 interpolation(u, x_grid,z_grid,xp,yp,zp, n_order, interpolated_u);
	 
	 ny=(int)(((p_y[i]-(sigma/2.0))/del_y));   
	 int ny_part=locate(y_air_grid,Max_dns_grid,p_y[i]);
	 avg_vx_air_at_py= splint(y_air_grid,avg_vx_air,y2_air_velo,ny_part,p_y[i]); // interpolation of average velocity
	 
	


	 ubase = (1.0 - square(abs(p_y[i]-(b_+a_)/2.0)/((b_-a_)/2.0))); /// chennel flow
	 
	 vx_dash_air=ubase+interpolated_u[0]-avg_vx_air_at_py;

	 vy_dash_air=interpolated_u[1];  // nondimensionalized with centerline air velocity of prabolic profile
	 vz_dash_air=interpolated_u[2];


	 accx_dash[i]= vx_dash_air/tau_vp;
	 accy_dash[i]= vy_dash_air/tau_vp;
	 accz_dash[i]= vz_dash_air/tau_vp;


	 //}
	   }

 for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
	 Nofp_from_air_[ny]=0;
    accx_dash_from_air_2[ny]=0.0;
    accy_dash_from_air_2[ny]=0.0;
    accz_dash_from_air_2[ny]=0.0;

    accx_dash_from_air_3[ny]=0.0;
    accy_dash_from_air_3[ny]=0.0;
    accz_dash_from_air_3[ny]=0.0;

    accx_dash_from_air_4[ny]=0.0;
    accy_dash_from_air_4[ny]=0.0;
    accz_dash_from_air_4[ny]=0.0;
     
        }


 for(i=0;i<n_atom;i++)
	{
	  ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
	  
	  Nofp_from_air_[ny]=Nofp_from_air_[ny]+1;

	   accx_dash_from_air_2[ny]=accx_dash_from_air_2[ny]+accx_dash[i]*accx_dash[i];
	   accy_dash_from_air_2[ny]=accy_dash_from_air_2[ny]+accy_dash[i]*accy_dash[i];
	   accz_dash_from_air_2[ny]=accz_dash_from_air_2[ny]+accz_dash[i]*accz_dash[i];

	   accx_dash_from_air_3[ny]=accx_dash_from_air_3[ny]+accx_dash[i]*accx_dash[i]*accx_dash[i];
	   accy_dash_from_air_3[ny]=accy_dash_from_air_3[ny]+accy_dash[i]*accy_dash[i]*accy_dash[i];
	   accz_dash_from_air_3[ny]=accz_dash_from_air_3[ny]+accz_dash[i]*accz_dash[i]*accz_dash[i];

	   accx_dash_from_air_4[ny]=accx_dash_from_air_4[ny]+accx_dash[i]*accx_dash[i]*accx_dash[i]*accx_dash[i];
	   accy_dash_from_air_4[ny]=accy_dash_from_air_4[ny]+accy_dash[i]*accy_dash[i]*accy_dash[i]*accy_dash[i];
	   accz_dash_from_air_4[ny]=accz_dash_from_air_4[ny]+accz_dash[i]*accz_dash[i]*accz_dash[i]*accz_dash[i];

	     }
     
 for(ny=0;ny<Ndel_y;ny++) 
	{
        if(Nofp_from_air_[ny]>0.0)
	  {
	    accx_dash_from_air_2[ny]=accx_dash_from_air_2[ny]/Nofp_from_air_[ny];
	    accy_dash_from_air_2[ny]=accy_dash_from_air_2[ny]/Nofp_from_air_[ny];
	    accz_dash_from_air_2[ny]=accz_dash_from_air_2[ny]/Nofp_from_air_[ny];

	    accx_dash_from_air_3[ny]=accx_dash_from_air_3[ny]/Nofp_from_air_[ny];
	    accy_dash_from_air_3[ny]=accy_dash_from_air_3[ny]/Nofp_from_air_[ny];
	    accz_dash_from_air_3[ny]=accz_dash_from_air_3[ny]/Nofp_from_air_[ny];

	    accx_dash_from_air_4[ny]=accx_dash_from_air_4[ny]/Nofp_from_air_[ny];
	    accy_dash_from_air_4[ny]=accy_dash_from_air_4[ny]/Nofp_from_air_[ny];
	    accz_dash_from_air_4[ny]=accz_dash_from_air_4[ny]/Nofp_from_air_[ny];

	  }

	avg_accx_dash_from_air_2[ny]=avg_accx_dash_from_air_2[ny]+accx_dash_from_air_2[ny];
	avg_accy_dash_from_air_2[ny]=avg_accy_dash_from_air_2[ny]+accy_dash_from_air_2[ny];
	avg_accz_dash_from_air_2[ny]=avg_accz_dash_from_air_2[ny]+accz_dash_from_air_2[ny];

	avg_accx_dash_from_air_3[ny]=avg_accx_dash_from_air_3[ny]+accx_dash_from_air_3[ny];
	avg_accy_dash_from_air_3[ny]=avg_accy_dash_from_air_3[ny]+accy_dash_from_air_3[ny];
	avg_accz_dash_from_air_3[ny]=avg_accz_dash_from_air_3[ny]+accz_dash_from_air_3[ny];

	avg_accx_dash_from_air_4[ny]=avg_accx_dash_from_air_4[ny]+accx_dash_from_air_4[ny];
	avg_accy_dash_from_air_4[ny]=avg_accy_dash_from_air_4[ny]+accy_dash_from_air_4[ny];
	avg_accz_dash_from_air_4[ny]=avg_accz_dash_from_air_4[ny]+accz_dash_from_air_4[ny];

	}

 /*-------------calculation starts for distribution function----------------*/
       


          
          offset_accx_from_air = 1.5+ accx_max_from_air/del_accx_from_air;
	  offset_accy_from_air = 1.5+ accy_max_from_air/del_accy_from_air;
	  offset_accz_from_air = 1.5+ accz_max_from_air/del_accz_from_air;
	  offset_accxy_from_air= 1.5+ accxy_max_from_air/del_accxy_from_air;

	 
     for(jj=0;jj<Max_dim_position;jj++)for(i=0;i<=Ndel_accx_from_air;i++) Nof_accx_inst[jj][i]=0;
     for(jj=0;jj<Max_dim_position;jj++)for(i=0;i<=Ndel_accy_from_air;i++) Nof_accy_inst[jj][i]=0;       
     for(jj=0;jj<Max_dim_position;jj++)for(i=0;i<=Ndel_accz_from_air;i++) Nof_accz_inst[jj][i]=0;
     for(jj=0;jj<Max_dim_position;jj++)for(i=0;i<=Ndel_accxy_from_air;i++)Nof_accxy_inst[jj][i]=0;         

    for(jj=0;jj<Max_dim_position;jj++) atom_acc_inst_dist[jj]=0;	 

     for(i=0;i<n_atom;i++)
           {
	     position_index=0;    //initialization of position_index

	     if(p_y[i]>y_min_for_distribution_1 && p_y[i]<y_max_for_distribution_1)position_index=1;
	     if(p_y[i]>y_min_for_distribution_2 && p_y[i]<y_max_for_distribution_2)position_index=2;
	     if(p_y[i]>y_min_for_distribution_3 && p_y[i]<y_max_for_distribution_3)position_index=3;
	     

	     if(position_index>0)
	 
	       {  
       if((accx_dash[i]>=-accx_max_from_air && accx_dash[i]<=accx_max_from_air) && (accy_dash[i]>=-accy_max_from_air && accy_dash[i]<=accy_max_from_air)
&& (accz_dash[i]>=-accz_max_from_air && accz_dash[i]<=accz_max_from_air)&& 
(accx_dash[i]*accy_dash[i]>=-accxy_max_from_air && accx_dash[i]*accy_dash[i]<=accxy_max_from_air))  // to avoid exceeding arrey subscript accidentally
   {
		 N_shell_accx=(int)(accx_dash[i]/del_accx_from_air+offset_accx_from_air);


                 Nof_accx_from_air[position_index][N_shell_accx]= Nof_accx_from_air[position_index][N_shell_accx] +1 ;
		 Nof_accx_inst[position_index][N_shell_accx]= Nof_accx_inst[position_index][N_shell_accx] +1 ;



		 N_shell_accy=(int)(accy_dash[i]/del_accy_from_air + offset_accy_from_air);
                 Nof_accy_from_air[position_index][N_shell_accy]= Nof_accy_from_air[position_index][N_shell_accy] +1 ; 
		 Nof_accy_inst[position_index][N_shell_accy]= Nof_accy_inst[position_index][N_shell_accy] +1 ; 

		 N_shell_accz=(int)(accz_dash[i]/del_accz_from_air + offset_accz_from_air);
                 Nof_accz_from_air[position_index][N_shell_accz]= Nof_accz_from_air[position_index][N_shell_accz] +1 ;
		 Nof_accz_inst[position_index][N_shell_accz]= Nof_accz_inst[position_index][N_shell_accz] +1 ;
	 
		 N_shell_accxy=(int)(accx_dash[i]*accy_dash[i]/del_accxy_from_air + offset_accxy_from_air);

		 
                 Nof_accxy_from_air[position_index][N_shell_accxy]= Nof_accxy_from_air[position_index][N_shell_accxy] +1 ; 
		 Nof_accxy_inst[position_index][N_shell_accxy]= Nof_accxy_inst[position_index][N_shell_accxy] +1 ; 


		 atom_acc_inst_dist[position_index]= atom_acc_inst_dist[position_index]+1;
		 tot_atom_acc_dist_from_air[position_index]= tot_atom_acc_dist_from_air[position_index] +1;


   }
	       }

	   }



   for(position_index=1;position_index<=max_position_index;position_index++)    
 for( kk=1;kk<=Ndel_accx_from_air;kk++) 
   { 
     fof_accx_from_air_sq[position_index][kk]=fof_accx_from_air_sq[position_index][kk]+(double)(Nof_accx_inst[position_index][kk]/
(atom_acc_inst_dist[position_index]*del_accx_from_air))*(double)(Nof_accx_inst[position_index][kk]/(atom_acc_inst_dist[position_index]*del_accx_from_air));
     // if(kk==173){cout<<"++++++++++++++++++++++++   "<<fof_accx_sq[kk]<<"    "<<" "<<atom_acc_inst_dist<<" "<<" "<<del_accx<<" "<<Nof_accx[kk]<<endl;
     //  cin.get();}

       }

 for(position_index=1;position_index<=max_position_index;position_index++) 
 for( kk=1;kk<=Ndel_accy_from_air;kk++) 
   { 
     fof_accy_from_air_sq[position_index][kk]=fof_accy_from_air_sq[position_index][kk]+(double)(Nof_accy_inst[position_index][kk]/
(atom_acc_inst_dist[position_index]*del_accy_from_air))*(double)(Nof_accy_inst[position_index][kk]/(atom_acc_inst_dist[position_index]*del_accy_from_air));
       }

 for(position_index=1;position_index<=max_position_index;position_index++) 
 for( kk=1;kk<=Ndel_accz_from_air;kk++) 
   { 
     fof_accz_from_air_sq[position_index][kk]=fof_accz_from_air_sq[position_index][kk]+(double)(Nof_accz_inst[position_index][kk]/
(atom_acc_inst_dist[position_index]*del_accz_from_air))*(double)(Nof_accz_inst[position_index][kk]/(atom_acc_inst_dist[position_index]*del_accz_from_air));
       }

 for(position_index=1;position_index<=max_position_index;position_index++) 
 for( kk=1;kk<=Ndel_accxy_from_air;kk++) 
   { 
     fof_accxy_from_air_sq[position_index][kk]=fof_accxy_from_air_sq[position_index][kk]+(double)(Nof_accxy_inst[position_index][kk]/
(atom_acc_inst_dist[position_index]*del_accxy_from_air))*(double)(Nof_accxy_inst[position_index][kk]/(atom_acc_inst_dist[position_index]*del_accxy_from_air));
       }

               delete [] accx_dash ;
	       delete [] accy_dash ;
	       delete [] accz_dash ;

	 }


 void acc_distri_func_from_air() //rutine to calculate acc dstribution (PDF) from air velocity fluctuation and free flight particle-
                                 //-velocity distribution
    {
  int k,ny;
  double fnorm,y_position;
  double accx, fof_accx, accy, fof_accy,accz,fof_accz,accxy, fof_accxy;
  double var_fof_accx, var_fof_accy,var_fof_accz,var_fof_accxy;

              FILE *fp30;
	      FILE *fp301;
	      FILE *fp302;


	      FILE *fp31;
	      FILE *fp311;
	      FILE *fp312;


	      FILE *fp32;
	      FILE *fp321;
	      FILE *fp322;

	      FILE *fp33;
	      FILE *fp331;
	      FILE *fp332;

	      FILE *fp333;
	      FILE *fp334;
	      FILE *fp335;
	      
 fp30=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_from_air_center_x.txt","w");
 fp301=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_from_air_middle_x.txt","w");
 fp302=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_from_air_wall_x.txt","w");



 fp31=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_from_air_center_y.txt","w");
 fp311=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_from_air_middle_y.txt","w");
 fp312=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_from_air_wall_y.txt","w");


 fp32=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_from_air_center_z.txt","w");
 fp321=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_from_air_middle_z.txt","w");
 fp322=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_from_air_wall_z.txt","w");

 fp33=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_from_air_center_xy.txt","w");
 fp331=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_from_air_middle_xy.txt","w");
 fp332=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_dist_from_air_wall_xy.txt","w");

	   
  fnorm =(double)(tot_atom_acc_dist_from_air[1]);       


              for(k=1;k<=Ndel_accx_from_air;k++)           // k should start from 1 not from 0  
                 {
		   accx= -accx_max_from_air + del_accx_from_air*(k-1);

		   fof_accx= (double)( Nof_accx_from_air[1][k])/((double)(tot_atom_acc_dist_from_air[1])*del_accx_from_air);
		   var_fof_accx=(fof_accx_from_air_sq[1][k]/(double)n_sample_acc_dist_from_air)- (fof_accx* fof_accx);
		   fprintf(fp30,"%20.14lf %20.14lf %20.14lf\n",accx,fof_accx,var_fof_accx);

		   fof_accx= (double)( Nof_accx_from_air[2][k])/((double)(tot_atom_acc_dist_from_air[2])*del_accx_from_air);
		   var_fof_accx=(fof_accx_from_air_sq[2][k]/(double)n_sample_acc_dist_from_air)- (fof_accx* fof_accx);
		   fprintf(fp301,"%20.14lf %20.14lf %20.14lf\n",accx,fof_accx,var_fof_accx);

		  fof_accx= (double)( Nof_accx_from_air[3][k])/((double)(tot_atom_acc_dist_from_air[3])*del_accx_from_air);
		   var_fof_accx=(fof_accx_from_air_sq[3][k]/(double)n_sample_acc_dist_from_air)- (fof_accx* fof_accx);
		   fprintf(fp302,"%20.14lf %20.14lf %20.14lf\n",accx,fof_accx,var_fof_accx); 
		 }
        
	       for(k=1;k<=Ndel_accy_from_air;k++)                   // k should start from 1 not from 0  
		 {
		   accy= -accy_max_from_air+ del_accy_from_air*(k-1);

		   fof_accy= (double)(Nof_accy_from_air[1][k])/((double)(tot_atom_acc_dist_from_air[1])*del_accy_from_air);
		   var_fof_accy= (fof_accy_from_air_sq[1][k]/(double)n_sample_acc_dist_from_air)- (fof_accy* fof_accy);
		   fprintf(fp31,"%20.14lf %20.14lf  %20.14lf\n",accy,fof_accy,var_fof_accy);

		   fof_accy= (double)(Nof_accy_from_air[2][k])/((double)(tot_atom_acc_dist_from_air[2])*del_accy_from_air);
		   var_fof_accy= (fof_accy_from_air_sq[2][k]/(double)n_sample_acc_dist_from_air)- (fof_accy* fof_accy);
		   fprintf(fp311,"%20.14lf %20.14lf  %20.14lf\n",accy,fof_accy,var_fof_accy);

		   fof_accy= (double)(Nof_accy_from_air[3][k])/((double)(tot_atom_acc_dist_from_air[3])*del_accy_from_air);
		   var_fof_accy= (fof_accy_from_air_sq[3][k]/(double)n_sample_acc_dist_from_air)- (fof_accy* fof_accy);
		   fprintf(fp312,"%20.14lf %20.14lf  %20.14lf\n",accy,fof_accy,var_fof_accy);


		 }

	       for(k=1;k<=Ndel_accz_from_air;k++)
		 {
		   accz= -accz_max_from_air+ del_accz_from_air*(k-1);
		   
		   fof_accz= (double)(Nof_accz_from_air[1][k])/((double)(tot_atom_acc_dist_from_air[1])*del_accz_from_air);
		   var_fof_accz= (fof_accz_from_air_sq[1][k]/(double)n_sample_acc_dist_from_air)- (fof_accz* fof_accz);
		   fprintf(fp32 ,"%20.14lf %20.14lf  %20.14lf\n",accz,fof_accz, var_fof_accz);

		   fof_accz= (double)(Nof_accz_from_air[2][k])/((double)(tot_atom_acc_dist_from_air[2])*del_accz_from_air);
		   var_fof_accz= (fof_accz_from_air_sq[2][k]/(double)n_sample_acc_dist_from_air)- (fof_accz* fof_accz);
		   fprintf(fp321 ,"%20.14lf %20.14lf  %20.14lf\n",accz,fof_accz, var_fof_accz);

		   fof_accz= (double)(Nof_accz_from_air[3][k])/((double)(tot_atom_acc_dist_from_air[3])*del_accz_from_air);
		   var_fof_accz= (fof_accz_from_air_sq[3][k]/(double)n_sample_acc_dist_from_air)- (fof_accz* fof_accz);
		   fprintf(fp322 ,"%20.14lf %20.14lf  %20.14lf\n",accz,fof_accz, var_fof_accz);
		   
		 }
	       
	       for(k=1;k<=Ndel_accxy_from_air;k++)
		 {
		   accxy= -accxy_max_from_air+ del_accxy_from_air*(k-1);


		   fof_accxy= (double)(Nof_accxy_from_air[1][k])/((double)(tot_atom_acc_dist_from_air[1])*del_accxy_from_air);
		   var_fof_accxy= (fof_accxy_from_air_sq[1][k]/(double)n_sample_acc_dist_from_air)- (fof_accxy* fof_accxy);
		   fprintf(fp33 ,"%20.14lf %20.14lf \n",accxy,fof_accxy);


		   fof_accxy= (double)(Nof_accxy_from_air[2][k])/((double)(tot_atom_acc_dist_from_air[2])*del_accxy_from_air);
		   var_fof_accxy= (fof_accxy_from_air_sq[2][k]/(double)n_sample_acc_dist_from_air)- (fof_accxy* fof_accxy);
		   fprintf(fp331 ,"%20.14lf %20.14lf \n",accxy,fof_accxy);

		   
		   fof_accxy= (double)(Nof_accxy_from_air[3][k])/((double)(tot_atom_acc_dist_from_air[3])*del_accxy_from_air);
		   var_fof_accxy= (fof_accxy_from_air_sq[3][k]/(double)n_sample_acc_dist_from_air)- (fof_accxy* fof_accxy);
		   fprintf(fp332 ,"%20.14lf %20.14lf  \n",accxy,fof_accxy);

	 
		 }

 for(ny=0;ny<Ndel_y;ny++) 
   {
	avg_accx_dash_from_air_2[ny]=avg_accx_dash_from_air_2[ny]/n_sample_acc_dist_from_air;
	avg_accy_dash_from_air_2[ny]=avg_accy_dash_from_air_2[ny]/n_sample_acc_dist_from_air;
	avg_accz_dash_from_air_2[ny]=avg_accz_dash_from_air_2[ny]/n_sample_acc_dist_from_air;

	avg_accx_dash_from_air_3[ny]=avg_accx_dash_from_air_3[ny]/n_sample_acc_dist_from_air;
	avg_accy_dash_from_air_3[ny]=avg_accy_dash_from_air_3[ny]/n_sample_acc_dist_from_air;
	avg_accz_dash_from_air_3[ny]=avg_accz_dash_from_air_3[ny]/n_sample_acc_dist_from_air;

	avg_accx_dash_from_air_4[ny]=avg_accx_dash_from_air_4[ny]/n_sample_acc_dist_from_air;
	avg_accy_dash_from_air_4[ny]=avg_accy_dash_from_air_4[ny]/n_sample_acc_dist_from_air;
	avg_accz_dash_from_air_4[ny]=avg_accz_dash_from_air_4[ny]/n_sample_acc_dist_from_air;

   }


fp333=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_moment_from_air_2.txt","w");
 fprintf(fp333,"%s Grid_index     Grid_point       moment_2_x     moment_2_y       moment_2_z       \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp333,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,     avg_accx_dash_from_air_2[ny] ,  avg_accy_dash_from_air_2[ny] ,  avg_accz_dash_from_air_2[ny]  );
   }
 fclose(fp333);


fp334=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_moment_from_air_3.txt","w");
 fprintf(fp334,"%s Grid_index     Grid_point       moment_3_x     moment_3_y       moment_3_z       \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp334,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,     avg_accx_dash_from_air_3[ny] ,  avg_accy_dash_from_air_3[ny] ,  avg_accz_dash_from_air_3[ny]  );
   }
 fclose(fp334);

fp335=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_moment_from_air_4.txt","w");
 fprintf(fp335,"%s Grid_index     Grid_point       moment_4_x     moment_4_y       moment_4_z       \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp335,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,     avg_accx_dash_from_air_4[ny] ,  avg_accy_dash_from_air_4[ny] ,  avg_accz_dash_from_air_4[ny]  );
   }
 fclose(fp335);



	       return;
    }


void angular_acc_dist_from_air(const FlowField& omega,const Vector& x_grid,const Vector& z_grid)
         {
        int ny,kk,jj; 
	int i,N_shell_angaccx, N_shell_angaccy,N_shell_angaccz,N_shell_angaccxy;

       int atom_angacc_inst_dist[Max_dim_position];

       double interpolated_omega[3];
       int n_order=5;
      // double a_=u.a();
      // double b_=u.b();
       //double ubase;
       double angvelx_dash_air,angvely_dash_air, angvelz_dash_air;
       //double  avg_angvelx_air_at_py;
       double xp,yp,zp;
       double  offset_angaccx_from_air, offset_angaccy_from_air , offset_angaccz_from_air ,offset_angaccxy_from_air;
       

       n_sample_angacc_dist_from_air=n_sample_angacc_dist_from_air+1;
       
       angaccx_dash=new double[n_atom];
       angaccy_dash=new double[n_atom];
       angaccz_dash=new double[n_atom];

       
        for(i=0;i<n_atom;i++)
           {
	     //position_index=0;    //initialization of position_index

	     //if(p_y[i]>y_min_for_distribution_1 && p_y[i]<y_max_for_distribution_1)position_index=1;
	     //if(p_y[i]>y_min_for_distribution_2 && p_y[i]<y_max_for_distribution_2)position_index=2;
	     //if(p_y[i]>y_min_for_distribution_3 && p_y[i]<y_max_for_distribution_3)position_index=3;
	     

	     //if(position_index>0)
	     //{  
		  xp=p_x[i];
		  yp=p_y[i];
		  zp=p_z[i];
		 
    	 interpolation(omega, x_grid,z_grid,xp,yp,zp, n_order, interpolated_omega);
	 
	 ny=(int)(((p_y[i]-(sigma/2.0))/del_y));   
	 //int ny_part=locate(y_air_grid,Max_dns_grid,p_y[i]);
	 //avg_angvelx_air_at_py= splint(y_air_grid,avg_angvelx_air,y2_air_angvelo,ny_part,p_y[i]); // interpolation of average angular velocity
	 
	


	// ubase = (1.0 - square(abs(p_y[i]-(b_+a_)/2.0)/((b_-a_)/2.0))); /// chennel flow
	 
	 angvelx_dash_air=0.5*interpolated_omega[0];//-avg_angvelx_air_at_py;

	 angvely_dash_air=0.5*interpolated_omega[1];  // nondimensionalized with centerline air angvelelocity of prabolic profile
	 angvelz_dash_air=0.5*interpolated_omega[2];


	 angaccx_dash[i]= angvelx_dash_air/tau_vp;
	 angaccy_dash[i]= angvely_dash_air/tau_vp;
	 angaccz_dash[i]= angvelz_dash_air/tau_vp;

	   }

 for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
	 Nofp_from_air_[ny]=0;
    angaccx_dash_from_air_2[ny]=0.0;
    angaccy_dash_from_air_2[ny]=0.0;
    angaccz_dash_from_air_2[ny]=0.0;

    angaccx_dash_from_air_3[ny]=0.0;
    angaccy_dash_from_air_3[ny]=0.0;
    angaccz_dash_from_air_3[ny]=0.0;

    angaccx_dash_from_air_4[ny]=0.0;
    angaccy_dash_from_air_4[ny]=0.0;
    angaccz_dash_from_air_4[ny]=0.0;
     
        }


 for(i=0;i<n_atom;i++)
	{
	  ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
	  
	  Nofp_from_air_[ny]=Nofp_from_air_[ny]+1;

	   angaccx_dash_from_air_2[ny]=angaccx_dash_from_air_2[ny]+angaccx_dash[i]*angaccx_dash[i];
	   angaccy_dash_from_air_2[ny]=angaccy_dash_from_air_2[ny]+angaccy_dash[i]*angaccy_dash[i];
	   angaccz_dash_from_air_2[ny]=angaccz_dash_from_air_2[ny]+angaccz_dash[i]*angaccz_dash[i];

	   angaccx_dash_from_air_3[ny]=angaccx_dash_from_air_3[ny]+angaccx_dash[i]*angaccx_dash[i]*angaccx_dash[i];
	   angaccy_dash_from_air_3[ny]=angaccy_dash_from_air_3[ny]+angaccy_dash[i]*angaccy_dash[i]*angaccy_dash[i];
	   angaccz_dash_from_air_3[ny]=angaccz_dash_from_air_3[ny]+angaccz_dash[i]*angaccz_dash[i]*angaccz_dash[i];

	   angaccx_dash_from_air_4[ny]=angaccx_dash_from_air_4[ny]+angaccx_dash[i]*angaccx_dash[i]*angaccx_dash[i]*angaccx_dash[i];
	   angaccy_dash_from_air_4[ny]=angaccy_dash_from_air_4[ny]+angaccy_dash[i]*angaccy_dash[i]*angaccy_dash[i]*angaccy_dash[i];
	   angaccz_dash_from_air_4[ny]=angaccz_dash_from_air_4[ny]+angaccz_dash[i]*angaccz_dash[i]*angaccz_dash[i]*angaccz_dash[i];

	     }
     
 for(ny=0;ny<Ndel_y;ny++) 
	{
        if(Nofp_from_air_[ny]>0.0)
	  {
	    angaccx_dash_from_air_2[ny]=angaccx_dash_from_air_2[ny]/Nofp_from_air_[ny];
	    angaccy_dash_from_air_2[ny]=angaccy_dash_from_air_2[ny]/Nofp_from_air_[ny];
	    angaccz_dash_from_air_2[ny]=angaccz_dash_from_air_2[ny]/Nofp_from_air_[ny];

	    angaccx_dash_from_air_3[ny]=angaccx_dash_from_air_3[ny]/Nofp_from_air_[ny];
	    angaccy_dash_from_air_3[ny]=angaccy_dash_from_air_3[ny]/Nofp_from_air_[ny];
	    angaccz_dash_from_air_3[ny]=angaccz_dash_from_air_3[ny]/Nofp_from_air_[ny];

	    angaccx_dash_from_air_4[ny]=angaccx_dash_from_air_4[ny]/Nofp_from_air_[ny];
	    angaccy_dash_from_air_4[ny]=angaccy_dash_from_air_4[ny]/Nofp_from_air_[ny];
	    angaccz_dash_from_air_4[ny]=angaccz_dash_from_air_4[ny]/Nofp_from_air_[ny];

	  }

	avg_angaccx_dash_from_air_2[ny]=avg_angaccx_dash_from_air_2[ny]+angaccx_dash_from_air_2[ny];
	avg_angaccy_dash_from_air_2[ny]=avg_angaccy_dash_from_air_2[ny]+angaccy_dash_from_air_2[ny];
	avg_angaccz_dash_from_air_2[ny]=avg_angaccz_dash_from_air_2[ny]+angaccz_dash_from_air_2[ny];

	avg_angaccx_dash_from_air_3[ny]=avg_angaccx_dash_from_air_3[ny]+angaccx_dash_from_air_3[ny];
	avg_angaccy_dash_from_air_3[ny]=avg_angaccy_dash_from_air_3[ny]+angaccy_dash_from_air_3[ny];
	avg_angaccz_dash_from_air_3[ny]=avg_angaccz_dash_from_air_3[ny]+angaccz_dash_from_air_3[ny];

	avg_angaccx_dash_from_air_4[ny]=avg_angaccx_dash_from_air_4[ny]+angaccx_dash_from_air_4[ny];
	avg_angaccy_dash_from_air_4[ny]=avg_angaccy_dash_from_air_4[ny]+angaccy_dash_from_air_4[ny];
	avg_angaccz_dash_from_air_4[ny]=avg_angaccz_dash_from_air_4[ny]+angaccz_dash_from_air_4[ny];

	}

 /*-------------calculation starts for distribution function----------------*/
       


          
          offset_angaccx_from_air = 1.5+ angaccx_max_from_air/del_angaccx_from_air;
	  offset_angaccy_from_air = 1.5+ angaccy_max_from_air/del_angaccy_from_air;
	  offset_angaccz_from_air = 1.5+ angaccz_max_from_air/del_angaccz_from_air;
	  offset_angaccxy_from_air= 1.5+ angaccxy_max_from_air/del_angaccxy_from_air;

	 
     for(jj=0;jj<Max_dim_position;jj++)for(i=0;i<=Ndel_angaccx_from_air;i++) Nof_angaccx_inst[jj][i]=0;
     for(jj=0;jj<Max_dim_position;jj++)for(i=0;i<=Ndel_angaccy_from_air;i++) Nof_angaccy_inst[jj][i]=0;       
     for(jj=0;jj<Max_dim_position;jj++)for(i=0;i<=Ndel_angaccz_from_air;i++) Nof_angaccz_inst[jj][i]=0;
     for(jj=0;jj<Max_dim_position;jj++)for(i=0;i<=Ndel_angaccxy_from_air;i++)Nof_angaccxy_inst[jj][i]=0;         

    for(jj=0;jj<Max_dim_position;jj++) atom_angacc_inst_dist[jj]=0;	 

     for(i=0;i<n_atom;i++)
           {
	     position_index=0;    //initialization of position_index

	     if(p_y[i]>y_min_for_distribution_1 && p_y[i]<y_max_for_distribution_1)position_index=1;
	     if(p_y[i]>y_min_for_distribution_2 && p_y[i]<y_max_for_distribution_2)position_index=2;
	     if(p_y[i]>y_min_for_distribution_3 && p_y[i]<y_max_for_distribution_3)position_index=3;
	     

	     if(position_index>0)
	 
	       {  
       if((angaccx_dash[i]>=-angaccx_max_from_air && angaccx_dash[i]<=angaccx_max_from_air) && (angaccy_dash[i]>=-angaccy_max_from_air && angaccy_dash[i]<=angaccy_max_from_air)
&& (angaccz_dash[i]>=-angaccz_max_from_air && angaccz_dash[i]<=angaccz_max_from_air)&& 
(angaccx_dash[i]*angaccy_dash[i]>=-angaccxy_max_from_air && angaccx_dash[i]*angaccy_dash[i]<=angaccxy_max_from_air))  // to avoid exceeding arrey subscript angaccidentally
   {
		 N_shell_angaccx=(int)(angaccx_dash[i]/del_angaccx_from_air+offset_angaccx_from_air);


                 Nof_angaccx_from_air[position_index][N_shell_angaccx]= Nof_angaccx_from_air[position_index][N_shell_angaccx] +1 ;
		 Nof_angaccx_inst[position_index][N_shell_angaccx]= Nof_angaccx_inst[position_index][N_shell_angaccx] +1 ;



		 N_shell_angaccy=(int)(angaccy_dash[i]/del_angaccy_from_air + offset_angaccy_from_air);
                 Nof_angaccy_from_air[position_index][N_shell_angaccy]= Nof_angaccy_from_air[position_index][N_shell_angaccy] +1 ; 
		 Nof_angaccy_inst[position_index][N_shell_angaccy]= Nof_angaccy_inst[position_index][N_shell_angaccy] +1 ; 

		 N_shell_angaccz=(int)(angaccz_dash[i]/del_angaccz_from_air + offset_angaccz_from_air);
                 Nof_angaccz_from_air[position_index][N_shell_angaccz]= Nof_angaccz_from_air[position_index][N_shell_angaccz] +1 ;
		 Nof_angaccz_inst[position_index][N_shell_angaccz]= Nof_angaccz_inst[position_index][N_shell_angaccz] +1 ;
	 
		 N_shell_angaccxy=(int)(angaccx_dash[i]*angaccy_dash[i]/del_angaccxy_from_air + offset_angaccxy_from_air);

		 
                 Nof_angaccxy_from_air[position_index][N_shell_angaccxy]= Nof_angaccxy_from_air[position_index][N_shell_angaccxy] +1 ; 
		 Nof_angaccxy_inst[position_index][N_shell_angaccxy]= Nof_angaccxy_inst[position_index][N_shell_angaccxy] +1 ; 


		 atom_angacc_inst_dist[position_index]= atom_angacc_inst_dist[position_index]+1;
		 tot_atom_angacc_dist_from_air[position_index]= tot_atom_angacc_dist_from_air[position_index] +1;


   }
	       }

	   }



   for(position_index=1;position_index<=max_position_index;position_index++)    
 for( kk=1;kk<=Ndel_angaccx_from_air;kk++) 
   { 
     fof_angaccx_from_air_sq[position_index][kk]=fof_angaccx_from_air_sq[position_index][kk]+(double)(Nof_angaccx_inst[position_index][kk]/
(atom_angacc_inst_dist[position_index]*del_angaccx_from_air))*(double)(Nof_angaccx_inst[position_index][kk]/(atom_angacc_inst_dist[position_index]*del_angaccx_from_air));
     // if(kk==173){cout<<"++++++++++++++++++++++++   "<<fof_angaccx_sq[kk]<<"    "<<" "<<atom_angacc_inst_dist<<" "<<" "<<del_angaccx<<" "<<Nof_angaccx[kk]<<endl;
     //  cin.get();}

       }

 for(position_index=1;position_index<=max_position_index;position_index++) 
 for( kk=1;kk<=Ndel_angaccy_from_air;kk++) 
   { 
     fof_angaccy_from_air_sq[position_index][kk]=fof_angaccy_from_air_sq[position_index][kk]+(double)(Nof_angaccy_inst[position_index][kk]/
(atom_angacc_inst_dist[position_index]*del_angaccy_from_air))*(double)(Nof_angaccy_inst[position_index][kk]/(atom_angacc_inst_dist[position_index]*del_angaccy_from_air));
       }

 for(position_index=1;position_index<=max_position_index;position_index++) 
 for( kk=1;kk<=Ndel_angaccz_from_air;kk++) 
   { 
     fof_angaccz_from_air_sq[position_index][kk]=fof_angaccz_from_air_sq[position_index][kk]+(double)(Nof_angaccz_inst[position_index][kk]/
(atom_angacc_inst_dist[position_index]*del_angaccz_from_air))*(double)(Nof_angaccz_inst[position_index][kk]/(atom_angacc_inst_dist[position_index]*del_angaccz_from_air));
       }

 for(position_index=1;position_index<=max_position_index;position_index++) 
 for( kk=1;kk<=Ndel_angaccxy_from_air;kk++) 
   { 
     fof_angaccxy_from_air_sq[position_index][kk]=fof_angaccxy_from_air_sq[position_index][kk]+(double)(Nof_angaccxy_inst[position_index][kk]/
(atom_angacc_inst_dist[position_index]*del_angaccxy_from_air))*(double)(Nof_angaccxy_inst[position_index][kk]/(atom_angacc_inst_dist[position_index]*del_angaccxy_from_air));
       }

               delete [] angaccx_dash ;
	       delete [] angaccy_dash ;
	       delete [] angaccz_dash ;

	 }


 void angular_acc_distri_func_from_air() //rutine to calculate acc dstribution (PDF) from air velocity fluctuation and free flight particle-
                                 //-velocity distribution
    {
  int k,ny;
  double fnorm,y_position;
  double angaccx, fof_angaccx, angaccy, fof_angaccy,angaccz,fof_angaccz,angaccxy, fof_angaccxy;
  double var_fof_angaccx, var_fof_angaccy,var_fof_angaccz,var_fof_angaccxy;

              FILE *fp1030;
	      FILE *fp10301;
	      FILE *fp10302;


	      FILE *fp1031;
	      FILE *fp10311;
	      FILE *fp10312;


	      FILE *fp1032;
	      FILE *fp10321;
	      FILE *fp10322;

	      FILE *fp1033;
	      FILE *fp10331;
	      FILE *fp10332;

	      FILE *fp10333;
	      FILE *fp10334;
	      FILE *fp10335;
	      
 fp1030=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_from_air_center_x.txt","w");
 fp10301=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_from_air_middle_x.txt","w");
 fp10302=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_from_air_wall_x.txt","w");



 fp1031=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_from_air_center_y.txt","w");
 fp10311=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_from_air_middle_y.txt","w");
 fp10312=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_from_air_wall_y.txt","w");


 fp1032=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_from_air_center_z.txt","w");
 fp10321=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_from_air_middle_z.txt","w");
 fp10322=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_from_air_wall_z.txt","w");

 fp1033=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_from_air_center_xy.txt","w");
 fp10331=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_from_air_middle_xy.txt","w");
 fp10332=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_dist_from_air_wall_xy.txt","w");

	   
  fnorm =(double)(tot_atom_angacc_dist_from_air[1]);       


              for(k=1;k<=Ndel_angaccx_from_air;k++)           // k should start from 1 not from 0  
                 {
		   angaccx= -angaccx_max_from_air + del_angaccx_from_air*(k-1);

		   fof_angaccx= (double)( Nof_angaccx_from_air[1][k])/((double)(tot_atom_angacc_dist_from_air[1])*del_angaccx_from_air);
		   var_fof_angaccx=(fof_angaccx_from_air_sq[1][k]/(double)n_sample_angacc_dist_from_air)- (fof_angaccx* fof_angaccx);
		   fprintf(fp1030,"%20.14lf %20.14lf %20.14lf\n",angaccx,fof_angaccx,var_fof_angaccx);

		   fof_angaccx= (double)( Nof_angaccx_from_air[2][k])/((double)(tot_atom_angacc_dist_from_air[2])*del_angaccx_from_air);
		   var_fof_angaccx=(fof_angaccx_from_air_sq[2][k]/(double)n_sample_angacc_dist_from_air)- (fof_angaccx* fof_angaccx);
		   fprintf(fp10301,"%20.14lf %20.14lf %20.14lf\n",angaccx,fof_angaccx,var_fof_angaccx);

		  fof_angaccx= (double)( Nof_angaccx_from_air[3][k])/((double)(tot_atom_angacc_dist_from_air[3])*del_angaccx_from_air);
		   var_fof_angaccx=(fof_angaccx_from_air_sq[3][k]/(double)n_sample_angacc_dist_from_air)- (fof_angaccx* fof_angaccx);
		   fprintf(fp10302,"%20.14lf %20.14lf %20.14lf\n",angaccx,fof_angaccx,var_fof_angaccx); 
		 }
        
	       for(k=1;k<=Ndel_angaccy_from_air;k++)                   // k should start from 1 not from 0  
		 {
		   angaccy= -angaccy_max_from_air+ del_angaccy_from_air*(k-1);

		   fof_angaccy= (double)(Nof_angaccy_from_air[1][k])/((double)(tot_atom_angacc_dist_from_air[1])*del_angaccy_from_air);
		   var_fof_angaccy= (fof_angaccy_from_air_sq[1][k]/(double)n_sample_angacc_dist_from_air)- (fof_angaccy* fof_angaccy);
		   fprintf(fp1031,"%20.14lf %20.14lf  %20.14lf\n",angaccy,fof_angaccy,var_fof_angaccy);

		   fof_angaccy= (double)(Nof_angaccy_from_air[2][k])/((double)(tot_atom_angacc_dist_from_air[2])*del_angaccy_from_air);
		   var_fof_angaccy= (fof_angaccy_from_air_sq[2][k]/(double)n_sample_angacc_dist_from_air)- (fof_angaccy* fof_angaccy);
		   fprintf(fp10311,"%20.14lf %20.14lf  %20.14lf\n",angaccy,fof_angaccy,var_fof_angaccy);

		   fof_angaccy= (double)(Nof_angaccy_from_air[3][k])/((double)(tot_atom_angacc_dist_from_air[3])*del_angaccy_from_air);
		   var_fof_angaccy= (fof_angaccy_from_air_sq[3][k]/(double)n_sample_angacc_dist_from_air)- (fof_angaccy* fof_angaccy);
		   fprintf(fp10312,"%20.14lf %20.14lf  %20.14lf\n",angaccy,fof_angaccy,var_fof_angaccy);


		 }

	       for(k=1;k<=Ndel_angaccz_from_air;k++)
		 {
		   angaccz= -angaccz_max_from_air+ del_angaccz_from_air*(k-1);
		   
		   fof_angaccz= (double)(Nof_angaccz_from_air[1][k])/((double)(tot_atom_angacc_dist_from_air[1])*del_angaccz_from_air);
		   var_fof_angaccz= (fof_angaccz_from_air_sq[1][k]/(double)n_sample_angacc_dist_from_air)- (fof_angaccz* fof_angaccz);
		   fprintf(fp1032 ,"%20.14lf %20.14lf  %20.14lf\n",angaccz,fof_angaccz, var_fof_angaccz);

		   fof_angaccz= (double)(Nof_angaccz_from_air[2][k])/((double)(tot_atom_angacc_dist_from_air[2])*del_angaccz_from_air);
		   var_fof_angaccz= (fof_angaccz_from_air_sq[2][k]/(double)n_sample_angacc_dist_from_air)- (fof_angaccz* fof_angaccz);
		   fprintf(fp10321 ,"%20.14lf %20.14lf  %20.14lf\n",angaccz,fof_angaccz, var_fof_angaccz);

		   fof_angaccz= (double)(Nof_angaccz_from_air[3][k])/((double)(tot_atom_angacc_dist_from_air[3])*del_angaccz_from_air);
		   var_fof_angaccz= (fof_angaccz_from_air_sq[3][k]/(double)n_sample_angacc_dist_from_air)- (fof_angaccz* fof_angaccz);
		   fprintf(fp10322 ,"%20.14lf %20.14lf  %20.14lf\n",angaccz,fof_angaccz, var_fof_angaccz);
		   
		 }
	       
	       for(k=1;k<=Ndel_angaccxy_from_air;k++)
		 {
		   angaccxy= -angaccxy_max_from_air+ del_angaccxy_from_air*(k-1);


		   fof_angaccxy= (double)(Nof_angaccxy_from_air[1][k])/((double)(tot_atom_angacc_dist_from_air[1])*del_angaccxy_from_air);
		   var_fof_angaccxy= (fof_angaccxy_from_air_sq[1][k]/(double)n_sample_angacc_dist_from_air)- (fof_angaccxy* fof_angaccxy);
		   fprintf(fp1033 ,"%20.14lf %20.14lf \n",angaccxy,fof_angaccxy);


		   fof_angaccxy= (double)(Nof_angaccxy_from_air[2][k])/((double)(tot_atom_angacc_dist_from_air[2])*del_angaccxy_from_air);
		   var_fof_angaccxy= (fof_angaccxy_from_air_sq[2][k]/(double)n_sample_angacc_dist_from_air)- (fof_angaccxy* fof_angaccxy);
		   fprintf(fp10331 ,"%20.14lf %20.14lf \n",angaccxy,fof_angaccxy);

		   
		   fof_angaccxy= (double)(Nof_angaccxy_from_air[3][k])/((double)(tot_atom_angacc_dist_from_air[3])*del_angaccxy_from_air);
		   var_fof_angaccxy= (fof_angaccxy_from_air_sq[3][k]/(double)n_sample_angacc_dist_from_air)- (fof_angaccxy* fof_angaccxy);
		   fprintf(fp10332 ,"%20.14lf %20.14lf  \n",angaccxy,fof_angaccxy);

	 
		 }

 for(ny=0;ny<Ndel_y;ny++) 
   {
	avg_angaccx_dash_from_air_2[ny]=avg_angaccx_dash_from_air_2[ny]/n_sample_angacc_dist_from_air;
	avg_angaccy_dash_from_air_2[ny]=avg_angaccy_dash_from_air_2[ny]/n_sample_angacc_dist_from_air;
	avg_angaccz_dash_from_air_2[ny]=avg_angaccz_dash_from_air_2[ny]/n_sample_angacc_dist_from_air;

	avg_angaccx_dash_from_air_3[ny]=avg_angaccx_dash_from_air_3[ny]/n_sample_angacc_dist_from_air;
	avg_angaccy_dash_from_air_3[ny]=avg_angaccy_dash_from_air_3[ny]/n_sample_angacc_dist_from_air;
	avg_angaccz_dash_from_air_3[ny]=avg_angaccz_dash_from_air_3[ny]/n_sample_angacc_dist_from_air;

	avg_angaccx_dash_from_air_4[ny]=avg_angaccx_dash_from_air_4[ny]/n_sample_angacc_dist_from_air;
	avg_angaccy_dash_from_air_4[ny]=avg_angaccy_dash_from_air_4[ny]/n_sample_angacc_dist_from_air;
	avg_angaccz_dash_from_air_4[ny]=avg_angaccz_dash_from_air_4[ny]/n_sample_angacc_dist_from_air;

   }


fp10333=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_moment_from_air_2.txt","w");
 fprintf(fp10333,"%s Grid_index     Grid_point       moment_2_x     moment_2_y       moment_2_z       \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp10333,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,     avg_angaccx_dash_from_air_2[ny] ,  avg_angaccy_dash_from_air_2[ny] ,  avg_angaccz_dash_from_air_2[ny]  );
   }
 fclose(fp10333);


fp10334=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_moment_from_air_3.txt","w");
 fprintf(fp10334,"%s Grid_index     Grid_point       moment_3_x     moment_3_y       moment_3_z       \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp10334,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,     avg_angaccx_dash_from_air_3[ny] ,  avg_angaccy_dash_from_air_3[ny] ,  avg_angaccz_dash_from_air_3[ny]  );
   }
 fclose(fp10334);

fp10335=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_moment_from_air_4.txt","w");
 fprintf(fp10335,"%s Grid_index     Grid_point       moment_4_x     moment_4_y       moment_4_z       \n ","%"); 
 for(ny=0;ny<Ndel_y;ny++)
   {
     y_position=(ny+0.5)*del_y+(0.5*sigma);   // ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
  fprintf(fp10335,"%8d %20.14lf %.14e %.14e %.14e  \n",
	  ny,    y_position,     avg_angaccx_dash_from_air_4[ny] ,  avg_angaccy_dash_from_air_4[ny] ,  avg_angaccz_dash_from_air_4[ny]  );
   }
 fclose(fp10335);



	       return;
    }


/*--------------  ---old routine-------------------------

//---------------routines for calculating particle acceleration correlation------

       void set_acc_corre()
        {

	  int i;
	  int k=0;

	  part_acc_corre=new int[Nmax_part_acc_corre];
	  accx_dash_initial=new double[Nmax_part_acc_corre];
	  accy_dash_initial=new double[Nmax_part_acc_corre];
	  accz_dash_initial=new double[Nmax_part_acc_corre];
	  
	 norm_x=new double[Nmax_part_acc_corre];
	 norm_y=new double[Nmax_part_acc_corre];
	 norm_z=new double[Nmax_part_acc_corre]; 
	 norm_xy=new double[Nmax_part_acc_corre];
	 norm_yx=new double[Nmax_part_acc_corre];


	  // ------select the particles to consider for acceleration correlation---

	  for(i=0; i<n_atom;i++){

	    if(p_y[i]>y_min_for_distribution && p_y[i]<y_max_for_distribution && k<Nmax_part_acc_corre)
	          {
		    part_acc_corre[k]=i;   //-- particle to be considered for acc- correlation calculation
		  k=k+1;
		  Nof_part_corre=k;
		    }
	  }

	  if(Nof_part_corre==0)
	    { // --if by chanse there is no atom in selected region
	      //----(contd) for safety factor it is given
 
	      Nof_part_corre=1;
	      part_acc_corre[0]=n_atom/2;
		}

	  acc_initial_time=run_time;

	  for(i=0;i<Nof_part_corre;i++)
	    {
	      int ii=part_acc_corre[i];
	      accx_dash_initial[i]=accx_dash[ii];
	      accy_dash_initial[i]=accy_dash[ii];
	      accz_dash_initial[i]=accz_dash[ii];
	     }

   fp18=fopen("/home/psg/dns_out/channel/with_pp_col/volfrac_0.0004/run_density_rough_2000/corre_part_history.txt","a");
   fprintf(fp18 ,"%20.14lf %20.14lf %20.14lf %20.14lf  \n",run_time,p_x[part_acc_corre[0]],p_y[part_acc_corre[0]],
	   p_z[part_acc_corre[0]]);
   fclose(fp18);
	  return;


	}


    void func_part_acc_corre()
     {
       FILE *fp15;
       FILE *fp16;
       FILE *fp17;

       int i;
       double *acc_corre_x,*acc_corre_y,*acc_corre_z,*acc_corre_xy,*acc_corre_yx;
       double tau_acc;
       double avg_acc_corre_x,avg_acc_corre_y,avg_acc_corre_z,avg_acc_corre_xy,avg_acc_corre_yx;

	 acc_corre_x=new double[Nmax_part_acc_corre];
         acc_corre_y=new double[Nmax_part_acc_corre];
         acc_corre_z=new double[Nmax_part_acc_corre];
         acc_corre_xy=new double[Nmax_part_acc_corre];
         acc_corre_yx=new double[Nmax_part_acc_corre];

	 


	 cout<<"N_count_acc_corre=   "<<N_count_acc_corre<<endl;
	 //cin.get();

       tau_acc=run_time-acc_initial_time;

       for(i=0;i<Nof_part_corre;i++)
	 {
	   int ii=part_acc_corre[i];
	   acc_corre_x[i] =accx_dash_initial[i]*accx_dash[ii];
	   acc_corre_y[i] =accy_dash_initial[i]*accy_dash[ii];
	   acc_corre_z[i] =accz_dash_initial[i]*accz_dash[ii];
	   acc_corre_xy[i]=accx_dash_initial[i]*accy_dash[ii];
	   acc_corre_yx[i]=accy_dash_initial[i]*accx_dash[ii];

	   //cout<<"----acc corre---   "<<ii<<"  "<<Nof_part_corre<<"    "<<acc_corre_x[i]<<endl;
	 }

      

  //---- doing averaging over all the particles, considered for calculation--

       avg_acc_corre_x=0.0;
       avg_acc_corre_y=0.0;
       avg_acc_corre_z=0.0;
       avg_acc_corre_xy=0.0;
       avg_acc_corre_yx=0.0;

       for(i=0;i<Nof_part_corre;i++)
	 {
       avg_acc_corre_x= avg_acc_corre_x+acc_corre_x[i];
       avg_acc_corre_y= avg_acc_corre_y+acc_corre_y[i];
       avg_acc_corre_z= avg_acc_corre_z+acc_corre_z[i];
       avg_acc_corre_xy=avg_acc_corre_xy+acc_corre_xy[i];
       avg_acc_corre_yx=avg_acc_corre_yx+acc_corre_yx[i];

	 }
       
       avg_acc_corre_x= avg_acc_corre_x/Nof_part_corre;
       avg_acc_corre_y= avg_acc_corre_y/Nof_part_corre;
       avg_acc_corre_z= avg_acc_corre_z/Nof_part_corre;
       avg_acc_corre_xy=avg_acc_corre_xy/Nof_part_corre;
       avg_acc_corre_yx=avg_acc_corre_yx/Nof_part_corre;

       //--------------- storing the Normalization factor----------

       if(N_count_acc_corre==0)
	 {
        for(i=0;i<Nof_part_corre;i++)
	  {
	    norm_x[i]=acc_corre_x[i];
	    norm_y[i]=acc_corre_y[i];
	    norm_z[i]=acc_corre_z[i];
	    norm_xy[i]=acc_corre_xy[i];
	    norm_yx[i]=acc_corre_yx[i];

	    //cout<<"---------   "<<norm_x[i]<<"   "<<acc_corre_x[i]<<endl;

	  }
	 }

       if(N_count_acc_corre==0)
	 {
	  avg_norm_x=avg_acc_corre_x;
	  avg_norm_y=avg_acc_corre_y;
	  avg_norm_z=avg_acc_corre_z;
	  avg_norm_xy=avg_acc_corre_xy;
	  avg_norm_yx=avg_acc_corre_yx;

	 }


       //------------Doing the Normalization of the correlation coefficient-----

       for(i=0;i<Nof_part_corre;i++)
	 {
	   acc_corre_x[i]=acc_corre_x[i]/norm_x[i];
	   acc_corre_y[i]=acc_corre_y[i]/norm_y[i];
	   acc_corre_z[i]=acc_corre_z[i]/norm_z[i];
	   acc_corre_xy[i]=acc_corre_xy[i]/norm_xy[i];
	   acc_corre_yx[i]=acc_corre_yx[i]/norm_yx[i];
	 }


       avg_acc_corre_x=avg_acc_corre_x/avg_norm_x;
       avg_acc_corre_y=avg_acc_corre_y/avg_norm_y;
       avg_acc_corre_z=avg_acc_corre_z/avg_norm_z;
       avg_acc_corre_xy=avg_acc_corre_xy/avg_norm_xy;
       avg_acc_corre_yx=avg_acc_corre_yx/avg_norm_yx;


       //cout<<"----acc corre---   "<<"  "<<Nof_part_corre<<"    "<<acc_corre_x[0]<<"  "<<norm_x[0]<<endl;
	      
   fp15=fopen("/home/psg/dns_out/channel/with_pp_col/volfrac_0.0004/run_density_rough_2000/acc_corre.txt","a");
   fprintf(fp15 ,"%20.14lf %20.14lf  %20.14lf  %20.14lf  %20.14lf  %20.14lf \n",tau_acc,acc_corre_x[0],acc_corre_y[0],acc_corre_z[0],
	   acc_corre_xy[0],acc_corre_yx[0]);

   fclose(fp15);


   fp16=fopen("/home/psg/dns_out/channel/with_pp_col/volfrac_0.0004/run_density_rough_2000/acc_corre_avg.txt","a");

   fprintf(fp16 ,"%20.14lf %20.14lf %20.14lf %20.14lf %20.14lf %20.14lf \n",tau_acc, avg_acc_corre_x, avg_acc_corre_y, avg_acc_corre_z, 
	   avg_acc_corre_xy, avg_acc_corre_yx);

   fclose(fp16);


   fp17=fopen("/home/psg/dns_out/channel/with_pp_col/volfrac_0.0004/run_density_rough_2000/acc_details.txt","a");
   fprintf(fp17 ,"%20.14lf %20.14lf %20.14lf %20.14lf  %20.14lf\n",run_time,accx_dash[part_acc_corre[5]],accy_dash[part_acc_corre[5]],
	   accz_dash[part_acc_corre[5]],p_y[part_acc_corre[5]]);
   fclose(fp17);

   delete [] acc_corre_x;
   delete [] acc_corre_y;
   delete [] acc_corre_z;
   delete [] acc_corre_xy;
   delete [] acc_corre_yx;

     }




*/



//---------------routines for calculating particle acceleration correlation------

     void call_part_acc_corre()
     {
       int ny;      
       int i;
       int *Nofp;
       double *accx_inst,*accy_inst,*accz_inst;
       double  avg_vx_air_at_py, avg_vx_part_at_py;

       if(index_return_from_part_acc_corre==1)
	 {
          return;
	 }

       Nofp=new int[Ndel_y];

       accx_inst=new double[Ndel_y];   
       accy_inst=new double[Ndel_y];
       accz_inst=new double[Ndel_y];

       accx_dash=new double[n_atom];
       accy_dash=new double[n_atom];
       accz_dash=new double[n_atom];
       
       if(index_average_from_input==1 && N_count_acc_corre==0)intake_average(); //--- colling is required 
                                                           //--when input average values are supplied.
       
   for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
    Nofp[ny]=0;
    accx_inst[ny]=0.0;
    accy_inst[ny]=0.0;
    accz_inst[ny]=0.0;
       }

   if(index_average_from_input==1)
     {
       for(i=0;i<n_atom;i++)
	 {
       ny=(int)(((p_y[i]-(sigma/2.0))/del_y));     
       if(p_y[i]>y_min_for_corre && p_y[i]<y_max_for_corre )//--- only the particles within limit of y_min to y_max have been considered.
	  {
	   int ny_part=locate(y_air_grid,Max_dns_grid,p_y[i]);
	   avg_vx_air_at_py= splint(y_air_grid,avg_vx_air,y2_air_velo,ny_part,p_y[i]); 

	   int ny1_part=locate(y_part_grid,Ndel_y,p_y[i]);
	   avg_vx_part_at_py= splint(y_part_grid,avg_vx_part,y2_part_velo,ny1_part,p_y[i]);

	   accx_dash[i]=acc_x[i]- (avg_vx_air_at_py-avg_vx_part_at_py)/tau_vp; //--taking care of the sign
	   accy_dash[i]=acc_y[i];
	   accz_dash[i]=acc_z[i];
	  }


	 }


     }  //-end of the if loop

   else{

 for(i=0;i<n_atom;i++)// summing up all the acceleration in each the y-grid for all the particle 
    {
      ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

      // cout<<"in property"<<p_y[i]<<" "<<ny<<endl;


              Nofp[ny]= Nofp[ny]+ 1;

              accx_inst[ny]= accx_inst[ny]+acc_x[i];
              accy_inst[ny]= accy_inst[ny]+acc_y[i];  
              accz_inst[ny]= accz_inst[ny]+acc_z[i];
      }

 for(ny=0;ny<Ndel_y;ny++) // doing the particle average in each y-grid-in each sample frame
       {
    if(Nofp[ny]>0.0)
      {
     accx_inst[ny]= accx_inst[ny]/Nofp[ny];  // accx_inst[ny] is basically particle averaged inst in time.
     accy_inst[ny]= accy_inst[ny]/Nofp[ny];
     accz_inst[ny]= accz_inst[ny]/Nofp[ny];
      }
    
       }



 for(i=0;i<n_atom;i++)
   {
 ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

 accx_dash[i]=acc_x[i]-accx_inst[ny];                // calculating fluctuation for each particle
 accy_dash[i]=acc_y[i]-accy_inst[ny];
 accz_dash[i]=acc_z[i]-accz_inst[ny];
 
 
   }
          } //---------------end of the else loop


 if(N_count_acc_corre==0)set_acc_corre();
  
	     func_part_acc_corre();              //------ for calculationg particle acceleration corre.
	     N_count_acc_corre=N_count_acc_corre+1;

  

	       delete [] Nofp;

	       delete []accx_inst ; 
	       delete []accy_inst ;
	       delete []accz_inst ;

	       delete [] accx_dash ;
	       delete [] accy_dash ;
	       delete [] accz_dash ;



     }

void call_part_angular_acc_corre()
     {
       int ny;      
       int i;
       int *Nofp;
       double *angaccx_inst,*angaccy_inst,*angaccz_inst;
       double  /*avg_angvelx_air_at_py,*/ avg_angvelx_part_at_py;

       if(index_return_from_part_angacc_corre==1)
	 {
          return;
	 }

       Nofp=new int[Ndel_y];

       angaccx_inst=new double[Ndel_y];   
       angaccy_inst=new double[Ndel_y];
       angaccz_inst=new double[Ndel_y];

       angaccx_dash=new double[n_atom];
       angaccy_dash=new double[n_atom];
       angaccz_dash=new double[n_atom];
       
       if(index_average_from_input==1 && N_count_angacc_corre==0)intake_average(); //--- colling is required 
                                                           //--when input average values are supplied.
       
   for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
    Nofp[ny]=0;
    angaccx_inst[ny]=0.0;
    angaccy_inst[ny]=0.0;
    angaccz_inst[ny]=0.0;
       }

   if(index_average_from_input==1)
     {
       for(i=0;i<n_atom;i++)
	 {
       ny=(int)(((p_y[i]-(sigma/2.0))/del_y));     
       if(p_y[i]>y_min_for_corre && p_y[i]<y_max_for_corre )//--- only the particles within limit of y_min to y_max have been considered.
	  {
	  // int ny_part=locate(y_air_grid,Max_dns_grid,p_y[i]);
	  // avg_angvelx_air_at_py= splint(y_air_grid,avg_angvelx_air,y2_air_angvelo,ny_part,p_y[i]); 

	   int ny1_part=locate(y_part_grid,Ndel_y,p_y[i]);
	   avg_angvelx_part_at_py= splint(y_part_grid,avg_angvelx_part,y2_part_velo,ny1_part,p_y[i]);

	   angaccx_dash[i]=angaccx[i] + (avg_angvelx_part_at_py)/tau_vp; //--taking care of the sign
	   angaccy_dash[i]=angaccy[i];
	   angaccz_dash[i]=angaccz[i];
	  }


	 }


     }  //-end of the if loop

   else{

 for(i=0;i<n_atom;i++)// summing up all the acceleration in each the y-grid for all the particle 
    {
      ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

      // cout<<"in property"<<p_y[i]<<" "<<ny<<endl;


              Nofp[ny]= Nofp[ny]+ 1;

              angaccx_inst[ny]= angaccx_inst[ny]+angaccx[i];
              angaccy_inst[ny]= angaccy_inst[ny]+angaccy[i];  
              angaccz_inst[ny]= angaccz_inst[ny]+angaccz[i];
      }

 for(ny=0;ny<Ndel_y;ny++) // doing the particle average in each y-grid-in each sample frame
       {
    if(Nofp[ny]>0.0)
      {
     angaccx_inst[ny]= angaccx_inst[ny]/Nofp[ny];  // angaccx_inst[ny] is basically particle averaged inst in time.
     angaccy_inst[ny]= angaccy_inst[ny]/Nofp[ny];
     angaccz_inst[ny]= angaccz_inst[ny]/Nofp[ny];
      }
    
       }



 for(i=0;i<n_atom;i++)
   {
 ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

 angaccx_dash[i]=angaccx[i]-angaccx_inst[ny];                // calculating fluctuation for each particle
 angaccy_dash[i]=angaccy[i]-angaccy_inst[ny];
 angaccz_dash[i]=angaccz[i]-angaccz_inst[ny];
 
 
   }
          } //---------------end of the else loop


 if(N_count_angacc_corre==0)set_angular_acc_corre();
  
	     func_part_angular_acc_corre();              //------ for calculationg particle angacceleration corre.
	     N_count_angacc_corre=N_count_angacc_corre+1;

  

	       delete [] Nofp;

	       delete []angaccx_inst ; 
	       delete []angaccy_inst ;
	       delete []angaccz_inst ;

	       delete [] angaccx_dash ;
	       delete [] angaccy_dash ;
	       delete [] angaccz_dash ;



     }

void intake_average()
{ int i;
 double rough; 
 string rough_line;
  y2_part_velo= new double[Ndel_y];
  y2_air_velo = new double[Max_dns_grid];
  y_part_grid = new double[Ndel_y];
  y_air_grid  = new double[Max_dns_grid];
  avg_vx_air  = new double[Max_dns_grid];


ifstream infile_avg_air_velo("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/Umean_0.asc");

 if (!infile_avg_air_velo.good()) {
   cerr << "can't open file Umean_0.asc for air average velocity" <<endl;
	
    exit(1);
 }



 getline(infile_avg_air_velo,rough_line);

 for(i=Max_dns_grid-1;i>=0;i--)
   {
     infile_avg_air_velo>>avg_vx_air[i];
    
     // cout<<y_air_grid[i]<<" "<<avg_vx_air[i]<<endl;
       }

 infile_avg_air_velo.close();

ifstream infile_air_grid("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/y_0.asc");

 if (!infile_air_grid.good()) {
   cerr << "can't open file y_0.asc "<<endl; 
       
    exit(1);
 }

 getline(infile_air_grid,rough_line);

 for(i=Max_dns_grid-1;i>=0;i--)
   {
         infile_air_grid>>y_air_grid[i];
	 //cout<<y_air_grid[i]<<" "<<avg_vx_air[i]<<endl;
       }

 infile_air_grid.close();



ifstream infile_avg_part_velo("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/velocity_0.txt");

 if (!infile_avg_part_velo.good()) {
   cerr << "can't open file velocity_0.txt for particle average velocity " <<endl;
	 
    exit(1);
 }

 getline(infile_avg_part_velo,rough_line);
 // cout<<"-------------------  "<<rough_line<<endl; 
 getline(infile_avg_part_velo,rough_line);
 //cout<<rough_line<<endl; 

 for(i=0;i<Ndel_y;i++)
   {
     infile_avg_part_velo>>rough>>y_part_grid[i]>>avg_vx_part[i]>>rough>>rough>>rough>>rough;
     
        cout<<y_part_grid[i]<<" "<<avg_vx_part[i]<<endl;
     }

 infile_avg_part_velo.close();
 

 spline(y_air_grid,avg_vx_air,Max_dns_grid,y2_air_velo);
 spline(y_part_grid,avg_vx_part,Ndel_y,y2_part_velo); 



ifstream infile_avg_part_angvelo("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angvelocity_0.txt");

 if (!infile_avg_part_angvelo.good()) {
   cerr << "can't open file angvelocity_0.txt for particle average angular velocity " <<endl;
	 
    exit(1);
 }

 getline(infile_avg_part_angvelo,rough_line);
 // cout<<"-------------------  "<<rough_line<<endl; 
 getline(infile_avg_part_angvelo,rough_line);
 //cout<<rough_line<<endl; 

 for(i=0;i<Ndel_y;i++)
   {
     infile_avg_part_angvelo>>rough>>y_part_grid[i]>>avg_angvelx_part[i]>>rough>>rough>>rough>>rough;
     
        cout<<y_part_grid[i]<<" "<<avg_angvelx_part[i]<<endl;
     }

 infile_avg_part_angvelo.close();
 

 //spline(y_air_grid,avg_angvelx_air,Max_dns_grid,y2_air_velo);
 spline(y_part_grid,avg_angvelx_part,Ndel_y,y2_part_velo);  
  return;
}



      

void set_acc_corre()
        {

	  int i;
	  //int k=0;

	  del_y_acc_corre=b_ly/(double)n_del_y_acc_corre;

	  
	  part_acc_corre=new int[Nmax_part_acc_corre];
	  accx_dash_initial=new double[Nmax_part_acc_corre];
	  accy_dash_initial=new double[Nmax_part_acc_corre];
	  accz_dash_initial=new double[Nmax_part_acc_corre];
	  


	 avg_norm_x=new double[n_del_y_acc_corre];
	 avg_norm_y=new double[n_del_y_acc_corre];
	 avg_norm_z=new double[n_del_y_acc_corre];
	 avg_norm_xy=new double[n_del_y_acc_corre];
	 avg_norm_yx=new double[n_del_y_acc_corre];


	 acc_initial_time=run_time;

	 for(i=0;i<n_atom;i++)part_acc_corre[i]=1;

	 for(i=0;i<n_atom;i++)
	    {
	     
	      index_part_y_post[i]=(int)(p_y[i]/del_y_acc_corre);
	      
	      accx_dash_initial[i]=accx_dash[i];
	      accy_dash_initial[i]=accy_dash[i];
	      accz_dash_initial[i]=accz_dash[i];
	     }


	 /* old before 23.02.2008-----to calculate acc corre only at center

	  // ------select the particles to consider for acceleration correlation---

	  for(i=0; i<n_atom;i++){

	    if(p_y[i]>y_min_for_corre && p_y[i]<y_max_for_corre && k<Nmax_part_acc_corre)
	          {
		    part_acc_corre[k]=i;   //-- particle to be considered for acc- correlation calculation
		  k=k+1;
		  Nof_part_corre=k;
		    }
	  }

	  if(Nof_part_corre==0)
	    { // --if by chanse there is no atom in selected region
	      //----(contd) for safety factor it is given
 
	      Nof_part_corre=1;
	      part_acc_corre[0]=n_atom/2;
		}

	  acc_initial_time=run_time;

	  for(i=0;i<Nof_part_corre;i++)
	    {
	      int ii=part_acc_corre[i];
	      accx_dash_initial[i]=accx_dash[ii];
	      accy_dash_initial[i]=accy_dash[ii];
	      accz_dash_initial[i]=accz_dash[ii];
	     }
	 */



   fp18=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/corre_part_history.txt","a");
   fprintf(fp18 ,"%20.14lf %20.14lf %20.14lf %20.14lf  \n",run_time,p_x[n_atom/2],p_y[n_atom/2], p_z[n_atom/2]);
   fclose(fp18);
	  return;


	}


    void func_part_acc_corre()
     {
       
       FILE *fp16;
      

       FILE *fp161;
       FILE *fp162;
       FILE *fp163;
       FILE *fp164;


	 


       int i,k;
       double *acc_corre_x,*acc_corre_y,*acc_corre_z,*acc_corre_xy,*acc_corre_yx;
       double tau_acc;
       double *avg_acc_corre_x,*avg_acc_corre_y,*avg_acc_corre_z,*avg_acc_corre_xy,*avg_acc_corre_yx;
       int *total_part_count_acc_corre;
       
	 acc_corre_x=new double[Nmax_part_acc_corre];
         acc_corre_y=new double[Nmax_part_acc_corre];
         acc_corre_z=new double[Nmax_part_acc_corre];
         acc_corre_xy=new double[Nmax_part_acc_corre];
         acc_corre_yx=new double[Nmax_part_acc_corre];

	
 
	 avg_acc_corre_x=new double[n_del_y_acc_corre];
	 avg_acc_corre_y=new double[n_del_y_acc_corre];
	 avg_acc_corre_z=new double[n_del_y_acc_corre];
	 avg_acc_corre_xy=new double[n_del_y_acc_corre];
	 avg_acc_corre_yx=new double[n_del_y_acc_corre];

	 total_part_count_acc_corre=new int[n_del_y_acc_corre];
	 
	 for(i=0;i<n_del_y_acc_corre;i++)
	   total_part_count_acc_corre[i]=0;


	 // cout<<"N_count_acc_corre=   "<<N_count_acc_corre<<endl;
	 //cin.get();

       tau_acc=run_time-acc_initial_time;
        
       
       for(i=0;i<n_atom;i++)
	 {
           int ii=(int)(p_y[i]/del_y_acc_corre); //index_part_y_post[i]
	  
	if( ii==index_part_y_post[i]&& part_acc_corre[i]>=0)
	 {
	   
	   acc_corre_x[i] =accx_dash_initial[i]*accx_dash[i];
	   acc_corre_y[i] =accy_dash_initial[i]*accy_dash[i];
	   acc_corre_z[i] =accz_dash_initial[i]*accz_dash[i];
	   acc_corre_xy[i]=accx_dash_initial[i]*accy_dash[i];
	   acc_corre_yx[i]=accy_dash_initial[i]*accx_dash[i];
	   
	   total_part_count_acc_corre[ii]=total_part_count_acc_corre[ii]+1;

	   //cout<<"----acc corre---   "<<ii<<"  "<<Nof_part_corre<<"    "<<acc_corre_x[i]<<endl;
	 }
	else
	  {
	    part_acc_corre[i]=-1;  // removing this particle from correlation list
	    acc_corre_x[i]=0.0;   // these are put to zero, to be added to the sum, with no contribution.
	    acc_corre_y[i]=0.0;
	    acc_corre_z[i]=0.0;
	    acc_corre_xy[i]=0.0;
	    acc_corre_yx[i]=0.0; 

	  }
	 }




  //---- doing averaging over all the particles, considered for calculation--

       for(k=0;k<n_del_y_acc_corre;k++)
	 {
       avg_acc_corre_x[k]=0.0;
       avg_acc_corre_y[k]=0.0;
       avg_acc_corre_z[k]=0.0;
       avg_acc_corre_xy[k]=0.0;
       avg_acc_corre_yx[k]=0.0;

	 }

       for(i=0;i<n_atom;i++)
	 {
	   k=(int)(p_y[i]/del_y_acc_corre);

       avg_acc_corre_x[k]= avg_acc_corre_x[k]+acc_corre_x[i];
       avg_acc_corre_y[k]= avg_acc_corre_y[k]+acc_corre_y[i];
       avg_acc_corre_z[k]= avg_acc_corre_z[k]+acc_corre_z[i];
       avg_acc_corre_xy[k]=avg_acc_corre_xy[k]+acc_corre_xy[i];
       avg_acc_corre_yx[k]=avg_acc_corre_yx[k]+acc_corre_yx[i];

	 }

        for(k=0;k<n_del_y_acc_corre;k++)
	  {
	    if(total_part_count_acc_corre[k]>0)
	      {
	      avg_acc_corre_x[k]= avg_acc_corre_x[k]/total_part_count_acc_corre[k];
	      avg_acc_corre_y[k]= avg_acc_corre_y[k]/total_part_count_acc_corre[k];
	      avg_acc_corre_z[k]= avg_acc_corre_z[k]/total_part_count_acc_corre[k];
	      avg_acc_corre_xy[k]=avg_acc_corre_xy[k]/total_part_count_acc_corre[k];
	      avg_acc_corre_yx[k]=avg_acc_corre_yx[k]/total_part_count_acc_corre[k];
	      }
	  }
       //--------------- storing the Normalization factor----------

     

       if(N_count_acc_corre==0)
         for(k=0;k<n_del_y_acc_corre;k++)
	   {
	    avg_norm_x[k]=avg_acc_corre_x[k];
	    avg_norm_y[k]=avg_acc_corre_y[k];
	    avg_norm_z[k]=avg_acc_corre_z[k];
	    avg_norm_xy[k]=avg_acc_corre_xy[k];
	    avg_norm_yx[k]=avg_acc_corre_yx[k];

	   }


       //------------Doing the Normalization of the correlation coefficient-----

      
       for(k=0;k<n_del_y_acc_corre;k++)
	   {
       avg_acc_corre_x[k]=avg_acc_corre_x[k]/avg_norm_x[k];
       avg_acc_corre_y[k]=avg_acc_corre_y[k]/avg_norm_y[k];
       avg_acc_corre_z[k]=avg_acc_corre_z[k]/avg_norm_z[k];
       avg_acc_corre_xy[k]=avg_acc_corre_xy[k]/avg_norm_xy[k];
       avg_acc_corre_yx[k]=avg_acc_corre_yx[k]/avg_norm_yx[k];
	   }
       

       for(k=0;k<n_del_y_acc_corre;k++)
	   if(total_part_count_acc_corre[k]<min_Nop_for_corre)
	 {
	   index_return_from_part_acc_corre=1;
	   return;
	 }
       //cout<<"----acc corre---   "<<"  "<<Nof_part_corre<<"    "<<acc_corre_x[0]<<"  "<<norm_x[0]<<endl;
	      
  

   fp16=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_corre_avg_x.txt","a");
   if(N_count_vel_corre==0)
     { 
       fprintf(fp16,"%s  avg_acc_norm_x[ at diff y]\n","%");
       fprintf(fp16,"%s ","%");
	 for(i=0;i<n_del_y_acc_corre;i++)
	   fprintf(fp16, "%20.14lf", avg_norm_x[i]);
       fprintf(fp16,"\n");

    }
   fprintf(fp16 ,"%20.14lf",tau_acc);
 for(i=0;i<n_del_y_acc_corre;i++)
   fprintf(fp16 ,"%20.14lf", avg_acc_corre_x[i]);
    fprintf(fp16 ,"\n");
   fclose(fp16);



 fp161=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_corre_avg_y.txt","a");
   if(N_count_vel_corre==0)
     { 
       fprintf(fp161,"%s  avg_acc_norm_y[ at diff y]\n","%");
       fprintf(fp161,"%s ","%");
	 for(i=0;i<n_del_y_acc_corre;i++)
	   fprintf(fp161, "%20.14lf", avg_norm_y[i]);
       fprintf(fp161,"\n\n");

    }
   fprintf(fp161 ,"%20.14lf",tau_acc);
 for(i=0;i<n_del_y_acc_corre;i++)
   fprintf(fp161 ,"%20.14lf", avg_acc_corre_y[i]);
    fprintf(fp161 ,"\n");
   fclose(fp161);



 fp162=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_corre_avg_z.txt","a");
   if(N_count_vel_corre==0)
     { 
       fprintf(fp162,"%s  avg_acc_norm_z[ at diff y]\n","%");
       fprintf(fp162,"%s ","%");
	 for(i=0;i<n_del_y_acc_corre;i++)
	   fprintf(fp162, "%20.14lf", avg_norm_z[i]);
       fprintf(fp162,"\n");
    }
   fprintf(fp162 ,"%20.14lf",tau_acc);
 for(i=0;i<n_del_y_acc_corre;i++)
   fprintf(fp162 ,"%20.14lf", avg_acc_corre_z[i]);
    fprintf(fp162 ,"\n");
   fclose(fp162);



 fp163=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_corre_avg_xy.txt","a");
   if(N_count_vel_corre==0)
     { 
       fprintf(fp163,"%s  avg_acc_norm_xy[ at diff y]\n","%");
       fprintf(fp163,"%s ","%");
	 for(i=0;i<n_del_y_acc_corre;i++)
	   fprintf(fp163, "%20.14lf", avg_norm_xy[i]);
       fprintf(fp163,"\n");

    }
   fprintf(fp163 ,"%20.14lf",tau_acc);
 for(i=0;i<n_del_y_acc_corre;i++)
   fprintf(fp163 ,"%20.14lf", avg_acc_corre_xy[i]);
    fprintf(fp163 ,"\n");
   fclose(fp163);


 fp164=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/acc_corre_avg_yx.txt","a");
   if(N_count_vel_corre==0)
     { 
       fprintf(fp164,"%s  avg_acc_norm_x[ at diff y]\n","%");
       fprintf(fp164,"%s ","%");
	 for(i=0;i<n_del_y_acc_corre;i++)
	   fprintf(fp164, "%20.14lf", avg_norm_yx[i]);
       fprintf(fp164,"\n");

    }
   fprintf(fp164 ,"%20.14lf",tau_acc);
 for(i=0;i<n_del_y_acc_corre;i++)
   fprintf(fp164 ,"%20.14lf", avg_acc_corre_yx[i]);
    fprintf(fp164 ,"\n");
   fclose(fp164);




   delete [] acc_corre_x;
   delete [] acc_corre_y;
   delete [] acc_corre_z;
   delete [] acc_corre_xy;
   delete [] acc_corre_yx;


   delete[] avg_acc_corre_x;
   delete[] avg_acc_corre_y;
   delete[] avg_acc_corre_z;
   delete[] avg_acc_corre_xy;
   delete[] avg_acc_corre_yx;
   
     }

 void set_angular_acc_corre()
        {

	  int i;
	  //int k=0;

	  del_y_angacc_corre=b_ly/(double)n_del_y_angacc_corre;

	  
	  part_angacc_corre=new int[Nmax_part_angacc_corre];
	  angaccx_dash_initial=new double[Nmax_part_angacc_corre];
	  angaccy_dash_initial=new double[Nmax_part_angacc_corre];
	  angaccz_dash_initial=new double[Nmax_part_angacc_corre];
	  


	 angavg_norm_x=new double[n_del_y_angacc_corre];
	 angavg_norm_y=new double[n_del_y_angacc_corre];
	 angavg_norm_z=new double[n_del_y_angacc_corre];
	 angavg_norm_xy=new double[n_del_y_angacc_corre];
	 angavg_norm_yx=new double[n_del_y_angacc_corre];


	 angacc_initial_time=run_time;

	 for(i=0;i<n_atom;i++)part_angacc_corre[i]=1;

	 for(i=0;i<n_atom;i++)
	    {
	     
	      index_part_y_post[i]=(int)(p_y[i]/del_y_angacc_corre);
	      
	      angaccx_dash_initial[i]=angaccx_dash[i];
	      angaccy_dash_initial[i]=angaccy_dash[i];
	      angaccz_dash_initial[i]=angaccz_dash[i];
	     }


	 /* old before 23.02.2008-----to calculate angacc corre only at center

	  // ------select the particles to consider for angacceleration correlation---

	  for(i=0; i<n_atom;i++){

	    if(p_y[i]>y_min_for_corre && p_y[i]<y_max_for_corre && k<Nmax_part_angacc_corre)
	          {
		    part_angacc_corre[k]=i;   //-- particle to be considered for angacc- correlation calculation
		  k=k+1;
		  Nof_part_corre=k;
		    }
	  }

	  if(Nof_part_corre==0)
	    { // --if by chanse there is no atom in selected region
	      //----(contd) for safety factor it is given
 
	      Nof_part_corre=1;
	      part_angacc_corre[0]=n_atom/2;
		}

	  angacc_initial_time=run_time;

	  for(i=0;i<Nof_part_corre;i++)
	    {
	      int ii=part_angacc_corre[i];
	      angaccx_dash_initial[i]=angaccx_dash[ii];
	      angaccy_dash_initial[i]=angaccy_dash[ii];
	      angaccz_dash_initial[i]=angaccz_dash[ii];
	     }
	 */



   fp1018=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/corre_part_history.txt","a");
   fprintf(fp1018 ,"%20.14lf %20.14lf %20.14lf %20.14lf  \n",run_time,p_x[n_atom/2],p_y[n_atom/2], p_z[n_atom/2]);
   fclose(fp1018);
	  return;


	}


    void func_part_angular_acc_corre()
     {
       
       FILE *fp1016;
      

       FILE *fp10161;
       FILE *fp10162;
       FILE *fp10163;
       FILE *fp10164;


	 


       int i,k;
       double *angacc_corre_x,*angacc_corre_y,*angacc_corre_z,*angacc_corre_xy,*angacc_corre_yx;
       double tau_angacc;
       double *avg_angacc_corre_x,*avg_angacc_corre_y,*avg_angacc_corre_z,*avg_angacc_corre_xy,*avg_angacc_corre_yx;
       int *total_part_count_angacc_corre;
       
	 angacc_corre_x=new double[Nmax_part_angacc_corre];
         angacc_corre_y=new double[Nmax_part_angacc_corre];
         angacc_corre_z=new double[Nmax_part_angacc_corre];
         angacc_corre_xy=new double[Nmax_part_angacc_corre];
         angacc_corre_yx=new double[Nmax_part_angacc_corre];

	
 
	 avg_angacc_corre_x=new double[n_del_y_angacc_corre];
	 avg_angacc_corre_y=new double[n_del_y_angacc_corre];
	 avg_angacc_corre_z=new double[n_del_y_angacc_corre];
	 avg_angacc_corre_xy=new double[n_del_y_angacc_corre];
	 avg_angacc_corre_yx=new double[n_del_y_angacc_corre];

	 total_part_count_angacc_corre=new int[n_del_y_angacc_corre];
	 
	 for(i=0;i<n_del_y_angacc_corre;i++)
	   total_part_count_angacc_corre[i]=0;


	 // cout<<"N_count_angacc_corre=   "<<N_count_angacc_corre<<endl;
	 //cin.get();

       tau_angacc=run_time-angacc_initial_time;
        
       
       for(i=0;i<n_atom;i++)
	 {
           int ii=(int)(p_y[i]/del_y_angacc_corre); //index_part_y_post[i]
	  
	if( ii==index_part_y_post[i]&& part_angacc_corre[i]>=0)
	 {
	   
	   angacc_corre_x[i] =angaccx_dash_initial[i]*angaccx_dash[i];
	   angacc_corre_y[i] =angaccy_dash_initial[i]*angaccy_dash[i];
	   angacc_corre_z[i] =angaccz_dash_initial[i]*angaccz_dash[i];
	   angacc_corre_xy[i]=angaccx_dash_initial[i]*angaccy_dash[i];
	   angacc_corre_yx[i]=angaccy_dash_initial[i]*angaccx_dash[i];
	   
	   total_part_count_angacc_corre[ii]=total_part_count_angacc_corre[ii]+1;

	   //cout<<"----angacc corre---   "<<ii<<"  "<<Nof_part_corre<<"    "<<angacc_corre_x[i]<<endl;
	 }
	else
	  {
	    part_angacc_corre[i]=-1;  // removing this particle from correlation list
	    angacc_corre_x[i]=0.0;   // these are put to zero, to be added to the sum, with no contribution.
	    angacc_corre_y[i]=0.0;
	    angacc_corre_z[i]=0.0;
	    angacc_corre_xy[i]=0.0;
	    angacc_corre_yx[i]=0.0; 

	  }
	 }




  //---- doing averaging over all the particles, considered for calculation--

       for(k=0;k<n_del_y_angacc_corre;k++)
	 {
       avg_angacc_corre_x[k]=0.0;
       avg_angacc_corre_y[k]=0.0;
       avg_angacc_corre_z[k]=0.0;
       avg_angacc_corre_xy[k]=0.0;
       avg_angacc_corre_yx[k]=0.0;

	 }

       for(i=0;i<n_atom;i++)
	 {
	   k=(int)(p_y[i]/del_y_angacc_corre);

       avg_angacc_corre_x[k]= avg_angacc_corre_x[k]+angacc_corre_x[i];
       avg_angacc_corre_y[k]= avg_angacc_corre_y[k]+angacc_corre_y[i];
       avg_angacc_corre_z[k]= avg_angacc_corre_z[k]+angacc_corre_z[i];
       avg_angacc_corre_xy[k]=avg_angacc_corre_xy[k]+angacc_corre_xy[i];
       avg_angacc_corre_yx[k]=avg_angacc_corre_yx[k]+angacc_corre_yx[i];

	 }

        for(k=0;k<n_del_y_angacc_corre;k++)
	  {
	    if(total_part_count_angacc_corre[k]>0)
	      {
	      avg_angacc_corre_x[k]= avg_angacc_corre_x[k]/total_part_count_angacc_corre[k];
	      avg_angacc_corre_y[k]= avg_angacc_corre_y[k]/total_part_count_angacc_corre[k];
	      avg_angacc_corre_z[k]= avg_angacc_corre_z[k]/total_part_count_angacc_corre[k];
	      avg_angacc_corre_xy[k]=avg_angacc_corre_xy[k]/total_part_count_angacc_corre[k];
	      avg_angacc_corre_yx[k]=avg_angacc_corre_yx[k]/total_part_count_angacc_corre[k];
	      }
	  }
       //--------------- storing the Normalization factor----------

     

       if(N_count_angacc_corre==0)
         for(k=0;k<n_del_y_angacc_corre;k++)
	   {
	    angavg_norm_x[k]=avg_angacc_corre_x[k];
	    angavg_norm_y[k]=avg_angacc_corre_y[k];
	    angavg_norm_z[k]=avg_angacc_corre_z[k];
	    angavg_norm_xy[k]=avg_angacc_corre_xy[k];
	    angavg_norm_yx[k]=avg_angacc_corre_yx[k];

	   }


       //------------Doing the Normalization of the correlation coefficient-----

      
       for(k=0;k<n_del_y_angacc_corre;k++)
	   {
       avg_angacc_corre_x[k]=avg_angacc_corre_x[k]/angavg_norm_x[k];
       avg_angacc_corre_y[k]=avg_angacc_corre_y[k]/angavg_norm_y[k];
       avg_angacc_corre_z[k]=avg_angacc_corre_z[k]/angavg_norm_z[k];
       avg_angacc_corre_xy[k]=avg_angacc_corre_xy[k]/angavg_norm_xy[k];
       avg_angacc_corre_yx[k]=avg_angacc_corre_yx[k]/angavg_norm_yx[k];
	   }
       

       for(k=0;k<n_del_y_angacc_corre;k++)
	   if(total_part_count_angacc_corre[k]<min_Nop_for_corre)
	 {
	   index_return_from_part_angacc_corre=1;
	   return;
	 }
       //cout<<"----angacc corre---   "<<"  "<<Nof_part_corre<<"    "<<angacc_corre_x[0]<<"  "<<norm_x[0]<<endl;
	      
  

   fp1016=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_corre_avg_x.txt","a");
   if(N_count_vel_corre==0)
     { 
       fprintf(fp1016,"%s  avg_angacc_norm_x[ at diff y]\n","%");
       fprintf(fp1016,"%s ","%");
	 for(i=0;i<n_del_y_angacc_corre;i++)
	   fprintf(fp1016, "%20.14lf", angavg_norm_x[i]);
       fprintf(fp1016,"\n");

    }
   fprintf(fp1016 ,"%20.14lf",tau_angacc);
 for(i=0;i<n_del_y_angacc_corre;i++)
   fprintf(fp1016 ,"%20.14lf", avg_angacc_corre_x[i]);
    fprintf(fp1016 ,"\n");
   fclose(fp1016);



 fp10161=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_corre_avg_y.txt","a");
   if(N_count_vel_corre==0)
     { 
       fprintf(fp10161,"%s  avg_angacc_norm_y[ at diff y]\n","%");
       fprintf(fp10161,"%s ","%");
	 for(i=0;i<n_del_y_angacc_corre;i++)
	   fprintf(fp10161, "%20.14lf", angavg_norm_y[i]);
       fprintf(fp10161,"\n\n");

    }
   fprintf(fp10161 ,"%20.14lf",tau_angacc);
 for(i=0;i<n_del_y_angacc_corre;i++)
   fprintf(fp10161 ,"%20.14lf", avg_angacc_corre_y[i]);
    fprintf(fp10161 ,"\n");
   fclose(fp10161);



 fp10162=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_corre_avg_z.txt","a");
   if(N_count_vel_corre==0)
     { 
       fprintf(fp10162,"%s  avg_angacc_norm_z[ at diff y]\n","%");
       fprintf(fp10162,"%s ","%");
	 for(i=0;i<n_del_y_angacc_corre;i++)
	   fprintf(fp10162, "%20.14lf", angavg_norm_z[i]);
       fprintf(fp10162,"\n");
    }
   fprintf(fp10162 ,"%20.14lf",tau_angacc);
 for(i=0;i<n_del_y_angacc_corre;i++)
   fprintf(fp10162 ,"%20.14lf", avg_angacc_corre_z[i]);
    fprintf(fp10162 ,"\n");
   fclose(fp10162);



 fp10163=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_corre_avg_xy.txt","a");
   if(N_count_vel_corre==0)
     { 
       fprintf(fp10163,"%s  avg_angacc_norm_xy[ at diff y]\n","%");
       fprintf(fp10163,"%s ","%");
	 for(i=0;i<n_del_y_angacc_corre;i++)
	   fprintf(fp10163, "%20.14lf", angavg_norm_xy[i]);
       fprintf(fp10163,"\n");

    }
   fprintf(fp10163 ,"%20.14lf",tau_angacc);
 for(i=0;i<n_del_y_angacc_corre;i++)
   fprintf(fp10163 ,"%20.14lf", avg_angacc_corre_xy[i]);
    fprintf(fp10163 ,"\n");
   fclose(fp10163);


 fp10164=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angacc_corre_avg_yx.txt","a");
   if(N_count_vel_corre==0)
     { 
       fprintf(fp10164,"%s  avg_angacc_norm_x[ at diff y]\n","%");
       fprintf(fp10164,"%s ","%");
	 for(i=0;i<n_del_y_angacc_corre;i++)
	   fprintf(fp10164, "%20.14lf", angavg_norm_yx[i]);
       fprintf(fp10164,"\n");

    }
   fprintf(fp10164 ,"%20.14lf",tau_angacc);
 for(i=0;i<n_del_y_angacc_corre;i++)
   fprintf(fp10164 ,"%20.14lf", avg_angacc_corre_yx[i]);
    fprintf(fp10164 ,"\n");
   fclose(fp10164);




   delete [] angacc_corre_x;
   delete [] angacc_corre_y;
   delete [] angacc_corre_z;
   delete [] angacc_corre_xy;
   delete [] angacc_corre_yx;


   delete[] avg_angacc_corre_x;
   delete[] avg_angacc_corre_y;
   delete[] avg_angacc_corre_z;
   delete[] avg_angacc_corre_xy;
   delete[] avg_angacc_corre_yx;
   
     }


//------------------------- Fnuction for air velocity correlation-- following the particle path--------

 void call_fluid_vel_corre(const FlowField& u,const Vector& x_grid,const Vector& z_grid)
{
  int i;
  double interpolated_u[3];
  int n_order=5;
  double a_=u.a();
  double b_=u.b();
  double ubase;

  double  avg_vx_air_at_py;
  double xp,yp,zp;


  if(index_return_from_air_vel_corre==1)
	 {
          return;
	 }



  for(i=0;i<n_atom;i++)
           {
	     
		  xp=p_x[i];
		  yp=p_y[i];
		  zp=p_z[i];
		 
    	 interpolation(u, x_grid,z_grid,xp,yp,zp, n_order, interpolated_u);
	 
	 
	 int ny_part=locate(y_air_grid,Max_dns_grid,p_y[i]);
	 avg_vx_air_at_py= splint(y_air_grid,avg_vx_air,y2_air_velo,ny_part,p_y[i]); // interpolation of average velocity
	 
	 // ubase =(p_y[i]-(a_+b_)/2.0)/((b_-a_)/2.0);    // for couette flow 
	 ubase = (1.0 - square(abs(p_y[i]-(b_+a_)/2.0)/((b_-a_)/2.0))); /// chennel flow
	 
	 vx_dash_air_vel_corre[i]=ubase+interpolated_u[0]-avg_vx_air_at_py;
	 vy_dash_air_vel_corre[i]=interpolated_u[1];  // nondimensionalized with centerline air velocity of prabolic profile
	 vz_dash_air_vel_corre[i]=interpolated_u[2];
	   }

	 if(N_count_fluid_vel_corre==0)set_fluid_vel_corre();



	 func_fluid_vel_corre();              //------ for calculationg particle accelerationfluid velocity corre.
	     N_count_fluid_vel_corre=N_count_fluid_vel_corre+1; 

}

void set_fluid_vel_corre()
{
  int i;
	  //int k=0;

	  del_y_air_vel_corre=b_ly/(double)n_del_y_air_vel_corre;

	  air_vel_initial_time=run_time;

	  //air_vel_corre_for_part=new int[Nmax_part_air_vel_corre];

	  for(i=0;i<n_atom;i++) consider_part_air_vel_corre[i]=1;


	  vx_dash_air_initial=new double[Nmax_part_air_vel_corre];
	  vy_dash_air_initial=new double[Nmax_part_air_vel_corre];
	  vz_dash_air_initial=new double[Nmax_part_air_vel_corre];
	  
	  avg_norm_air_x=new double[n_del_y_air_vel_corre];
	  avg_norm_air_y=new double[n_del_y_air_vel_corre];
	  avg_norm_air_z=new double[n_del_y_air_vel_corre];
	  avg_norm_air_xy=new double[n_del_y_air_vel_corre];
	  avg_norm_air_yx=new double[n_del_y_air_vel_corre];
	  
	  for(i=0;i<n_atom;i++)
	    {
	     
	      part_y_pos_for_air_vel_corre[i]=(int)(p_y[i]/del_y_air_vel_corre);
	      
	      vx_dash_air_initial[i]=vx_dash_air_vel_corre[i];
	      vy_dash_air_initial[i]=vy_dash_air_vel_corre[i];
	      vz_dash_air_initial[i]=vz_dash_air_vel_corre[i];
	     }

}




void func_fluid_vel_corre()// calculation of the fluid velocity correlation at the particle position
{

  FILE *fp260;
  FILE *fp261;
  FILE *fp262;
  FILE *fp263;
  FILE *fp264;

       int i,k;
       double *vel_corre_air_x,*vel_corre_air_y,*vel_corre_air_z,*vel_corre_air_xy,*vel_corre_air_yx;
       double tau_vel_air;
       double *avg_vel_corre_air_x,*avg_vel_corre_air_y,*avg_vel_corre_air_z,*avg_vel_corre_air_xy,*avg_vel_corre_air_yx;
       int *total_part_count_air_vel_corre;



       vel_corre_air_x=new double[Nmax_part_air_vel_corre];
       vel_corre_air_y=new double[Nmax_part_air_vel_corre]; 
       vel_corre_air_z=new double[Nmax_part_air_vel_corre];
       vel_corre_air_xy=new double[Nmax_part_air_vel_corre];
       vel_corre_air_yx=new double[Nmax_part_air_vel_corre];

       avg_vel_corre_air_x=new double[n_del_y_air_vel_corre];
       avg_vel_corre_air_y=new double[n_del_y_air_vel_corre];
       avg_vel_corre_air_z=new double[n_del_y_air_vel_corre];
       avg_vel_corre_air_xy=new double[n_del_y_air_vel_corre];
       avg_vel_corre_air_yx=new double[n_del_y_air_vel_corre];

       total_part_count_air_vel_corre=new int[n_del_y_air_vel_corre];

       for(i=0;i<n_del_y_air_vel_corre;i++)
	 total_part_count_air_vel_corre[i]=0;

       tau_vel_air=run_time-air_vel_initial_time;



       for(i=0;i<n_atom;i++)
	 {
           int ii=(int)(p_y[i]/del_y_air_vel_corre); //index_part_y_post[i]
	   
	if( ii==part_y_pos_for_air_vel_corre[i] && consider_part_air_vel_corre[i]>=0)
	 {
      vel_corre_air_x[i] =vx_dash_air_initial[i]*vx_dash_air_vel_corre[i];
      vel_corre_air_y[i] =vy_dash_air_initial[i]*vy_dash_air_vel_corre[i];
      vel_corre_air_z[i] =vz_dash_air_initial[i]*vz_dash_air_vel_corre[i];
      vel_corre_air_xy[i]=vx_dash_air_initial[i]*vy_dash_air_vel_corre[i];
      vel_corre_air_yx[i]=vy_dash_air_initial[i]*vx_dash_air_vel_corre[i];

      total_part_count_air_vel_corre[ii]= total_part_count_air_vel_corre[ii]+1;
	 
	 }
	else
	  {
	    consider_part_air_vel_corre[i]=-1;
	    
	    vel_corre_air_x[i] =0.0;
	    vel_corre_air_y[i] =0.0;
	    vel_corre_air_z[i] =0.0;
	    vel_corre_air_xy[i]=0.0;
	    vel_corre_air_yx[i]=0.0;
	  }
	 }

 //---- doing averaging over all the particle-positions, considered for calculation--

       for(k=0;k<n_del_y_air_vel_corre;k++)
	 {
       avg_vel_corre_air_x[k]=0.0;
       avg_vel_corre_air_y[k]=0.0;
       avg_vel_corre_air_z[k]=0.0;
       avg_vel_corre_air_xy[k]=0.0;
       avg_vel_corre_air_yx[k]=0.0;

	 }

       for(i=0;i<n_atom;i++)
	 {
	   k=(int)(p_y[i]/del_y_air_vel_corre);
	   
	   avg_vel_corre_air_x[k]= avg_vel_corre_air_x[k]+vel_corre_air_x[i];
	    avg_vel_corre_air_y[k]= avg_vel_corre_air_y[k]+vel_corre_air_y[i];
	     avg_vel_corre_air_z[k]= avg_vel_corre_air_z[k]+vel_corre_air_z[i];
	      avg_vel_corre_air_xy[k]= avg_vel_corre_air_xy[k]+vel_corre_air_xy[i];
	       avg_vel_corre_air_yx[k]= avg_vel_corre_air_yx[k]+vel_corre_air_yx[i];
	 }



       for(k=0;k<n_del_y_air_vel_corre;k++)
	  {
	    if(total_part_count_air_vel_corre[k]>0)
	      {
	      avg_vel_corre_air_x[k]= avg_vel_corre_air_x[k]/total_part_count_air_vel_corre[k];
	      avg_vel_corre_air_y[k]= avg_vel_corre_air_y[k]/total_part_count_air_vel_corre[k];
	      avg_vel_corre_air_z[k]= avg_vel_corre_air_z[k]/total_part_count_air_vel_corre[k];
	      avg_vel_corre_air_xy[k]= avg_vel_corre_air_xy[k]/total_part_count_air_vel_corre[k];
	      avg_vel_corre_air_yx[k]= avg_vel_corre_air_yx[k]/total_part_count_air_vel_corre[k];

	      }
	  }



       if(N_count_fluid_vel_corre==0)
         for(k=0;k<n_del_y_air_vel_corre;k++)
	   {
	     avg_norm_air_x[k]=avg_vel_corre_air_x[k];
	     avg_norm_air_y[k]=avg_vel_corre_air_y[k];
	     avg_norm_air_z[k]=avg_vel_corre_air_z[k];
	     avg_norm_air_xy[k]=avg_vel_corre_air_xy[k];
	     avg_norm_air_yx[k]=avg_vel_corre_air_yx[k];
	   }



       //------------Doing the Normalization of the correlation coefficient-----

      
       for(k=0;k<n_del_y_acc_corre;k++)
	   {
	     avg_vel_corre_air_x[k]=avg_vel_corre_air_x[k]/avg_norm_air_x[k];
	     avg_vel_corre_air_y[k]=avg_vel_corre_air_y[k]/avg_norm_air_y[k];
	     avg_vel_corre_air_z[k]=avg_vel_corre_air_z[k]/avg_norm_air_z[k];
	     avg_vel_corre_air_xy[k]=avg_vel_corre_air_xy[k]/avg_norm_air_xy[k];	 
	     avg_vel_corre_air_yx[k]=avg_vel_corre_air_yx[k]/avg_norm_air_yx[k];

	   }
       for(k=0;k<n_del_y_air_vel_corre;k++)
	   if(total_part_count_air_vel_corre[k]<min_Nop_for_corre)
	 {
	   index_return_from_air_vel_corre=1;
	   return;
	 }

 fp260=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_corre_avg_x.txt","a");
   if(N_count_fluid_vel_corre==0)
     { 
       fprintf(fp260,"%s  avg_air_vel_corre_norm_x[ at diff y]\n","%");
       fprintf(fp260,"%s ","%");
	 for(i=0;i<n_del_y_air_vel_corre;i++)
	   fprintf(fp260, "%20.14lf", avg_norm_air_x[i]);
       fprintf(fp260,"\n\n");

     }
   fprintf(fp260 ,"%20.14lf",tau_vel_air);
  for(i=0;i<n_del_y_air_vel_corre;i++)
   fprintf(fp260 ,"%20.14lf", avg_vel_corre_air_x[i]);
    fprintf(fp260 ,"\n");
   fclose(fp260);


 fp261=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_corre_avg_y.txt","a");
   if(N_count_fluid_vel_corre==0)
     { 
       fprintf(fp261,"%s  avg_air_vel_corre_norm_y[ at diff y]\n","%");
       fprintf(fp261,"%s ","%");
	 for(i=0;i<n_del_y_air_vel_corre;i++)
	   fprintf(fp261, "%20.14lf", avg_norm_air_y[i]);
       fprintf(fp261,"\n\n");

    }
   fprintf(fp261 ,"%20.14lf",tau_vel_air);
  for(i=0;i<n_del_y_air_vel_corre;i++)
   fprintf(fp261 ,"%20.14lf", avg_vel_corre_air_y[i]);
    fprintf(fp261 ,"\n");
   fclose(fp261);





 fp262=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_corre_avg_z.txt","a");
   if(N_count_fluid_vel_corre==0)
     { 
       fprintf(fp262,"%s  avg_air_vel_corre_norm_z[ at diff y]\n","%");
       fprintf(fp262,"%s ","%");
	 for(i=0;i<n_del_y_air_vel_corre;i++)
	   fprintf(fp262, "%20.14lf", avg_norm_air_z[i]);
       fprintf(fp262,"\n\n");

    }
   fprintf(fp262 ,"%20.14lf",tau_vel_air);
  for(i=0;i<n_del_y_air_vel_corre;i++)
   fprintf(fp262 ,"%20.14lf", avg_vel_corre_air_z[i]);
    fprintf(fp262 ,"\n");
   fclose(fp262);

   


 fp263=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_corre_avg_xy.txt","a");
   if(N_count_fluid_vel_corre==0)
     { 
       fprintf(fp263,"%s  avg_air_vel_corre_norm_xy[ at diff y]\n","%");
       fprintf(fp263,"%s ","%");
	 for(i=0;i<n_del_y_air_vel_corre;i++)
	   fprintf(fp263, "%20.14lf", avg_norm_air_xy[i]);
       fprintf(fp263,"\n\n");

    }
   fprintf(fp263 ,"%20.14lf",tau_vel_air);
  for(i=0;i<n_del_y_air_vel_corre;i++)
   fprintf(fp263 ,"%20.14lf", avg_vel_corre_air_xy[i]);
    fprintf(fp263 ,"\n");
   fclose(fp263);

   



 fp264=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_vel_corre_avg_yx.txt","a");
   if(N_count_fluid_vel_corre==0)
     { 
       fprintf(fp264,"%s  avg_air_vel_corre_norm_yx[ at diff y]\n","%");
       fprintf(fp264,"%s ","%");
	 for(i=0;i<n_del_y_air_vel_corre;i++)
	   fprintf(fp264, "%20.14lf", avg_norm_air_yx[i]);
       fprintf(fp264,"\n\n");

    }
   fprintf(fp264 ,"%20.14lf",tau_vel_air);
  for(i=0;i<n_del_y_air_vel_corre;i++)
   fprintf(fp264 ,"%20.14lf", avg_vel_corre_air_yx[i]);
    fprintf(fp264 ,"\n");
   fclose(fp264);



   delete [] vel_corre_air_x;
   delete [] vel_corre_air_y;
   delete [] vel_corre_air_z;
   delete [] vel_corre_air_xy;
   delete [] vel_corre_air_yx;

   delete []  avg_vel_corre_air_x;
   delete []  avg_vel_corre_air_y;
   delete []  avg_vel_corre_air_z; 
   delete []  avg_vel_corre_air_xy; 
   delete []  avg_vel_corre_air_yx; 
  


}

void call_fluid_angular_vel_corre(const FlowField& omega,const Vector& x_grid,const Vector& z_grid)
{
  int i;
  double interpolated_omega[3];
  int n_order=5;
  //double a_=omega.a();
  //double b_=omega.b();
  //double ubase;

  //double  avg_angvelx_air_at_py;
  double xp,yp,zp;


  if(index_return_from_air_angvel_corre==1)
	 {
          return;
	 }



  for(i=0;i<n_atom;i++)
           {
	     
		  xp=p_x[i];
		  yp=p_y[i];
		  zp=p_z[i];
		 
    	 interpolation(omega, x_grid,z_grid,xp,yp,zp, n_order, interpolated_omega);
	 
	 
	// int ny_part=locate(y_air_grid,Max_dns_grid,p_y[i]);
	// avg_angvelx_air_at_py= splint(y_air_grid,avg_angvelx_air,y2_air_angvelo,ny_part,p_y[i]); // interpolation of average angular velocity
	 
	 // ubase =(p_y[i]-(a_+b_)/2.0)/((b_-a_)/2.0);    // for couette flow 
	 //ubase = (1.0 - square(abs(p_y[i]-(b_+a_)/2.0)/((b_-a_)/2.0))); /// chennel flow
	 
	 angvelx_dash_air_angvel_corre[i]=0.5*interpolated_omega[0];//-avg_angvelx_air_at_py;
	 angvely_dash_air_angvel_corre[i]=0.5*interpolated_omega[1];  // nondimensionalized with centerline air angular velocity of prabolic profile
	 angvelz_dash_air_angvel_corre[i]=0.5*interpolated_omega[2];
	   }

	 if(N_count_fluid_angvel_corre==0)set_fluid_angular_vel_corre();



	 func_fluid_angular_vel_corre();              //------ for calculationg particle accelerationfluid angular velocity corre.
	     N_count_fluid_angvel_corre=N_count_fluid_angvel_corre+1; 

}

void set_fluid_angular_vel_corre()
{
  int i;
	  //int k=0;

	  del_y_air_angvel_corre=b_ly/(double)n_del_y_air_angvel_corre;

	  air_angvel_initial_time=run_time;

	  //air_angvel_corre_for_part=new int[Nmax_part_air_angvel_corre];

	  for(i=0;i<n_atom;i++) consider_part_air_angvel_corre[i]=1;


	  angvelx_dash_air_initial=new double[Nmax_part_air_angvel_corre];
	  angvely_dash_air_initial=new double[Nmax_part_air_angvel_corre];
	  angvelz_dash_air_initial=new double[Nmax_part_air_angvel_corre];
	  
	  angavg_norm_air_x=new double[n_del_y_air_angvel_corre];
	  angavg_norm_air_y=new double[n_del_y_air_angvel_corre];
	  angavg_norm_air_z=new double[n_del_y_air_angvel_corre];
	  angavg_norm_air_xy=new double[n_del_y_air_angvel_corre];
	  angavg_norm_air_yx=new double[n_del_y_air_angvel_corre];
	  
	  for(i=0;i<n_atom;i++)
	    {
	     
	      part_y_pos_for_air_angvel_corre[i]=(int)(p_y[i]/del_y_air_angvel_corre);
	      
	      angvelx_dash_air_initial[i]=angvelx_dash_air_angvel_corre[i];
	      angvely_dash_air_initial[i]=angvely_dash_air_angvel_corre[i];
	      angvelz_dash_air_initial[i]=angvelz_dash_air_angvel_corre[i];
	     }

}




void func_fluid_angular_vel_corre()// calculation of the fluid angular velocity correlation at the particle position
{

  FILE *fp10260;
  FILE *fp10261;
  FILE *fp10262;
  FILE *fp10263;
  FILE *fp10264;

       int i,k;
       double *angvel_corre_air_x,*angvel_corre_air_y,*angvel_corre_air_z,*angvel_corre_air_xy,*angvel_corre_air_yx;
       double tau_vel_air;
       double *avg_angvel_corre_air_x,*avg_angvel_corre_air_y,*avg_angvel_corre_air_z,*avg_angvel_corre_air_xy,*avg_angvel_corre_air_yx;
       int *total_part_count_air_angvel_corre;



       angvel_corre_air_x=new double[Nmax_part_air_angvel_corre];
       angvel_corre_air_y=new double[Nmax_part_air_angvel_corre]; 
       angvel_corre_air_z=new double[Nmax_part_air_angvel_corre];
       angvel_corre_air_xy=new double[Nmax_part_air_angvel_corre];
       angvel_corre_air_yx=new double[Nmax_part_air_angvel_corre];

       avg_angvel_corre_air_x=new double[n_del_y_air_angvel_corre];
       avg_angvel_corre_air_y=new double[n_del_y_air_angvel_corre];
       avg_angvel_corre_air_z=new double[n_del_y_air_angvel_corre];
       avg_angvel_corre_air_xy=new double[n_del_y_air_angvel_corre];
       avg_angvel_corre_air_yx=new double[n_del_y_air_angvel_corre];

       total_part_count_air_angvel_corre=new int[n_del_y_air_angvel_corre];

       for(i=0;i<n_del_y_air_angvel_corre;i++)
	 total_part_count_air_angvel_corre[i]=0;

       tau_vel_air=run_time-air_vel_initial_time;



       for(i=0;i<n_atom;i++)
	 {
           int ii=(int)(p_y[i]/del_y_air_angvel_corre); //index_part_y_post[i]
	   
	if( ii==part_y_pos_for_air_angvel_corre[i] && consider_part_air_angvel_corre[i]>=0)
	 {
      angvel_corre_air_x[i] =angvelx_dash_air_initial[i]*angvelx_dash_air_angvel_corre[i];
      angvel_corre_air_y[i] =angvely_dash_air_initial[i]*angvely_dash_air_angvel_corre[i];
      angvel_corre_air_z[i] =angvelz_dash_air_initial[i]*angvelz_dash_air_angvel_corre[i];
      angvel_corre_air_xy[i]=angvelx_dash_air_initial[i]*angvely_dash_air_angvel_corre[i];
      angvel_corre_air_yx[i]=angvely_dash_air_initial[i]*angvelx_dash_air_angvel_corre[i];

      total_part_count_air_angvel_corre[ii]= total_part_count_air_angvel_corre[ii]+1;
	 
	 }
	else
	  {
	    consider_part_air_angvel_corre[i]=-1;
	    
	    angvel_corre_air_x[i] =0.0;
	    angvel_corre_air_y[i] =0.0;
	    angvel_corre_air_z[i] =0.0;
	    angvel_corre_air_xy[i]=0.0;
	    angvel_corre_air_yx[i]=0.0;
	  }
	 }

 //---- doing averaging over all the particle-positions, considered for calculation--

       for(k=0;k<n_del_y_air_angvel_corre;k++)
	 {
       avg_angvel_corre_air_x[k]=0.0;
       avg_angvel_corre_air_y[k]=0.0;
       avg_angvel_corre_air_z[k]=0.0;
       avg_angvel_corre_air_xy[k]=0.0;
       avg_angvel_corre_air_yx[k]=0.0;

	 }

       for(i=0;i<n_atom;i++)
	 {
	   k=(int)(p_y[i]/del_y_air_angvel_corre);
	   
	   avg_angvel_corre_air_x[k]= avg_angvel_corre_air_x[k]+angvel_corre_air_x[i];
	    avg_angvel_corre_air_y[k]= avg_angvel_corre_air_y[k]+angvel_corre_air_y[i];
	     avg_angvel_corre_air_z[k]= avg_angvel_corre_air_z[k]+angvel_corre_air_z[i];
	      avg_angvel_corre_air_xy[k]= avg_angvel_corre_air_xy[k]+angvel_corre_air_xy[i];
	       avg_angvel_corre_air_yx[k]= avg_angvel_corre_air_yx[k]+angvel_corre_air_yx[i];
	 }



       for(k=0;k<n_del_y_air_angvel_corre;k++)
	  {
	    if(total_part_count_air_angvel_corre[k]>0)
	      {
	      avg_angvel_corre_air_x[k]= avg_angvel_corre_air_x[k]/total_part_count_air_angvel_corre[k];
	      avg_angvel_corre_air_y[k]= avg_angvel_corre_air_y[k]/total_part_count_air_angvel_corre[k];
	      avg_angvel_corre_air_z[k]= avg_angvel_corre_air_z[k]/total_part_count_air_angvel_corre[k];
	      avg_angvel_corre_air_xy[k]= avg_angvel_corre_air_xy[k]/total_part_count_air_angvel_corre[k];
	      avg_angvel_corre_air_yx[k]= avg_angvel_corre_air_yx[k]/total_part_count_air_angvel_corre[k];

	      }
	  }



       if(N_count_fluid_angvel_corre==0)
         for(k=0;k<n_del_y_air_angvel_corre;k++)
	   {
	     angavg_norm_air_x[k]=avg_angvel_corre_air_x[k];
	     angavg_norm_air_y[k]=avg_angvel_corre_air_y[k];
	     angavg_norm_air_z[k]=avg_angvel_corre_air_z[k];
	     angavg_norm_air_xy[k]=avg_angvel_corre_air_xy[k];
	     angavg_norm_air_yx[k]=avg_angvel_corre_air_yx[k];
	   }



       //------------Doing the Normalization of the correlation coefficient-----

      
       for(k=0;k<n_del_y_acc_corre;k++)
	   {
	     avg_angvel_corre_air_x[k]=avg_angvel_corre_air_x[k]/angavg_norm_air_x[k];
	     avg_angvel_corre_air_y[k]=avg_angvel_corre_air_y[k]/angavg_norm_air_y[k];
	     avg_angvel_corre_air_z[k]=avg_angvel_corre_air_z[k]/angavg_norm_air_z[k];
	     avg_angvel_corre_air_xy[k]=avg_angvel_corre_air_xy[k]/angavg_norm_air_xy[k];	 
	     avg_angvel_corre_air_yx[k]=avg_angvel_corre_air_yx[k]/angavg_norm_air_yx[k];

	   }
       for(k=0;k<n_del_y_air_angvel_corre;k++)
	   if(total_part_count_air_angvel_corre[k]<min_Nop_for_corre)
	 {
	   index_return_from_air_angvel_corre=1;
	   return;
	 }

 fp10260=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_angvel_corre_avg_x.txt","a");
   if(N_count_fluid_angvel_corre==0)
     { 
       fprintf(fp10260,"%s  avg_air_angvel_corre_norm_x[ at diff y]\n","%");
       fprintf(fp10260,"%s ","%");
	 for(i=0;i<n_del_y_air_angvel_corre;i++)
	   fprintf(fp10260, "%20.14lf", angavg_norm_air_x[i]);
       fprintf(fp10260,"\n\n");

     }
   fprintf(fp10260 ,"%20.14lf",tau_vel_air);
  for(i=0;i<n_del_y_air_angvel_corre;i++)
   fprintf(fp10260 ,"%20.14lf", avg_angvel_corre_air_x[i]);
    fprintf(fp10260 ,"\n");
   fclose(fp10260);


 fp10261=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_angvel_corre_avg_y.txt","a");
   if(N_count_fluid_angvel_corre==0)
     { 
       fprintf(fp10261,"%s  avg_air_angvel_corre_norm_y[ at diff y]\n","%");
       fprintf(fp10261,"%s ","%");
	 for(i=0;i<n_del_y_air_angvel_corre;i++)
	   fprintf(fp10261, "%20.14lf", angavg_norm_air_y[i]);
       fprintf(fp10261,"\n\n");

    }
   fprintf(fp10261 ,"%20.14lf",tau_vel_air);
  for(i=0;i<n_del_y_air_angvel_corre;i++)
   fprintf(fp10261 ,"%20.14lf", avg_angvel_corre_air_y[i]);
    fprintf(fp10261 ,"\n");
   fclose(fp10261);





 fp10262=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_angvel_corre_avg_z.txt","a");
   if(N_count_fluid_angvel_corre==0)
     { 
       fprintf(fp10262,"%s  avg_air_angvel_corre_norm_z[ at diff y]\n","%");
       fprintf(fp10262,"%s ","%");
	 for(i=0;i<n_del_y_air_angvel_corre;i++)
	   fprintf(fp10262, "%20.14lf", angavg_norm_air_z[i]);
       fprintf(fp10262,"\n\n");

    }
   fprintf(fp10262 ,"%20.14lf",tau_vel_air);
  for(i=0;i<n_del_y_air_angvel_corre;i++)
   fprintf(fp10262 ,"%20.14lf", avg_angvel_corre_air_z[i]);
    fprintf(fp10262 ,"\n");
   fclose(fp10262);

   


 fp10263=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_angvel_corre_avg_xy.txt","a");
   if(N_count_fluid_angvel_corre==0)
     { 
       fprintf(fp10263,"%s  avg_air_angvel_corre_norm_xy[ at diff y]\n","%");
       fprintf(fp10263,"%s ","%");
	 for(i=0;i<n_del_y_air_angvel_corre;i++)
	   fprintf(fp10263, "%20.14lf", angavg_norm_air_xy[i]);
       fprintf(fp10263,"\n\n");

    }
   fprintf(fp10263 ,"%20.14lf",tau_vel_air);
  for(i=0;i<n_del_y_air_angvel_corre;i++)
   fprintf(fp10263 ,"%20.14lf", avg_angvel_corre_air_xy[i]);
    fprintf(fp10263 ,"\n");
   fclose(fp10263);

   



 fp10264=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/air_angvel_corre_avg_yx.txt","a");
   if(N_count_fluid_angvel_corre==0)
     { 
       fprintf(fp10264,"%s  avg_air_angvel_corre_norm_yx[ at diff y]\n","%");
       fprintf(fp10264,"%s ","%");
	 for(i=0;i<n_del_y_air_angvel_corre;i++)
	   fprintf(fp10264, "%20.14lf", angavg_norm_air_yx[i]);
       fprintf(fp10264,"\n\n");

    }
   fprintf(fp10264 ,"%20.14lf",tau_vel_air);
  for(i=0;i<n_del_y_air_angvel_corre;i++)
   fprintf(fp10264 ,"%20.14lf", avg_angvel_corre_air_yx[i]);
    fprintf(fp10264 ,"\n");
   fclose(fp10264);



   delete [] angvel_corre_air_x;
   delete [] angvel_corre_air_y;
   delete [] angvel_corre_air_z;
   delete [] angvel_corre_air_xy;
   delete [] angvel_corre_air_yx;

   delete []  avg_angvel_corre_air_x;
   delete []  avg_angvel_corre_air_y;
   delete []  avg_angvel_corre_air_z; 
   delete []  avg_angvel_corre_air_xy; 
   delete []  avg_angvel_corre_air_yx; 
  


}



//------------------------- Fnuction for particle velocity correlation------------------------------


void call_part_velo_corre()
{
  int ny;
  int i;
double  avg_vx_part_at_py;

 if(index_return_from_corre==1)
	 {
          return;
	 }

 for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
    Nofp_inst[ny]=0;
    vx_inst[ny]=0.0;
    vy_inst[ny]=0.0;
    vz_inst[ny]=0.0;

       }

 if(index_average_from_input==1)
     {
       for(i=0;i<n_atom;i++)
	 {
	    ny=(int)(((p_y[i]-(sigma/2.0))/del_y));     
       if(p_y[i]>y_min_for_corre && p_y[i]<y_max_for_corre )   
	 {

       int ny1_part=locate(y_part_grid,Ndel_y,p_y[i]);
       avg_vx_part_at_py= splint(y_part_grid,avg_vx_part,y2_part_velo,ny1_part,p_y[i]);

       pvx_dash[i]=p_vx[i]-avg_vx_part_at_py;  // calculating fluctuation for each particle
       pvy_dash[i]=p_vy[i];
       pvz_dash[i]=p_vz[i];
	 }
	 }

     }


 else{
 for(i=0;i<n_atom;i++)// summing up all the velocities in each the y-grid for all the particle 
    {
      ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

      // cout<<"in property"<<p_y[i]<<" "<<ny<<endl;


              Nofp_inst[ny]= Nofp_inst[ny]+ 1;

              vx_inst[ny]= vx_inst[ny]+p_vx[i];
              vy_inst[ny]= vy_inst[ny]+p_vy[i];  
              vz_inst[ny]= vz_inst[ny]+p_vz[i];
      }

 for(ny=0;ny<Ndel_y;ny++) // doing the particle average in each y-grid-in each sample frame
       {
    if(Nofp_inst[ny]>0.0)
      {
     vx_inst[ny]= vx_inst[ny]/Nofp_inst[ny];  // vx_inst[ny] is basically particle averaged inst in time.
     vy_inst[ny]= vy_inst[ny]/Nofp_inst[ny];
     vz_inst[ny]= vz_inst[ny]/Nofp_inst[ny];
      }
       }

 for(i=0;i<n_atom;i++)
   {
     ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
     pvx_dash[i]=p_vx[i]-vx_inst[ny];  // calculating fluctuation for each particle
     pvy_dash[i]=p_vy[i]-vy_inst[ny];
     pvz_dash[i]=p_vz[i]-vz_inst[ny];
   }

 } //end of the else loop

   if(N_count_vel_corre==0)set_vel_corre();
     func_part_vel_corre();              //------ for calculationg particle acceleration corre.
     N_count_vel_corre=N_count_vel_corre+1;

}

 void set_vel_corre()
 {
  int i;
  int kk=0;

  part_vel_corre=new int[Nmax_part_vel_corre];
  pvx_dash_initial=new double[Nmax_part_vel_corre];
  pvy_dash_initial=new double[Nmax_part_vel_corre];
  pvz_dash_initial=new double[Nmax_part_vel_corre];
  
  


  // ------select the particles to consider for acceleration correlation---

	  for(i=0; i<n_atom;i++){

	    if(p_y[i]>y_min_for_corre && p_y[i]<y_max_for_corre && kk<Nmax_part_vel_corre)
	          {
		  part_vel_corre[kk]=i;   //-- particle to be considered for acc- correlation calculation
		  kk=kk+1;
		  Nof_part_vel_corre=kk;
		    }
	  }

   if(Nof_part_vel_corre==0)
	    { // --if by chanse there is no atom in selected region
	      //----(contd) for safety factor it is given
 
	      Nof_part_vel_corre=1;
	      part_vel_corre[0]=n_atom/2;
		}

	  vel_corre_initial_time=run_time;

	  for(i=0;i<Nof_part_vel_corre;i++)
	    {
	      int ii=part_vel_corre[i];
	      pvx_dash_initial[i]=pvx_dash[ii];
	      pvy_dash_initial[i]=pvy_dash[ii];
	      pvz_dash_initial[i]=pvz_dash[ii];
	    }


	  return;
}

  void func_part_vel_corre()
  {

  int i;
  double *vel_corre_x,*vel_corre_y,*vel_corre_z,*vel_corre_xy,*vel_corre_yx;
  double avg_vel_corre_x,avg_vel_corre_y,avg_vel_corre_z,avg_vel_corre_xy,avg_vel_corre_yx; 
  double tau_vel_corre;
  int total_part_count_vel_corre;

         vel_corre_x=new double[Nmax_part_vel_corre];
          vel_corre_y=new double[Nmax_part_vel_corre];
	   vel_corre_z=new double[Nmax_part_vel_corre];
	    vel_corre_xy=new double[Nmax_part_vel_corre];
	     vel_corre_yx=new double[Nmax_part_vel_corre];

	     tau_vel_corre=run_time- vel_corre_initial_time;


	     total_part_count_vel_corre=0;

      for(i=0;i<Nof_part_vel_corre;i++)
	       {
		 int ii=part_vel_corre[i];
        if(p_y[ii]>y_min_for_corre && p_y[ii]<y_max_for_corre && ii>=0)
	   {
	   
	   vel_corre_x[i] =pvx_dash_initial[i]*pvx_dash[ii];
	   vel_corre_y[i] =pvy_dash_initial[i]*pvy_dash[ii];
	   vel_corre_z[i] =pvz_dash_initial[i]*pvz_dash[ii];
	   vel_corre_xy[i]=pvx_dash_initial[i]*pvy_dash[ii];
	   vel_corre_yx[i]=pvy_dash_initial[i]*pvx_dash[ii];
	   
	   total_part_count_vel_corre=total_part_count_vel_corre+1;

	   //cout<<"----acc corre---   "<<ii<<"  "<<Nof_part_corre<<"    "<<acc_corre_x[i]<<endl;
	   }
	else
	  {
	    part_vel_corre[i]=-1;  // removing this particle from correlation list
	    vel_corre_x[i]=0.0;   // these are put to zero, to be added to the sum, with no contribution.
	    vel_corre_y[i]=0.0;
	    vel_corre_z[i]=0.0;
	    vel_corre_xy[i]=0.0;
	    vel_corre_yx[i]=0.0; 

	  }
	       }


       avg_vel_corre_x=0.0;
       avg_vel_corre_y=0.0;
       avg_vel_corre_z=0.0;
       avg_vel_corre_xy=0.0;
       avg_vel_corre_yx=0.0;

      
       for(i=0;i<Nof_part_vel_corre;i++)
	 { 
	         
       avg_vel_corre_x= avg_vel_corre_x + vel_corre_x[i];
       avg_vel_corre_y= avg_vel_corre_y + vel_corre_y[i];
       avg_vel_corre_z= avg_vel_corre_z + vel_corre_z[i];
       avg_vel_corre_xy=avg_vel_corre_xy + vel_corre_xy[i];
       avg_vel_corre_yx=avg_vel_corre_yx + vel_corre_yx[i];

	 }

       if(total_part_count_vel_corre>0)
	 {
     avg_vel_corre_x= avg_vel_corre_x/total_part_count_vel_corre;
     avg_vel_corre_y= avg_vel_corre_y/total_part_count_vel_corre;
     avg_vel_corre_z= avg_vel_corre_z/total_part_count_vel_corre;
     avg_vel_corre_xy= avg_vel_corre_xy/total_part_count_vel_corre;
     avg_vel_corre_yx= avg_vel_corre_yx/total_part_count_vel_corre;
	 }
 //--------------- storing the Normalization factor----------

     if(N_count_vel_corre==0)
	 {
	  avg_vel_norm_x=avg_vel_corre_x;
	  avg_vel_norm_y=avg_vel_corre_y;
	  avg_vel_norm_z=avg_vel_corre_z;
	  avg_vel_norm_xy=avg_vel_corre_xy;
	  avg_vel_norm_yx=avg_vel_corre_yx;
	 }

     
//------------Doing the Normalization of the correlation coefficient-----

     avg_vel_corre_x=avg_vel_corre_x/avg_vel_norm_x;
     avg_vel_corre_y=avg_vel_corre_y/avg_vel_norm_y;
     avg_vel_corre_z=avg_vel_corre_z/avg_vel_norm_z;
     avg_vel_corre_xy=avg_vel_corre_xy/avg_vel_norm_xy;
     avg_vel_corre_yx=avg_vel_corre_yx/avg_vel_norm_yx;

     if(total_part_count_vel_corre<min_Nop_for_corre)
       {
	 index_return_from_corre=1;
          return;
       }

fp19=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/vel_corre_avg.txt","a");

 if(N_count_vel_corre==0) fprintf(fp19,"%s  avg_vel_norm_x , y, z, xy, yx are = %20.14lf %20.14lf %20.14lf %20.14lf %20.14lf\n","%",
				  avg_vel_norm_x,avg_vel_norm_y,avg_vel_norm_z,avg_vel_norm_xy,avg_vel_norm_yx);

   fprintf(fp19 ,"%20.14lf %20.14lf %20.14lf %20.14lf %20.14lf %20.14lf \n",tau_vel_corre, avg_vel_corre_x, 
avg_vel_corre_y, avg_vel_corre_z,  avg_vel_corre_xy, avg_vel_corre_yx);

   fclose(fp19);


   delete [] vel_corre_x;
   delete [] vel_corre_y;
   delete [] vel_corre_z;
   delete [] vel_corre_xy;
   delete [] vel_corre_yx;
   return;

}

//------------------------- Fnuction for particle angular velocity correlation------------------------------


void call_part_angular_velo_corre()
{
  int ny;
  int i;
double  avg_angvelx_part_at_py;

 if(angular_index_return_from_corre==1)
	 {
          return;
	 }

 for(ny=0;ny<Ndel_y;ny++) // initiallizing the intantaneous velocities
       {
    Nofp_inst[ny]=0;
    angvelx_inst[ny]=0.0;
    angvely_inst[ny]=0.0;
    angvelz_inst[ny]=0.0;

       }

 if(index_average_from_input==1)
     {
       for(i=0;i<n_atom;i++)
	 {
	    ny=(int)(((p_y[i]-(sigma/2.0))/del_y));     
       if(p_y[i]>y_min_for_corre && p_y[i]<y_max_for_corre )   
	 {

       int ny1_part=locate(y_part_grid,Ndel_y,p_y[i]);
       avg_angvelx_part_at_py= splint(y_part_grid,avg_angvelx_part,y2_part_velo,ny1_part,p_y[i]);

       pangvelx_dash[i]=p_angvelx[i]-avg_angvelx_part_at_py;  // calculating fluctuation for each particle
       pangvely_dash[i]=p_angvely[i];
       pangvelz_dash[i]=p_angvelz[i];
	 }
	 }

     }


 else{
 for(i=0;i<n_atom;i++)// summing up all the velocities in each the y-grid for all the particle 
    {
      ny=(int)(((p_y[i]-(sigma/2.0))/del_y));

      // cout<<"in property"<<p_y[i]<<" "<<ny<<endl;


              Nofp_inst[ny]= Nofp_inst[ny]+ 1;

              angvelx_inst[ny]= angvelx_inst[ny]+p_angvelx[i];
              angvely_inst[ny]= angvely_inst[ny]+p_angvely[i];  
              angvelz_inst[ny]= angvelz_inst[ny]+p_angvelz[i];
      }

 for(ny=0;ny<Ndel_y;ny++) // doing the particle average in each y-grid-in each sample frame
       {
    if(Nofp_inst[ny]>0.0)
      {
     angvelx_inst[ny]= angvelx_inst[ny]/Nofp_inst[ny];  // angvelx_inst[ny] is basically particle averaged inst in time.
     angvely_inst[ny]= angvely_inst[ny]/Nofp_inst[ny];
     angvelz_inst[ny]= angvelz_inst[ny]/Nofp_inst[ny];
      }
       }

 for(i=0;i<n_atom;i++)
   {
     ny=(int)(((p_y[i]-(sigma/2.0))/del_y));
     pangvelx_dash[i]=p_angvelx[i]-angvelx_inst[ny];  // calculating fluctuation for each particle
     pangvely_dash[i]=p_angvely[i]-angvely_inst[ny];
     pangvelz_dash[i]=p_angvelz[i]-angvelz_inst[ny];
   }

 } //end of the else loop

   if(N_count_angvel_corre==0)set_angular_vel_corre();
     func_part_angular_vel_corre();              //------ for calculationg particle acceleration corre.
     N_count_angvel_corre=N_count_angvel_corre+1;

}

 void set_angular_vel_corre()
 {
  int i;
  int kk=0;

  part_angvel_corre=new int[Nmax_part_angvel_corre];
  pangvelx_dash_initial=new double[Nmax_part_angvel_corre];
  pangvely_dash_initial=new double[Nmax_part_angvel_corre];
  pangvelz_dash_initial=new double[Nmax_part_angvel_corre];
  
  


  // ------select the particles to consider for acceleration correlation---

	  for(i=0; i<n_atom;i++){

	    if(p_y[i]>y_min_for_corre && p_y[i]<y_max_for_corre && kk<Nmax_part_angvel_corre)
	          {
		  part_angvel_corre[kk]=i;   //-- particle to be considered for acc- correlation calculation
		  kk=kk+1;
		  Nof_part_angvel_corre=kk;
		    }
	  }

   if(Nof_part_angvel_corre==0)
	    { // --if by chanse there is no atom in selected region
	      //----(contd) for safety factor it is given
 
	      Nof_part_angvel_corre=1;
	      part_angvel_corre[0]=n_atom/2;
		}

	  angvel_corre_initial_time=run_time;

	  for(i=0;i<Nof_part_angvel_corre;i++)
	    {
	      int ii=part_angvel_corre[i];
	      pangvelx_dash_initial[i]=pangvelx_dash[ii];
	      pangvely_dash_initial[i]=pangvely_dash[ii];
	      pangvelz_dash_initial[i]=pangvelz_dash[ii];
	    }


	  return;
}

  void func_part_angular_vel_corre()
  {

  int i;
  double *angvel_corre_x,*angvel_corre_y,*angvel_corre_z,*angvel_corre_xy,*angvel_corre_yx;
  double avg_angvel_corre_x,avg_angvel_corre_y,avg_angvel_corre_z,avg_angvel_corre_xy,avg_angvel_corre_yx; 
  double tau_vel_corre;
  int total_part_count_angvel_corre;

         angvel_corre_x=new double[Nmax_part_angvel_corre];
          angvel_corre_y=new double[Nmax_part_angvel_corre];
	   angvel_corre_z=new double[Nmax_part_angvel_corre];
	    angvel_corre_xy=new double[Nmax_part_angvel_corre];
	     angvel_corre_yx=new double[Nmax_part_angvel_corre];

	     tau_vel_corre=run_time- vel_corre_initial_time;


	     total_part_count_angvel_corre=0;

      for(i=0;i<Nof_part_angvel_corre;i++)
	       {
		 int ii=part_angvel_corre[i];
        if(p_y[ii]>y_min_for_corre && p_y[ii]<y_max_for_corre && ii>=0)
	   {
	   
	   angvel_corre_x[i] =pvx_dash_initial[i]*pvx_dash[ii];
	   angvel_corre_y[i] =pvy_dash_initial[i]*pvy_dash[ii];
	   angvel_corre_z[i] =pvz_dash_initial[i]*pvz_dash[ii];
	   angvel_corre_xy[i]=pvx_dash_initial[i]*pvy_dash[ii];
	   angvel_corre_yx[i]=pvy_dash_initial[i]*pvx_dash[ii];
	   
	   total_part_count_angvel_corre=total_part_count_angvel_corre+1;

	   //cout<<"----acc corre---   "<<ii<<"  "<<Nof_part_corre<<"    "<<acc_corre_x[i]<<endl;
	   }
	else
	  {
	    part_angvel_corre[i]=-1;  // removing this particle from correlation list
	    angvel_corre_x[i]=0.0;   // these are put to zero, to be added to the sum, with no contribution.
	    angvel_corre_y[i]=0.0;
	    angvel_corre_z[i]=0.0;
	    angvel_corre_xy[i]=0.0;
	    angvel_corre_yx[i]=0.0; 

	  }
	       }


       avg_angvel_corre_x=0.0;
       avg_angvel_corre_y=0.0;
       avg_angvel_corre_z=0.0;
       avg_angvel_corre_xy=0.0;
       avg_angvel_corre_yx=0.0;

      
       for(i=0;i<Nof_part_angvel_corre;i++)
	 { 
	         
       avg_angvel_corre_x= avg_angvel_corre_x + angvel_corre_x[i];
       avg_angvel_corre_y= avg_angvel_corre_y + angvel_corre_y[i];
       avg_angvel_corre_z= avg_angvel_corre_z + angvel_corre_z[i];
       avg_angvel_corre_xy=avg_angvel_corre_xy + angvel_corre_xy[i];
       avg_angvel_corre_yx=avg_angvel_corre_yx + angvel_corre_yx[i];

	 }

       if(total_part_count_angvel_corre>0)
	 {
     avg_angvel_corre_x= avg_angvel_corre_x/total_part_count_angvel_corre;
     avg_angvel_corre_y= avg_angvel_corre_y/total_part_count_angvel_corre;
     avg_angvel_corre_z= avg_angvel_corre_z/total_part_count_angvel_corre;
     avg_angvel_corre_xy= avg_angvel_corre_xy/total_part_count_angvel_corre;
     avg_angvel_corre_yx= avg_angvel_corre_yx/total_part_count_angvel_corre;
	 }
 //--------------- storing the Normalization factor----------

     if(N_count_angvel_corre==0)
	 {
	  avg_angvel_norm_x=avg_angvel_corre_x;
	  avg_angvel_norm_y=avg_angvel_corre_y;
	  avg_angvel_norm_z=avg_angvel_corre_z;
	  avg_angvel_norm_xy=avg_angvel_corre_xy;
	  avg_angvel_norm_yx=avg_angvel_corre_yx;
	 }

     
//------------Doing the Normalization of the correlation coefficient-----

     avg_angvel_corre_x=avg_angvel_corre_x/avg_angvel_norm_x;
     avg_angvel_corre_y=avg_angvel_corre_y/avg_angvel_norm_y;
     avg_angvel_corre_z=avg_angvel_corre_z/avg_angvel_norm_z;
     avg_angvel_corre_xy=avg_angvel_corre_xy/avg_angvel_norm_xy;
     avg_angvel_corre_yx=avg_angvel_corre_yx/avg_angvel_norm_yx;

     if(total_part_count_angvel_corre<min_Nop_for_corre)
       {
	 angular_index_return_from_corre=1;
          return;
       }

fp1019=fopen("/localscratch/cheantya/multipole/density_2000/volfrac_0.001/angvel_corre_avg.txt","a");

 if(N_count_angvel_corre==0) fprintf(fp1019,"%s  avg_angvel_norm_x , y, z, xy, yx are = %20.14lf %20.14lf %20.14lf %20.14lf %20.14lf\n","%",
				  avg_angvel_norm_x,avg_angvel_norm_y,avg_angvel_norm_z,avg_angvel_norm_xy,avg_angvel_norm_yx);

   fprintf(fp1019 ,"%20.14lf %20.14lf %20.14lf %20.14lf %20.14lf %20.14lf \n",tau_vel_corre, avg_angvel_corre_x, 
avg_angvel_corre_y, avg_angvel_corre_z,  avg_angvel_corre_xy, avg_angvel_corre_yx);

   fclose(fp1019);


   delete [] angvel_corre_x;
   delete [] angvel_corre_y;
   delete [] angvel_corre_z;
   delete [] angvel_corre_xy;
   delete [] angvel_corre_yx;
   return;

}

/*---------functions to get interpolated air and particle average velocities at
- particle position to get particle average acceleration which is used to -
-- calculate the particle acceleration fluctuation-----*/



void spline(double *x,double *y, int n, double *y2)
{
  int i,k;
double p,qn,sig,un,*u;

 u=new double[n-1];
 y2[0]=0;

 for(i=1;i<=n-2;i++)
   {
     sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
     p=sig*y2[i-1]+2.0;
     y2[i]=(sig-1.0)/p;
     u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
     u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
   }
 qn=un=0.0;
 y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
for(k=n-2;k>=0;k--)
  y2[k]=y2[k]*y2[k+1]+u[k];

delete []u;
}



int locate(const double xx[], int n, double x)
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


double splint(double *xa, double *ya, double *y2a, int ny_part, double x)

{
  int klo,khi;
  double h,b,a;
  double y;
 

  klo=ny_part;
  khi=ny_part+1;

   h=xa[khi]-xa[klo];
  if(h==0)cout<< "error in routine splint"<<endl;
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;

  y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  //cout<<"-----   "<<y<<endl;
  return y;
}





//------------------END OF THE PROGRAMME--------------------------------------








#include "channelflow/flowfield.h"
#include "channelflow/def"
using namespace channelflow;

extern double eta, sigma,sigma_sq,t_equb,time_max,del_ts,t_sample,radius;
extern double et,en,ep;
extern double b_lx,b_ly,b_lz;
extern double del_v,p_v_max,del_vx, p_vx_max,del_vy,p_vy_max,del_vz,p_vz_max;
extern double del_accx,del_accy,del_accz,del_accxy,accx_max,accy_max,accz_max,accxy_max;

extern double del_angvel,p_angvel_max,del_angvelx, p_angvelx_max,del_angvely,p_angvely_max,del_angvelz,p_angvelz_max;
extern double del_angaccx,del_angaccy,del_angaccz,del_angaccxy,angaccx_max,angaccy_max,angaccz_max,angaccxy_max;              //angular

extern double accx_max_from_air,accy_max_from_air,accz_max_from_air,accxy_max_from_air;
extern  double del_accx_from_air,del_accy_from_air,del_accz_from_air,del_accxy_from_air;

extern double angaccx_max_from_air,angaccy_max_from_air,angaccz_max_from_air,angaccxy_max_from_air;
extern  double del_angaccx_from_air,del_angaccy_from_air,del_angaccz_from_air,del_angaccxy_from_air;                          //angular

extern double max_air_velx,max_air_vely,max_air_velz;
extern double del_air_velx_dist,del_air_vely_dist,del_air_velz_dist;

extern double max_air_angvelx,max_air_angvely,max_air_angvelz;
extern double del_air_angvelx_dist,del_air_angvely_dist,del_air_angvelz_dist;                                                 //angular

double  y_min_for_distribution,y_max_for_distribution;


extern double vx_max_air,vy_air,vz_air;
extern double angvelx_max_air,angvely_air,angvelz_air;                                                                       //angular
extern double den_gas,vis_gas,den_part;
extern double g_x,g_y,g_z;
extern double tau_vp,del_tb,del_tb_actual;
extern double tau_x,tau_y,tau_z;
extern double length_scale,time_scale,acc_scale,vel_scale;
extern double diff_x[diff_nmax],diff_y[diff_nmax],diff_z[diff_nmax],diff_xy[diff_nmax],y_ch[diff_nmax];
extern double const_diff_x,const_diff_y,const_diff_z,const_diff_xy;
extern double del_x,del_y,del_z,del_r;

extern int  N_diff;
extern int max_coll,k_write,time_step_write,eqb_coll;
extern int N_rdels,Ndel_v,Ndel_vx,Ndel_vy,Ndel_vz;
extern int Ndel_accx, Ndel_accy,Ndel_accz,Ndel_accxy;
extern int Ndel_accx_from_air,Ndel_accy_from_air, Ndel_accz_from_air, Ndel_accxy_from_air;

extern int Ndel_air_velx,Ndel_air_vely,Ndel_air_velz;

extern int Ndel_angvel,Ndel_angvelx,Ndel_angvely,Ndel_angvelz;
extern int Ndel_angaccx, Ndel_angaccy,Ndel_angaccz,Ndel_angaccxy;
extern int Ndel_angaccx_from_air,Ndel_angaccy_from_air, Ndel_angaccz_from_air, Ndel_angaccxy_from_air;                                        //angular

extern int Ndel_air_angvelx,Ndel_air_angvely,Ndel_air_angvelz;

extern int N_diff;
extern int Ndel_x,Ndel_y,Ndel_z;

extern double p_x[n_atom+2],p_y[n_atom+2],p_z[n_atom+2];
//extern double dist;
extern int n_unit;

extern double  p_vx[n_atom+2],p_vy[n_atom+2],p_vz[n_atom+2],p_v[n_atom+2];
extern double  p_angvelx[n_atom+2],p_angvely[n_atom+2],p_angvelz[n_atom+2],p_angvel[n_atom+2];
extern double part_mass, moment_inertia,tot_KE;
extern double part_moment_inertia, tot_angularKE;

extern void   read_parameters();
extern void   initial_position();
extern void   initial_velocity(const FlowField& u,const Vector& x_grid, const Vector&  z_grid);
extern void   initial_position_velocity_angularvelocity_restart(const string& );
extern void   initial_angular_velocity(const FlowField& omega,const Vector& x_grid, const Vector&  z_grid);

extern double  Reynolds;
extern int n_order; // --5 th order lagrangian interpolation
extern void interpolation(const FlowField& , const Vector& ,const Vector& ,double ,double , double , int, double[]);



extern int    data_saving_index;
extern int    restart_for_part_index;

extern double del_data_saving_time,starting_time,data_saving_time;


extern int Nmax_part_acc_corre;
extern int N_count_acc_corre;
extern int Nmax_part_vel_corre;
extern int N_count_vel_corre;

extern int Nmax_part_angacc_corre;
extern int N_count_angacc_corre;
extern int Nmax_part_angvel_corre;
extern int N_count_angvel_corre;

extern int index_average_from_input, Max_dns_grid;
extern int min_Nop_for_corre;
extern int index_return_from_corre;
extern int index_return_from_part_acc_corre;
extern double y_min_for_corre,y_max_for_corre;

extern int angular_index_return_from_corre;
extern int index_return_from_part_angacc_corre;                    //angular

extern int index_stress_varience;   // -- put 1 if want to calculate the varience of stress else 0
extern int index_fourth_moment ;
extern int index_acc_dist_from_air;
extern int index_air_vel_fluc_dist;
extern int  n_del_y_acc_corre,n_del_y_angacc_corre;

extern int n_del_y_air_vel_corre;
       int index_air_vel_corre;
extern int index_return_from_air_vel_corre;
extern int N_count_fluid_vel_corre;
extern int Nmax_part_air_vel_corre;

extern int n_del_y_air_angvel_corre;
extern int index_return_from_air_angvel_corre;
extern int N_count_fluid_angvel_corre;                              //angular
extern int Nmax_part_air_angvel_corre;

 double  y_min_for_distribution_1,y_max_for_distribution_1,
  y_min_for_distribution_2,y_max_for_distribution_2,
  y_min_for_distribution_3,y_max_for_distribution_3;
 extern double pi;

extern int max_position_index;

void  read_parameters()
{
  //FILE *fp1;
  // double y_dns[100],part_dia;
 double  density,Sigma_scale;
 double pi=3.14159265;
  
   
 // eta =0.0000944;
// eta =0.0007;
//eta = 0.0004;
 //eta = 0.0004;
// eta =0.007;
   eta = 0.001;
  // eta = 0.04;
  // eta = 0.07;
  //eta =0.00001888;
  /* for(i=0;i<331;++i)
      sigma=0.00006;
   for(i=331;i<942;++i)
      sigma=0.00008;
   for(i=942;i<1909;++i)
      sigma=0.0001;
   for(i=1909;i<3221;++i)
      sigma=0.00012;
   for(i=3221;i<4779;++i)
      sigma=0.00015;
   for(i=4779;i<6091;++i)
      sigma=0.00018;
   for(i=6091;i<7057;++i)
      sigma=0.0002;
   for(i=7057;i<7668;++i)
      sigma=0.00022;
   for(i=7668;i<n_atom;++i)
      sigma=0.00024;
     
     sigma_mean=0.00015;*/
  sigma=0.000039;       // insert it in meter unit



  Sigma_scale=1.e0;

 
  /*-------------Air related---------------*/
  
    den_gas=1.1790*pow(Sigma_scale,-1);      // if sigma_scale=1, then unit is MKS
    vis_gas=1.75e-5*Sigma_scale;            //
     
   den_part=2000.0*pow(Sigma_scale,-1);//  glass particle
//den_part=8000.0*pow(Sigma_scale,-1);//  steel
    //den_part=400.0*pow(Sigma_scale,-1);// lycopodium


  density=6.0*eta/pi;
  max_coll=2000;
  k_write=20;
  eqb_coll=1000;
  time_step_write=10;

  /*------------ different time for simulation----------------*/
   t_equb=6000.0; 
   time_max=6001.1;



  /*-----------intermediate data saving for air and particle-----*/

   data_saving_index=1;    /*--- it is 1 if we want to save----*/
   restart_for_part_index=1;

   starting_time=6000.0;//2100.0;
   
   del_data_saving_time=50.0  ;
   data_saving_time=starting_time+del_data_saving_time ;


  
   /*---------- Sampling time related------------------------*/

  del_ts=0.01;            //interval of sampling ///initially it was .02
  t_sample=t_equb;       //  starting of the sampling time


  
  /*------------particle elastricity related--------------*/

  et=1.0;
  en=1.0;
  ep=1.0;
  
 
  double a=8.0*pi,b=4.0*pi/3.0;/*------ratio of the length in x,y, and z direction--*/

  b_ly=pow(((4.0*pi*sigma/*sigma_mean*/*sigma/*sigma_mean*/*sigma/*sigma_mean*/*n_atom)/(6.0*a*b*eta)),(1.0/3.0));
 // b_ly=0.002196;
  b_lx=a*b_ly/2.0;
  b_lz=b*b_ly/2.0;

   
  cout<<"  particle dia = "<<sigma <<"      particle  density =  "<<den_part<<endl;
  cout<<"b_ly="<<b_ly<<" "<<"b_ly/sigma= "<<b_ly/sigma/*sigma_mean*/<<endl;
 
  vx_max_air=(Reynolds*vis_gas)/(0.5*b_ly*den_gas);  


   length_scale=0.5*b_ly;                             //length scale
   acc_scale= (vx_max_air* vx_max_air)/length_scale;  //scaling of acceleration
   time_scale=(length_scale/vx_max_air);              //time scale  
   vel_scale=vx_max_air;                              //velocity scale


   cout<<"n_atom= "<<n_atom<<endl;
   cout<<" volume fraction   " <<eta <<endl;

   cout<<"scalse- length- time- velo- acc "<<length_scale<<" "<<time_scale<<" "<< vel_scale<<" "<<acc_scale<<endl;
   cout<<"shear_rate=  "<<vx_max_air/length_scale<<endl;
   cout<<endl;
   cout<<" t_sampling =  "<<t_equb <<" to  "<<time_max<<endl;

 tau_vp=(den_part*sigma*sigma)/(18.0*vis_gas); //tau_vp is dimentional,
   cout<<"tau_v (S)" << tau_vp<<endl;


   b_ly=b_ly/length_scale;
   b_lx=b_lx/length_scale;
   b_lz=b_lz/length_scale;

  sigma=sigma/length_scale;

  
  sigma_sq=sigma*sigma;
  radius = sigma/2;   
      
  /* del_r=0.025D0   //!Sampling increment and delta r for g(r)
     del_r=del_r*sigma  //!Scale Rdel to units of box length */


  N_rdels=500;    
  del_r=0.50*b_lx/500.0;
  cout<<"del_r= "<<del_r<<endl;
  //N_rdels = (int)(0.50*b_lx/del_r -1.0);

  /*-------Parameters for sampling the velocity distribution----*/

  del_v=0.008;  
  p_v_max=2.0;
  Ndel_v=(int)(2.0*p_v_max/del_v+1.0);
      
  del_vx=0.008;
  p_vx_max=2.00;
  Ndel_vx=(int)(2.0*p_vx_max/del_vx+1.0);

  del_vy=.002;
  p_vy_max=0.6;
  Ndel_vy=(int)(2.0*p_vy_max/del_vy+1.0);

  del_vz=0.002;
  p_vz_max=0.6;
  Ndel_vz=(int)(2.0*p_vz_max/del_vz+1.0);


/*--------------parameters for sampling air velocity(fluctuation) distribution function----------*/

  del_air_velx_dist=.004  ;
  max_air_velx=1.0;
  Ndel_air_velx=(int)(2.0* max_air_velx/del_air_velx_dist+1.0);
  del_air_velx_dist=(2.0*max_air_velx/(Ndel_air_velx-1));

  del_air_vely_dist=.0025  ;
  max_air_vely=.5;
  Ndel_air_vely=(int)(2.0* max_air_vely/del_air_vely_dist+1.0);
  del_air_vely_dist=(2.0*max_air_vely/(Ndel_air_vely-1));

  del_air_velz_dist=.0025  ;
  max_air_velz=.5;
  Ndel_air_velz=(int)(2.0* max_air_velz/del_air_velz_dist+1.0);
  del_air_velz_dist=(2.0*max_air_velz/(Ndel_air_velz-1));

/*-------Parameters for sampling the angular velocity distribution----*/

  del_angvel=0.008;  
  p_angvel_max=2.0;
  Ndel_angvel=(int)(2.0*p_angvel_max/del_angvel+1.0);
      
  del_angvelx=0.008;
  p_angvelx_max=2.00;
  Ndel_angvelx=(int)(2.0*p_angvelx_max/del_angvelx+1.0);

  del_angvely=.002;
  p_angvely_max=0.6;
  Ndel_angvely=(int)(2.0*p_angvely_max/del_angvely+1.0);

  del_angvelz=0.002;
  p_angvelz_max=0.6;
  Ndel_angvelz=(int)(2.0*p_angvelz_max/del_angvelz+1.0);


/*--------------parameters for sampling air velocity(fluctuation) distribution function----------*/

  del_air_angvelx_dist=.004  ;
  max_air_angvelx=1.0;
  Ndel_air_angvelx=(int)(2.0* max_air_angvelx/del_air_angvelx_dist+1.0);
  del_air_angvelx_dist=(2.0*max_air_angvelx/(Ndel_air_angvelx-1));

  del_air_angvely_dist=.0025  ;
  max_air_angvely=.5;
  Ndel_air_angvely=(int)(2.0* max_air_angvely/del_air_angvely_dist+1.0);
  del_air_angvely_dist=(2.0*max_air_angvely/(Ndel_air_angvely-1));

  del_air_angvelz_dist=.0025  ;
  max_air_angvelz=.5;
  Ndel_air_angvelz=(int)(2.0* max_air_angvelz/del_air_angvelz_dist+1.0);
  del_air_angvelz_dist=(2.0*max_air_angvelz/(Ndel_air_angvelz-1));


  /*-------Parameters for sampling  accleration distribution function ---*/


 del_accx=4.0e-5;
 accx_max=0.03;
 Ndel_accx=(int)(2.0*accx_max/del_accx+1.0);
 del_accx=(2.0*accx_max/(Ndel_accx-1));   // --recalculation of del_accx


 del_accy=2.e-5;
 accy_max=0.015;
 Ndel_accy=(int)(2.0*accy_max/del_accy+1.0);
 del_accy=(2.0*accy_max/(Ndel_accy-1));


 del_accz=2.5e-5;
 accz_max=0.02;
 Ndel_accz=(int)(2.0*accz_max/del_accz+1.0);;
 del_accz=(2.0*accz_max/(Ndel_accz-1));

 del_accxy= 4.0e-7;
 accxy_max=.0003;
 Ndel_accxy= (int)(2.0*accxy_max/del_accxy+1.0);;
 del_accxy=(2.0*accxy_max/(Ndel_accxy-1));
            
           
 /*-------Parameters for sampling  accleration distribution function from air vel fluctuation---*/
 del_accx_from_air=4.0e-5;
 accx_max_from_air=0.025;
 Ndel_accx_from_air=(int)(2.0*accx_max_from_air/del_accx_from_air+1.0);
 del_accx_from_air=(2.0*accx_max_from_air/(Ndel_accx_from_air-1));   // --recalculation of del_accx


 del_accy_from_air=2.e-5;
 accy_max_from_air=0.015;
 Ndel_accy_from_air=(int)(2.0*accy_max_from_air/del_accy_from_air+1.0);
 del_accy_from_air=(2.0*accy_max_from_air/(Ndel_accy_from_air-1));   // --recalculation of del_accx

 del_accz_from_air=2.5e-5;
 accz_max_from_air=0.02;
 Ndel_accz_from_air=(int)(2.0*accz_max_from_air/del_accz_from_air+1.0);
 del_accz_from_air=(2.0*accz_max_from_air/(Ndel_accz_from_air-1));   // --recalculation of del_accx

 del_accxy_from_air=2*4.0e-7;
 accxy_max_from_air=0.0003;
 Ndel_accxy_from_air=(int)(2.0*accxy_max_from_air/del_accxy_from_air+1.0);
 del_accxy_from_air=(2.0*accxy_max_from_air/(Ndel_accxy_from_air-1));   // --recalculation of del_accx

 cout<<"+++++++++++  "<<del_accxy_from_air<<endl;

 /*-------Parameters for sampling  angular accleration distribution function ---*/


 del_angaccx=4.0e-5;
 angaccx_max=0.03;
 Ndel_angaccx=(int)(2.0*angaccx_max/del_angaccx+1.0);
 del_angaccx=(2.0*angaccx_max/(Ndel_angaccx-1));   // --recalculation of del_accx


 del_angaccy=2.e-5;
 angaccy_max=0.015;
 Ndel_angaccy=(int)(2.0*angaccy_max/del_angaccy+1.0);
 del_angaccy=(2.0*angaccy_max/(Ndel_angaccy-1));


 del_angaccz=2.5e-5;
 angaccz_max=0.02;
 Ndel_angaccz=(int)(2.0*angaccz_max/del_angaccz+1.0);;
 del_angaccz=(2.0*angaccz_max/(Ndel_angaccz-1));

 del_angaccxy= 4.0e-7;
 angaccxy_max=.0003;
 Ndel_angaccxy= (int)(2.0*angaccxy_max/del_angaccxy+1.0);;
 del_angaccxy=(2.0*angaccxy_max/(Ndel_angaccxy-1));
            
           
 /*-------Parameters for sampling  angular accleration distribution function from air vel fluctuation---*/
 del_angaccx_from_air=4.0e-5;
 angaccx_max_from_air=0.025;
 Ndel_angaccx_from_air=(int)(2.0*angaccx_max_from_air/del_angaccx_from_air+1.0);
 del_angaccx_from_air=(2.0*angaccx_max_from_air/(Ndel_angaccx_from_air-1));   // --recalculation of del_accx


 del_angaccy_from_air=2.e-5;
 angaccy_max_from_air=0.015;
 Ndel_angaccy_from_air=(int)(2.0*angaccy_max_from_air/del_angaccy_from_air+1.0);
 del_angaccy_from_air=(2.0*angaccy_max_from_air/(Ndel_angaccy_from_air-1));   // --recalculation of del_accx

 del_angaccz_from_air=2.5e-5;
 angaccz_max_from_air=0.02;
 Ndel_angaccz_from_air=(int)(2.0*angaccz_max_from_air/del_angaccz_from_air+1.0);
 del_angaccz_from_air=(2.0*angaccz_max_from_air/(Ndel_angaccz_from_air-1));   // --recalculation of del_accx

 del_angaccxy_from_air=2*4.0e-7;
 angaccxy_max_from_air=0.0003;
 Ndel_angaccxy_from_air=(int)(2.0*angaccxy_max_from_air/del_angaccxy_from_air+1.0);
 del_angaccxy_from_air=(2.0*angaccxy_max_from_air/(Ndel_angaccxy_from_air-1));   // --recalculation of del_accx

 cout<<"+++++++++++  "<<del_angaccxy_from_air<<endl;


  max_position_index=3;
 
 y_min_for_distribution=0.8;
 y_max_for_distribution=1.2;


 y_min_for_distribution_1=.9;
 y_max_for_distribution_1=1.1;

 y_min_for_distribution_2=1.6;
 y_max_for_distribution_2=1.7;

  y_min_for_distribution_3=1.85;
  y_max_for_distribution_3=1.95;





 /*----------parameter for calculating particle acceleration correlation---*/

  Nmax_part_acc_corre=n_atom;
  Nmax_part_vel_corre=1500;
  Nmax_part_air_vel_corre=n_atom;

  Nmax_part_angacc_corre=n_atom;
  Nmax_part_angvel_corre=1500;
  Nmax_part_air_angvel_corre=n_atom;

  min_Nop_for_corre=100; // minimum no of particle to do averaging for corre.

  index_return_from_corre=0;
  index_return_from_air_vel_corre=0;
  index_return_from_part_acc_corre=0;

  angular_index_return_from_corre=0;
  index_return_from_air_angvel_corre=0;
  index_return_from_part_angacc_corre=0;

  N_count_acc_corre=0;     // just to initialize
  N_count_vel_corre=0;
  N_count_fluid_vel_corre=0;

  N_count_angacc_corre=0;     // just to initialize
  N_count_angvel_corre=0;
  N_count_fluid_angvel_corre=0;


  index_average_from_input=0; // if you want to load average air and particle velocity for part_acc_corre and
                              // particle velocity correlation,  
                              //-- put it 1 otherwise it is 0, then it will do only particle
                              // -- average.

  index_stress_varience=0;   // -- put 1 if want to calculate the varience of stress else 0
  index_fourth_moment=0 ;    //--- put 1 if you want to calculate fourth moment by 
                             //--loading avg particle velocity else 0
  index_acc_dist_from_air=0;
  index_air_vel_fluc_dist=0;
  index_air_vel_corre=0;


  Max_dns_grid=65; // this is equal to Ny used in main()

  y_min_for_corre=0.8;   /* these are omly for part vel corre.*/
   y_max_for_corre=1.2;  

  // y_min_for_corre=1.6;
  // y_max_for_corre=1.7; 

  //y_min_for_corre=1.85;
  // y_max_for_corre=1.95; 

  n_del_y_acc_corre=10;
  n_del_y_air_vel_corre=10;

  n_del_y_angacc_corre=10;
  n_del_y_air_angvel_corre=10;

  /*-----------Gravity---------*/

    double g=9.8;  /*------ in m/s^2 unit--- */
    g_x=g/acc_scale;
    g_y=0.0;
    g_z=0.0;
      
  /* ---------time related parameter--------*/
  
 
  
  del_tb=0.025*time_scale;   // dns integration interval,see the rescaling below  //earlier 0.025
  del_tb_actual=0.03*time_scale; // not required
  tau_x=8.e-3;
  tau_y=8.e-3;
  tau_z=8.e-3 ; 

  /*----------times in non dimensional form---rescaling----*/
  tau_vp=tau_vp/time_scale;
  del_tb=del_tb/time_scale;
  tau_x= tau_x/time_scale;
  tau_y= tau_y/time_scale;
  tau_z= tau_z/time_scale;


//-------Brownian Noise related -----------** not required


 const_diff_x=0.0025/(tau_vp*tau_vp);
 const_diff_y=3.3784e-4/(tau_vp*tau_vp);
 const_diff_z=8.3333e-4/(tau_vp*tau_vp);
 const_diff_xy=-7.4325e-04/(tau_vp*tau_vp);
       

 /*----------Grid related---------------*/


  
  Ndel_x=100;
  Ndel_z=100;
  Ndel_y=100;
  del_x=(b_lx)/(double)(Ndel_x);
  del_y=(b_ly-sigma)/(double)(Ndel_y); //readjustment.
  del_z=(b_lz)/(double)(Ndel_z);
  return;
}

/*--------------------------------------------------
 declearing the initial position of the particle
---------------------------------------------------*/
   void initial_position()
     {
       int i,j;  //ij,m, k,kct;
       //double t;
       double p_x_ij,p_y_ij,p_z_ij,r_ij;

       /* t=pow(((double)n_atom/4.0),((double)1.0/3.0));
       if ((t-(int)t)>0.9999999999)t=(int)t+1.0;
      
       printf("%30.20f\n",t);


       n_unit=(int)t;
       cout<< n_unit<<" "<< t<<" "<<(int)t << endl;
       dist=b_lx/(double)n_unit;

cout<<"nunit"<<" "<< n_unit<<" "<<(pow(((double)n_atom/4.0),(1./3.)))<<endl;
// cin.get();

       p_x[0]=0.0;
       p_y[0]=0.0+1.50*sigma/2.0;
       p_z[0]=0.0 ;    //!Assign position to the four spheres of the unit cell
        
	  p_x[1]=0.0;
	  p_y[1]=0.50*dist+1.50*sigma/2.0;
	  p_z[1]=0.50*dist;
        
	  p_x[2]=0.50*dist; 
	  p_y[2]=0.0+1.50*sigma/2.0;
	  p_z[2]=0.50*dist;

	  p_x[3]=0.50*dist ;
	  p_y[3]=0.50*dist+1.5*sigma/2.0;
	  p_z[3]=0.0;
	  
	  m=0;
	  kct=0;
         
	for(i=0;i<n_unit;i++)       //! Replicate unit cell n_unit times in each direction
	   for( j=0;j<n_unit;j++)
	     for(k=0;k<n_unit;k++){
	       for(ij=0;ij<4;ij++){
		 
            if(kct >= n_atom) 
	      return;              //! If number of positions assigned-
           
	    p_x[ij+m]=p_x[ij]+dist*(k);  // !- equals number of spheres, get out
	    p_y[ij+m]=p_y[ij]+dist*(j);
            p_z[ij+m]=p_z[ij]+dist*(i);
            kct=kct+1 ;                         
            
	    //cout<<"going on"<<ij+m<<ENDL;
           
            if((p_y[ij+m]+sqrt(sigma_sq)/2.0)>b_ly)
	      {
	      cout<<"PARTICLE IS GOING OUT OF THE WALL"<<p_y[ij+m]<<endl;
	    exit(1);//------------------------------check
	      }
	       }
	       m=m+4;
	        }

       */

       /*---------------- To give the particle at random position---------------*/

	p_y[n_atom]=0.0;           // position of the left wall
	p_y[n_atom+1]=b_ly;         // position of the right wall
	
	i=0;
	while(i<n_atom)
	  {
	  repeat:
	    p_x[i]=(b_lx-2.0*sigma)*(double)rand()/RAND_MAX + sigma;
	    p_y[i]=(b_ly-2.0*sigma)*(double)rand()/RAND_MAX + sigma;
	    p_z[i]=(b_lz-2.0*sigma)*(double)rand()/RAND_MAX + sigma;
	    
	    if(i>1)
	    for(j=0;j<i;j++)
	      {
		p_x_ij=p_x[i]-p_x[j];
		p_y_ij=p_y[i]-p_y[j];
		p_z_ij=p_z[i]-p_z[j];
		r_ij=sqrt(p_x_ij*p_x_ij+p_y_ij*p_y_ij+p_z_ij*p_z_ij);
		
		if((r_ij-sigma)<0.0)goto repeat;
	      }

	    if((p_y[i]+sqrt(sigma_sq)/2.0)>b_ly)
	      {
	      cout<<"PARTICLE IS GOING OUT OF THE WALL"<<i<<endl;
	    exit(1);//------------------------------check
	      }

	    i=i+1;
		  }

		


	   return;
     }

/*------------------------------------------------------
//                    initial velocities
//-----------------------------------------------------*/


void initial_velocity(const FlowField& u,const Vector& x_grid, const Vector&  z_grid)


      {
       int i;                           
       //double  sum_pvx=0.0,sum_pvy=0.0,sum_pvz=0.0;
     //float xx,yy,zz;
     //double fatom;
     double interpolated_u[3];
     int n_order=5;

     //srand(SEED);

           double a_=u.a();
	   double b_=u.b();
	   double ubase_air;

 //-----loop for interpolated initial velocity-------

     for(i=0;i<n_atom;i++)
       {
interpolation(u, x_grid,z_grid,p_x[i],p_y[i],p_z[i], n_order, interpolated_u);

ubase_air = (1.0 - square(abs(p_y[i]-(b_+a_)/2.0)/((b_-a_)/2.0)));
 //ubase_air = p_y[i]-(a_+b_)/2.0;
p_vx[i]=ubase_air+interpolated_u[0];
p_vy[i]=interpolated_u[1];
p_vz[i]=interpolated_u[2];


       }

 for(i=n_atom;i<n_atom+2;++i) //here i=n_atom means the left wall
       {

       p_vx[i]=0;
       p_vy[i]=0;
       p_vz[i]=0;
     }

tot_KE=0.0;

 for(i=0;i<n_atom;i++)
       {
        
	 tot_KE=tot_KE + p_vx[i]* p_vx[i]+p_vy[i]* p_vy[i]+p_vz[i]* p_vz[i];
     }


 /* -----loop for random initial velocity-------*/

     /*    for(i=0;i<n_atom;i++)
       { //xx=5.0*((double)rand()/RAND_MAX);  // vx is given from 0-5

	 xx=(2.0*((double)rand()/RAND_MAX)-1.0);
         yy=(2.0*(double)rand()/RAND_MAX-1.0);  //vy is -1 to 1

         zz=(2.0*(double)rand()/RAND_MAX-1.0);

     //cout<<"xx="<<xx<<" "<<yy<<" "<<zz<<endl;
     p_vx[i]=xx;
     p_vy[i]=yy;
     p_vz[i]=zz;
     
     sum_pvx=sum_pvx+p_vx[i];
     sum_pvy=sum_pvy+p_vy[i];
     sum_pvz=sum_pvz+p_vz[i];
      }
     
     for(i=n_atom;i<n_atom+2;++i) //here i=n_atom means the left wall
       {
       p_vx[i]=0;
       p_vy[i]=0;
       p_vz[i]=0;
     }

     tot_KE=0.0;
     fatom=(double)(n_atom);
     for(i=0;i<n_atom;i++)
       {
         p_vx[i]= p_vx[i]-sum_pvx/fatom;  //Adjust linear velocities so that
	 p_vy[i]= p_vy[i]-sum_pvy/fatom;   //total linear momentum is zero.
	 p_vz[i]= p_vz[i]-sum_pvz/fatom;
           
	 tot_KE=tot_KE + p_vx[i]* p_vx[i]+p_vy[i]* p_vy[i]+p_vz[i]* p_vz[i];
     }
     */







     tot_KE=0.50*part_mass*tot_KE;

     cout<<"returning from initial vel"<<endl;
     return;
     }


  
/*------------------------------------------------------
//                    initial angular velocities
//-----------------------------------------------------*/


void initial_angular_velocity(const FlowField& omega, const Vector& x_grid, const Vector&  z_grid)


      {
       int i;                           
       double interpolated_omega[3];
     int n_order=5;

    
          // double a_=u.a();
	  // double b_=u.b();
	  // double ubase_air;

 //-----loop for interpolated initial velocity-------

       

     for(i=0;i<n_atom;i++)
       {
 interpolation(omega, x_grid,z_grid,p_x[i],p_y[i],p_z[i], n_order, interpolated_omega);

//ubase_air = (1.0 - square(abs(p_y[i]-(b_+a_)/2.0)/((b_-a_)/2.0)));
 //ubase_air = p_y[i]-(a_+b_)/2.0;
p_angvelx[i]=0.5*interpolated_omega[0];
p_angvely[i]=0.5*interpolated_omega[1];
p_angvelz[i]=0.5*interpolated_omega[2];

}

 for(i=n_atom;i<n_atom+2;++i) //here i=n_atom means the left wall
       {

       p_angvelx[i]=0;
       p_angvely[i]=0;
       p_angvelz[i]=0;
     }

tot_angularKE=0.0;

 for(i=0;i<n_atom;i++)
       {
        
	 tot_angularKE=tot_angularKE + p_angvelx[i]* p_angvelx[i]+p_angvely[i]* p_angvely[i]+p_angvelz[i]* p_angvelz[i];
     }



     tot_angularKE=0.50*moment_inertia*tot_angularKE;

     cout<<"returning from initial angular vel"<<endl;
     return;
     }

void initial_position_velocity_angularvelocity_restart(const string& filebase)
  {
    int ik,junk_index;
  string infile1(filebase);
  infile1 += ".asc";
  ifstream is1(infile1.c_str(),ios::in );

 if (!is1.good()) {
    cerr << "can't open file " 
	 << infile1 << endl;
    exit(1);
 }

 for(ik=0;ik<n_atom;ik++)
   {
     is1>>junk_index>>p_x[ik]>>p_y[ik]>>p_z[ik]>>p_vx[ik]>>p_vy[ik]>>p_vz[ik]>>p_angvelx[ik]>>p_angvely[ik]>>p_angvelz[ik];
  }
  is1.close();

   p_y[n_atom]=0.0;           // position of the left wall
   p_y[n_atom+1]=b_ly; // position of the right wall

 for(ik=n_atom;ik<n_atom+2;++ik) //here i=n_atom means the left wall
       {

       p_vx[ik]=0;
       p_vy[ik]=0;
       p_vz[ik]=0;

       p_angvelx[ik]=0;
       p_angvely[ik]=0;
       p_angvelz[ik]=0;
       
     }



   return;
  }



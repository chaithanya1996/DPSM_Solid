/* -------------------------------------------------

Date : 14 October 2018
Program Name : Solid_transducer_u.cpp

Description : Calculation the Displacement and the Stress produced in solids due to the Transducer placed on the source location

Written By : 
              __      __     _____ 
     /\       \ \    / /    / ____|
    /  \       \ \  / /    | (___  
   / /\ \       \ \/ /      \___ \ 
  / ____ \       \  /       ____) |
 /_/    \_\       \/       |_____/ 
                                   

-------------------------------------------------- */


#include<iostream>
#include "dpsm_interface.hpp"

using std::cout;
using std::endl;
using std::cin;
using std::pow;
using std::vector;

// Wrtiting Log Files

using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;
#include <cstdlib>

// This function is not working in the header file 
dpsm::solid_interface<double> combine_interface(dpsm::solid_interface<double> interface1,dpsm::solid_interface<double> interface2){
  mat joined_mats = join_cols(interface1.interface_source_1,interface2.interface_source_1);
  dpsm::solid_interface <double> joint_interface_solid(joined_mats,interface1.normal_to_solid_interface);
  return(joint_interface_solid);
};

template <typename TYPE>
inline TYPE calc_lamda(TYPE E, TYPE V ){
  return(E * V / (1+V) / (1-2*V));
}

template <typename TYPE>
inline TYPE calc_mu (TYPE E, TYPE V ){
  return(E / 2 /(1+V));
}

int main(){
  
   // Defining Parameters (work in progress)
  double mm = 0.001;
  double kHz = 1000;
  double km = 1000;
  // Cubic centimenter
  double cubic_centimeter = pow(10,-6);
  double gms = pow(10,-3);
  
  

  // Wave Speed in Water
  double c_w = 340;
  double c_al = 1.48 * km;
  // P-Wave Speed
  double  c_al_pwave = 6.5 * km;
  // S-Wave Speed
  double  c_al_swave = 3.13 * km;
  
  // Density of Water
  double rho_w = 1 * gms / cubic_centimeter;
  double rho_al = 2.7 ;

  // Aluminium Properties
  // double nu_al = 0.349;
  // double ral_al_crit_angle = 30.4196 * M_PI / 180 ;


  // double f = 20 * kHz; double  w = 2*M_PI*f; double k = w/c;
  // double a0 = 100.0 * mm, v0 = 1;
  // cx_double img_start (0,1.0);
  // double k_s = 0.01, k_p = 0.03, omega = w;



  /* --------------------------------
        Transducer Properties
   ----------------------------------*/

    // Transducer Properties
  double trans_freq = 1000 * kHz;
  double omega_trans = 2 * M_PI * trans_freq ;
  double Transducer_diameter = 1 * mm ;

  // Transducer Properties in Aluminium

  double k_s_aluminium = c_al_swave / omega_trans ;
  double k_p_aluminium = c_al_pwave / omega_trans ;

  // calculating the r_s for given Transducer and Medium 

  double r_s_tran = r_s_calculator<double>(trans_freq,c_al);
  int No_of_transducer_points = Transducer_diameter / r_s_tran;
  //int No_of_transducer_points = 3;
  
  rowvec normal_to_location = {1,0,0};
  rowvec start_line_source = {0,0,0};
  rowvec end_line_source = {0,0.002,0};
  // mat source_point_locations = line_generator(start_line_source,end_line_source,No_of_transducer_points); // transducer Location in first array
  mat source_point_locations,source_point_locations_normals;
  vec origin_of_circle = {0,0,0};
  // Cicular Source Generration

  std::tie(source_point_locations,source_point_locations_normals) = circle_maker(r_s_tran,Transducer_diameter, origin_of_circle , normal_to_location);

  /* Design Of Experiment 
     Phase - 1 : Apply Boundary Condiitons and Accuire Stregnth of sources determined
     Phase - 2 : Get Actual Strengths as an image
     Phase - 3 : Plot and Profit !
   */

  // Start Phase  1
    /* --------------------------------
      Alumnium  Material Properties
   ----------------------------------*/
  double GPa = pow(10,9);
  double E_al = 69 * GPa, V_al = 0.334;
  double E_st = 200 * GPa, V_st = 0.305;
  
  double lamda_al  = calc_lamda(E_al,V_al);
  double mu_al  = calc_mu(E_al,V_al);
  
  double lamda_st  = calc_lamda(E_st,V_st);
  double mu_st  = calc_mu(E_st,V_st);

  
  dpsm::solid_material<double> aluminium(lamda_al, mu_al);
  dpsm::solid_material<double> steel(lamda_st, mu_st);


  
  
  // place of Observation 

  rowvec x_space_ply_1 = regspace<rowvec>(0,0.001,0.05);
  rowvec x_space_ply_2 = regspace<rowvec>(0.05,0.001,0.1);
  rowvec y_space_ply_1 = regspace<rowvec>(0,0.001,0.2);
  rowvec y_space_ply_2 = regspace<rowvec>(0,0.001,0.2);

  mat observation_target = grid_target_generator(x_space_ply_1,y_space_ply_1); // Observation Point in Ply 1

  observation_target.save("OBSER.csv",csv_ascii);

  cout << observation_target.n_rows << " ------" << observation_target.n_cols << "---------" << observation_target.n_elem << endl;

  // Interface Vector 

  rowvec interface_start = {0.05,0,0};
  rowvec interface_end = {0.05,0.2,0};
  mat interface_line = line_generator(interface_start,interface_end,No_of_transducer_points);

  // cout << interface_line << endl;
  
  rowvec normal_to_interface = {1,0,0};
  rowvec normal_to_interface_reverse = {-1,0,0};

  rowvec interface_st_end_start = {0.1,0,0};
  rowvec interface_st_end_end = {0.1,0.2,0};
  mat interface_line_st = line_generator(interface_start,interface_end,No_of_transducer_points);

 

  //cout << interface_line_st << endl;
  
  // List of Interfaces Used in secondary Ply

  dpsm::solid_interface<double> interface_al(interface_line,normal_to_interface_reverse);
  dpsm::solid_interface<double> interface_al_st(interface_line,normal_to_interface);
  dpsm::solid_interface<double> interface_st_end(interface_line_st,normal_to_interface);

  // Creating the Transducer Type For Geometry preparation

  dpsm::transducer<double> ultrasonice_trasducer(source_point_locations,normal_to_location);
  
  // Preparing all Things for Dpsm Geometry thing
  vector<dpsm::solid_material<double>> material_array_exp(2,dpsm::solid_material<double>());
  vector<dpsm::solid_interface<double>> Interface_array_exp(2,dpsm::solid_interface<double>());
  vector<dpsm::transducer<double>> transducer_array_exp(2,dpsm::transducer<double>());


  cout << "Black Sheep 3" << endl; 
  // Assign The vector of Things
  
  material_array_exp[0] = aluminium;
  material_array_exp[1] = steel;

  // Assigning the interface
  
  Interface_array_exp[0] = interface_al;
  Interface_array_exp[1] = combine_interface(interface_al_st,interface_st_end);
  
  // Assginig the Transducer
  // We simply leavbe the other as blank in the array thus we nly define one transducer 
  transducer_array_exp[0] = ultrasonice_trasducer;
    
  //grid generation
  cout << "Black Sheep 4" << endl; 
  //mat obs_plane_location_ply_1 = grid_target_generator(x_space_ply_1,y_space_ply_1);  //   logically the program is complete

  dpsm::geometry<double> TOTAL_LAMINATE = dpsm::geometry<double>(2,material_array_exp,transducer_array_exp,Interface_array_exp,r_s_tran);
  cout << "Black Sheep 5" << endl; 
  std::vector<dpsm::solid_ply<double>> T_lam_ply_vec = TOTAL_LAMINATE.get_the_material_grid();
  // Completed the Generation of Plane

  cout << "Black Sheep 6" << endl;

  // mat Transducer_sterss_matrix(source_point_locations.n_rows,3,fill::zeros);
  // for (int i = 0; i <Transducer_sterss_matrix.n_rows ; ++i) {
  //   Transducer_sterss_matrix(i,0) = 1;
  // }

  
  cx_3d<double> Transducer_stress_Cx_3d(source_point_locations.n_rows,cx_mat(3,3,fill::zeros));
  for (int i= 0; i < source_point_locations.n_rows; ++i) {
    Transducer_stress_Cx_3d[i](0,0) = 1 * pow(10,5); 
  }

  // Path - /home/chaithanya/Documents/DPSM/Package/Results

  
  // Commencing the operation of Migties

   // get_strength_hetro(mat source_point_list,mat passive_source_list ,double r_s, rowvec direction_cosine,cx_3d stress_matrix,double k_s,double k_p,double rho , double omega,double mu,double lamda)
  // cout << "-------------------------" << endl;
  // cout << T_lam_ply_vec[0].active_sources << endl;
  // cout << T_lam_ply_vec[0].passve_sources << endl;
  // cout << r_s_tran << endl;
  // cout << normal_to_location << endl;
  // cout << Transducer_stress_Cx_3d.size()  << endl;
  // cout << T_lam_ply_vec[0].active_sources.n_rows << endl;
  cx_mat  Result = get_strength_hetro(T_lam_ply_vec[0].active_sources, T_lam_ply_vec[0].passve_sources,r_s_tran,normal_to_location,Transducer_stress_Cx_3d,k_s_aluminium,k_p_aluminium,rho_al,omega_trans,mu_al,lamda_al);

  cout << "Black sheep 11" << endl ;
  cx_3d<double> Stress_resultants_Result =stress_from_points(T_lam_ply_vec[0].active_sources, Result, observation_target,k_s_aluminium,k_p_aluminium,rho_al,omega_trans,lamda_al,mu_al);
  cout << "Black sheep 12" << endl ;


  Result.save("point_strength.dat",csv_ascii);


  std::string PATH_TO_SAVE_RESULTS = "/home/chaithanya/Documents/DPSM/Package/Results/";

  int Sucess = save_cx_3d(Stress_resultants_Result,PATH_TO_SAVE_RESULTS);

  

  
  return(0);
}

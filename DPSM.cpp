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

template<typename T>
dpsm::solid_interface<T> combine_interface(dpsm::solid_interface<T> interface1,dpsm::solid_interface<T> interface2){
  Mat<T> joined_mats = join_cols(interface1.interface_source_1,interface2.interface_source_1);
  dpsm::solid_interface <T> joint_interface_solid(joined_mats,interface1.normal_to_solid_interface);
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
 
  using P_DTYPE = float;  // Defining the Precision of the Computation
  
   // Defining Parameters (work in progress)
  P_DTYPE mm = 0.001;
  P_DTYPE kHz = 1000;
  P_DTYPE km = 1000;
  // Cubic centimenter
  P_DTYPE cubic_centimeter = pow(10,-6);
  P_DTYPE gms = pow(10,-3);
  
  

  // Wave Speed in Water
  P_DTYPE c_w = 340;
  P_DTYPE c_al = 1.48 * km;
  // P-Wave Speed
  P_DTYPE  c_al_pwave = 6.5 * km;
  // S-Wave Speed
  P_DTYPE  c_al_swave = 3.13 * km;
  
  // Density of Water
  P_DTYPE rho_w = 1 * gms / cubic_centimeter;
  P_DTYPE rho_al = 2.7 ;

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
  P_DTYPE trans_freq = 200 * kHz;
  P_DTYPE omega_trans = 2 * M_PI * trans_freq ;
  P_DTYPE Transducer_diameter = 20 * mm ;

  // Transducer Properties in Aluminium

  P_DTYPE k_s_aluminium = c_al_swave / omega_trans ;
  P_DTYPE k_p_aluminium = c_al_pwave / omega_trans ;

  // calculating the r_s for given Transducer and Medium 

  P_DTYPE r_s_tran = r_s_calculator<P_DTYPE>(trans_freq,c_al);
  //P_DTYPE r_s_tran = 0.0002;
  cout << " The Distance Between Transducer Source Points --" << r_s_tran << endl; 
  
  Row<P_DTYPE> normal_to_location = {1,0,0};
  Row<P_DTYPE> normal_to_location_reverse = {-1,0,0};

  
  // Row<P_DTYPE> source_origin = {0,-0.01,-0.01};
  // Row<P_DTYPE> source_y = {0,0.01,-0.01};
  // Row<P_DTYPE> source_x = {0,-0.01,0.01};

  // int  source_x_div = int(vec_mag<P_DTYPE>(source_x) / r_s_tran) ;
  // int  source_y_div = int(vec_mag<P_DTYPE>(source_y) / r_s_tran);

  // cout << " No of Divisions in Direction 1 := " << vec_mag<P_DTYPE>(source_x) << "--" << source_x_div << endl ;
  // cout << " No of Divisions in Direction 2 := " << vec_mag<P_DTYPE>(source_y) << "--" <<  source_y_div << endl ;
  
  //  // Cicular Source Generration
  // Mat<P_DTYPE> source_point_locations = rectangle_generator<P_DTYPE>(source_x,source_y,source_origin,source_x_div,source_y_div);

  Row<P_DTYPE> source_st = {0,-0.2,0};
  Row<P_DTYPE> source_end = {0,0.2,0};
  int source_divs = vec_mag<P_DTYPE>(source_st - source_end)/r_s_tran;
  cout << "No of Sources" <<  source_divs << endl;
  Mat<P_DTYPE> source_point_locations = line_generator<P_DTYPE>(source_st,source_end,source_divs-1);
  
  
  /* Design Of Experiment 
     Phase - 1 : Apply Boundary Condiitons and Accuire Stregnth of sources determined
     Phase - 2 : Get Actual Strengths as an image
     Phase - 3 : Plot and Profit !
   */

  // Start Phase  1
    /* --------------------------------
      Alumnium  Material Properties
   ----------------------------------*/
  
  P_DTYPE GPa = pow(10,9);
  P_DTYPE E_al = 69 * GPa, V_al = 0.334;
  P_DTYPE E_st = 200 * GPa, V_st = 0.305;
  
  P_DTYPE lamda_al  = calc_lamda<P_DTYPE>(E_al,V_al);
  P_DTYPE mu_al  = calc_mu<P_DTYPE>(E_al,V_al);
  
  P_DTYPE lamda_st  = calc_lamda(E_st,V_st);
  P_DTYPE mu_st  = calc_mu(E_st,V_st);

  
  dpsm::solid_material<P_DTYPE> aluminium(lamda_al, mu_al);
  dpsm::solid_material<P_DTYPE> steel(lamda_st, mu_st);
  
  cout << "Materials Are Defined" << endl;
  // place of Observation 

  Row<P_DTYPE> x_space_ply_1 = regspace<Row<P_DTYPE>>(0.0,0.001,0.05);
  Row<P_DTYPE> x_space_ply_2 = regspace<Row<P_DTYPE>>(0.05,0.001,0.1);
  Row<P_DTYPE> y_space_ply_1 = regspace<Row<P_DTYPE>>(0.0,0.001,0.2);
  Row<P_DTYPE> y_space_ply_2 = regspace<Row<P_DTYPE>>(0.0,0.001,0.2);

  //Mat<P_DTYPE> observation_target_1 = grid_target_generator<P_DTYPE>(x_space_ply_1,y_space_ply_1); // Observation Point in Ply 1
  //Mat<P_DTYPE> observation_target_2 = grid_target_generator<P_DTYPE>(x_space_ply_2,y_space_ply_2); // Observation Point in Ply 1

  //observation_target_1.save("OBSER_1.csv",csv_ascii);
  //observation_target_2.save("OBSER_2.csv",csv_ascii);
  
  //cout << observation_target_1.n_rows << " ------" << observation_target_1.n_cols << "---------" << observation_target_1.n_elem << endl;
  //cout << observation_target_2.n_rows << " ------" << observation_target_2.n_cols << "---------" << observation_target_2.n_elem << endl;


  /* Defining And Assigning the Interfaces */
  
  // Interface Vector 

  // Row<P_DTYPE> al_st_interface_origin = {0.05,-0.01,-0.01};
  // Row<P_DTYPE> al_st_interface_y = {0.05,0.01,-0.01};
  // Row<P_DTYPE> al_st_interface_x = {0.05,-0.01,0.01};
  // int  al_st_x_div = int(vec_mag<P_DTYPE>(al_st_interface_x) / r_s_tran) ;
  // int  al_st_y_div = int(vec_mag<P_DTYPE>(al_st_interface_y) / r_s_tran) ;
  // Mat<P_DTYPE> al_st_interface = rectangle_generator<P_DTYPE>(al_st_interface_x,al_st_interface_y,al_st_interface_origin,al_st_x_div,al_st_y_div);

  Row<P_DTYPE> al_st_int_st = {0.05,-0.2,0};
  Row<P_DTYPE> al_st_int_end = {0.05,0.2,0};
  int al_st_divs = vec_mag<P_DTYPE>(al_st_int_st - al_st_int_end)/r_s_tran;
  Mat<P_DTYPE> al_st_interface = line_generator<P_DTYPE>(source_st,source_end,al_st_divs-1);
  
  al_st_interface.save("al_st_interface.csv",csv_ascii);


  // STEEL interface END
  
  Row<P_DTYPE> normal_to_interface = {1,0,0};
  Row<P_DTYPE> normal_to_interface_reverse = {-1,0,0};

  // Row<P_DTYPE> st_end_interface_origin = {0.1,-0.01,-0.01};
  // Row<P_DTYPE> st_end_interface_y = {0.1,0.01,-0.01};
  // Row<P_DTYPE> st_end_interface_x = {0.1,-0.01,0.01};
  // int  st_end_x_div = int(vec_mag<P_DTYPE>(st_end_interface_x) / r_s_tran) ;
  // int  st_end_y_div = int(vec_mag<P_DTYPE>(st_end_interface_y) / r_s_tran) ;
  // Mat<P_DTYPE> st_end_interface = rectangle_generator<P_DTYPE>(st_end_interface_x,st_end_interface_y,st_end_interface_origin,st_end_x_div,st_end_y_div);

  Row<P_DTYPE> st_end_int_st = {0.1,-0.2,0};
  Row<P_DTYPE> st_end_int_end = {0.1,0.2,0};
  int st_end_divs = vec_mag<P_DTYPE>(al_st_int_st - al_st_int_end)/r_s_tran;
  Mat<P_DTYPE> st_end_interface = line_generator<P_DTYPE>(source_st,source_end,st_end_divs-1);
  
  st_end_interface.save("st_end_interface.csv",csv_ascii); 
  // List of Interfaces Used in secondary Ply

  /*---------------------------------------
    We are using Physical Interfaces While Creating the interfaces
   -----------------------------------------*/
  dpsm::solid_interface<P_DTYPE> interface_al(al_st_interface,normal_to_interface_reverse);
  dpsm::solid_interface<P_DTYPE> interface_al_st(al_st_interface,normal_to_interface);
  dpsm::solid_interface<P_DTYPE> interface_st_end(st_end_interface,normal_to_interface);

  // Creating the Transducer Type For Geometry preparation

  dpsm::transducer<P_DTYPE> ultrasonice_trasducer(source_point_locations,normal_to_location_reverse);
  
  // Preparing all Things for Dpsm Geometry thing

  cout << "Started Creating the Containers for Holding the Geometrical Contents" << endl ;
  
  vector<dpsm::solid_material<P_DTYPE>> material_array_exp(2,dpsm::solid_material<P_DTYPE>());
  vector<dpsm::solid_interface<P_DTYPE>> Interface_array_exp(2,dpsm::solid_interface<P_DTYPE>());
  vector<dpsm::transducer<P_DTYPE>> transducer_array_exp(2,dpsm::transducer<P_DTYPE>());


  cout << "Created the Containters For Holding Geometrical Entities" << endl; 
  // Assign The vector of Things
  
  material_array_exp[0] = aluminium;
  material_array_exp[1] = steel;

  // Assigning the interface
  
  Interface_array_exp[0] = interface_al;
  Interface_array_exp[1] = combine_interface<P_DTYPE>(interface_al_st,interface_st_end);
  
  // Assginig the Transducer
  // We simply leavbe the other as blank in the array thus we nly define one transducer 
  transducer_array_exp[0] = ultrasonice_trasducer;
    
  //grid generation
  cout << "Assigned the Entities in Respective Positions" << endl; 
  //mat obs_plane_location_ply_1 = grid_target_generator(x_space_ply_1,y_space_ply_1);  //   logically the program is complete

  dpsm::geometry<P_DTYPE> TOTAL_LAMINATE(2,material_array_exp,transducer_array_exp,Interface_array_exp,r_s_tran);
  
  cout << "Created The Geometry With The Propery Matreial" << endl; 
  std::vector<dpsm::solid_ply<P_DTYPE>> T_lam_ply_vec = TOTAL_LAMINATE.get_the_material_grid();
  /* Note 
   The get_material_grid funtion puts the percieved or imaginary location of active passive sources and thus doesnt contain the physical location 
   */
  // Completed the Generation of Plane

  cout << "Generated the Compuation Grid With Proper Boundaries" << endl;

  // mat Transducer_sterss_matrix(source_point_locations.n_rows,3,fill::zeros);
  // for (int i = 0; i <Transducer_sterss_matrix.n_rows ; ++i) {
  //   Transducer_sterss_matrix(i,0) = 1;
  // }

  
  cx_3d<P_DTYPE> Transducer_stress_Cx_3d(source_point_locations.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));
  for (int i= 0; i < source_point_locations.n_rows; ++i) {
    Transducer_stress_Cx_3d[i](0,0) = 1 * pow(10,8); 
  }
  //cout <<  Transducer_stress_Cx_3d.size() << endl;
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
  Mat<complex<P_DTYPE>>  Result = get_strength_hetro(T_lam_ply_vec[0].active_sources, T_lam_ply_vec[0].passve_sources,r_s_tran,normal_to_location,Transducer_stress_Cx_3d,k_s_aluminium,k_p_aluminium,rho_al,omega_trans,mu_al,lamda_al);

   cout << "Black sheep 11" << endl ;
  // cx_3d<double> Stress_resultants_Result =stress_from_points(T_lam_ply_vec[0].active_sources, Result, observation_target,k_s_aluminium,k_p_aluminium,rho_al,omega_trans,lamda_al,mu_al);
  // cout << "Black sheep 12" << endl ;


  Result.save("point_strength.dat",csv_ascii);


  std::string PATH_TO_SAVE_RESULTS = "/home/chaithanya/Documents/DPSM_Solid/Results/";

 
  // int Sucess = save_cx_3d(Stress_resultants_Result,PATH_TO_SAVE_RESULTS);

  
  /* -----------------------------------------------------

     ------------------------------------------------------ */
  cx_3d<P_DTYPE> sigma_resultant_interface =stress_from_points<P_DTYPE>(T_lam_ply_vec[0].active_sources, Result, al_st_interface ,k_s_aluminium,k_p_aluminium,rho_al,omega_trans,lamda_al,mu_al);

   cout << "Black sheep 12" << endl ;
 
   // cout << T_lam_ply_vec[1].active_sources.n_rows << endl;
   // cout << "Black sheep 13" << endl ;
   // cout << T_lam_ply_vec[1].passve_sources.n_rows << endl;

   // cout << al_st_interface << endl;
   // cout << "-------------------------------" << endl;
   // cout << T_lam_ply_vec[1].passve_sources << endl;
   // cout << "-------------------------------" << endl;
   Row<P_DTYPE> dummmy_normal = {-1,0,0};
   int Sucess = save_cx_3d<complex<P_DTYPE>>(sigma_resultant_interface,PATH_TO_SAVE_RESULTS);
   Mat<complex<P_DTYPE>> Result_2 = get_strength_hetro(al_st_interface, T_lam_ply_vec[1].passve_sources,r_s_tran,dummmy_normal,sigma_resultant_interface,k_s_aluminium,k_p_aluminium,rho_al,omega_trans,mu_al,lamda_al);
   
  return(0);
}

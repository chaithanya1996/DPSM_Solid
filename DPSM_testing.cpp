/* -------------------------------------------------

Date : 4th April 2019
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
#include "dpsm_core.hpp"
#include <omp.h>

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


template <typename T>
cx_3d<T> join_cx_3d(cx_3d<T> former, cx_3d<T> latter){

  // vector<Mat<complex<T>>> combined;
  cx_3d<T> combined;
  combined.reserve( former.size() + latter.size() ); // preallocate memory
  combined.insert( combined.end(), former.begin(), former.end() );
  combined.insert( combined.end(), latter.begin(), latter.end() );

 return(combined);
}


template <typename TYPE>
inline TYPE calc_lamda(TYPE E, TYPE V ){
  return(E * V / (1+V) / (1-2*V));
}

template <typename TYPE>
inline TYPE calc_mu (TYPE E, TYPE V ){
  return(E / 2 /(1+V));
}

// template <typename T>
// T combine_vector 

int main(){

  // Controlling Parameters For OpenMP
  omp_set_num_threads(8);


   
  using P_DTYPE = double;  // Defining the Precision of the Computation
  
   // Defining Parameters (work in progress)
  P_DTYPE mm = 0.001;
  P_DTYPE cm = 0.01;
  P_DTYPE kHz = 1000;
  P_DTYPE km = 1000;
  // Cubic centimenter
  P_DTYPE cubic_centimeter = pow(10,-6);
  P_DTYPE gms = pow(10,-3);
  
  

  // Wave Speed in Water
  P_DTYPE c_al = 1.48 * km;
  // P-Wave Speed
  P_DTYPE  c_al_pwave = 6.42 * km;
  // S-Wave Speed
  P_DTYPE  c_al_swave = 3.04 * km;


  
  // Density of Water
  P_DTYPE rho_al = 2.7 * gms / cubic_centimeter;


    // Start Phase  1
    /* --------------------------------
      Alumnium  Material Properties
   ----------------------------------*/
  
  P_DTYPE GPa = pow(10,9);
  P_DTYPE E_al = 69 * GPa, V_al = 0.334;
  P_DTYPE E_st = 200 * GPa, V_st = 0.305;
  
  P_DTYPE lamda_al  = 61.17 * GPa;
  P_DTYPE mu_al  = 26.45 * GPa;
  
  P_DTYPE lamda_st  = calc_lamda(E_st,V_st);
  P_DTYPE mu_st  = calc_mu(E_st,V_st);

   /* --------------------------------
        Transducer Properties
   ----------------------------------*/

    // Transducer Properties
  P_DTYPE trans_freq =  50 * kHz;
  P_DTYPE omega_trans = 2 * M_PI * trans_freq ;

  // Transducer Properties in Aluminium

  P_DTYPE k_s_aluminium = omega_trans / c_al_swave  ;
  P_DTYPE k_p_aluminium = omega_trans / c_al_pwave  ;


  // calculating the r_s for given Transducer and Medium 
  int r_s_trans_factor = 5;
  P_DTYPE r_s_tran = r_s_calculator<P_DTYPE>(trans_freq,c_al,r_s_trans_factor);
  cout << " The Distance Between Transducer Source Points --" << r_s_tran << endl;


  
  // Defining The Edges
  Row<P_DTYPE> S_1_origin = {0,0,0};
  Row<P_DTYPE> S_1_vec_x = {20*cm,0,0};
  Row<P_DTYPE> S_1_vec_y = {0,20*cm,0};
  
  int S_1_X_DIVS = vec_mag<P_DTYPE>(S_1_vec_x - S_1_origin)/r_s_tran;
  int S_1_Y_DIVS = vec_mag<P_DTYPE>(S_1_vec_y - S_1_origin)/r_s_tran;
  Row<P_DTYPE> S_1_normal = {0,0,1};
  Row<P_DTYPE> S_1_normal_reverse = {0,0,-1};

  Row<P_DTYPE> S_2_origin = {0,0,5*cm};
  Row<P_DTYPE> S_2_vec_x = {20*cm,0,5*cm};
  Row<P_DTYPE> S_2_vec_y = {0,20*cm,5*cm};
  
  int S_2_X_DIVS = vec_mag<P_DTYPE>(S_2_vec_x - S_2_origin)/r_s_tran;
  int S_2_Y_DIVS = vec_mag<P_DTYPE>(S_2_vec_y - S_2_origin)/r_s_tran;
  Row<P_DTYPE> S_2_normal = {0,0,-1};
  Row<P_DTYPE> S_2_normal_reverse = {0,0,1};

  
  Row<P_DTYPE> S_3_origin = {0,0,0};
  Row<P_DTYPE> S_3_vec_x = {0,0,5*cm};
  Row<P_DTYPE> S_3_vec_y = {0,20*cm,0};
  
  int S_3_X_DIVS = vec_mag<P_DTYPE>(S_3_vec_x - S_3_origin)/r_s_tran;
  int S_3_Y_DIVS = vec_mag<P_DTYPE>(S_3_vec_y - S_3_origin)/r_s_tran;
  Row<P_DTYPE> S_3_normal = {1,0,0};
  Row<P_DTYPE> S_3_normal_reverse = {-1,0,0};
			     
  Row<P_DTYPE> S_4_origin = {20*cm,0,0};
  Row<P_DTYPE> S_4_vec_x = {20*cm,0,5*cm};
  Row<P_DTYPE> S_4_vec_y = {20*cm,20*cm,0};
  
  int S_4_X_DIVS = vec_mag<P_DTYPE>(S_4_vec_x - S_4_origin)/r_s_tran;
  int S_4_Y_DIVS = vec_mag<P_DTYPE>(S_4_vec_y - S_4_origin)/r_s_tran;
  Row<P_DTYPE> S_4_normal = {-1,0,0};
  Row<P_DTYPE> S_4_normal_reverse = {1,0,0};
  
  Row<P_DTYPE> S_5_origin = {0,20*cm,0};
  Row<P_DTYPE> S_5_vec_x = {20*cm,20*cm,0};
  Row<P_DTYPE> S_5_vec_y = {0,20*cm,5*cm};
  
  int S_5_X_DIVS = vec_mag<P_DTYPE>(S_5_vec_x - S_5_origin)/r_s_tran;
  int S_5_Y_DIVS = vec_mag<P_DTYPE>(S_5_vec_y - S_5_origin)/r_s_tran;
  Row<P_DTYPE> S_5_normal = {0,-1,0};
  Row<P_DTYPE> S_5_normal_reverse = {0,1,0};
  
  Row<P_DTYPE> S_6_origin = {0,0,0};
  Row<P_DTYPE> S_6_vec_x = {20*cm,0,0};
  Row<P_DTYPE> S_6_vec_y = {0,0,5*cm};
  
  int S_6_X_DIVS = vec_mag<P_DTYPE>(S_6_vec_x - S_6_origin)/r_s_tran;
  int S_6_Y_DIVS = vec_mag<P_DTYPE>(S_6_vec_y - S_6_origin)/r_s_tran;
  Row<P_DTYPE> S_6_normal = {0,1,0};
  Row<P_DTYPE> S_6_normal_reverse = {0,-1,0};


  
  // Generating the Surfaces
  
  Mat<P_DTYPE> SURFACE_1_MAT =  rectangle_generator<P_DTYPE>( S_1_vec_x, S_1_vec_y, S_1_origin ,S_1_X_DIVS, S_1_Y_DIVS);
  Mat<P_DTYPE> SURFACE_2_MAT =  rectangle_generator<P_DTYPE>( S_2_vec_x, S_2_vec_y, S_2_origin ,S_2_X_DIVS, S_2_Y_DIVS);
  Mat<P_DTYPE> SURFACE_3_MAT =  rectangle_generator<P_DTYPE>( S_3_vec_x, S_3_vec_y, S_3_origin ,S_3_X_DIVS, S_3_Y_DIVS);
  Mat<P_DTYPE> SURFACE_4_MAT =  rectangle_generator<P_DTYPE>( S_4_vec_x, S_4_vec_y, S_4_origin ,S_4_X_DIVS, S_4_Y_DIVS);
  Mat<P_DTYPE> SURFACE_5_MAT =  rectangle_generator<P_DTYPE>( S_5_vec_x, S_5_vec_y, S_5_origin ,S_5_X_DIVS, S_5_Y_DIVS);
  Mat<P_DTYPE> SURFACE_6_MAT =  rectangle_generator<P_DTYPE>( S_6_vec_x, S_6_vec_y, S_6_origin ,S_6_X_DIVS, S_6_Y_DIVS);

  // Exact Locations

  Mat<P_DTYPE> SURFACE_1_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(SURFACE_1_MAT, S_1_normal_reverse, r_s_tran);
  Mat<P_DTYPE> SURFACE_2_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(SURFACE_2_MAT, S_2_normal_reverse, r_s_tran);
  Mat<P_DTYPE> SURFACE_3_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(SURFACE_3_MAT, S_3_normal_reverse, r_s_tran);
  Mat<P_DTYPE> SURFACE_4_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(SURFACE_4_MAT, S_4_normal_reverse, r_s_tran);
  Mat<P_DTYPE> SURFACE_5_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(SURFACE_5_MAT, S_5_normal_reverse, r_s_tran);
  Mat<P_DTYPE> SURFACE_6_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(SURFACE_6_MAT, S_6_normal_reverse, r_s_tran);

  // Defining the Stress Coniditons
  
  Mat<complex<P_DTYPE>> SURFACE_1_STR_BC(SURFACE_1_MAT.n_rows,3,fill::zeros);
  Mat<complex<P_DTYPE>> SURFACE_2_STR_BC(SURFACE_2_MAT.n_rows,3,fill::zeros);
  Mat<complex<P_DTYPE>> SURFACE_3_STR_BC(SURFACE_3_MAT.n_rows,3,fill::zeros);
  Mat<complex<P_DTYPE>> SURFACE_4_STR_BC(SURFACE_4_MAT.n_rows,3,fill::zeros);
  Mat<complex<P_DTYPE>> SURFACE_5_STR_BC(SURFACE_5_MAT.n_rows,3,fill::zeros);
  Mat<complex<P_DTYPE>> SURFACE_6_STR_BC(SURFACE_6_MAT.n_rows,3,fill::zeros);

  // SURFACE_1_STR_BC = SURFACE_1_STR_BC * r_s_tran * 0.00001;
  // SURFACE_2_STR_BC = SURFACE_2_STR_BC * r_s_tran * 0.00001;  
  for (int i = 0; i < SURFACE_1_STR_BC.n_rows; ++i) {
    SURFACE_1_STR_BC(i,0) = SURFACE_1_MAT(i,0) *  r_s_tran* 0.001 + SURFACE_1_MAT(i,1) *  r_s_tran* 0.001 ;
    SURFACE_1_STR_BC(i,1) = SURFACE_1_MAT(i,0) *  r_s_tran* 0.001 + SURFACE_1_MAT(i,1) *  r_s_tran* 0.001 ;
  }

  for (int i = 0; i < SURFACE_2_STR_BC.n_rows; ++i) {
    SURFACE_2_STR_BC(i,0) = SURFACE_2_MAT(i,0) *  r_s_tran* 0.002 + SURFACE_2_MAT(i,1) *  r_s_tran* 0.002 ;
    SURFACE_2_STR_BC(i,2) = SURFACE_2_MAT(i,0) *  r_s_tran* 0.002 + SURFACE_2_MAT(i,1) *  r_s_tran* 0.002 ;
  }


  cout << " Start Stiching" << endl;
  
  Mat<P_DTYPE> ACTIVE_SOURCES,ACTIVE_SOURCES_DPSM;


  ACTIVE_SOURCES = join_cols(SURFACE_1_MAT,SURFACE_2_MAT);
  ACTIVE_SOURCES = join_cols(ACTIVE_SOURCES,SURFACE_3_MAT);
  ACTIVE_SOURCES = join_cols(ACTIVE_SOURCES,SURFACE_4_MAT);
  ACTIVE_SOURCES = join_cols(ACTIVE_SOURCES,SURFACE_5_MAT);
  ACTIVE_SOURCES = join_cols(ACTIVE_SOURCES,SURFACE_6_MAT);


  ACTIVE_SOURCES_DPSM = join_cols(SURFACE_1_MAT_DPSM_SOURCE,SURFACE_2_MAT_DPSM_SOURCE);
  ACTIVE_SOURCES_DPSM = join_cols(ACTIVE_SOURCES_DPSM,SURFACE_3_MAT_DPSM_SOURCE);
  ACTIVE_SOURCES_DPSM = join_cols(ACTIVE_SOURCES_DPSM,SURFACE_4_MAT_DPSM_SOURCE);
  ACTIVE_SOURCES_DPSM = join_cols(ACTIVE_SOURCES_DPSM,SURFACE_5_MAT_DPSM_SOURCE);
  ACTIVE_SOURCES_DPSM = join_cols(ACTIVE_SOURCES_DPSM,SURFACE_6_MAT_DPSM_SOURCE);

  cout << " Joining the Stress Conditions" << endl;


  Mat<complex<P_DTYPE>> ACTIVE_STRESS_BC;
  ACTIVE_STRESS_BC = join_cols(SURFACE_1_STR_BC,SURFACE_2_STR_BC);
  ACTIVE_STRESS_BC = join_cols(ACTIVE_STRESS_BC,SURFACE_3_STR_BC);
  ACTIVE_STRESS_BC = join_cols(ACTIVE_STRESS_BC,SURFACE_4_STR_BC);
  ACTIVE_STRESS_BC = join_cols(ACTIVE_STRESS_BC,SURFACE_5_STR_BC);
  ACTIVE_STRESS_BC = join_cols(ACTIVE_STRESS_BC,SURFACE_6_STR_BC);
  
  

  // complex<P_DTYPE> img_start(0,1.0);
  // complex<P_DTYPE> e_p = exp(img_start * k_p_aluminium * vec_mag<P_DTYPE>(target_point-source_point))/vec_mag<P_DTYPE>(target_point-source_point);
  // //cout <<  G_ps_mat_for_point(target_point,source_point,k_s_aluminium ,k_p_aluminium,rho_al ,omega_trans) << endl;
  // complex<P_DTYPE> e_s = exp(img_start * k_s_aluminium * vec_mag<P_DTYPE>(target_point-source_point))/vec_mag<P_DTYPE>(target_point-source_point);

  // mat A = randu<mat>(1,1);
  // mat B = randu<mat>(1,1);

  // mat K =  kron(A,B);
  // cout << e_p << "-------" << e_s << endl;
  // cout << exp(img_start) << endl;
  // cout << K << endl;

  //cx_4d<P_DTYPE> Result = G_p_diff_ijk(source_point,Target_plane,k_s_aluminium ,k_p_aluminium,rho_al ,omega_trans);
  // Mat<complex<P_DTYPE>> Result = EQN_assembler(Source_Mat,Target_plane ,k_s_aluminium ,k_p_aluminium,rho_al ,omega_trans,mu_al,lamda_al);

  Mat<P_DTYPE> PASS_SOURCES(0,0,fill::zeros);
  
  //Mat<complex<P_DTYPE>> Result = EQN_assembler_DISP(Source_Mat,Target_plane ,k_s_aluminium ,k_p_aluminium,rho_al ,omega_trans);
  Mat<complex<P_DTYPE>> Result =  solve_dpsm_disp<P_DTYPE>(ACTIVE_SOURCES_DPSM,PASS_SOURCES,ACTIVE_STRESS_BC,ACTIVE_SOURCES,k_s_aluminium,k_p_aluminium,rho_al,omega_trans);
  Result.save("D_BUG_Result.csv",csv_ascii);
  // Area of Observation
  
  Row<P_DTYPE> S_OBS_origin = {10*cm,0,0};
  Row<P_DTYPE> S_OBS_vec_x = {10*cm,20*cm,0};
  Row<P_DTYPE> S_OBS_vec_y = {10*cm,0,5*cm};
  
  int S_OBS_X_DIVS = vec_mag<P_DTYPE>(S_OBS_vec_x - S_OBS_origin)/r_s_tran;
  int S_OBS_Y_DIVS = vec_mag<P_DTYPE>(S_OBS_vec_y - S_OBS_origin)/r_s_tran;
  Row<P_DTYPE> S_OBS_normal = {1,0,0};
  Row<P_DTYPE> S_OBS_normal_reverse = {-1,0,0};


  
  // Generating the Surfaces
  
  Mat<P_DTYPE> SURFACE_OBS_MAT =  rectangle_generator<P_DTYPE>( S_OBS_vec_x, S_OBS_vec_y, S_OBS_origin ,S_OBS_X_DIVS * 20, S_OBS_Y_DIVS*10);
  SURFACE_OBS_MAT.save("OBS_SURFACE.csv",csv_ascii);
  Mat<complex<P_DTYPE>> DISP_OBS_CAL = disp_calc_3d_ver<P_DTYPE>(ACTIVE_SOURCES_DPSM, SURFACE_OBS_MAT,Result, k_s_aluminium,k_p_aluminium,rho_al,omega_trans);
  Mat<P_DTYPE> DISP_OBS_CAL_ABS = abs(DISP_OBS_CAL);
  DISP_OBS_CAL_ABS.save("OBS_DISP.csv",csv_ascii);

  cx_3d<P_DTYPE> DISP_OBS_CAL_stress = stress_calc_3d_ver<P_DTYPE>(ACTIVE_SOURCES_DPSM, SURFACE_OBS_MAT,Result, k_s_aluminium,k_p_aluminium,rho_al,omega_trans,mu_al,lamda_al);
  vector<Mat<P_DTYPE>> DISP_OBS_CAL_stress_stress(DISP_OBS_CAL_stress.size(),Mat<P_DTYPE>(3,3));
  for(size_t i = 0; i < DISP_OBS_CAL_stress.size(); i++)
  {
    DISP_OBS_CAL_stress_stress[i] = abs(DISP_OBS_CAL_stress[i]);
  }
  
  save_cx_3d(DISP_OBS_CAL_stress_stress, "OBS_STRESS.csv");
  
  Mat<complex<P_DTYPE>> DISP__CALC_ON_ACTIVE = disp_calc_3d_ver<P_DTYPE>(ACTIVE_SOURCES_DPSM, ACTIVE_SOURCES,Result, k_s_aluminium,k_p_aluminium,rho_al,omega_trans);
  Mat<P_DTYPE> ABS_DISP_ACTIVE =  abs(DISP__CALC_ON_ACTIVE);
  ABS_DISP_ACTIVE.save("DBUG_ACTIVE_DISP.csv",csv_ascii);
  //vector<P_DTYPE> DISP_OBS_CAL_stress_stress = {0,0,0};
  // cout << Result << endl;
 }

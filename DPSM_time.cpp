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
#include "dpsm_interface.hpp"
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

template<typename T>
dpsm::solid_interface<T> combine_interface(dpsm::solid_interface<T> interface1,dpsm::solid_interface<T> interface2){
  Mat<T> joined_mats = join_cols(interface1.interface_source_1,interface2.interface_source_1);
  dpsm::solid_interface <T> joint_interface_solid(joined_mats,interface1.normal_to_solid_interface);
  return(joint_interface_solid);
};

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


   
  using P_DTYPE = float;  // Defining the Precision of the Computation
  
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
  
  P_DTYPE lamda_al  = calc_lamda<P_DTYPE>(E_al,V_al);
  P_DTYPE mu_al  = calc_mu<P_DTYPE>(E_al,V_al);
  
  P_DTYPE lamda_st  = calc_lamda(E_st,V_st);
  P_DTYPE mu_st  = calc_mu(E_st,V_st);

   /* --------------------------------
        Transducer Properties
   ----------------------------------*/

    // Transducer Properties
  P_DTYPE trans_freq =  3000 * kHz;
  P_DTYPE omega_trans = 2 * M_PI * trans_freq ;

  // Transducer Properties in Aluminium

  P_DTYPE k_s_aluminium = omega_trans / c_al_swave  ;
  
  P_DTYPE k_p_aluminium = omega_trans / c_al_pwave  ;


  // calculating the r_s for given Transducer and Medium 
  int r_s_trans_factor = 2;
  P_DTYPE r_s_tran = r_s_calculator<P_DTYPE>(trans_freq,c_al,r_s_trans_factor);
  cout << " The Distance Between Transducer Source Points --" << r_s_tran << endl;

  
  // Defining The Edges

  
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
  
  // Mat<P_DTYPE> SURFACE_1_MAT =  rectangle_generator<P_DTYPE>( S_1_vec_x, S_1_vec_y, S_1_origin ,S_1_X_DIVS, S_1_Y_DIVS);
  // Mat<P_DTYPE> SURFACE_2_MAT =  rectangle_generator<P_DTYPE>( S_2_vec_x, S_2_vec_y, S_2_origin ,S_2_X_DIVS, S_2_Y_DIVS);
  Mat<P_DTYPE> SURFACE_3_MAT =  rectangle_generator<P_DTYPE>( S_3_vec_x, S_3_vec_y, S_3_origin ,S_3_X_DIVS, S_3_Y_DIVS);
  Mat<P_DTYPE> SURFACE_4_MAT =  rectangle_generator<P_DTYPE>( S_4_vec_x, S_4_vec_y, S_4_origin ,S_4_X_DIVS, S_4_Y_DIVS);
  Mat<P_DTYPE> SURFACE_5_MAT =  rectangle_generator<P_DTYPE>( S_5_vec_x, S_5_vec_y, S_5_origin ,S_5_X_DIVS, S_5_Y_DIVS);
  Mat<P_DTYPE> SURFACE_6_MAT =  rectangle_generator<P_DTYPE>( S_6_vec_x, S_6_vec_y, S_6_origin ,S_6_X_DIVS, S_6_Y_DIVS);

  // Exact Locations

  // Mat<P_DTYPE> SURFACE_1_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(SURFACE_1_MAT, S_1_normal_reverse, r_s_tran);
  // Mat<P_DTYPE> SURFACE_2_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(SURFACE_2_MAT, S_2_normal_reverse, r_s_tran);
  Mat<P_DTYPE> SURFACE_3_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(SURFACE_3_MAT, S_3_normal_reverse, r_s_tran);
  Mat<P_DTYPE> SURFACE_4_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(SURFACE_4_MAT, S_4_normal_reverse, r_s_tran);
  Mat<P_DTYPE> SURFACE_5_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(SURFACE_5_MAT, S_5_normal_reverse, r_s_tran);
  Mat<P_DTYPE> SURFACE_6_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(SURFACE_6_MAT, S_6_normal_reverse, r_s_tran);

  // Defining the Stress Coniditons
  
  // Mat<P_DTYPE> SURFACE_1_STR_BC(size(SURFACE_1_MAT_DPSM_SOURCE),fill::zeros);
  // Mat<P_DTYPE> SURFACE_2_STR_BC(size(SURFACE_2_MAT_DPSM_SOURCE),fill::zeros);
  cx_3d<P_DTYPE> SURFACE_3_STR_BC(SURFACE_3_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));
  cx_3d<P_DTYPE> SURFACE_4_STR_BC(SURFACE_4_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));
  cx_3d<P_DTYPE> SURFACE_5_STR_BC(SURFACE_5_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));
  cx_3d<P_DTYPE> SURFACE_6_STR_BC(SURFACE_6_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));




  
    /*==================================================*/
   /*===============Transducer Modelling===============*/
  /*==================================================*/

  P_DTYPE STRESS_VAL_X = pow(10,6);
  P_DTYPE STRESS_VAL_Y = pow(10,6);
  
  Row<P_DTYPE> S_1_TD_P1_origin = {8*cm,8*cm,0};
  Row<P_DTYPE> S_1_TD_P1_vec_x = {10*cm,8*cm,0};
  Row<P_DTYPE> S_1_TD_P1_vec_y = {8*cm,10*cm,0};
  
  int S_1_TD_P1_X_DIVS = vec_mag<P_DTYPE>(S_1_TD_P1_vec_x - S_1_TD_P1_origin)/r_s_tran;
  int S_1_TD_P1_Y_DIVS = vec_mag<P_DTYPE>(S_1_TD_P1_vec_y - S_1_TD_P1_origin)/r_s_tran;
			     
  Row<P_DTYPE> S_1_TD_P1_normal = {0,0,1};
  Row<P_DTYPE> S_1_TD_P1_normal_reverse = {0,0,-1};

  cout << vec_mag<P_DTYPE>(S_1_TD_P1_vec_x - S_1_TD_P1_origin) /r_s_tran << endl;
  cout << vec_mag<P_DTYPE>(S_1_TD_P1_vec_y - S_1_TD_P1_origin)/r_s_tran << endl;
  
  cout << "S_1_TD_P1_MAT--" << S_1_TD_P1_X_DIVS << "--" << S_1_TD_P1_Y_DIVS << endl;

  Mat<P_DTYPE> S_1_TD_P1_MAT =  rectangle_generator<P_DTYPE>( S_1_TD_P1_vec_x , S_1_TD_P1_vec_y, S_1_TD_P1_origin ,S_1_TD_P1_X_DIVS, S_1_TD_P1_Y_DIVS);
  Mat<P_DTYPE> S_1_TD_P1_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(S_1_TD_P1_MAT,  S_1_TD_P1_normal_reverse, r_s_tran);

  
  cx_3d<P_DTYPE> S_1_TD_P1_BC(S_1_TD_P1_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));

  // The Stress Condiitons Needs To Applied
  #pragma omp parallel for
  for (int i= 0; i < S_1_TD_P1_BC.size(); ++i) {
    S_1_TD_P1_BC[i](0,0)=-STRESS_VAL_X;
    S_1_TD_P1_BC[i](1,1)=-STRESS_VAL_Y;
  }


  
  Row<P_DTYPE> S_1_TD_P2_origin = {10*cm,8*cm,0};
  Row<P_DTYPE> S_1_TD_P2_vec_x = {12*cm,8*cm,0};
  Row<P_DTYPE> S_1_TD_P2_vec_y = {10*cm,10*cm,0};
  
  int S_1_TD_P2_X_DIVS = vec_mag<P_DTYPE>(S_1_TD_P2_vec_x - S_1_TD_P2_origin)/r_s_tran;
  int S_1_TD_P2_Y_DIVS = vec_mag<P_DTYPE>(S_1_TD_P2_vec_y - S_1_TD_P2_origin)/r_s_tran;
  Row<P_DTYPE> S_1_TD_P2_normal = {0,0,1};
  Row<P_DTYPE> S_1_TD_P2_normal_reverse = {0,0,-1};

  Mat<P_DTYPE> S_1_TD_P2_MAT =  rectangle_generator<P_DTYPE>( S_1_TD_P2_vec_x , S_1_TD_P2_vec_y, S_1_TD_P2_origin ,S_1_TD_P2_X_DIVS, S_1_TD_P2_Y_DIVS);
  Mat<P_DTYPE> S_1_TD_P2_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(S_1_TD_P2_MAT,  S_1_TD_P2_normal_reverse, r_s_tran);

  cout << "S_1_TD_P1_MAT--" << S_1_TD_P2_X_DIVS << "--" << S_1_TD_P2_Y_DIVS << endl;
    
  cx_3d<P_DTYPE> S_1_TD_P2_BC(S_1_TD_P2_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));

  // The Stress Condiitons Needs To Applied
  #pragma omp parallel for
  for (int i= 0; i < S_1_TD_P2_BC.size(); ++i) {
    S_1_TD_P2_BC[i](0,0)=+STRESS_VAL_X;
    S_1_TD_P2_BC[i](1,1)=-STRESS_VAL_Y;
  }

  // Generating the Stress Conditons For Active Surface
  
  Row<P_DTYPE> S_1_TD_P3_origin = {8*cm,10*cm,0};
  Row<P_DTYPE> S_1_TD_P3_vec_x = {10*cm,10*cm,0};
  Row<P_DTYPE> S_1_TD_P3_vec_y = {8*cm,12*cm,0};
  
  int S_1_TD_P3_X_DIVS = vec_mag<P_DTYPE>(S_1_TD_P3_vec_x - S_1_TD_P3_origin)/r_s_tran;
  int S_1_TD_P3_Y_DIVS = vec_mag<P_DTYPE>(S_1_TD_P3_vec_y - S_1_TD_P3_origin)/r_s_tran;
  Row<P_DTYPE> S_1_TD_P3_normal = {0,0,1};
  Row<P_DTYPE> S_1_TD_P3_normal_reverse = {0,0,-1};

  Mat<P_DTYPE> S_1_TD_P3_MAT =  rectangle_generator<P_DTYPE>( S_1_TD_P3_vec_x , S_1_TD_P3_vec_y, S_1_TD_P3_origin ,S_1_TD_P3_X_DIVS, S_1_TD_P3_Y_DIVS);
  Mat<P_DTYPE> S_1_TD_P3_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(S_1_TD_P3_MAT,  S_1_TD_P3_normal_reverse, r_s_tran);
  cx_3d<P_DTYPE> S_1_TD_P3_BC(S_1_TD_P3_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));

  // The Stress Condiitons Needs To Applied
#pragma omp parallel for
  for (int i= 0; i < S_1_TD_P3_BC.size(); ++i) {
    S_1_TD_P3_BC[i](0,0)=-STRESS_VAL_X;
    S_1_TD_P3_BC[i](1,1)=+STRESS_VAL_Y;
  }
    // Generating the Stress Conditons For Active Surface
  
  Row<P_DTYPE> S_1_TD_P4_origin = {10*cm,10*cm,0};
  Row<P_DTYPE> S_1_TD_P4_vec_x = {12*cm,10*cm,0};
  Row<P_DTYPE> S_1_TD_P4_vec_y = {12*cm,12*cm,0};
  
  int S_1_TD_P4_X_DIVS = vec_mag<P_DTYPE>(S_1_TD_P4_vec_x - S_1_TD_P4_origin)/r_s_tran;
  int S_1_TD_P4_Y_DIVS = vec_mag<P_DTYPE>(S_1_TD_P4_vec_y - S_1_TD_P4_origin)/r_s_tran;
  Row<P_DTYPE> S_1_TD_P4_normal = {0,0,1};
  Row<P_DTYPE> S_1_TD_P4_normal_reverse = {0,0,-1};

  Mat<P_DTYPE> S_1_TD_P4_MAT =  rectangle_generator<P_DTYPE>( S_1_TD_P4_vec_x , S_1_TD_P4_vec_y, S_1_TD_P4_origin ,S_1_TD_P4_X_DIVS, S_1_TD_P4_Y_DIVS);
  Mat<P_DTYPE> S_1_TD_P4_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(S_1_TD_P4_MAT,  S_1_TD_P4_normal_reverse, r_s_tran);

  cx_3d<P_DTYPE> S_1_TD_P4_BC(S_1_TD_P4_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));

  // The Stress Condiitons Needs To Applied
#pragma omp parallel for
  for (int i= 0; i < S_1_TD_P4_BC.size(); ++i) {
    S_1_TD_P4_BC[i](0,0)=+STRESS_VAL_X;
    S_1_TD_P4_BC[i](1,1)=+STRESS_VAL_Y;
  }

  /*=====================================================================*/
  /* First Active Surface is Created  Commencing Creation on the Sencond*/
  /*===================================================================*/

  
  Row<P_DTYPE> S_2_TD_P1_origin = {8*cm,8*cm,5*cm};
  Row<P_DTYPE> S_2_TD_P1_vec_x = {10*cm,8*cm,5*cm};
  Row<P_DTYPE> S_2_TD_P1_vec_y = {8*cm,10*cm,5*cm};
  
  int S_2_TD_P1_X_DIVS = vec_mag<P_DTYPE>(S_2_TD_P1_vec_x - S_2_TD_P1_origin)/r_s_tran;
  int S_2_TD_P1_Y_DIVS = vec_mag<P_DTYPE>(S_2_TD_P1_vec_y - S_2_TD_P1_origin)/r_s_tran;
  Row<P_DTYPE> S_2_TD_P1_normal = {0,0,1};
  Row<P_DTYPE> S_2_TD_P1_normal_reverse = {0,0,-1};

  Mat<P_DTYPE> S_2_TD_P1_MAT =  rectangle_generator<P_DTYPE>( S_2_TD_P1_vec_x , S_2_TD_P1_vec_y, S_2_TD_P1_origin ,S_2_TD_P1_X_DIVS, S_2_TD_P1_Y_DIVS);
  Mat<P_DTYPE> S_2_TD_P1_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(S_2_TD_P1_MAT,  S_2_TD_P1_normal, r_s_tran);

  
  cx_3d<P_DTYPE> S_2_TD_P1_BC(S_2_TD_P1_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));

  // The Stress Condiitons Needs To Applied
  #pragma omp parallel for
  for (int i= 0; i < S_2_TD_P1_BC.size(); ++i) {
    S_2_TD_P1_BC[i](0,0)=-STRESS_VAL_X;
    S_2_TD_P1_BC[i](1,1)=-STRESS_VAL_Y;
  }


  
  Row<P_DTYPE> S_2_TD_P2_origin = {10*cm,8*cm,5*cm};
  Row<P_DTYPE> S_2_TD_P2_vec_x = {12*cm,8*cm,5*cm};
  Row<P_DTYPE> S_2_TD_P2_vec_y = {10*cm,10*cm,5*cm};

  
  int S_2_TD_P2_X_DIVS = vec_mag<P_DTYPE>(S_2_TD_P2_vec_x - S_2_TD_P2_origin)/r_s_tran;
  int S_2_TD_P2_Y_DIVS = vec_mag<P_DTYPE>(S_2_TD_P2_vec_y - S_2_TD_P2_origin)/r_s_tran;
  Row<P_DTYPE> S_2_TD_P2_normal = {0,0,1};
  Row<P_DTYPE> S_2_TD_P2_normal_reverse = {0,0,-1};

  Mat<P_DTYPE> S_2_TD_P2_MAT =  rectangle_generator<P_DTYPE>( S_2_TD_P2_vec_x , S_2_TD_P2_vec_y, S_2_TD_P2_origin ,S_2_TD_P2_X_DIVS, S_2_TD_P2_Y_DIVS);
  Mat<P_DTYPE> S_2_TD_P2_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(S_2_TD_P2_MAT,  S_2_TD_P2_normal, r_s_tran);

    
  cx_3d<P_DTYPE> S_2_TD_P2_BC(S_2_TD_P2_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));

  // The Stress Condiitons Needs To Applied
  #pragma omp parallel for
  for (int i= 0; i < S_2_TD_P2_BC.size(); ++i) {
    S_2_TD_P2_BC[i](0,0)=+STRESS_VAL_X;
    S_2_TD_P2_BC[i](1,1)=-STRESS_VAL_Y;
  }

  // Generating the Stress Conditons For Active Surface
  
  Row<P_DTYPE> S_2_TD_P3_origin = {8*cm,10*cm,5*cm};
  Row<P_DTYPE> S_2_TD_P3_vec_x = {10*cm,10*cm,5*cm};
  Row<P_DTYPE> S_2_TD_P3_vec_y = {8*cm,12*cm,5*cm};
  
  int S_2_TD_P3_X_DIVS = vec_mag<P_DTYPE>(S_2_TD_P3_vec_x - S_2_TD_P3_origin)/r_s_tran;
  int S_2_TD_P3_Y_DIVS = vec_mag<P_DTYPE>(S_2_TD_P3_vec_y - S_2_TD_P3_origin)/r_s_tran;
  Row<P_DTYPE> S_2_TD_P3_normal = {0,0,1};
  Row<P_DTYPE> S_2_TD_P3_normal_reverse = {0,0,-1};

  Mat<P_DTYPE> S_2_TD_P3_MAT =  rectangle_generator<P_DTYPE>( S_2_TD_P3_vec_x , S_2_TD_P3_vec_y, S_2_TD_P3_origin ,S_2_TD_P3_X_DIVS, S_2_TD_P3_Y_DIVS);
  Mat<P_DTYPE> S_2_TD_P3_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(S_2_TD_P3_MAT,  S_2_TD_P3_normal, r_s_tran);
  cx_3d<P_DTYPE> S_2_TD_P3_BC(S_2_TD_P3_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));

  // The Stress Condiitons Needs To Applied
#pragma omp parallel for
  for (int i= 0; i < S_2_TD_P3_BC.size(); ++i) {
    S_2_TD_P3_BC[i](0,0)=-STRESS_VAL_X;
    S_2_TD_P3_BC[i](1,1)=+STRESS_VAL_Y;
  }
    // Generating the Stress Conditons For Active Surface
  
  Row<P_DTYPE> S_2_TD_P4_origin = {10*cm,10*cm,5*cm};
  Row<P_DTYPE> S_2_TD_P4_vec_x = {12*cm,10*cm,5*cm};
  Row<P_DTYPE> S_2_TD_P4_vec_y = {12*cm,12*cm,5*cm};

  
  int S_2_TD_P4_X_DIVS = vec_mag<P_DTYPE>(S_2_TD_P4_vec_x - S_2_TD_P4_origin)/r_s_tran;
  int S_2_TD_P4_Y_DIVS = vec_mag<P_DTYPE>(S_2_TD_P4_vec_y - S_2_TD_P4_origin)/r_s_tran;
  Row<P_DTYPE> S_2_TD_P4_normal = {0,0,1};
  Row<P_DTYPE> S_2_TD_P4_normal_reverse = {0,0,-1};

  Mat<P_DTYPE> S_2_TD_P4_MAT =  rectangle_generator<P_DTYPE>( S_2_TD_P4_vec_x , S_2_TD_P4_vec_y, S_2_TD_P4_origin ,S_2_TD_P4_X_DIVS, S_2_TD_P4_Y_DIVS);
  Mat<P_DTYPE> S_2_TD_P4_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(S_2_TD_P4_MAT,  S_2_TD_P4_normal, r_s_tran);

  cx_3d<P_DTYPE> S_2_TD_P4_BC(S_2_TD_P4_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));

  // The Stress Condiitons Needs To Applied
#pragma omp parallel for
  for (int i= 0; i < S_2_TD_P4_BC.size(); ++i) {
    S_2_TD_P4_BC[i](0,0)=+STRESS_VAL_X;
    S_2_TD_P4_BC[i](1,1)=+STRESS_VAL_Y;
  }

  /*====================================================*/


  
  Row<P_DTYPE> S_1_P1_origin = {0,0,0};
  Row<P_DTYPE> S_1_P1_vec_x = {20*cm,0,0};
  Row<P_DTYPE> S_1_P1_vec_y = {0,8*cm,0};
  
  int S_1_P1_X_DIVS = vec_mag<P_DTYPE>(S_1_P1_vec_x - S_1_P1_origin)/r_s_tran;
  int S_1_P1_Y_DIVS = vec_mag<P_DTYPE>(S_1_P1_vec_y - S_1_P1_origin)/r_s_tran;
  Row<P_DTYPE> S_1_P1_normal = {0,0,1};
  Row<P_DTYPE> S_1_P1_normal_reverse = {0,0,-1};

  Mat<P_DTYPE> S_1_P1_MAT =  rectangle_generator<P_DTYPE>( S_1_P1_vec_x , S_1_P1_vec_y ,S_1_P1_origin , S_1_P1_X_DIVS , S_1_P1_Y_DIVS );
  Mat<P_DTYPE>  S_1_P1_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(S_1_P1_MAT,  S_1_P1_normal_reverse , r_s_tran);
  cx_3d<P_DTYPE> S_1_P1_BC(S_1_P1_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));
  
  
  Row<P_DTYPE> S_1_P2_origin = {0,12*cm,0};
  Row<P_DTYPE> S_1_P2_vec_x = {20*cm,12*cm,0};
  Row<P_DTYPE> S_1_P2_vec_y = {0,20*cm,0};
  
  int S_1_P2_X_DIVS = vec_mag<P_DTYPE>(S_1_P2_vec_x - S_1_P2_origin)/r_s_tran;
  int S_1_P2_Y_DIVS = vec_mag<P_DTYPE>(S_1_P2_vec_y - S_1_P2_origin)/r_s_tran;
  Row<P_DTYPE> S_1_P2_normal = {0,0,-1};
  Row<P_DTYPE> S_1_P2_normal_reverse = {0,0,1};

  Mat<P_DTYPE> S_1_P2_MAT =  rectangle_generator<P_DTYPE>( S_1_P2_vec_x , S_1_P2_vec_y ,S_1_P2_origin , S_1_P2_X_DIVS , S_1_P2_Y_DIVS );
  Mat<P_DTYPE>  S_1_P2_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(S_1_P2_MAT,  S_1_P2_normal_reverse , r_s_tran);
  cx_3d<P_DTYPE> S_1_P2_BC(S_1_P2_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));
  

  Row<P_DTYPE> S_1_P3_origin = {0,8*cm,0};
  Row<P_DTYPE> S_1_P3_vec_x = {8*cm,8*cm,0};
  Row<P_DTYPE> S_1_P3_vec_y = {0,12*cm,0};
  
  int S_1_P3_X_DIVS = vec_mag<P_DTYPE>(S_1_P3_vec_x - S_1_P3_origin)/r_s_tran;
  int S_1_P3_Y_DIVS = vec_mag<P_DTYPE>(S_1_P3_vec_y - S_1_P3_origin)/r_s_tran;
  Row<P_DTYPE> S_1_P3_normal = {0,0,-1};
  Row<P_DTYPE> S_1_P3_normal_reverse = {0,0,1};

  Mat<P_DTYPE> S_1_P3_MAT =  rectangle_generator<P_DTYPE>( S_1_P3_vec_x , S_1_P3_vec_y ,S_1_P3_origin , S_1_P3_X_DIVS , S_1_P3_Y_DIVS );
  Mat<P_DTYPE>  S_1_P3_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(S_1_P3_MAT,  S_1_P3_normal_reverse , r_s_tran);
  cx_3d<P_DTYPE> S_1_P3_BC(S_1_P3_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));
  


  Row<P_DTYPE> S_1_P4_origin = {12*cm,8*cm,0};
  Row<P_DTYPE> S_1_P4_vec_x = {20*cm,8*cm,0};
  Row<P_DTYPE> S_1_P4_vec_y = {12*cm,12*cm,0};
  
  int S_1_P4_X_DIVS = vec_mag<P_DTYPE>(S_1_P4_vec_x - S_1_P4_origin)/r_s_tran;
  int S_1_P4_Y_DIVS = vec_mag<P_DTYPE>(S_1_P4_vec_y - S_1_P4_origin)/r_s_tran;
  Row<P_DTYPE> S_1_P4_normal = {0,0,-1};
  Row<P_DTYPE> S_1_P4_normal_reverse = {0,0,1};

  Mat<P_DTYPE> S_1_P4_MAT =  rectangle_generator<P_DTYPE>( S_1_P4_vec_x , S_1_P4_vec_y ,S_1_P4_origin , S_1_P4_X_DIVS , S_1_P4_Y_DIVS );
  Mat<P_DTYPE>  S_1_P4_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(S_1_P4_MAT,  S_1_P4_normal_reverse , r_s_tran);
  cx_3d<P_DTYPE> S_1_P4_BC(S_1_P4_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));
  


  // Defining the Next Surface

  Row<P_DTYPE> S_2_P1_origin = {0,0,5*cm};
  Row<P_DTYPE> S_2_P1_vec_x = {20*cm,0,5*cm};
  Row<P_DTYPE> S_2_P1_vec_y = {0,8*cm,5*cm};
  
  int S_2_P1_X_DIVS = vec_mag<P_DTYPE>(S_2_P1_vec_x - S_2_P1_origin)/r_s_tran;
  int S_2_P1_Y_DIVS = vec_mag<P_DTYPE>(S_2_P1_vec_y - S_2_P1_origin)/r_s_tran;
  Row<P_DTYPE> S_2_P1_normal = {0,0,1};
  Row<P_DTYPE> S_2_P1_normal_reverse = {0,0,-1};

  Mat<P_DTYPE> S_2_P1_MAT =  rectangle_generator<P_DTYPE>( S_2_P1_vec_x , S_2_P1_vec_y ,S_2_P1_origin , S_2_P1_X_DIVS , S_2_P1_Y_DIVS );
  Mat<P_DTYPE>  S_2_P1_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(S_2_P1_MAT,  S_2_P1_normal_reverse , r_s_tran);
  cx_3d<P_DTYPE> S_2_P1_BC(S_2_P1_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));
  

  
  Row<P_DTYPE> S_2_P2_origin = {0,12*cm,5*cm};
  Row<P_DTYPE> S_2_P2_vec_x = {20*cm,12*cm,5*cm};
  Row<P_DTYPE> S_2_P2_vec_y = {0,20*cm,5*cm};

  
  int S_2_P2_X_DIVS = vec_mag<P_DTYPE>(S_2_P2_vec_x - S_2_P2_origin)/r_s_tran;
  int S_2_P2_Y_DIVS = vec_mag<P_DTYPE>(S_2_P2_vec_y - S_2_P2_origin)/r_s_tran;
  Row<P_DTYPE> S_2_P2_normal = {0,0,-1};
  Row<P_DTYPE> S_2_P2_normal_reverse = {0,0,1};

  Mat<P_DTYPE> S_2_P2_MAT =  rectangle_generator<P_DTYPE>( S_2_P2_vec_x , S_2_P2_vec_y ,S_2_P2_origin , S_2_P2_X_DIVS , S_2_P2_Y_DIVS );
  Mat<P_DTYPE>  S_2_P2_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(S_2_P2_MAT,  S_2_P2_normal_reverse , r_s_tran);
  cx_3d<P_DTYPE> S_2_P2_BC(S_2_P2_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));
  


  Row<P_DTYPE> S_2_P3_origin = {0,8*cm,5*cm};
  Row<P_DTYPE> S_2_P3_vec_x = {8*cm,8*cm,5*cm};
  Row<P_DTYPE> S_2_P3_vec_y = {0,12*cm,5*cm};

  
  int S_2_P3_X_DIVS = vec_mag<P_DTYPE>(S_2_P3_vec_x - S_2_P3_origin)/r_s_tran;
  int S_2_P3_Y_DIVS = vec_mag<P_DTYPE>(S_2_P3_vec_y - S_2_P3_origin)/r_s_tran;
  Row<P_DTYPE> S_2_P3_normal = {0,0,-1};
  Row<P_DTYPE> S_2_P3_normal_reverse = {0,0,1};


  Mat<P_DTYPE> S_2_P3_MAT =  rectangle_generator<P_DTYPE>( S_2_P3_vec_x , S_2_P3_vec_y ,S_2_P3_origin , S_2_P3_X_DIVS , S_2_P3_Y_DIVS );
  Mat<P_DTYPE>  S_2_P3_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(S_2_P3_MAT,  S_2_P3_normal_reverse , r_s_tran);
  cx_3d<P_DTYPE> S_2_P3_BC(S_2_P3_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));
  

  
  Row<P_DTYPE> S_2_P4_origin = {12*cm,8*cm,5*cm};
  Row<P_DTYPE> S_2_P4_vec_x = {20*cm,8*cm,5*cm};
  Row<P_DTYPE> S_2_P4_vec_y = {12*cm,12*cm,5*cm};

  
  int S_2_P4_X_DIVS = vec_mag<P_DTYPE>(S_2_P4_vec_x - S_2_P4_origin)/r_s_tran;
  int S_2_P4_Y_DIVS = vec_mag<P_DTYPE>(S_2_P4_vec_y - S_2_P4_origin)/r_s_tran;
  Row<P_DTYPE> S_2_P4_normal = {0,0,-1};
  Row<P_DTYPE> S_2_P4_normal_reverse = {0,0,1};

  Mat<P_DTYPE> S_2_P4_MAT =  rectangle_generator<P_DTYPE>( S_2_P4_vec_x , S_2_P4_vec_y ,S_2_P4_origin , S_2_P4_X_DIVS , S_2_P4_Y_DIVS );
  Mat<P_DTYPE>  S_2_P4_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(S_2_P4_MAT,  S_2_P4_normal_reverse , r_s_tran);
  cx_3d<P_DTYPE> S_2_P4_BC(S_2_P4_MAT_DPSM_SOURCE.n_rows,Mat<complex<P_DTYPE>>(3,3,fill::zeros));
  
  // Combining the Whole Matrices and Cx_3d

  // Now We need to Combine All Active Sources

  Mat<P_DTYPE> ACTIVE_SOURCES,ACTIVE_SOURCES_DPSM;
  Mat<P_DTYPE> SURFACE_1_MAT,SURFACE_2_MAT;

  cout << " Joining the Geometries" << endl;
  /*Joining The Surface Parts */
  
  SURFACE_1_MAT = join_cols(S_1_TD_P1_MAT, S_1_TD_P2_MAT);
  SURFACE_1_MAT = join_cols(SURFACE_1_MAT, S_1_TD_P3_MAT);
  SURFACE_1_MAT = join_cols(SURFACE_1_MAT, S_1_TD_P4_MAT);

  
  // SURFACE_1_MAT = join_cols(SURFACE_1_MAT, S_1_P1_MAT);
  // SURFACE_1_MAT = join_cols(SURFACE_1_MAT, S_1_P2_MAT);
  // SURFACE_1_MAT = join_cols(SURFACE_1_MAT, S_1_P3_MAT);
  // SURFACE_1_MAT = join_cols(SURFACE_1_MAT, S_1_P4_MAT);

  /* Similarly For Surface 2 */
  SURFACE_2_MAT = join_cols(S_2_TD_P1_MAT, S_2_TD_P2_MAT);
  SURFACE_2_MAT = join_cols(SURFACE_2_MAT, S_2_TD_P3_MAT);
  SURFACE_2_MAT = join_cols(SURFACE_2_MAT, S_2_TD_P4_MAT);

  // SURFACE_2_MAT = join_cols(SURFACE_2_MAT, S_2_P1_MAT);
  // SURFACE_2_MAT = join_cols(SURFACE_2_MAT, S_2_P2_MAT);
  // SURFACE_2_MAT = join_cols(SURFACE_2_MAT, S_2_P3_MAT);
  // SURFACE_2_MAT = join_cols(SURFACE_2_MAT, S_2_P4_MAT);
 // Combining All the Surfaces

  ACTIVE_SOURCES = join_cols(SURFACE_1_MAT,SURFACE_2_MAT);
  // ACTIVE_SOURCES = join_cols(ACTIVE_SOURCES,SURFACE_3_MAT);
  // ACTIVE_SOURCES = join_cols(ACTIVE_SOURCES,SURFACE_4_MAT);
  // ACTIVE_SOURCES = join_cols(ACTIVE_SOURCES,SURFACE_5_MAT);
  // ACTIVE_SOURCES = join_cols(ACTIVE_SOURCES,SURFACE_6_MAT);

    /*Joining The Surface Parts of DPSM*/
  Row<P_DTYPE>S_1_normal = {0,0,1};
  Row<P_DTYPE>S_1_normal_reverse = {0,0,-1};
  Row<P_DTYPE>S_2_normal = {0,0,-1};
  Row<P_DTYPE>S_2_normal_reverse = {0,0,1};


  cout << "Creating the DPSM for Split Geometries" << endl;
  Mat<P_DTYPE> SURFACE_1_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(SURFACE_1_MAT, S_1_normal_reverse, r_s_tran);
  Mat<P_DTYPE> SURFACE_2_MAT_DPSM_SOURCE = source_point_placer<P_DTYPE>(SURFACE_2_MAT, S_2_normal_reverse, r_s_tran);

  // Combining the DPSM
  
  ACTIVE_SOURCES_DPSM = join_cols(SURFACE_1_MAT_DPSM_SOURCE,SURFACE_2_MAT_DPSM_SOURCE);
  // ACTIVE_SOURCES_DPSM = join_cols(ACTIVE_SOURCES_DPSM,SURFACE_3_MAT_DPSM_SOURCE);
  // ACTIVE_SOURCES_DPSM = join_cols(ACTIVE_SOURCES_DPSM,SURFACE_4_MAT_DPSM_SOURCE);
  // ACTIVE_SOURCES_DPSM = join_cols(ACTIVE_SOURCES_DPSM,SURFACE_5_MAT_DPSM_SOURCE);
  // ACTIVE_SOURCES_DPSM = join_cols(ACTIVE_SOURCES_DPSM,SURFACE_6_MAT_DPSM_SOURCE);

  cout << " Joining the Stress Conditions" << endl;


  
    /*Joining The Surface Parts */
  cx_3d<P_DTYPE> SURFACE_1_STR_BC,SURFACE_2_STR_BC;

  
  SURFACE_1_STR_BC = join_cx_3d(S_1_TD_P1_BC, S_1_TD_P2_BC);
  SURFACE_1_STR_BC = join_cx_3d(SURFACE_1_STR_BC, S_1_TD_P3_BC);
  SURFACE_1_STR_BC = join_cx_3d(SURFACE_1_STR_BC, S_1_TD_P4_BC);

  // SURFACE_1_STR_BC = join_cx_3d(SURFACE_1_STR_BC, S_1_P1_BC);
  // SURFACE_1_STR_BC = join_cx_3d(SURFACE_1_STR_BC, S_1_P2_BC);
  // SURFACE_1_STR_BC = join_cx_3d(SURFACE_1_STR_BC, S_1_P3_BC);
  // SURFACE_1_STR_BC = join_cx_3d(SURFACE_1_STR_BC, S_1_P4_BC);



  // Combining All the Surface DPSM Sources

  
  /* Similarly For Surface 2 */
  SURFACE_2_STR_BC = join_cx_3d(S_2_TD_P1_BC, S_2_TD_P2_BC);
  SURFACE_2_STR_BC = join_cx_3d(SURFACE_2_STR_BC, S_2_TD_P3_BC);
  SURFACE_2_STR_BC = join_cx_3d(SURFACE_2_STR_BC, S_2_TD_P4_BC);

  // SURFACE_2_STR_BC = join_cx_3d(SURFACE_2_STR_BC, S_2_P1_BC);
  // SURFACE_2_STR_BC = join_cx_3d(SURFACE_2_STR_BC, S_2_P2_BC);
  // SURFACE_2_STR_BC = join_cx_3d(SURFACE_2_STR_BC, S_2_P3_BC);
  // SURFACE_2_STR_BC = join_cx_3d(SURFACE_2_STR_BC, S_2_P4_BC);

  cx_3d<P_DTYPE> ACTIVE_STRESS_BC;
  ACTIVE_STRESS_BC = join_cx_3d<P_DTYPE>(SURFACE_1_STR_BC,SURFACE_2_STR_BC);
  // ACTIVE_STRESS_BC = join_cx_3d<P_DTYPE>(ACTIVE_STRESS_BC,SURFACE_3_STR_BC);
  // ACTIVE_STRESS_BC = join_cx_3d<P_DTYPE>(ACTIVE_STRESS_BC,SURFACE_4_STR_BC);
  // ACTIVE_STRESS_BC = join_cx_3d<P_DTYPE>(ACTIVE_STRESS_BC,SURFACE_5_STR_BC);
  // ACTIVE_STRESS_BC = join_cx_3d<P_DTYPE>(ACTIVE_STRESS_BC,SURFACE_6_STR_BC);
  
  // process_mem_usage(vm, rss);
  // cout << "VM: " << vm << "; RSS: " << rss << endl;
  // Mat<complex<P_DTYPE>> Result =  get_strength_hetro(Mat<T> source_point_list,Mat<T> passive_source_list ,T r_s, Row<T> direction_cosine,cx_3d<T> stress_matrix,T k_s,T k_p,T rho , T omega,T mu,T lamda,bool have_active_sources)
  ACTIVE_SOURCES.save("DBUG_TIME_ACTIVE.csv",csv_ascii);
  ACTIVE_SOURCES_DPSM.save("DBUG_TIME_ACTIVE_DPSM.csv",csv_ascii);
  cout << ACTIVE_STRESS_BC.size() << "-------" << ACTIVE_SOURCES_DPSM.n_rows << endl;
  // Mat<P_DTYPE> PASS_SOURCES(0,0,fill::zeros);
  
  // Mat<complex<P_DTYPE>> Result =  solve_dpsm_str<P_DTYPE>(ACTIVE_SOURCES_DPSM,PASS_SOURCES,ACTIVE_STRESS_BC,ACTIVE_SOURCES,k_s_aluminium,k_p_aluminium,rho_al,omega_trans,mu_al,lamda_al);
  // Result.save("D_BUG_Result.csv",csv_ascii);

  // Mat<P_DTYPE> Observation_point = {{10*cm,10*cm,0},{5*cm,5*cm,0}};

  // cx_3d<P_DTYPE> DBUG_Stress_Result = stress_from_points<P_DTYPE>(ACTIVE_SOURCES_DPSM ,Result , Observation_point ,k_s_aluminium,k_p_aluminium,rho_al,omega_trans,lamda_al,mu_al);
  // int Sucess = save_cx_3d<complex<P_DTYPE>>(DBUG_Stress_Result,"./");
  cout << "Sucess Reading " << endl;	
 }

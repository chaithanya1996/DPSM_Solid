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

  cout << "lamda_al--" << lamda_al << " mu_al--" << mu_al <<endl;
   /* --------------------------------
        Transducer Properties
   ----------------------------------*/

    // Transducer Properties
  P_DTYPE trans_freq =  3 * kHz;
  P_DTYPE omega_trans = 2 * M_PI * trans_freq ;

  cout << "trans_freq--" << trans_freq << endl;
  cout << "omega_trans--" << omega_trans << endl;
  
  // Transducer Properties in Aluminium

  P_DTYPE k_s_aluminium = omega_trans / c_al_swave  ;
  P_DTYPE k_p_aluminium = omega_trans / c_al_pwave  ;

  cout << "k_s_aluminium--" << k_s_aluminium << " k_p_aluminium--" << k_p_aluminium << endl;

  // calculating the r_s for given Transducer and Medium 
  int r_s_trans_factor = 5;
  P_DTYPE r_s_tran = r_s_calculator<P_DTYPE>(trans_freq,c_al,r_s_trans_factor);
  cout << " The Distance Between Transducer Source Points --" << r_s_tran << endl;

  Row<P_DTYPE> source_point = {1,1,1};
  Row<P_DTYPE> target_point = {2,2,2};
  Mat<P_DTYPE> Target_plane = {{2,2,2}};

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

  cx_4d<P_DTYPE> Result = G_p_diff_ijk(source_point,Target_plane,k_s_aluminium ,k_p_aluminium,rho_al ,omega_trans);
  cout << Result[0] << endl;
 }

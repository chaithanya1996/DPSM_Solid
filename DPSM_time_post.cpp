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



template <typename TYPE>
inline TYPE calc_lamda(TYPE E, TYPE V ){
  return(E * V / (1+V) / (1-2*V));
}

template <typename TYPE>
inline TYPE calc_mu (TYPE E, TYPE V ){
  return(E / 2 /(1+V));
}

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
  
  P_DTYPE lamda_al  = calc_lamda<P_DTYPE>(E_al,V_al);
  P_DTYPE mu_al  = calc_mu<P_DTYPE>(E_al,V_al);
  
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
  Row<P_DTYPE> source_point_r = {1,1,1};
  Row<P_DTYPE> target_point_r = {2,2,2};
  Cube<complex<P_DTYPE>> G_CUBE = G_diff_ijk_point(source_point_r,target_point_r,k_s_aluminium, k_p_aluminium,rho_al ,omega_trans);
  cout <<  G_diff_ijk_point(source_point_r,target_point_r,k_s_aluminium, k_p_aluminium,rho_al ,omega_trans) << endl;
  cout << "--------------------------------------------" << endl;
  for (int j = 0; j < 3; ++j) {
    cout <<  CUBE_EXTRACT(G_CUBE,j) << endl;
  }
  cout << "--------------------------------------------" << endl;
  cout <<  stress_coff_point<P_DTYPE>(G_CUBE, lamda_al,mu_al) << endl;

  Cube<complex<P_DTYPE>> G_CUB(3,3,3,fill::randu);
  cout << "--------------------------------------------" << endl;
  cout << G_CUB << endl;

  cout << "--------------------------------------------" << endl;
   for (int j = 0; j < 3; ++j) {
    cout <<  CUBE_EXTRACT(G_CUB,j) << endl;
  }
}  

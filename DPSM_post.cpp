
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

int main(int argc, char** argv){
  
   
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
  
  Row<P_DTYPE> S_OBS_origin = {10*cm,0,0};
  Row<P_DTYPE> S_OBS_vec_x = {10*cm,20*cm,0};
  Row<P_DTYPE> S_OBS_vec_y = {10*cm,0,5*cm};
  
  int S_OBS_X_DIVS = vec_mag<P_DTYPE>(S_OBS_vec_x - S_OBS_origin)/r_s_tran;
  int S_OBS_Y_DIVS = vec_mag<P_DTYPE>(S_OBS_vec_y - S_OBS_origin)/r_s_tran;
  Row<P_DTYPE> S_OBS_normal = {1,0,0};
  Row<P_DTYPE> S_OBS_normal_reverse = {-1,0,0};

  int XDIVS = S_OBS_X_DIVS * 20;
  int YDIVS = S_OBS_Y_DIVS*10;

  // Generating the Surfaces
  Mat<P_DTYPE> SURFACE_OBS_MAT =  rectangle_generator<P_DTYPE>( S_OBS_vec_x, S_OBS_vec_y, S_OBS_origin ,XDIVS, YDIVS);

  Mat<P_DTYPE> LOAD_MAT;
  LOAD_MAT.load(argv[1],csv_ascii);
  cout << argv[0] << endl;
  cout << LOAD_MAT.n_rows << endl;
  cout << LOAD_MAT.n_cols << endl;
  
  Col<P_DTYPE> COL_1 = LOAD_MAT.col(0);
  Col<P_DTYPE> COL_2 = LOAD_MAT.col(1);
  Col<P_DTYPE> COL_3 = LOAD_MAT.col(2);
  
  Mat<P_DTYPE> MESH_GRID_X(YDIVS,XDIVS ,fill::zeros);
  Mat<P_DTYPE> MESH_GRID_Y(YDIVS,XDIVS ,fill::zeros);
  Mat<P_DTYPE> MESH_GRID_Z(YDIVS,XDIVS ,fill::zeros);
  
  Mat<P_DTYPE> MESH_GRID_DISP_1(YDIVS,XDIVS,fill::zeros);
  Mat<P_DTYPE> MESH_GRID_DISP_2(YDIVS,XDIVS,fill::zeros);
  Mat<P_DTYPE> MESH_GRID_DISP_3(YDIVS,XDIVS,fill::zeros);
  
  for (int i= 0; i < YDIVS; ++i) {
    cout << i << endl;
    MESH_GRID_DISP_1.row(i) = trans(COL_1(span(i*XDIVS,(i+1) * XDIVS - 1)));
    MESH_GRID_DISP_2.row(i) = trans(COL_2(span(i*XDIVS,(i+1) * XDIVS - 1)));
    MESH_GRID_DISP_3.row(i) = trans(COL_3(span(i*XDIVS,(i+1) * XDIVS - 1)));
    MESH_GRID_X.row(i) = trans(SURFACE_OBS_MAT(span(i*XDIVS,(i+1) * XDIVS - 1),0));
    MESH_GRID_Y.row(i) = trans(SURFACE_OBS_MAT(span(i*XDIVS,(i+1) * XDIVS -1 ),1));
    MESH_GRID_Z.row(i) = trans( SURFACE_OBS_MAT(span(i*XDIVS,(i+1) * XDIVS -1 ),2));
  }
  MESH_GRID_DISP_1.save("MESH_GRID_DISP_1.csv",csv_ascii);
  MESH_GRID_DISP_2.save("MESH_GRID_DISP_2.csv",csv_ascii);
  MESH_GRID_DISP_3.save("MESH_GRID_DISP_3.csv",csv_ascii);
  MESH_GRID_X.save("MESH_GRID_X.csv",csv_ascii);
  MESH_GRID_Y.save("MESH_GRID_Y.csv",csv_ascii);
  MESH_GRID_Z.save("MESH_GRID_Z.csv",csv_ascii);
}
   
 

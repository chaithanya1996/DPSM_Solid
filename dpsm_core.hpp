#include "dpsm_helpers.hpp"


template <typename TYPE>TYPE r_s_calculator(TYPE freq, TYPE c);

// Kron Delta _funtion

template<typename T>
T kron_delta(int i, int j);

// Templated G_ij Basic Calulation

template <typename T>
Row<complex<T>> G_ij(int i_sub, int j_sub, Row<T> affect_point , Row<T> source_point,T k_s, T k_p ,T rho, T omega);

// Green Function Calculation Source to Target ( One Point - One Point)
template <typename T>
Cube<complex<T>> G_ps_mat_for_point(Row<T> affect_point, Row<T> source_point,T k_s , T k_p, T rho , T omega);

// Second Version for facilitate Displacemnt Caclualtuon
template <typename T>
Mat<complex<T>> G_mat_for_point(Row<T> affect_point , Row<T> source_point,T k_s ,T k_p,T rho ,T omega);

// Genrate G from Source Point to Entire Target point Matric  
template <typename T>
cx_3d<T> G_gen(Row<T> source_point_vector,Mat<T> target_point ,T k_s ,T k_p,T rho ,T omega);

// Complete G_generator For Source To Target

template<typename T>
cx_4d_mat<T> G_4d_full(Mat<T> source_point_mat, Mat<T> target_point_mat,T k_s , T k_p,T rho , T omega);

// ------------------------------------------------------------------//
// -----------------G Diff calculation procedures-------------------//
// ----------------------------------------------------------------//

template <typename T>
complex<T> r_diff(T radius_mag,T k,T R_val);

template <typename T>
Mat<T> r_diff_mat_gen(Row<T> affect_point , Row<T> source_point,T k_s , T k_p);
  
template <typename T>
T eoiidi(T i_val ,Row<T> R_vector,T r_mag);

template <typename T>
T eoijdk(T i_val ,T j_val ,Row<T> R_vector,T r_mag);

template <typename T>
T eoiidj(T i_val ,T j_val ,Row<T> R_vector,T r_mag);

template <typename T>
T eoijdi(T i_val ,T j_val,Row<T> R_vector,T r_mag);

template <typename T>
T eo_d_universal(Row<T> i_arr,Row<T> R_vec, T r_mag);


  
Mat<int> Combination_gen();
double G_ijk_helper(rowvec ijk_thing);
cx_double G_p_ijk (rowvec ijk_row_passed, rowvec affect_point, rowvec source_point,double k_s , double k_p,double rho , double omega);
cx_4d G_p_diff_ijk(rowvec source_point_vector,mat observing_plane_cube,double k_s , double k_p,double rho , double omega);
cx_5d G_ijk_full_matrix(mat source_point_mat,mat target_point_mat,double k_s , double k_p,double rho , double omega);
cx_3d strain_cal_for_one_source(cx_4d G_ijk_cube, cx_rowvec strength_cube);
cx_3d strain_cal(cx_5d G_ijk_source,cx_mat src_str);
cx_cube transpose_cube_in_dir_2(cx_cube cube_to_trans);
cx_cube stress_coff_point(cx_cube G_ijk,double lamda,double mu);
cx_5d stress_coff_calc(const cx_5d& G_ijk_source,double lamda,double mu);
cx_3d stress_cal_one_source(cx_4d G_ijk_cube, cx_rowvec strength_cube);
cx_3d stress_cal(const cx_5d &G_ijk_source,cx_mat src_str,double lamda,double mu);
cx_3d stress_from_points(mat source_points, cx_mat src_str, mat target_points,double k_s , double k_p,double rho , double omega , double lamda,double mu);
cx_colvec source_strength_derivation(cx_5d G_ijk_cube_of_sources,cx_3d stress_matrix);
cx_mat get_strength(mat source_point_list,double r_s, rowvec direction_cosine,cx_3d stress_matrix,double k_s,double k_p,double rho , double omega);
cx_mat get_strength_hetro(mat source_point_list,mat passive_source_list ,double r_s, rowvec direction_cosine,cx_3d stress_matrix,double k_s,double k_p,double rho , double omega,double mu,double lamda);

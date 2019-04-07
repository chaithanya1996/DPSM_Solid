#include "dpsm_helpers.hpp"


template <typename TYPE>TYPE r_s_calculator(TYPE freq, TYPE c,int safety_factor);

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
Mat<complex<T>> r_diff_mat_gen(Row<T> affect_point , Row<T> source_point,T k_s , T k_p);
  
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


template <typename T>
T G_ijk_helper(Row<T> ijk_thing);

template <typename T> 
complex<T> G_p_ijk (Row<T> ijk_row_passed, Row<T> affect_point, Row<T> source_point,T k_s , T k_p,T rho , T omega);

template <typename T>
cx_4d<T> G_p_diff_ijk(Row<T> source_point_vector,Mat<T> observing_plane_cube,T k_s , T k_p,T rho , T omega);

template <typename T>
cx_5d<T> G_ijk_full_matrix(Mat<T> source_point_mat,Mat<T> target_point_mat,T k_s , T k_p,T rho ,T omega);

  
// ------------------------------------------------------------------//
// -----------------Strain Calculation Procedur---------------------//
// ----------------------------------------------------------------//


template <typename T>
cx_3d<T> strain_cal_for_one_source(cx_4d<T> G_ijk_cube, Row<complex<T>> strength_cube);

template <typename T>
cx_3d<T> strain_cal(cx_5d<T> G_ijk_source,Mat<complex<T>> src_str);

// ------------------------------------------------------------------//
// -----------------Stress Calculation Procedure---------------------//
// ----------------------------------------------------------------//

template<typename T>
Cube<T> transpose_cube_in_dir_2(Cube<T> cube_to_trans);
  
template <typename T>
Cube<T> stress_coff_point(Cube<T> G_ijk,T lamda,T mu);

template <typename T>
cx_5d<T> stress_coff_calc(const cx_5d<T> G_ijk_source,T lamda,T mu);

template <typename T>
cx_3d<T> stress_cal_one_source(cx_4d<T> G_ijk_cube, Row<complex<T>> strength_cube);

template <typename T>
cx_3d<T> stress_cal(const cx_5d<T> G_ijk_source,Mat<complex<T>> src_str,T lamda,T mu);


template <typename T>
cx_3d<T> stress_from_points(Mat<T> source_points, Mat<complex<T>> src_str, Mat<T> target_points,T k_s , T k_p,T rho ,T omega , T lamda,T mu);

template <typename T>
Col<complex<T>> source_strength_derivation(cx_5d<T> G_ijk_cube_of_sources,cx_3d<T> stress_matrix);


template<typename T>
Mat<complex<T>> get_strength(Mat<T> source_point_list,T r_s,Row<T> direction_cosine,cx_3d<T> stress_matrix,T k_s,T k_p,T rho , T omega);

template <typename T>
Mat<complex<T>> get_strength_hetro(Mat<T> source_point_list,Mat<T> passive_source_list ,T r_s, Row<T> direction_cosine,cx_3d<T> stress_matrix,T k_s,T k_p,T rho , T omega,T mu,T lamda,bool have_active_sources);

template <typename T>
Mat<complex<T>> solve_dpsm_str (Mat<T> ACTIVE_SOURCES_DPSM_POINT, Mat<T> PASSIVE_SOURCES_DPSM_POINT, cx_3d<T> STRESS_CX_3D_MATRIX , Mat<T> Points_of_enforcement, T k_s,T k_p,T rho , T omega,T mu,T lamda);

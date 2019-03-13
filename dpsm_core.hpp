#include "dpsm_helpers.hpp"

template <typename TYPE>TYPE r_s_calculator(TYPE freq, TYPE c){
  int safety_factor = 2;
  TYPE r_s;
  r_s = c / freq / 2 / M_PI / safety_factor;
  return(r_s);
};


double kron_delta(int i, int j);
cx_rowvec G_ij(int i_sub, int j_sub, rowvec affect_point , rowvec source_point,double k_s , double k_p ,double rho , double omega);
cx_cube G_ps_mat_for_point(rowvec affect_point, rowvec source_point,double k_s , double k_p, double rho , double omega);
cx_mat G_mat_for_point(rowvec affect_point , rowvec source_point,double k_s , double k_p,double rho , double omega);
cx_3d G_gen(rowvec source_point_vector,mat target_point ,double k_s , double k_p,double rho , double omega);
cx_4d_mat G_4d_full(mat source_point_mat, mat target_point_mat,double k_s , double k_p,double rho , double omega);
std::complex<double> r_diff(double radius_mag,double k,double R_val);
cx_mat r_diff_mat_gen(rowvec affect_point , rowvec source_point,double k_s , double k_p);
double eoiidi(double i_val ,rowvec R_vector,double r_mag);
double eoijdk(double i_val ,double j_val ,rowvec R_vector,double r_mag);
double eoiidj(double i_val ,double j_val ,rowvec R_vector,double r_mag);
double eoijdi(double i_val ,double j_val,rowvec R_vector,double r_mag);
double eo_d_universal(rowvec i_arr,rowvec R_vec, double r_mag);
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

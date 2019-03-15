#include "dpsm_core.hpp"


template <typename TYPE>
TYPE r_s_calculator(TYPE freq, TYPE c){
  int safety_factor = 2;
  TYPE r_s;
  r_s = c / freq / 2 / M_PI / safety_factor;
  return(r_s);
};

// Explicit Instantiation

template float r_s_calculator(float,float);
template double r_s_calculator(double,double);

// ------------------------------------------------------------------//
// ---------------Green generator Catergory-------------------------//
// ----------------------------------------------------------------//


// Kron Delta _funtion

template<typename T>
T kron_delta(int i, int j){
  if(i == j){
    return(T (1));
  }
  else{
    return(T (0));
  }
}

// Explicit Instantiation

template float kron_delta(int,int);
template double kron_delta(int,int);


// ------------------------------------------------------------------//
// Templated Green Function 

template <typename T>
Row<complex<T>> G_ij(int i_sub, int j_sub, Row<T> affect_point , Row<T> source_point,T k_s, T k_p ,T rho, T omega){
  Row<T> pos_vector  = source_point - affect_point ;
  T radius_mag = sqrt(pow(pos_vector(0),2)+pow(pos_vector(1),2)+pow(pos_vector(2),2));
  
  // Defining The Required parameters
  
  std::complex<T> img_start (0,1.0);
  std::complex<T> e_p = exp(img_start * k_p * radius_mag)/radius_mag; 
  std::complex<T> e_s = exp(img_start * k_s * radius_mag)/radius_mag;
  std::complex<T> r_p = img_start * k_p / radius_mag - 1/pow(radius_mag,2);
  std::complex<T> r_s = img_start * k_s / radius_mag - 1/pow(radius_mag,2);

  
  Row<T> R_vec(3,fill::zeros);
  for(size_t i = 0;i<3;++i){
    R_vec(i)= affect_point(i) - source_point(i);
  }

  // Calculating G_sub

  Row<complex<T>> G_point(2,fill::zeros);
  
  G_point(0) = e_p *( pow(k_p,2) * R_vec(i_sub-1) *  R_vec(j_sub-1) + ( T (3) * R_vec(i_sub-1) *  R_vec(j_sub-1) -  kron_delta<T>(i_sub,j_sub) ) * r_p) / ( 4 * M_PI *  rho * pow(omega,2));
  G_point(1) = e_s *( pow(k_s,2) * ( kron_delta<T>(i_sub,j_sub) -  R_vec(i_sub-1) *  R_vec(j_sub-1) ) - ( T (3) * R_vec(i_sub-1) *  R_vec(j_sub-1) -  kron_delta<T>(i_sub,j_sub) ) * r_s)/ ( 4 * M_PI *  rho * pow(omega,2));

  return(G_point);
}


// Explicit Instantiation

template Row<complex<float>> G_ij(int,int,Row<float>,Row<float>,float,float,float,float);

template Row<complex<double>> G_ij(int,int,Row<double>,Row<double>,double,double,double,double);

// ------------------------------------------------------------------//



// Genrate the Transformation cube for the observing plane for each source point

template <typename T>
Cube<complex<T>> G_ps_mat_for_point(Row<T> affect_point, Row<T> source_point,T k_s , T k_p, T rho , T omega){
  Cube<complex<T>> G_ps_matrix(3,3,2,fill::zeros);
  for(size_t i=0;i<3;++i){
    for(size_t j=0;j<3;++j){
      G_ps_matrix.tube(i,j) = G_ij(i+1,j+1,affect_point,source_point,k_s ,k_p,rho,omega);
    }
  }
  return(G_ps_matrix);
}

// Explicit Instantiation
template  Cube<complex<float>> G_ps_mat_for_point(Row<float>,Row<float>,float,float,float,float);

template Cube<complex<double>> G_ps_mat_for_point(Row<double>,Row<double>,double,double,double,double);


// ------------------------------------------------------------------//


// Second Version for facilitate Displacemnt Caclualtuon
template <typename T>
Mat<complex<T>> G_mat_for_point(Row<T> affect_point , Row<T> source_point,T k_s ,T k_p,T rho ,T omega){
  Mat<complex<T>> G_ps_matrix(3,3,fill::zeros);
  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j){
      Row<complex<T>> place_h = G_ij(i+1,j+1,affect_point,source_point,k_s ,k_p,rho,omega);
      G_ps_matrix(i,j) = place_h(0) + place_h(1);
    }
  }
  return(G_ps_matrix);
}

template Mat<complex<float>> G_mat_for_point(Row<float>,Row<float>,float,float,float ,float);

template Mat<complex<double>> G_mat_for_point(Row<double>, Row<double>,double,double,double,double);
// ------------------------------------------------------------------//


// Cube Generation of G for a plane
template <typename T>
cx_3d<T> G_gen(Row<T> source_point_vector,Mat<T> target_point ,T k_s ,T k_p,T rho ,T omega){ 
  int tar_vec_len = target_point.n_rows;
  cx_3d<T> G_vec_source(tar_vec_len,cx_mat(3,3));
  for(int i =0;i<tar_vec_len;++i){
    G_vec_source[i] = G_mat_for_point<T>(target_point.row(i),source_point_vector,k_s,k_p,rho,omega);
  }
  //std::cout << "ol" << endl ;
  return(G_vec_source);
}

template cx_3d<float> G_gen(Row<float>,Mat<float> ,float ,float,float ,float);

template cx_3d<double> G_gen(Row<double>,Mat<double> ,double ,double,double ,double);


// ------------------------------------------------------------------//


// The Whole matrix generation
template<typename T>
cx_4d_mat<T> G_4d_full(Mat<T> source_point_mat, Mat<T> target_point_mat,T k_s , T k_p,T rho , T omega){
  int no_of_source_point = source_point_mat.n_rows;
  int no_of_target_point = target_point_mat.n_rows;
  cx_4d_mat<T> G_for_all_source_and_target(no_of_source_point,cx_3d(no_of_target_point,cx_mat(3,3,fill::zeros)));
  for (int i = 0; i < G_for_all_source_and_target.size(); ++i) {
    G_for_all_source_and_target[i] = G_gen<T>(source_point_mat.row(i),target_point_mat ,k_s , k_p,rho ,omega);
  }
  return(G_for_all_source_and_target);
}

template cx_4d_mat<float> G_4d_full(Mat<float>, Mat<float>,float, float,float , float);

template cx_4d_mat<double> G_4d_full(Mat<double>, Mat<double>,double, double,double , double);

// ------------------------------------------------------------------//
// -----------------G Diff calculation procedures-------------------//
// ----------------------------------------------------------------//
template <typename T>
complex<T> r_diff(T radius_mag,T k,T R_val){
  complex<T> img_start(0,1.0);
  complex<T> value_of_differential = 2 * R_val / pow(radius_mag,3) - img_start * k * R_val /  pow(radius_mag,2);
  return(value_of_differential);
}

template complex<float> r_diff(float,float,float);
template complex<double> r_diff(double,double,double);

// ----------------------------------------------------------------//
template <typename T>
Mat<T> r_diff_mat_gen(Row<T> affect_point , Row<T> source_point,T k_s , T k_p){
  Row<T> R_vec(3,fill::zeros);
  for(int i = 0;i<3;++i){
    R_vec(i)= affect_point(i) - source_point(i);
  }
  Row<T> r_diff_matrix(3,2,fill::zeros);

  // Now we are going to generate matrix in the same way formuals written in book
  // P - 1st column | S - Second Column
  // The rows are for the directional cosines

  T radius_mag = vec_mag<T>(source_point - affect_point);

  for(int i = 0 ; i < 3 ; ++i){ // i is for rows in matrix
    for(int j = 0 ; j < 2 ; ++j){ // J is for columns in matrix
      if(j == 0){
	r_diff_matrix(i,j) = r_diff<T>(radius_mag, k_p, R_vec(i));
      }else{
	r_diff_matrix(i,j) = r_diff<T>(radius_mag,k_s, R_vec(i));
      } 
    }
  }
  return(r_diff_matrix);
}

template Mat<float> r_diff_mat_gen(Row<float>,Row<float>,float, float);

template Mat<double> r_diff_mat_gen(Row<double>,Row<double>,double, double);

// ---------------------------------------------------------------------------------//
// small subset_functions


template <typename T>
T eoiidi(T i_val ,Row<T> R_vector,T r_mag){
  T return_value =   - 2 * pow(R_vector(i_val),3) / r_mag + 2 * R_vector(i_val) / r_mag;
  return(return_value);
}

template float eoiidi(float ,Row<float>,float);
template double eoiidi(double ,Row<double>,double);

template <typename T>
T eoijdk(T i_val ,T j_val ,Row<T> R_vector,T r_mag){
  T return_value =   -2 * R_vector(0) * R_vector(1) * R_vector(2) / r_mag ;
  return(return_value);
}

template float eoijdk(float ,float ,Row<float>,float);
template double eoijdk(double ,double ,Row<double>,double);
  

template <typename T>
T eoiidj(T i_val ,T j_val ,Row<T> R_vector,T r_mag){
  T return_value =   -2 * R_vector(0) * R_vector(0) * R_vector(2) / r_mag ;
  return(return_value);
}

template float eoiidj(float ,float ,Row<float> ,float);
template double eoiidj(double ,double ,Row<double> ,double);

template <typename T>
T eoijdi(T i_val ,T j_val,Row<T> R_vector,T r_mag){
  T return_value =   - 2 * pow(R_vector(i_val),2) * R_vector(j_val) / r_mag +  R_vector(j_val) / r_mag;
  return(return_value);
}

template float eoijdi(float ,float,Row<float>,float);
template double eoijdi(double ,double,Row<double>,double);

template <typename T>
T eo_d_universal(Row<T> i_arr,Row<T> R_vec, T r_mag){
  T eo_val = 0;
  if(i_arr(0) == i_arr(1)){
    if(i_arr(1) == i_arr(2)){
      return(eoiidi<T>(i_arr(0),R_vec,r_mag));
    }
    else{
      return(eoiidj<T>(i_arr(0),i_arr(2),R_vec,r_mag));
    }
  }else{
    if(i_arr(1) == i_arr(2)){
      return(eoijdi<T>(i_arr(0),i_arr(2),R_vec,r_mag));
    }
    else{
      return(eoijdk<T>(i_arr(0),i_arr(2),R_vec,r_mag));
    }
  }
}

template float eo_d_universal(Row<float> ,Row<float> , float);
template double eo_d_universal(Row<double> ,Row<double> , double);

// ********************************************************************************************************************************//
// Keeping in mind that the et is 3 times eo simply and we can use the above fucntions rather than writing the new fucntionsfor the
// Calculations


Mat<int> Combination_gen(){
  int n_comb_rows = 27;
  Mat<int> all_poss_comb(27,3,fill::zeros);
  for(int i = 0; i < 3;++i){
    for(int j = 0;j < 3;++j){
      for(int k = 0;k< 3;++k){
	      Row<int> generate_helper = {i,j,k};
        std::cout << i*9+j*3+k << endl;
	      all_poss_comb.row(i*9+j*3+k) = generate_helper;
      }
    }
  }
  return(all_poss_comb);
}


// ------------------------------------------------------------------//
// ---------------G p point calculation procedures------------------//
// ----------------------------------------------------------------//
template <typename T>
T G_ijk_helper(Row<T> ijk_thing){
  if(ijk_thing[0] == ijk_thing[1]){
    return -1;
  }
  else{
    return 0;
  }
}

template T G_ijk_helper(Row<T> ijk_thing);
template T G_ijk_helper(Row<T> ijk_thing);
  


cx_double G_p_ijk (rowvec ijk_row_passed, rowvec affect_point, rowvec source_point,double k_s , double k_p,double rho , double omega){

  cx_mat R_differ_mat = r_diff_mat_gen(affect_point , source_point, k_s , k_p);
  
  // Now we are going to generate matrix in the same way formuals written in book
  // P - 1st column | S - Second Column
  // The rows are for the directional cosines
     
  // Repetitive openations of calculations
  
  rowvec pos_vector  = source_point - affect_point ;
  double radius_mag = pow(pos_vector(0),2)+pow(pos_vector(1),2)+pow(pos_vector(2),2);
  
  // Defining The Required parameters
  
  std::complex<double> img_start (0,1.0);
  std::complex<double> e_p = exp(img_start * k_p * radius_mag)/radius_mag; 
  std::complex<double> e_s = exp(img_start * k_s * radius_mag)/radius_mag;
  std::complex<double> r_p = img_start * k_p / radius_mag - 1/pow(radius_mag,2);
  std::complex<double> r_s = img_start * k_s / radius_mag - 1/pow(radius_mag,2);
  rowvec R_vec(3,fill::zeros);
  
  for(int i = 0;i<3;++i){
    R_vec(i)= affect_point(i) - source_point(i);
  }
  rowvec ijk_row =  ijk_row_passed + 1; // For compatibility with the Lower level implementations in the Helper functions 

  // cx_mat G_p_s_matrix = G_ps_mat_for_point(affect_point,source_point,k_s , k_p);
  //std::cout << "iam here" << endl;  
  cx_rowvec G_place_holder =  G_ij(ijk_row[0],ijk_row[1],affect_point,source_point,k_s ,k_p,rho ,omega) ;
  
  ijk_row =  ijk_row_passed ;
  cx_double Gpiidi  = e_p * ( pow(k_p,2) *  eo_d_universal(ijk_row,R_vec,radius_mag) +  ( G_ijk_helper(ijk_row) +3 * R_vec[ijk_row[0]] * R_vec[ijk_row[1]]) * R_differ_mat(ijk_row[2],0) +  r_p * double (3) * eo_d_universal(ijk_row,R_vec,radius_mag)) / ( 4 * M_PI *  rho * pow(omega,2)) + G_place_holder[0] * (img_start * R_vec[ijk_row[2]] * e_p - R_vec[ijk_row[2]] * e_p / radius_mag);			 

  cx_double Gsiidi  = e_s * ( -pow(k_s,2) *  eo_d_universal(ijk_row,R_vec,radius_mag) -  ( G_ijk_helper(ijk_row) +3 * R_vec[ijk_row[0]] * R_vec[ijk_row[1]]) * R_differ_mat(ijk_row[2],0) +  r_s * double (3) * eo_d_universal(ijk_row,R_vec,radius_mag)) / ( 4 * M_PI *  rho * pow(omega,2))+ G_place_holder[1] * (img_start * R_vec[ijk_row[2]] * e_s - R_vec[ijk_row[2]] * e_s / radius_mag);			 
  
  return Gpiidi + Gsiidi ;
}




// ------------------------------------------------------------------//
// ---------------G_ijk Matrix Computation--------------------------//
// ----------------------------------------------------------------//


cx_4d G_p_diff_ijk(rowvec source_point_vector,mat observing_plane_cube,double k_s , double k_p,double rho , double omega){
  int obs_pln_rows = observing_plane_cube.n_rows;
  
  cx_4d G_diff_4dcube(obs_pln_rows,cx_cube(3,3,3,fill::zeros));

  for (int i = 0; i < observing_plane_cube.n_rows; ++i){
      // For Generating the Differential matrices
      for (double k = 0; k < 3; ++k){
        for (double l = 0; l < 3; ++l){
          for (double m = 0; m < 3; ++m){
            //std::cout << "----------degbug begin--------------" << endl ;
            //std::cout << "i,j" << i << j << endl ;
            rowvec ijk_row_gen_vector = {k,l,m};
            //std::cout << "k,l,m" << ijk_row_gen_vector << endl ;
            cx_double G_p_i_k_point_value = G_p_ijk(ijk_row_gen_vector,source_point_vector,observing_plane_cube.row(i),k_s,k_p,rho ,omega);
            //std::cout << " G_value "<< G_p_i_k_point_value << endl ;
            G_diff_4dcube[i](k,l,m) =   G_p_i_k_point_value;
            //std::cout << "----------degbug end--------------" << endl ;
          }
        }
      }
      // End the Differenetial Loops
      //std::cout << G_diff_5dcube[i][j] << endl ;
    }
  return G_diff_4dcube;
}

cx_5d G_ijk_full_matrix(mat source_point_mat,mat target_point_mat,double k_s , double k_p,double rho , double omega){
 double no_of_source_point = source_point_mat.n_rows;
 double no_of_target_point = target_point_mat.n_rows;
 cx_5d  G_for_all_source_and_target(no_of_source_point,cx_4d(no_of_target_point,cx_cube(3,3,3,fill::zeros)));
 #pragma omp parallel for
  for (int i = 0; i < G_for_all_source_and_target.size(); ++i) {
    G_for_all_source_and_target[i] = G_p_diff_ijk(source_point_mat.row(i),target_point_mat ,k_s , k_p,rho ,omega);
  }
  return(G_for_all_source_and_target);  
}

// ------------------------------------------------------------------//
// -----------------Strain Calculation Procedur---------------------//
// ----------------------------------------------------------------//

cx_3d strain_cal_for_one_source(cx_4d G_ijk_cube, cx_rowvec strength_cube){

  cx_colvec strength_cube_col = strength_cube.t();
  int n_row_cube = G_ijk_cube.size();

  cx_3d obs_plane_stress_comp(n_row_cube,cx_mat(3,3,fill::zeros));
  
  /*
  for(size_t i = 0; i < n_row_cube; i++){
    obs_plane_stress_comp[i] = cx_mat(3,3,fill::zeros);
  }
  */
  
  for(size_t i = 0; i < n_row_cube; i++){
    for(size_t k = 0; k < 3; k++){
      cx_cube mat_to_operate = G_ijk_cube[i].subcube( 0, k, 0, 2, k, 2 );
      cx_mat converted_cub_tomat = conv_cube_2_mat(mat_to_operate);
      //  obs_plane_stress_comp[i][j].col(k)= (converted_cub_tomat  + converted_cub_tomat.t()) * strength_cube_col ;
      obs_plane_stress_comp[i]= obs_plane_stress_comp[i] + (converted_cub_tomat  + converted_cub_tomat.t()) * strength_cube_col(k) ;
    }
  }
  return(obs_plane_stress_comp);
}

cx_3d strain_cal(cx_5d G_ijk_source,cx_mat src_str){
  double target_size = G_ijk_source.size();
  cx_3d strains_on_plane(target_size,cx_mat(3,3,fill::zeros));
  #pragma omp parallel for
  for(size_t i = 0; i < src_str.n_cols; i++){
    strains_on_plane = add_cx_3d(strains_on_plane,strain_cal_for_one_source(G_ijk_source[i],src_str.row(i)));
  }
  return(strains_on_plane);
}


// ------------------------------------------------------------------//
// -----------------Stress Calculation Procedure---------------------//
// ----------------------------------------------------------------//

cx_cube transpose_cube_in_dir_2(cx_cube cube_to_trans){
  cx_cube return_cx_cube;
  for (int i=0; i < 0; ++i) {
    return_cx_cube.slice(i) = cube_to_trans.slice(i).t();
  }
  return(return_cx_cube);
}


cx_cube stress_coff_point(cx_cube G_ijk,double lamda,double mu){

  cx_cube Stress_coff(3,3,3,fill::zeros);

  Stress_coff(0,0,0) = (2*mu+lamda)* G_ijk(0,0,0) + lamda * ( G_ijk(1,0,1) + G_ijk(2,0,2)) ;
  Stress_coff(0,0,1) = (2*mu+lamda)* G_ijk(0,1,0) + lamda * ( G_ijk(1,1,1) + G_ijk(2,1,2)) ;
  Stress_coff(0,0,2) = (2*mu+lamda)* G_ijk(0,2,0) + lamda * ( G_ijk(1,2,1) + G_ijk(2,2,2)) ;

  Stress_coff(1,1,0) = (2*mu+lamda)* G_ijk(1,0,1) + lamda * ( G_ijk(0,0,0) + G_ijk(2,0,2)) ;
  Stress_coff(1,1,1) = (2*mu+lamda)* G_ijk(1,1,1) + lamda * ( G_ijk(0,1,0) + G_ijk(2,1,2)) ;
  Stress_coff(1,1,2) = (2*mu+lamda)* G_ijk(1,2,1) + lamda * ( G_ijk(0,2,0) + G_ijk(2,2,2)) ; 

  Stress_coff(2,2,0) = (2*mu+lamda)* G_ijk(2,0,2) + lamda * ( G_ijk(0,0,0) + G_ijk(1,0,1)) ;
  Stress_coff(2,2,1) = (2*mu+lamda)* G_ijk(2,1,2) + lamda * ( G_ijk(0,1,0) + G_ijk(1,1,1)) ;
  Stress_coff(2,2,2) = (2*mu+lamda)* G_ijk(2,2,2) + lamda * ( G_ijk(0,2,0) + G_ijk(1,2,1)) ; 


  Stress_coff(2,0,0) = (mu)* G_ijk(2,0,2) + G_ijk(0,0,2) ;
  Stress_coff(2,0,1) = (mu)* G_ijk(2,1,2) + G_ijk(0,1,2) ;
  Stress_coff(2,0,2) = (mu)* G_ijk(2,2,2) + G_ijk(0,2,2) ;
  
  Stress_coff(0,1,0) = (mu)* G_ijk(2,0,2) + G_ijk(1,0,0) ;
  Stress_coff(0,1,1) = (mu)* G_ijk(2,1,2) + G_ijk(1,1,0) ;
  Stress_coff(0,1,2) = (mu)* G_ijk(2,2,2) + G_ijk(1,2,0) ;

  Stress_coff(1,2,0) = (mu)* G_ijk(1,0,2) + G_ijk(2,0,1) ;
  Stress_coff(1,2,1) = (mu)* G_ijk(1,1,2) + G_ijk(1,1,1) ;
  Stress_coff(1,2,2) = (mu)* G_ijk(1,2,2) + G_ijk(2,2,1) ;


  
  Stress_coff(0,2,0) = (mu)* G_ijk(2,0,2) + G_ijk(0,0,2) ;
  Stress_coff(0,2,1) = (mu)* G_ijk(2,1,2) + G_ijk(0,1,2) ;
  Stress_coff(0,2,2) = (mu)* G_ijk(2,2,2) + G_ijk(0,2,2) ;
  
  Stress_coff(1,0,0) = (mu)* G_ijk(2,0,2) + G_ijk(1,0,0) ;
  Stress_coff(1,0,1) = (mu)* G_ijk(2,1,2) + G_ijk(1,1,0) ;
  Stress_coff(1,0,2) = (mu)* G_ijk(2,2,2) + G_ijk(1,2,0) ;

  Stress_coff(2,1,0) = (mu)* G_ijk(1,0,2) + G_ijk(2,0,1) ;
  Stress_coff(2,1,1) = (mu)* G_ijk(1,1,2) + G_ijk(1,1,1) ;
  Stress_coff(2,1,2) = (mu)* G_ijk(1,2,2) + G_ijk(2,2,1) ;


  return(Stress_coff);
}

cx_5d stress_coff_calc(const cx_5d& G_ijk_source,double lamda,double mu){
  double size_1 = G_ijk_source.size();
  double size_2 = G_ijk_source[0].size();
  cx_5d Converterd_stress_coff(size_1,cx_4d(size_2,cx_cube(3,3,3,fill::zeros)));
  
  for (int i = 0; i < size_1; ++i) {
    for (int j = 0; j < size_2; ++j) {
      Converterd_stress_coff[i][j] = stress_coff_point(G_ijk_source[i][j],lamda,mu);
    }
  }
  return(Converterd_stress_coff);
}


cx_3d stress_cal_one_source(cx_4d G_ijk_cube, cx_rowvec strength_cube){

  cx_colvec strength_cube_col = strength_cube.t();
  int n_row_cube = G_ijk_cube.size();

  cx_3d obs_plane_stress_comp(n_row_cube,cx_mat(3,3,fill::zeros));
  
  
  for(size_t i = 0; i < n_row_cube; i++){
    for(size_t k = 0; k < 3; k++){
      cx_cube mat_to_operate = G_ijk_cube[i].subcube( 0, k, 0, 2, k, 2 );
      cx_mat converted_cub_tomat = conv_cube_2_mat(mat_to_operate);
      obs_plane_stress_comp[i]= obs_plane_stress_comp[i] + (converted_cub_tomat) * strength_cube_col(k) ;
    }
  }
  return(obs_plane_stress_comp);
}



cx_3d stress_cal(const cx_5d &G_ijk_source,cx_mat src_str,double lamda,double mu){ // the G_ijk_must be Converted Soueces here

  double target_size = G_ijk_source[0].size(); // Thought as Correct 
  cx_3d stress_on_plane(target_size,cx_mat(3,3,fill::zeros));
  
  // MultiThreading the Calculations Using OPenMP
  #pragma omp parallel for
  for(size_t i = 0; i < src_str.n_cols; i++){
    stress_on_plane = add_cx_3d(stress_on_plane,stress_cal_one_source(G_ijk_source[i],src_str.row(i)));
  }
  return(stress_on_plane);
}


// The Entire Fucntion is Redundant Waste Memory And thereof Needed Care ful Restructuing of Code To Avoid Memory Problem in Furture
// This Paticular Function is used with constat Pointer Reference to Minimize the Memory usage Peak During the Computation
cx_3d stress_from_points(mat source_points, cx_mat src_str, mat target_points,double k_s , double k_p,double rho , double omega , double lamda,double mu){
  cx_5d G_vals_for_given = G_ijk_full_matrix(source_points,target_points,k_s ,k_p,rho , omega); 
  cx_5d stress_for_given  =  stress_coff_calc(G_vals_for_given,lamda, mu); // Pass by Reference
  cx_3d str_vals = stress_cal(stress_for_given,src_str,lamda,mu); // Pass By Reference 
  return(str_vals);
}

/* Redundant in Occasion of Convering all observing plane to Vectors 


std::vector<cx_double> conv_observing_plane_mat_to_vector(cx_mat Observing_plane){
  int no_of_rows = Observing_plane.n_rows;
  int no_of_columns = Observing_plane.n_cols;

  std::vector<cx_double> in_vecc_format(no_of_rows * no_of_columns);
  
  for(int i = 0; i < no_of_rows; ++i){
    for( int j = 0; j < no_of_columns; ++j){
      in_vecc_format[i*no_of_columns+j] = Observing_plane[i,j];
    }
  }
  return(in_vecc_format);
}

*/


// Stress Vector Calculation 



// cx_cube stess_cal(cx_cube G_ijk_cube, cx_rowvec strength_cube){
//   int n_row_cube = G_ijk_cube.n_rows;
//   int n_col_cube = G_ijk_cube.n_cols;
//   cx_cube obs_plane_stress_comp(n_row_cube,n_col_cube,3,fill::zeros);
  
//   for(size_t i = 0; i < n_row_cube; i++){
//     for(size_t j = 0; i < n_col_cube; j++){

//       cx_rowvec temp_var1 = G_ijk_cube.tube(i,j);
//       obs_plane_stress_comp(i,j,0) = sum((temp_var1(span(0,2))  + temp_var1(span(8,10)) + temp_var1(span(18,20)))* strength_cube(0));
//       obs_plane_stress_comp(i,j,1) = sum((temp_var1(span(3,5))  + temp_var1(span(8,10)) + temp_var1(span(18,20)))* strength_cube(0));
//       obs_plane_stress_comp(i,j,2) = sum((temp_var1(span(6,8))  + temp_var1(span(8,10)) + temp_var1(span(18,20)))* strength_cube(0));
//       }
//     }
//   return(obs_plane_stress_comp);
// }


// ------------------------------------------------------------------//
// -----------------Point Strength Calculation ---------------------//
// ----------------------------------------------------------------//



cx_colvec source_strength_derivation(cx_5d G_ijk_cube_of_sources,cx_3d stress_matrix){
  
  int sources_len = stress_matrix.size();
  cout << " Source Lenngth -" << sources_len << endl;
  // cx_mat str_mat(sources_len*stress_matrix[0].n_rows,stress_matrix[0].n_cols,fill::zeros);
  cx_colvec str_mat(sources_len * stress_matrix[0].n_elem);

  
  int source_no = G_ijk_cube_of_sources.size() ;
  int target_no =  G_ijk_cube_of_sources[0].size();

  
  cx_cube G_full(target_no * 3,sources_len * 3,3,fill::zeros);

  //std::cout << "Creating G_full" << endl;
  // std::cout << G_full << endl;

#pragma omp parallel for
  for(size_t i = 0; i < sources_len; i++){
    for(size_t j =0; j < 3 ; ++j){
      str_mat.subvec(i*9+j*3,i*9+j*3+2) = stress_matrix[i].col(j);
    }
  }

  //std::cout << "Creating srmat" << endl;
  //std::cout << str_mat << endl;

#pragma omp parallel for
  for(size_t j = 0; j < target_no; j++){ 
    for(size_t i = 0; i < sources_len; i++){
      G_full(span(j*3,(j*3)+2),span(i*3,i*3+2),span(0,2)) =  G_ijk_cube_of_sources[i][j];
    }
  }
  
  //std::cout << "Calculated  G_full" << endl;
  //std::cout << G_full << endl;
  
  // Converting Above Cube tot the matrix

  cx_mat G_full_mat_form(G_full.n_rows * 3,G_full.n_cols,fill::zeros);

#pragma omp parallel for
  for(size_t i = 0; i < G_full.n_rows /3; i++){
    for(size_t j =0; j < 3 ; ++j){
      cx_mat temper = G_full(span(0,2),span(0,G_full.n_cols-1),span(j,j));
      //      std::cout << "adjsnfdb" <<i << j << endl ;
      G_full_mat_form(span(i*9+j*3,i*9+j*3+2), span(0,G_full.n_cols-1) ) = temper ;
    }
  }

  //std::cout << G_full_mat_form <<endl;
  cout << " Start_solving " << sources_len << endl;
  cx_colvec col_one = solve(G_full_mat_form,str_mat);
  cout << " End Solving -" << sources_len << endl;
  return(col_one);
}

// ------------------------------------------------------------------//
// High Level Function For point Strenght caculation----------------//
// ----------------------------------------------------------------//


// Need to Redesign The Hight level function for Hetro geneous and homogeneous material


cx_mat get_strength(mat source_point_list,double r_s, rowvec direction_cosine,cx_3d stress_matrix,double k_s,double k_p,double rho , double omega){
  
  mat target_point(size(source_point_list));
  
  for(int i = 0 ; i < source_point_list.n_rows;++i){
    target_point.row(i) = source_point_list.row(i) - r_s * direction_cosine;
  }
  
  double no_of_source_points = source_point_list.n_rows;
  cx_5d G_matrix_calculations(no_of_source_points,cx_4d(no_of_source_points,cx_cube(3,3,3,fill::zeros)));
  
  for(int i = 0; i < no_of_source_points; ++i){
    G_matrix_calculations[i] = G_p_diff_ijk(source_point_list.row(i), target_point,k_s,k_p,rho ,omega);
  }

  
  cx_colvec point_strenghts = source_strength_derivation(G_matrix_calculations,stress_matrix);
  
  
  double out_mat_rows =  point_strenghts.n_elem/3;
  cx_mat P_mat_form(out_mat_rows,3,fill::zeros);
  
  for(int i = 0;i<out_mat_rows;++i){
    for( int j = 0; j < 3 ;++j){
      P_mat_form(i,j) = point_strenghts(i*3+j);
    }
  }
  return(P_mat_form);
}

cx_mat get_strength_hetro(mat source_point_list,mat passive_source_list ,double r_s, rowvec direction_cosine,cx_3d stress_matrix,double k_s,double k_p,double rho , double omega,double mu,double lamda){

  // incoming Point Locations are Exact Location on which the Sources are placed
  
  int no_of_active_sources = source_point_list.n_rows;
  int no_of_passive_sources = passive_source_list.n_rows;

  cx_mat active_source_str(no_of_active_sources,3,fill::zeros);
  cx_mat passive_source_str(no_of_passive_sources,3,fill::zeros);

  mat total_sources = join_cols(source_point_list,passive_source_list);
  
  mat exact_source_point_mat_stress(no_of_active_sources,3);
  cout << "Black Sheep 7" << endl;
#pragma omp parallel for 
  for (int i = 0; i < no_of_active_sources; ++i) {
    exact_source_point_mat_stress.row(i) = source_point_list.row(i)  + direction_cosine * r_s ;
  }
    cout << "Black Sheep 8" << endl;
  cx_5d G_matrix_calculations(total_sources.n_rows,cx_4d(no_of_active_sources ,cx_cube(3,3,3,fill::zeros)));
#pragma omp parallel for 
  for(int i = 0; i < total_sources.n_rows; ++i){
    //std::cout << i << "\n";
    G_matrix_calculations[i] = G_p_diff_ijk(total_sources.row(i),exact_source_point_mat_stress,k_s,k_p,rho ,omega);
  }
  cout << "Black Sheep 9" << endl;
  cx_colvec solid_point_strength_colvec = source_strength_derivation(G_matrix_calculations,stress_matrix);
  cout << "Black Sheep 10" << endl;
  double out_mat_rows =  solid_point_strength_colvec.n_elem/3;
  cx_mat P_mat_form(out_mat_rows,3,fill::zeros);
  
  for(int i = 0;i<out_mat_rows;++i){
    for( int j = 0; j < 3 ;++j){
      P_mat_form(i,j) = solid_point_strength_colvec(i*3+j);
    }
  }
  
  return(P_mat_form);
}

// Initializing The Template Class


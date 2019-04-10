#include "dpsm_core.hpp"


template <typename TYPE>
TYPE r_s_calculator(TYPE freq, TYPE c,int safety_factor){
  TYPE r_s;
  r_s = c / freq / 2 / M_PI / safety_factor;
  return(r_s);
};

// Explicit Instantiation

template float r_s_calculator(float,float,int safety_factor);
template double r_s_calculator(double,double,int safety_factor);

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
  T radius_mag = T(sqrt(pow<T>(pos_vector(0),2)+pow<T>(pos_vector(1),2)+pow<T>(pos_vector(2),2)));
  
  // Defining The Required parameters
  
  complex<T> img_start(0,1.0);
  complex<T> e_p = exp(img_start * k_p * radius_mag)/radius_mag; 
  complex<T> e_s = exp(img_start * k_s * radius_mag)/radius_mag;
  complex<T> r_p = img_start * k_p / radius_mag - 1/T(pow<T>(radius_mag,2));
  complex<T> r_s = img_start * k_s / radius_mag - 1/T(pow<T>(radius_mag,2));

  //cout << "e_p--" << e_p << " e_s--" << e_s << " r_p--" << r_p << " r_s--" << r_s << endl;

  Row<T> R_vec(3,fill::zeros);
  for(size_t i = 0;i<3;++i){
    R_vec(i)= (affect_point(i) - source_point(i))/radius_mag;
  }

  //cout << "R_vec--" << R_vec << endl;
  
  // Calculating G_sub

  Row<complex<T>> G_point(2,fill::zeros);

  T multiplier =  4 * M_PI * rho * omega * omega;
  
  G_point(0) = e_p *( T(pow<T>(k_p,2)) * R_vec(i_sub-1) *  R_vec(j_sub-1) + ( T (3) * R_vec(i_sub-1) *  R_vec(j_sub-1) -  kron_delta<T>(i_sub,j_sub) ) * r_p);
  
  G_point(1) =e_s *( T(pow<T>(k_s,2)) * ( kron_delta<T>(i_sub,j_sub) -  R_vec(i_sub-1) *  R_vec(j_sub-1) ) - ( T (3) * R_vec(i_sub-1) *  R_vec(j_sub-1) -  kron_delta<T>(i_sub,j_sub)) * r_s);

  G_point = G_point / multiplier;

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
      G_ps_matrix.tube(i,j) = G_ij<T>(i+1,j+1,affect_point,source_point,k_s ,k_p,rho,omega);
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
      Row<complex<T>> place_h = G_ij<T>(i+1,j+1,affect_point,source_point,k_s ,k_p,rho,omega);
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
  cx_3d<T> G_vec_source(tar_vec_len,Mat<complex<T>>(3,3));
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
  cx_3d<T> Dummy_cx_3d(no_of_target_point,Mat<complex<T>>(3,3,fill::zeros));
  cx_4d_mat<T> G_for_all_source_and_target(no_of_source_point,Dummy_cx_3d);
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
  complex<T> value_of_differential =complex<T>(2 * R_val / T(pow<T>(radius_mag,3)) - img_start * k * R_val / T(pow<T>(radius_mag,2)));
  return(value_of_differential);
}

template complex<float> r_diff(float,float,float);
template complex<double> r_diff(double,double,double);

// ----------------------------------------------------------------//



template <typename T>
Mat<complex<T>> r_diff_mat_gen(Row<T> affect_point , Row<T> source_point,T k_s , T k_p){
  T radius_mag = vec_mag<T>(source_point - affect_point);
  Row<T> R_vec(3,fill::zeros);
  
  for(int i = 0;i<3;++i){
    R_vec(i)= (affect_point(i) - source_point(i))/radius_mag;
  }
  Mat<complex<T>> r_diff_matrix(3,2,fill::zeros);
  // Now we are going to generate matrix in the same way formuals written in book
  // P - 1st column | S - Second Column
  // The rows are for the directional cosines

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

template Mat<complex<float>> r_diff_mat_gen(Row<float>,Row<float>,float, float);
template Mat<complex<double>> r_diff_mat_gen(Row<double>,Row<double>,double, double);

// ---------------------------------------------------------------------------------//
// small subset_functions


template <typename T>
T eoiidi(T i_val ,Row<T> R_vector,T r_mag){
  T return_value =   - 2 * pow<T>(R_vector(i_val),3) / r_mag + 2 * R_vector(i_val) / r_mag;
  //cout << "Triggered -- iii" << endl; 
  return(return_value);
}

template float eoiidi(float ,Row<float>,float);
template double eoiidi(double ,Row<double>,double);

template <typename T>
T eoijdk(T i_val ,T j_val ,Row<T> R_vector,T r_mag){
  T return_value =   -2 * R_vector(0) * R_vector(1) * R_vector(2) / r_mag ;
  //cout << "Triggered -- ijk" << endl;
  return(return_value);
}

template float eoijdk(float ,float ,Row<float>,float);
template double eoijdk(double ,double ,Row<double>,double);
  

template <typename T>
T eoiidj(T i_val ,T j_val ,Row<T> R_vector,T r_mag){
  T return_value =   -2 * R_vector(i_val) * R_vector(i_val) * R_vector(j_val) / r_mag ;
  //cout << "Triggered -- iij -- i_val" << i_val << "  j_val--" << j_val<< endl;
  return(return_value);
}

template float eoiidj(float ,float ,Row<float> ,float);
template double eoiidj(double ,double ,Row<double> ,double);

template <typename T>
T eoijdi(T i_val ,T j_val,Row<T> R_vector,T r_mag){
  T return_value =   - 2 * pow<T>(R_vector(i_val),2) * R_vector(j_val) / r_mag +  R_vector(j_val) / r_mag;
  //cout << "Triggered -- iji -- i_val" << i_val << "  j_val--" << j_val  << endl;
  return(return_value);
}

template float eoijdi(float ,float,Row<float>,float);
template double eoijdi(double ,double,Row<double>,double);

template <typename T>
T eo_d_universal(Row<T> i_arr,Row<T> R_vec, T r_mag){
  T eo_val = 0;
  if(i_arr(0) == i_arr(1)){
    if(i_arr(1) == i_arr(2)){
      //cout << "Triggering " << i_arr << endl;
      return(eoiidi<T>(i_arr(0),R_vec,r_mag));
    }
    else{
      //cout << "Triggering " << i_arr << endl;
      return(eoiidj<T>(i_arr(0),i_arr(2),R_vec,r_mag));
    }
  }else{
    if(i_arr(0) == i_arr(2)){
      //cout << "Triggering " << i_arr << endl;
      return(eoijdi<T>(i_arr(0),i_arr(1),R_vec,r_mag));
    }
    if(i_arr(1) == i_arr(2)){
      //cout << "Triggering " << i_arr << endl;
      return(eoijdi<T>(i_arr(1),i_arr(0),R_vec,r_mag));
    }
    else{
      //cout << "Triggering " << i_arr << endl;
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

template float G_ijk_helper(Row<float>);
template double G_ijk_helper(Row<double>);
  



template <typename T> 
complex<T> G_p_ijk (Row<T> ijk_row_passed, Row<T> affect_point_this, Row<T> source_point_this,T k_s , T k_p,T rho , T omega){
  // template Mat<complex<float>> r_diff_mat_gen(Row<float>,Row<float>,float, float);
  Mat<complex<T>> R_differ_mat = r_diff_mat_gen<T>(affect_point_this,source_point_this, k_s,k_p);
  // Now we are going to generate matrix in the same way formuals written in book
  // P - 1st column | S - Second Column
  // The rows are for the directional cosines
     
  // Repetitive openations of calculations
  
  Row<T> pos_vector  = source_point_this - affect_point_this ;
  T radius_mag = pow<T>(pos_vector(0),2)+pow<T>(pos_vector(1),2)+pow<T>(pos_vector(2),2);
  
  // Defining The Required parameters
  
  complex<T> img_start (0,1.0);
  complex<T> e_p = exp(img_start * k_p * radius_mag)/radius_mag; 
  complex<T> e_s = exp(img_start * k_s * radius_mag)/radius_mag;
  complex<T> r_p = img_start * k_p / radius_mag - 1/T(pow<T>(radius_mag,2));
  complex<T> r_s = img_start * k_s / radius_mag - 1/T(pow<T>(radius_mag,2));
  Row<T> R_vec(3,fill::zeros);
  
  for(size_t i = 0;i<3;++i){
    R_vec(i)= (affect_point_this(i) - source_point_this(i)) / radius_mag;
  }
  Row<T> ijk_row =  ijk_row_passed + 1; // For compatibility with the Lower level implementations in the Helper functions 

  // cx_mat G_p_s_matrix = G_ps_mat_for_point(affect_point,source_point,k_s , k_p);
  //std::cout << "iam here" << endl;  
  Row<complex<T>> G_place_holder =  G_ij<T>(ijk_row[0],ijk_row[1],affect_point_this,source_point_this,k_s ,k_p,rho ,omega);
  
  ijk_row =  ijk_row_passed ;
  T multiplier = 4 * M_PI *  rho * pow<T>(omega,2); 
  complex<T> Gpiidi  = (e_p * ( pow<T>(k_p,2) *  eo_d_universal(ijk_row,R_vec,radius_mag) +
   ( G_ijk_helper(ijk_row) + 3 * R_vec[ijk_row[0]] * R_vec[ijk_row[1]]) * R_differ_mat(ijk_row[2],0) +  r_p * T (3) * eo_d_universal(ijk_row,R_vec,radius_mag)))/ multiplier +
    G_place_holder[0] * (img_start * k_p * R_vec[ijk_row[2]] * e_p - R_vec[ijk_row[2]] * e_p / radius_mag);			 

  complex<T> Gsiidi  = (e_s * ( -pow<T>(k_s,2) *  eo_d_universal(ijk_row,R_vec,radius_mag) -  ( G_ijk_helper(ijk_row) +3 * R_vec[ijk_row[0]] * R_vec[ijk_row[1]]) * R_differ_mat(ijk_row[2],0) +  r_s * T (3) * eo_d_universal(ijk_row,R_vec,radius_mag))) / multiplier +
    G_place_holder[1] * (img_start * k_s * R_vec[ijk_row[2]] * e_s - R_vec[ijk_row[2]] * e_s / radius_mag);			 
  
  return Gpiidi + Gsiidi ;
}

template complex<float> G_p_ijk (Row<float>, Row<float>, Row<float>,float, float,float, float);
template complex<double> G_p_ijk (Row<double>, Row<double>, Row<double>,double, double,double, double);


// ------------------------------------------------------------------//
// ---------------G_ijk Matrix Computation--------------------------//
// ----------------------------------------------------------------//

// template <typename T>
// cx_4d<T> G_p_diff_ijk(Row<T> source_point_vector,Mat<T> observing_plane_cube,T k_s , T k_p,T rho , T omega){
//   int obs_pln_rows = observing_plane_cube.n_rows;
  
//   cx_4d<T> G_diff_4dcube(obs_pln_rows,Cube<complex<T>>(3,3,3,fill::zeros));

//   for (int i = 0; i < observing_plane_cube.n_rows; ++i){
//       // For Generating the Differential matrices
//       for (int k = 0; k < 3; ++k){
//         for (int l = 0; l < 3; ++l){
//           for (int m = 0; m < 3; ++m){
//             //std::cout << "----------degbug begin--------------" << endl ;
//             //std::cout << "i,j" << i << j << endl ;
//             Row<T> ijk_row_gen_vector = {T(k),T(l),T(m)}; // Type Conversion into Appropiate Type for Templating the Operation
//             //std::cout << "k,l,m" << ijk_row_gen_vector << endl ;
//             complex<T> G_p_i_k_point_value = G_p_ijk<T>(ijk_row_gen_vector,source_point_vector,observing_plane_cube.row(i),k_s,k_p,rho ,omega);
//             //std::cout << " G_value "<< G_p_i_k_point_value << endl ;
//             G_diff_4dcube[i](k,l,m) =   G_p_i_k_point_value;
//             //std::cout << "----------degbug end--------------" << endl ;
//           }
//         }
//       }
//       // End the Differenetial Loops
//       //std::cout << G_diff_5dcube[i][j] << endl ;
//     }
//   return G_diff_4dcube;
// }

// template cx_4d<float> G_p_diff_ijk(Row<float>,Mat<float>,float,float,float, float);
// template cx_4d<double> G_p_diff_ijk(Row<double>,Mat<double>,double,double,double, double);


template <typename T>
Cube<complex<T>> G_diff_ijk_point(Row<T> source_point_row,Row<T> target_point_row,T k_s , T k_p,T rho , T omega){
  
  Cube<complex<T>> G_diff(3,3,3,fill::zeros);
      // For Generating the Differential matrices
      for (int k = 0; k < 3; ++k){
        for (int l = 0; l < 3; ++l){
          for (int m = 0; m < 3; ++m){
            Row<T> ijk_row_gen_vector = {T(k),T(l),T(m)}; // Type Conversion into Appropiate Type for Templating the Operation
            G_diff(k,l,m) =  G_p_ijk<T>(ijk_row_gen_vector, source_point_row,target_point_row,k_s,k_p,rho ,omega);
          }
        }
      }
      // End the Differenetial Loops
      //std::cout << G_diff_5dcube[i][j] << endl ;
  return G_diff;
}





// template <typename T>
// cx_5d<T> G_ijk_full_matrix(Mat<T> source_point_mat,Mat<T> target_point_mat,T k_s , T k_p,T rho ,T omega){
//  int no_of_source_point = source_point_mat.n_rows;
//  int no_of_target_point = target_point_mat.n_rows;
//  cx_5d<T>  G_for_all_source_and_target(no_of_source_point,cx_4d<T>(no_of_target_point,Cube<complex<T>>(3,3,3,fill::zeros)));
//  #pragma omp parallel for
//   for (int i = 0; i < G_for_all_source_and_target.size(); ++i) {
//     Row<T> temp_row_gre = source_point_mat.row(i);
//     G_for_all_source_and_target[i] = G_p_diff_ijk<T>(temp_row_gre,target_point_mat ,k_s , k_p,rho ,omega);
//   }
//   return(G_for_all_source_and_target);  
// }

// template cx_5d<float> G_ijk_full_matrix(Mat<float>,Mat<float>,float,float,float,float);

// template cx_5d<double> G_ijk_full_matrix(Mat<double>,Mat<double>,double,double,double,double);
  
// ------------------------------------------------------------------//
// -----------------Strain Calculation Procedur---------------------//
// ----------------------------------------------------------------//


// Strain Calculation for Onesource For Entire Observation Target
// template <typename T>
// cx_3d<T> strain_cal_for_one_source(cx_4d<T> G_ijk_cube, Row<complex<T>> strength_cube){

//   Col<complex<T>> strength_cube_col = strength_cube.t();
//   int n_row_cube = G_ijk_cube.size();

//   cx_3d<T> obs_plane_stress_comp(n_row_cube,Mat<complex<T>>(3,3,fill::zeros));
  
//   /*
//   for(size_t i = 0; i < n_row_cube; i++){
//     obs_plane_stress_comp[i] = cx_mat(3,3,fill::zeros);
//   }
//   */
  
//   for(size_t i = 0; i < n_row_cube; i++){
//     for(size_t k = 0; k < 3; k++){
//       Cube<complex<T>> mat_to_operate = G_ijk_cube[i].subcube( 0, k, 0, 2, k, 2 );
//       Mat<complex<T>> converted_cub_tomat = conv_cube_2_mat(mat_to_operate);
//       //  obs_plane_stress_comp[i][j].col(k)= (converted_cub_tomat  + converted_cub_tomat.t()) * strength_cube_col ;
//       obs_plane_stress_comp[i]= obs_plane_stress_comp[i] + (converted_cub_tomat  + converted_cub_tomat.t()) * strength_cube_col(k) ;
//     }
//   }
//   return(obs_plane_stress_comp);
// }

// template cx_3d<double> strain_cal_for_one_source(cx_4d<double>, Row<complex<double>>);
// template cx_3d<float> strain_cal_for_one_source(cx_4d<float>, Row<complex<float>>);

// template <typename T>
// cx_3d<T> strain_cal(cx_5d<T> G_ijk_source,Mat<complex<T>> src_str){
//   int target_size = G_ijk_source.size();
//   cx_3d<T> strains_on_plane(target_size,Mat<complex<T>>(3,3,fill::zeros));
//   #pragma omp parallel for
//   for(size_t i = 0; i < src_str.n_cols; i++){
//     strains_on_plane = add_cx_3d<T>(strains_on_plane,strain_cal_for_one_source<T>(G_ijk_source[i],src_str.row(i)));
//   }
//   return(strains_on_plane);
// }

// template cx_3d<float> strain_cal(cx_5d<float>,Mat<complex<float>> );
// template cx_3d<double> strain_cal(cx_5d<double>,Mat<complex<double>> );

// ------------------------------------------------------------------//
// -----------------Stress Calculation Procedure---------------------//
// ----------------------------------------------------------------//
// template<typename T>
// Cube<T> transpose_cube_in_dir_2(Cube<T> cube_to_trans){
//   Cube<T> return_cx_cube;
//   for (int i=0; i < 0; ++i) {
//     return_cx_cube.slice(i) = cube_to_trans.slice(i).t();
//   }
//   return(return_cx_cube);
// }

// template Cube<double> transpose_cube_in_dir_2(Cube<double>);
// template Cube<float> transpose_cube_in_dir_2(Cube<float>);


template <typename T>
Cube<complex<T>> stress_coff_point(Cube<complex<T>> G_ijk,T lamda,T mu){

  Cube<complex<T>> Stress_coff(3,3,3,fill::zeros);

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

template Cube<complex<float>> stress_coff_point(Cube<complex<float>>,float,float);
template Cube<complex<double>> stress_coff_point(Cube<complex<double>>,double,double);

template <typename T>
Cube<complex<T>> S_coff_source_target (Row<T> source_point_row,Row<T> target_point_row,T k_s , T k_p,T rho , T omega,T lamda,T mu){
  Cube<complex<T>> G = G_diff_ijk_point( source_point_row, target_point_row,k_s ,  k_p, rho ,  omega);
  Cube<complex<T>> S_COFF = stress_coff_point(G,lamda,mu);
  return S_COFF;
}

template  Cube<complex<float>> S_coff_source_target (Row<float>,Row<float> ,float , float,float , float,float,float);
template  Cube<complex<double>> S_coff_source_target (Row<double>,Row<double> ,double , double,double , double,double,double);

// template <typename T>
// cx_5d<T> stress_coff_calc(const cx_5d<T> G_ijk_source,T lamda,T mu){
//   int size_1 = G_ijk_source.size();
//   int size_2 = G_ijk_source[0].size();
//   cx_5d<T> Converterd_stress_coff(size_1,cx_4d<T>(size_2,Cube<complex<T>>(3,3,3,fill::zeros)));
  
//   for (int i = 0; i < size_1; ++i) {
//     for (int j = 0; j < size_2; ++j) {
//       Converterd_stress_coff[i][j] = stress_coff_point<T>(G_ijk_source[i][j],lamda,mu);
//     }
//   }
//   return(Converterd_stress_coff);
// }

// template cx_5d<double> stress_coff_calc(const cx_5d<double>,double,double);
// template cx_5d<float> stress_coff_calc(const cx_5d<float>,float,float);



// template <typename T>
// cx_3d<T> stress_cal_one_source(cx_4d<T> G_ijk_cube, Row<complex<T>> strength_cube){

//   Col<complex<T>> strength_cube_col = strength_cube.t();
//   int n_row_cube = G_ijk_cube.size();

//   cx_3d<T> obs_plane_stress_comp(n_row_cube,Mat<complex<T>>(3,3,fill::zeros));
  
//   for(size_t i = 0; i < n_row_cube; i++){
//     for(size_t k = 0; k < 3; k++){
//       Cube<complex<T>> mat_to_operate = G_ijk_cube[i].subcube( 0, k, 0, 2, k, 2 );
//       Mat<complex<T>> converted_cub_tomat = conv_cube_2_mat(mat_to_operate);
//       obs_plane_stress_comp[i]= obs_plane_stress_comp[i] + (converted_cub_tomat) * strength_cube_col(k) ;
//     }
//   }
//   return(obs_plane_stress_comp);
// }

// template cx_3d<float> stress_cal_one_source(cx_4d<float>, Row<complex<float>>);
// template cx_3d<double> stress_cal_one_source(cx_4d<double>, Row<complex<double>>);


// template <typename T>
// cx_3d<T> stress_cal(const cx_5d<T> G_ijk_source,Mat<complex<T>> src_str,T lamda,T mu){ // the G_ijk_must be Converted Soueces here

//   int target_size = G_ijk_source[0].size(); // Thought as Correct 
//   cx_3d<T> stress_on_plane(target_size,Mat<complex<T>>(3,3,fill::zeros));
  
//   // MultiThreading the Calculations Using OPenMP
//   #pragma omp parallel for
//   for(size_t i = 0; i < src_str.n_cols; i++){
//     cx_3d<T> temp_call = stress_cal_one_source<T>(G_ijk_source[i],src_str.row(i));
//     stress_on_plane = add_cx_3d<T>(stress_on_plane,temp_call);
//   }
//   return(stress_on_plane);
// }

// template cx_3d<double> stress_cal(const cx_5d<double>,Mat<complex<double>>,double,double);
// template cx_3d<float> stress_cal(const cx_5d<float>,Mat<complex<float>>,float,float);



// // The Entire Fucntion is Redundant Waste Memory And thereof Needed Care ful Restructuing of Code To Avoid Memory Problem in Furture
// // This Paticular Function is used with constat Pointer Reference to Minimize the Memory usage Peak During the Computation
// template <typename T>
// cx_3d<T> stress_from_points(Mat<T> source_points, Mat<complex<T>> src_str, Mat<T> target_points,T k_s , T k_p,T rho ,T omega , T lamda,T mu){
//   cx_5d<T> G_vals_for_given = G_ijk_full_matrix<T>(source_points,target_points,k_s ,k_p,rho , omega); 
//   cx_5d<T> stress_for_given  =  stress_coff_calc<T>(G_vals_for_given,lamda, mu); // Pass by Reference
//   cx_3d<T> str_vals = stress_cal<T>(stress_for_given,src_str,lamda,mu); // Pass By Reference 
//   return(str_vals);
// }

// template cx_3d<float> stress_from_points(Mat<float>,Mat<complex<float>>,Mat<float>,float,float,float,float,float,float);
// template cx_3d<double> stress_from_points(Mat<double>,Mat<complex<double>>,Mat<double>,double,double,double,double,double,double);

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


// template <typename T>
// Col<complex<T>> source_strength_derivation(cx_5d<T> G_ijk_cube_of_sources,cx_3d<T> stress_matrix){
  
//   int sources_len = stress_matrix.size();
//   int cx_source = G_ijk_cube_of_sources.size();
//   cout << " Source Strength Derivation --"<< cx_source <<"No_of Source Supplied -" << sources_len << endl;
//   // cx_mat str_mat(sources_len*stress_matrix[0].n_rows,stress_matrix[0].n_cols,fill::zeros);
//   Col<complex<T>> str_mat(sources_len * stress_matrix[0].n_elem);

  
//   int source_no = G_ijk_cube_of_sources.size() ;
//   int target_no =  G_ijk_cube_of_sources[0].size();

//   for (int i = 0; i < source_no; ++i) {
//     for (int j = 0; j <  target_no; ++j) {
//       if (G_ijk_cube_of_sources[i][j].has_nan()) {
// 	cout << "NaN Detected" << "SOURCE --" << i << " Target--" << j << endl;
//      } 
//     }
//   }
//   cout << "Nan Detection Complete" << endl;
  
//   Cube<complex<T>> G_full(target_no * 3,source_no * 3,3,fill::zeros);

//   //std::cout << "Creating G_full" << endl;
//   // std::cout << G_full << endl;

// #pragma omp parallel for
//   for(size_t i = 0; i < sources_len; i++){
//     for(size_t j =0; j < 3 ; ++j){
//       str_mat.subvec(i*9+j*3,i*9+j*3+2) = stress_matrix[i].col(j);
//     }
//   }

//   //std::cout << "Creating srmat" << endl;
//   //std::cout << str_mat << endl;

// #pragma omp parallel for
//   for(size_t j = 0; j < target_no; j++){ 
//     for(size_t i = 0; i < source_no; i++){
//       G_full(span(j*3,(j*3)+2),span(i*3,i*3+2),span(0,2)) =  G_ijk_cube_of_sources[i][j];
//     }
//   }
  
//   //std::cout << "Calculated  G_full" << endl;
//   //std::cout << G_full << endl;
  
//   // Converting Above Cube tot the matrix

//   Mat<complex<T>> G_full_mat_form(G_full.n_rows * 3,G_full.n_cols,fill::zeros);

// #pragma omp parallel for
//   for(size_t i = 0; i < G_full.n_rows /3; i++){
//     for(size_t j =0; j < 3 ; ++j){
//         Mat<complex<T>> temper = G_full(span(0,2),span(0,G_full.n_cols-1),span(j,j));
//       //      std::cout << "adjsnfdb" <<i << j << endl ;
//       G_full_mat_form(span(i*9+j*3,i*9+j*3+2), span(0,G_full.n_cols-1) ) = temper ;
//     }
//   }

//   //std::cout << G_full_mat_form <<endl;
//   cout << " Start_solving " << sources_len << endl;

  
//   G_full_mat_form.save("D_BUG_G_Full.csv",csv_ascii);
//   str_mat.save("D_BUG_str_mat.csv",csv_ascii);

//   cout << "G_FULL ROWS --" << G_full_mat_form.n_rows << "Cols--" << G_full_mat_form.n_cols << endl;
//   Col<complex<T>> col_one = solve(G_full_mat_form,str_mat);

  
//   cout << " End Solving -" << sources_len << endl;
//   return(col_one);
// }


// template Col<complex<float>> source_strength_derivation(cx_5d<float>,cx_3d<float>);
// template Col<complex<double>> source_strength_derivation(cx_5d<double>,cx_3d<double>);
  
// ------------------------------------------------------------------//
// High Level Function For point Strenght caculation----------------//
// ----------------------------------------------------------------//


// // Need to Redesign The Hight level function for Hetro geneous and homogeneous material

// template<typename T>
// Mat<complex<T>> get_strength(Mat<T> source_point_list,T r_s,Row<T> direction_cosine,cx_3d<T> stress_matrix,T k_s,T k_p,T rho , T omega){
  
//   Mat<T> target_point(size(source_point_list));
//   for(int i = 0 ; i < source_point_list.n_rows;++i){
//     target_point.row(i) = source_point_list.row(i) - r_s * direction_cosine;
//   }
  
//   int no_of_source_points = source_point_list.n_rows;
//   cx_5d<T> G_matrix_calculations(no_of_source_points,cx_4d<T>(no_of_source_points,Cube<complex<T>>(3,3,3,fill::zeros)));
  
//   for(int i = 0; i < no_of_source_points; ++i){
//     G_matrix_calculations[i] = G_p_diff_ijk<T>(source_point_list.row(i), target_point,k_s,k_p,rho ,omega);
//   }

  
//   Col<complex<T>> point_strenghts = source_strength_derivation(G_matrix_calculations,stress_matrix);
  
  
//   T out_mat_rows =  point_strenghts.n_elem/3;
//   Mat<complex<T>> P_mat_form(out_mat_rows,3,fill::zeros);
  
//   for(int i = 0;i<out_mat_rows;++i){
//     for( int j = 0; j < 3 ;++j){
//       P_mat_form(i,j) = point_strenghts(i*3+j);
//     }
//   }
//   return(P_mat_form);
// }

// template Mat<complex<float>> get_strength(Mat<float>,float,Row<float>,cx_3d<float>,float,float,float,float);
// template Mat<complex<double>> get_strength(Mat<double>,double,Row<double>,cx_3d<double>,double,double,double,double);


// template <typename T>
// Mat<complex<T>> get_strength_hetro(Mat<T> source_point_list,Mat<T> passive_source_list ,T r_s, Row<T> direction_cosine,cx_3d<T> stress_matrix,T k_s,T k_p,T rho , T omega,T mu,T lamda,bool have_active_sources){

  
//   // incoming Point Locations are Exact Location on which the Sources are placed
  
//   int no_of_active_sources = source_point_list.n_rows;
//   int no_of_passive_sources = passive_source_list.n_rows;

//   Mat<complex<T>> active_source_str(no_of_active_sources,3,fill::zeros);
//   Mat<complex<T>> passive_source_str(no_of_passive_sources,3,fill::zeros);

//   // Mat<T> total_sources = join_cols(source_point_list,passive_source_list);

//   Mat<T> total_sources = (have_active_sources) ? join_cols(source_point_list,passive_source_list) : passive_source_list;
  
//   Mat<T> exact_source_point_mat_stress(no_of_active_sources,3);

  
//   cout << "Black Sheep 7" << endl;
// #pragma omp parallel for 
//   for (int i = 0; i < no_of_active_sources; ++i) {
//     exact_source_point_mat_stress.row(i) = source_point_list.row(i)  + direction_cosine * r_s ;
//   }

  
//    // cout << source_point_list << endl;
//    // cout << "-------------------------------" << endl;
   
//    // cout << passive_source_list << endl;
//    // cout << "-------------------------------" << endl;
   
//   cout << "Printing INPUT" << endl;
//   if (!have_active_sources) {
//     total_sources.save("D_BUG_TOATL_SOURCES.csv",csv_ascii);
//     exact_source_point_mat_stress.save("D_BUG_Exact_SOURCES.csv",csv_ascii);
//   }
  
//   cx_5d<T> G_matrix_calculations(total_sources.n_rows,cx_4d<T>(no_of_active_sources ,Cube<complex<T>>(3,3,3,fill::zeros)));
// #pragma omp parallel for 
//   for(int i = 0; i < total_sources.n_rows; ++i){
//     //std::cout << i << "\n";
//     G_matrix_calculations[i] = G_p_diff_ijk<T>(total_sources.row(i),exact_source_point_mat_stress,k_s,k_p,rho ,omega);
//   }
  
//   G_matrix_calculations = stress_coff_calc<T>(G_matrix_calculations,mu,lamda);
  
//   cout << "Starting Solving Linear Equations" << endl;
//   Col<complex<T>> solid_point_strength_colvec = source_strength_derivation(G_matrix_calculations,stress_matrix);
//   cout << "Completed Solving Linear Equations" << endl;
  
//   int out_mat_rows =  solid_point_strength_colvec.n_elem/3;
//   Mat<complex<T>> P_mat_form(out_mat_rows,3,fill::zeros);
  
//   for(int i = 0;i<out_mat_rows;++i){
//     for( int j = 0; j < 3 ;++j){
//       P_mat_form(i,j) = solid_point_strength_colvec(i*3+j);
//     }
//   }
  
//   return(P_mat_form);
// }


// // Initializing The Template Class
// template Mat<complex<float>> get_strength_hetro(Mat<float>,Mat<float>,float,Row<float>,cx_3d<float>,float,float,float,float,float,float,bool);
// template Mat<complex<double>> get_strength_hetro(Mat<double>,Mat<double>,double,Row<double>,cx_3d<double>,double,double,double,double,double,double,bool);

// template <typename T>
// Mat<complex<T>> solve_dpsm_str (Mat<T> ACTIVE_SOURCES_DPSM_POINT, Mat<T> PASSIVE_SOURCES_DPSM_POINT, cx_3d<T> STRESS_CX_3D_MATRIX , Mat<T> Points_of_enforcement, T k_s,T k_p,T rho , T omega,T mu,T lamda){
  
//   int total_source_no = ACTIVE_SOURCES_DPSM_POINT.n_rows + PASSIVE_SOURCES_DPSM_POINT.n_rows;
//   int no_points_enforce = Points_of_enforcement.n_rows;
  
//   Mat<T> TOTAL_SOURCE_DPSM_MAT = join_cols(ACTIVE_SOURCES_DPSM_POINT,PASSIVE_SOURCES_DPSM_POINT);
  
//   cx_5d<T> G_matrix(total_source_no,cx_4d<T>(no_points_enforce ,Cube<complex<T>>(3,3,3,fill::zeros)));
// #pragma omp parallel for 
//   for(int i = 0; i < total_source_no; ++i){
//     //std::cout << i << "\n";
//     G_matrix[i] = G_p_diff_ijk<T>(TOTAL_SOURCE_DPSM_MAT.row(i),Points_of_enforcement,k_s,k_p,rho ,omega);
//   }

//   G_matrix = stress_coff_calc<T>(G_matrix,mu,lamda);
  
//   cout << "Starting Solving Linear Equations" << endl;
//   Col<complex<T>> solid_point_strength_colvec = source_strength_derivation(G_matrix,STRESS_CX_3D_MATRIX);
//   cout << "Completed Solving Linear Equations" << endl;

//   int out_mat_rows =  solid_point_strength_colvec.n_elem/3;
//   Mat<complex<T>> P_mat_form(out_mat_rows,3,fill::zeros);
  
//   for(int i = 0;i<out_mat_rows;++i){
//     for( int j = 0; j < 3 ;++j){
//       P_mat_form(i,j) = solid_point_strength_colvec(i*3+j);
//     }
//   }
  
//   return(P_mat_form);
// }

// template Mat<complex<float>> solve_dpsm_str (Mat<float>, Mat<float>, cx_3d<float>, Mat<float>, float,float,float,float,float,float);
// template Mat<complex<double>> solve_dpsm_str (Mat<double>, Mat<double>, cx_3d<double>, Mat<double>, double,double,double,double,double,double);

template <typename T> Mat<complex<T>>
EQN_ROW(const Mat<T> &ACTIVE_SOURCES_DPSM_POINT, Row<T> Points_of_enforcement,const T &k_s,const T &k_p,const T & rho , const T & omega,const T & mu,const T & lamda){
  int no_active_points = ACTIVE_SOURCES_DPSM_POINT.n_rows;
  Cube<complex<T>> S_Coff_for_target(3,3,3,fill::zeros);
  Mat<complex<T>> EQN_MAT(9,no_active_points*3,fill::zeros);
  for (int i= 0; i < no_active_points; ++i) {
    Cube<complex<T>> G_Cube = S_coff_source_target<T>(ACTIVE_SOURCES_DPSM_POINT.row(i), Points_of_enforcement, k_s ,  k_p, rho ,  omega, lamda, mu);
    for (int j = 0; j < 3; ++j) {
      EQN_MAT(span(j*3,j*3+2),span(i*3,i*3+2)) = G_Cube.slice(j);
    }
  }
  return(EQN_MAT);
}

template Mat<complex<double>> EQN_ROW(const Mat<double> &, Row<double>,const double &,const double &,const double &, const double &,const double & ,const double &);
template Mat<complex<float>> EQN_ROW(const Mat<float> &, Row<float>,const float &,const float &,const float &, const float &,const float & ,const float &);
  
template <typename T>
Mat<complex<T>> EQN_assembler (Mat<T> ACTIVE_SOURCES_DPSM_POINT, Mat<T> Points_of_enforcement, T k_s,T k_p,T rho , T omega,T mu,T lamda){
  
  int no_active_points = ACTIVE_SOURCES_DPSM_POINT.n_rows;
  int no_target_points = Points_of_enforcement.n_rows;
  Mat<complex<T>> A_COFF(no_target_points*9,no_active_points * 3,fill::zeros); 
  for (int i = 0; i < no_target_points; ++i) {
    A_COFF(span(i*9,i*9+8),span::all) = EQN_ROW<T>(ACTIVE_SOURCES_DPSM_POINT,Points_of_enforcement.row(i),k_s,k_p, rho ,omega,mu, lamda);
  }
  return(A_COFF);
}

template Mat<complex<float>> EQN_assembler (Mat<float>, Mat<float>, float,float,float, float,float,float);
template Mat<complex<double>> EQN_assembler (Mat<double>, Mat<double>, double,double,double, double,double,double);

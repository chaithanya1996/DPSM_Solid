
#include "dpsm_helpers.hpp"

template <typename TYPE>
TYPE r_s_calculator(TYPE freq, TYPE c,int safety_factor){
  TYPE r_s;
  r_s = c / freq / 2 / M_PI / safety_factor;
  return(r_s);
};

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

// ------------------------------------------------------------------//
// -----------------G Diff calculation procedures-------------------//
// ----------------------------------------------------------------//
template <typename T>
complex<T> r_diff(T radius_mag,T k,T R_val){
  complex<T> img_start(0,1.0);
  complex<T> value_of_differential =complex<T>(2 * R_val / T(pow<T>(radius_mag,3)) - img_start * k * R_val / T(pow<T>(radius_mag,2)));
  return(value_of_differential);
}


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


// ---------------------------------------------------------------------------------//
// small subset_functions


template <typename T>
T eoiidi(T i_val ,Row<T> R_vector,T r_mag){
  T return_value =   - 2 * pow<T>(R_vector(i_val),3) / r_mag + 2 * R_vector(i_val) / r_mag;
  //cout << "Triggered -- iii" << endl; 
  return(return_value);
}


template <typename T>
T eoijdk(T i_val ,T j_val ,Row<T> R_vector,T r_mag){
  T return_value =   -2 * R_vector(0) * R_vector(1) * R_vector(2) / r_mag ;
  //cout << "Triggered -- ijk" << endl;
  return(return_value);
}

  

template <typename T>
T eoiidj(T i_val ,T j_val ,Row<T> R_vector,T r_mag){
  T return_value =   -2 * R_vector(i_val) * R_vector(i_val) * R_vector(j_val) / r_mag ;
  //cout << "Triggered -- iij -- i_val" << i_val << "  j_val--" << j_val<< endl;
  return(return_value);
}


template <typename T>
T eoijdi(T i_val ,T j_val,Row<T> R_vector,T r_mag){
  T return_value =   - 2 * pow<T>(R_vector(i_val),2) * R_vector(j_val) / r_mag +  R_vector(j_val) / r_mag;
  //cout << "Triggered -- iji -- i_val" << i_val << "  j_val--" << j_val  << endl;
  return(return_value);
}



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
  complex<T> Gpiidi  = (e_p * ( k_p * k_p *  eo_d_universal(ijk_row,R_vec,radius_mag) +
   ( G_ijk_helper(ijk_row) + 3 * R_vec[ijk_row[0]] * R_vec[ijk_row[1]]) * R_differ_mat(ijk_row[2],0) +  r_p * T (3) * eo_d_universal(ijk_row,R_vec,radius_mag)))/ multiplier +
    G_place_holder(0) * (img_start * k_p * R_vec[ijk_row[2]] * e_p - R_vec[ijk_row[2]] * e_p / radius_mag);			 

  complex<T> Gsiidi  = (e_s * ( -k_s* k_s *  eo_d_universal(ijk_row,R_vec,radius_mag) -  ( G_ijk_helper(ijk_row) +3 * R_vec[ijk_row[0]] * R_vec[ijk_row[1]]) * R_differ_mat(ijk_row[2],0) +  r_s * T (3) * eo_d_universal(ijk_row,R_vec,radius_mag))) / multiplier +
    G_place_holder(1) * (img_start * k_s * R_vec[ijk_row[2]] * e_s - R_vec[ijk_row[2]] * e_s / radius_mag);			 
  
  return Gpiidi + Gsiidi ;
}


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


template <typename T>
Cube<complex<T>> S_coff_source_target (Row<T> source_point_row,Row<T> target_point_row,T k_s , T k_p,T rho , T omega,T lamda,T mu){
  Cube<complex<T>> G = G_diff_ijk_point( source_point_row, target_point_row,k_s ,  k_p, rho ,  omega);
  Cube<complex<T>> S_COFF = stress_coff_point(G,lamda,mu);
  return S_COFF;
}




template <typename T> Mat<complex<T>>
EQN_ROW(const Mat<T> &ACTIVE_SOURCES_DPSM_POINT, Row<T> Points_of_enforcement, T k_s, T k_p,T rho , T omega, T mu,T lamda){
  int no_active_points = ACTIVE_SOURCES_DPSM_POINT.n_rows;
  Cube<complex<T>> S_Coff_for_target(3,3,3,fill::zeros);
  Mat<complex<T>> EQN_MAT(9,no_active_points*3,fill::zeros);
  for (int i= 0; i < no_active_points; ++i) {
    Cube<complex<T>> G_Cube = S_coff_source_target<T>(ACTIVE_SOURCES_DPSM_POINT.row(i), Points_of_enforcement, k_s ,  k_p, rho ,  omega, lamda, mu);
    for (int j = 0; j < 3; ++j) {
      Mat<complex<T>> TEMP_MAT_EXTRACT = CUBE_EXTRACT(G_Cube,j);
      EQN_MAT(span(j*3,j*3+2),span(i*3,i*3+2)) = TEMP_MAT_EXTRACT;
      // EQN_MAT(span(j*3,j*3+2),span(i*3,i*3+2)) = G_Cube.slice(j);
    }
  }
  return(EQN_MAT);
}


template <typename T>
Mat<complex<T>> EQN_assembler (const Mat<T> &ACTIVE_SOURCES_DPSM_POINT, const Mat<T> &Points_of_enforcement, T k_s,T k_p,T rho , T omega,T mu,T lamda){
  
  int no_active_points = ACTIVE_SOURCES_DPSM_POINT.n_rows;
  int no_target_points = Points_of_enforcement.n_rows;
  Mat<complex<T>> A_COFF(no_target_points*9,no_active_points * 3,fill::zeros);
  #pragma omp parallel for
  for (int i = 0; i < no_target_points; ++i) {
    A_COFF(span(i*9,i*9+8),span::all) = EQN_ROW<T>(ACTIVE_SOURCES_DPSM_POINT,Points_of_enforcement.row(i),k_s,k_p, rho ,omega,mu, lamda);
  }
  return(A_COFF);
}

template Mat<complex<float>> EQN_assembler (const Mat<float>&, const Mat<float>&, float,float,float, float,float,float);
template Mat<complex<double>> EQN_assembler (const Mat<double> &,const Mat<double>&, double,double,double, double,double,double);


template <typename T>
Mat<complex<T>> solve_dpsm_str (const Mat<T> & ACTIVE_SOURCES_DPSM_POINT,const Mat<T> & PASSIVE_SOURCES_DPSM_POINT, cx_3d<T> STRESS_CX_3D_MATRIX , const Mat<T> & Points_of_enforcement, T k_s,T k_p,T rho , T omega,T mu,T lamda,bool DBUG_STATUS=false){
  
  Mat<T> TOTAL_SOURCES_DPSM_MAT = join_cols(ACTIVE_SOURCES_DPSM_POINT,PASSIVE_SOURCES_DPSM_POINT);
  Mat<complex<T>> EQN_MAT_assembled = EQN_assembler(TOTAL_SOURCES_DPSM_MAT,Points_of_enforcement,k_s,k_p, rho ,omega,mu, lamda);
  
  if (EQN_MAT_assembled.has_nan()){
    cout << " NAN Detedted" << endl;
  }
  Col<complex<T>> STRESS_Col_ENFORCER(STRESS_CX_3D_MATRIX.size() * 9,fill::zeros);

  
  #pragma omp parallel for
  for (int i = 0; i < STRESS_CX_3D_MATRIX.size() ; ++i) {
    for (int j = 0; j < 3; ++j) {
      STRESS_Col_ENFORCER(span(i*9+j*3,i*9+j*3+2)) = STRESS_CX_3D_MATRIX[i].col(j);
    }
  }
  if(DBUG_STATUS){
    cout << "Started Writing To Disk " << endl;
    EQN_MAT_assembled.save("DBUG_Giii_A_MATRIX.csv",csv_ascii);
    STRESS_Col_ENFORCER.save("DBUG_Giii_B_MATRIX.csv",csv_ascii);
    cout << "Completed Writing To Disk " << endl;
  }
  cout << "Starting Solving Linear Equations" << endl;
  Col<complex<T>> solid_point_strength_colvec = solve(EQN_MAT_assembled,STRESS_Col_ENFORCER);
  cout << "Completed Solving Linear Equations" << endl;
  
  cout << solid_point_strength_colvec.n_elem/ 3 << endl;
  
  Mat<complex<T>> SOLID_STR_MAT(solid_point_strength_colvec.n_elem/ 3,3,fill::zeros) ;
  cout << "Cowabanga" << endl;
  for (int i = 0; i < SOLID_STR_MAT.n_rows; ++i) {
    Row<complex<T>> temp_rowvec =trans(solid_point_strength_colvec(span(i*3,i*3+2)));
    SOLID_STR_MAT.row(i) = temp_rowvec;
  }

  cout << "Cowabanga_2" << endl;
  return(SOLID_STR_MAT);
}



template <typename T>
Mat<complex<T>> stress_calc_mat_ver(const Mat<T> &ACTIVE_SOURCES_DPSM_POINT, const Mat<T> &Points_of_enforcement, const Mat<complex<T>> &ACTIVE_STR , T k_s,T k_p,T rho , T omega,T mu,T lamda){
  Mat<complex<T>> EQN_MAT  = EQN_assembler(ACTIVE_SOURCES_DPSM_POINT,Points_of_enforcement,k_s,k_p, rho ,omega,mu, lamda);
  Mat<complex<T>> stress_vals = EQN_MAT * ACTIVE_STR;
  return(stress_vals);
}


template <typename T>
cx_3d<T> stress_calc_3d_ver(const Mat<T> &ACTIVE_SOURCES_DPSM_POINT, const Mat<T> &Points_of_enforcement, const Mat<complex<T>> &ACTIVE_STR , T k_s,T k_p,T rho , T omega,T mu,T lamda){
  Mat<complex<T>> EQN_MAT  = EQN_assembler(ACTIVE_SOURCES_DPSM_POINT,Points_of_enforcement,k_s,k_p, rho ,omega,mu, lamda);
  Col<complex<T>> ACTIVE_STR_COL(ACTIVE_STR.n_elem,fill::zeros);
  #pragma omp parallel for
  for (int i = 0; i < ACTIVE_STR.n_rows; ++i) {
    ACTIVE_STR_COL(span(i*3,i*3+2)) = trans(ACTIVE_STR.row(i));
  }
  Col<complex<T>> stress_vals = EQN_MAT * ACTIVE_STR_COL;
  cx_3d<T> stress_vals_cx(stress_vals.n_elem/9,Mat<complex<T>>(3,3,fill::zeros));
  #pragma omp parallel for
  for (int i = 0; i < stress_vals_cx.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      stress_vals_cx[i].col(j) = stress_vals(span(i*9+j*3,i*9+j*3+2));
    }
  }
  return(stress_vals_cx);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <typename T> Mat<complex<T>>
EQN_ROW_DISP(const Mat<T> &ACTIVE_SOURCES_DPSM_POINT, Row<T> Points_of_enforcement, T k_s, T k_p,T rho , T omega){
  int no_active_points = ACTIVE_SOURCES_DPSM_POINT.n_rows;
  Mat<complex<T>> G_COFF(3,3,fill::zeros);
  Mat<complex<T>> EQN_MAT(3,no_active_points*3,fill::zeros);
  for (int i= 0; i < no_active_points; ++i) {
    Mat<complex<T>> G_COFF_MAT = G_mat_for_point<T>(ACTIVE_SOURCES_DPSM_POINT.row(i), Points_of_enforcement, k_s ,  k_p, rho ,  omega);
    EQN_MAT(span::all,span(i*3,i*3+2)) = G_COFF_MAT;
  }
  return(EQN_MAT);
}

  
template <typename T>
Mat<complex<T>> EQN_assembler_DISP (const Mat<T> &ACTIVE_SOURCES_DPSM_POINT, const Mat<T> &Points_of_enforcement, T k_s,T k_p,T rho , T omega){
  
  int no_active_points = ACTIVE_SOURCES_DPSM_POINT.n_rows;
  int no_target_points = Points_of_enforcement.n_rows;
  Mat<complex<T>> A_COFF(no_target_points*3,no_active_points*3,fill::zeros);
  #pragma omp parallel for
  for (int i = 0; i < no_target_points; ++i) {
    A_COFF(span(i*3,i*3+2),span::all) = EQN_ROW_DISP<T>(ACTIVE_SOURCES_DPSM_POINT,Points_of_enforcement.row(i),k_s,k_p, rho ,omega);
  }
  return(A_COFF);
}


template <typename T>
Mat<complex<T>> solve_dpsm_disp (const Mat<T> & ACTIVE_SOURCES_DPSM_POINT,const Mat<T> & PASSIVE_SOURCES_DPSM_POINT, Mat<complex<T>> disp_mat , const Mat<T> & Points_of_enforcement, T k_s,T k_p,T rho , T omega,bool DBUG_STATUS=false){
  
  Mat<T> TOTAL_SOURCES_DPSM_MAT = join_cols(ACTIVE_SOURCES_DPSM_POINT,PASSIVE_SOURCES_DPSM_POINT);
  Mat<complex<T>> EQN_MAT_assembled = EQN_assembler_DISP(TOTAL_SOURCES_DPSM_MAT,Points_of_enforcement,k_s,k_p, rho ,omega);
  
  if (EQN_MAT_assembled.has_nan()){
    cout << " NAN Detedted" << endl;
  }
  Col<complex<T>> STRESS_Col_ENFORCER(disp_mat.n_elem,fill::zeros);
  
  #pragma omp parallel for
  for (int i = 0; i < disp_mat.n_rows ; ++i) {
    STRESS_Col_ENFORCER(span(i*3,i*3+2)) = trans(disp_mat.row(i));
  }

  if(DBUG_STATUS){
    cout << "Started Writing To Disk " << endl;
    EQN_MAT_assembled.save("DBUG_A_MATRIX.csv",csv_ascii);
    STRESS_Col_ENFORCER.save("DBUG_B_MATRIX.csv",csv_ascii);
    cout << "Completed Writing To Disk " << endl;
  }
  cout << "Starting Solving Linear Equations" << endl;
  Col<complex<T>> solid_point_strength_colvec = solve(EQN_MAT_assembled,STRESS_Col_ENFORCER);
  cout << "Completed Solving Linear Equations" << endl;
  
  cout << solid_point_strength_colvec.n_elem/ 3 << endl;
  
  Mat<complex<T>> SOLID_STR_MAT(solid_point_strength_colvec.n_elem/ 3,3,fill::zeros) ;
  cout << "Started Ctrating STR_MAT" << endl;
  
  for (int i = 0; i < SOLID_STR_MAT.n_rows; ++i) {
    Row<complex<T>> temp_rowvec =trans(solid_point_strength_colvec(span(i*3,i*3+2)));
    SOLID_STR_MAT.row(i) = temp_rowvec;
  }
  cout << "Started Ctrating STR_MA" << endl;
  return(SOLID_STR_MAT);
}

template <typename T>
Mat<complex<T>> disp_calc_3d_ver(const Mat<T> &ACTIVE_SOURCES_DPSM_POINT, const Mat<T> &Points_of_enforcement, const Mat<complex<T>> &ACTIVE_STR , T k_s,T k_p,T rho , T omega){
  Mat<complex<T>> EQN_MAT  = EQN_assembler_DISP(ACTIVE_SOURCES_DPSM_POINT,Points_of_enforcement,k_s,k_p, rho ,omega);
  Col<complex<T>> ACTIVE_STR_COL(ACTIVE_STR.n_elem,fill::zeros);
  #pragma omp parallel for
  for (int i = 0; i < ACTIVE_STR.n_rows; ++i) {
    ACTIVE_STR_COL(span(i*3,i*3+2)) = trans(ACTIVE_STR.row(i));
  }
  Col<complex<T>> disp_vals_cal = EQN_MAT * ACTIVE_STR_COL;
  Mat<complex<T>> Disp_vals_cx(disp_vals_cal.n_elem/3,3);
  #pragma omp parallel for
  for (int i = 0; i < Disp_vals_cx.n_rows; ++i) {
    Disp_vals_cx.row(i) = trans(disp_vals_cal(span(i*3,i*3+2)));
  }
  return(Disp_vals_cx);
}

#include "dpsm_helpers.hpp"



// ------------------------------------------------------------------//
// --------------Basic Arthemetic on Defined Fucntion---------------//
// ----------------------------------------------------------------//

cx_3d add_cx_3d (cx_3d former, cx_3d latter){
  // No Explicit Checking is Defined and Dev is asssumed to  Use fucntions sane manner
  double lend_cx_3d = former.size();
  cx_3d ret_data(lend_cx_3d,cx_mat(size(former[0])));
  for(size_t i = 0; i < lend_cx_3d; i++){
    ret_data[i] = former[i] + latter[i];
  }
  return(ret_data);
}


cx_mat conv_cube_2_mat(cx_cube to_be_conv){
  int n_row_cube = to_be_conv.n_rows;
  int n_slice_cube = to_be_conv.n_slices;
  cx_mat mat_to_return(n_row_cube,n_slice_cube,fill::zeros);
  
  for(size_t i = 0; i < n_row_cube; i++){
    cx_rowvec temp_row = to_be_conv.tube(i,0);
    mat_to_return.row(i) = temp_row;
  }
  return(mat_to_return);
}



// ------------------------------------------------------------------//
// --------------vector Magnitude generator fucntion----------------//
// ----------------------------------------------------------------//


double vec_mag(rowvec pos_vector){
  double dim_vec = pos_vector.n_elem, mag = 0;
  for(int i=0;i<dim_vec;++i){
    mag = mag + pow(pos_vector(i),2);
  }
  return(mag);
}


// ------------------------------------------------------------------//
// -----------------Mesh Genrator Cube function---------------------//
// ----------------------------------------------------------------//




cube grid_cube_generator(rowvec x, rowvec y){
  int x_axis_size = x.n_elem, y_axis_size = y.n_elem;
  cube Mesh_cube(y_axis_size,x_axis_size,3,fill::zeros);
  for(int i = 0;i < y_axis_size; ++i){
    for(int j = 0; j < x_axis_size; ++j){
      Mesh_cube(i,j,0) = x[j];
      Mesh_cube(i,j,1) = y[i];
    }
  }
  return(Mesh_cube);
}

mat grid_target_generator(rowvec x, rowvec y){
  int x_axis_size = x.n_elem, y_axis_size = y.n_elem;
  mat Mesh_cube(y_axis_size * x_axis_size,3,fill::zeros);
  for(int i = 0;i < x_axis_size ; ++i){
    for (int j = 0; j < y_axis_size ; ++j) {
      Mesh_cube(i*y_axis_size+j,0) = x[i];
      Mesh_cube(i*y_axis_size+j,1) = y[j];
    }
  }
  return(Mesh_cube);
}



// ------------------------------------------------------------------//
// --------------Source Modelling Functions-------------------------//
// ----------------------------------------------------------------//



// Line Maker Funcion

mat line_generator(rowvec start, rowvec end , int no_divisions = 20){
  rowvec pos_vector = end - start;
  double mag_of_vec = vec_mag(pos_vector);
  rowvec unit_pos_vec = pos_vector/mag_of_vec;
  double d_x = mag_of_vec/no_divisions;
  
  mat generated_cords(no_divisions+1,3,fill::zeros);
  for(int i = 0; i <= no_divisions; ++i){
    generated_cords.row(i) = start + i * d_x  * unit_pos_vec;
  }

  return(generated_cords);
}


// Circle maker Fucntion

std::tuple<mat,mat> circle_maker(double d_r,double radii, vec origin_circ , rowvec normal){

  rowvec n = normal/norm(normal);
  rowvec y_axis_vector = {0,1,0};
  rowvec l = cross(normal,y_axis_vector);
  l = l/norm(l);
  rowvec m = cross(n, l);

  // Calculating the Rotation matrix

  mat rot;

  rot.insert_rows(0,l);
  rot.insert_rows(1,m);
  rot.insert_rows(2,n);
  

  
  mat pos, norm ;
  int row_counter_pos = 0;
  int row_counter_norms = 0;

  vec loop_radians = regspace<vec>(d_r,d_r,radii);

  for(int i = 0; i<loop_radians.n_elem;++i){

    double a_n = std::ceil(2*M_PI*loop_radians[i]/d_r);
    vec no_division = regspace<vec>(0,1,a_n-1);

    for(int j = 0;j < no_division.n_elem ; ++j){

      vec dpos = {std::cos(2*M_PI*no_division[j]/a_n),std::sin(2*M_PI*no_division[j]/a_n),0} ;
      dpos = dpos * loop_radians[i];

      
      vec npos = origin_circ +  rot.t() * dpos ;
     
      rowvec npos_rowvec =  conv_to<rowvec>::from(npos);
      pos.insert_rows(row_counter_pos,npos_rowvec);
      row_counter_pos++;
      norm.insert_rows(row_counter_norms,normal);  // this shit is messes up need to Rethink it
      row_counter_norms++;
      
    }
  }
  return std::make_tuple(pos,norm);
}



mat rectangle_generator(rowvec x, rowvec y , rowvec o , double x_div = 1000 , double y_div = 1000){
  rowvec x_dir_rat = x - o ;
  rowvec x_dir_cos = x_dir_rat  / vec_mag(x_dir_rat);
  rowvec y_dir_rat = y - o ;
  rowvec y_dir_cos = y_dir_rat  / vec_mag(y_dir_rat);
  
  mat Mesh_cube(x_div * y_div ,3,fill::zeros);
  for(int i = 0;i < y_div; ++i){
    Mesh_cube(span(i*x_div,(i+1)*x_div-1),span::all) =  line_generator(o + i * y_dir_cos, o + i * y_dir_cos + x_div * x_dir_cos , x_div-1);

  }
  return(Mesh_cube);
}

mat source_point_placer(mat source_location_mat, rowvec normal_location, double radiii){
  int no_of_sources = (int)source_location_mat.n_rows;
  mat displaced_sources(no_of_sources,3,fill::zeros);
#pragma omp parallel for
  for (int i = 0; i < no_of_sources; ++i) { // openMP is Crying if we compare the double and int in a loop so keep that in mind
   displaced_sources.row(i) = source_location_mat.row(i) - radiii * normal_location;
  }
  
  return(displaced_sources);
}




// ------------------------------------------------------------------//
// --------------Stress Saving Function File------------------------//
// ----------------------------------------------------------------//
template <typename T>
int save_cx_3d(std::vector<Mat<T>> mat_to_be_saved, const std::string &PATH_TO_SAVE){
  
  int size_of_cx_3d = mat_to_be_saved.size();

  std::vector<ofstream> outdata(9);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0 ; i < 3; ++i) {
      std::string File_prefix = "Sigma-";
      std::string File_suffix_i = std::to_string(i);
      std::string File_suffix_j = std::to_string(j);
      std::string File_Format = ".csv";
      std::string File_Name = PATH_TO_SAVE + File_prefix + File_suffix_i + File_suffix_j + File_Format;
      outdata[i*3+j].open(File_Name);
      if( !outdata[j*3+i] ) { //file couldn't be opened
	cerr << "Error: file could not be opened" << endl;
	exit(1);
      }
    }
  }


// parallel Writing to The Disk
#pragma omp parallel
  {
    for (int k=0; k < 3; ++k) {
      for (int j=0; j < 3; ++j) {
#pragma omp single
	for (int i = 0; i < size_of_cx_3d; ++i) {
	  outdata[k*3+j] << std::abs(mat_to_be_saved[i](0,0)) << endl;  
	}
      }
    }
  }

  for (int i = 0; i < 9; ++i) {
    outdata[i].close();
  }
  return 0;
}

template int save_cx_3d<cx_double>(std::vector<Mat<cx_double>>, const std::string &);
template int save_cx_3d<double>(std::vector<Mat<double>>, const std::string &);
template int save_cx_3d<float>(std::vector<Mat<float>>, const std::string &);

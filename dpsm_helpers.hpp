#include "data_dec.hpp"


// ------------------------------------------------------------------//
//----- Note : You Need The Default argument values in Header only--//
// ----------------------------------------------------------------//



// ------------------------------------------------------------------//
// --------------Basic Arthemetic on Defined Fucntion---------------//
// ----------------------------------------------------------------//
template <typename T>
cx_3d<T> add_cx_3d (cx_3d<T> former, cx_3d<T> latter){
  // No Explicit Checking is Defined and Dev is asssumed to  Use fucntions sane manner
  int lend_cx_3d = former.size();
  cx_3d<T> ret_data(lend_cx_3d,Mat<std::complex<T>>(size(former[0])));
  for(int i = 0; i < lend_cx_3d; i++){
    ret_data[i] = former[i] + latter[i];
  }
  return(ret_data);
}



template <typename T>
Mat<std::complex<T>> conv_cube_2_mat(Cube<std::complex<T>> to_be_conv){
  int n_row_cube = to_be_conv.n_rows;
  int n_slice_cube = to_be_conv.n_slices;
  Mat<std::complex<T>> mat_to_return(n_row_cube,n_slice_cube,fill::zeros);
  
  for(size_t i = 0; i < n_row_cube; i++){
    Row<std::complex<T>> temp_row = to_be_conv.tube(i,0);
    mat_to_return.row(i) = temp_row;
  }
  return(mat_to_return);
}


// ------------------------------------------------------------------//
// --------------vector Magnitude generator fucntion----------------//
// ----------------------------------------------------------------//

template <typename T>
T vec_mag(Row<T> pos_vector){
  T  mag = 0;
  size_t dim_vec = pos_vector.n_elem;
  for(size_t i=0;i<dim_vec;++i){
    mag = mag + pow(pos_vector(i),2);
  }
  return(mag);
}



// ------------------------------------------------------------------//
// -----------------Mesh Genrator Cube function---------------------//
// ----------------------------------------------------------------//


template <typename T>
Mat<T> grid_target_generator(Row<T> x, Row<T> y){
  size_t x_axis_size = x.n_elem, y_axis_size = y.n_elem;
  Mat<T> Mesh_cube(y_axis_size * x_axis_size,3,fill::zeros);
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
template<typename T>
Mat<T> line_generator(Row<T> start, Row<T> end , size_t no_divisions){
  Row<T> pos_vector = end - start;
  T mag_of_vec = vec_mag<T>(pos_vector);
  Row<T> unit_pos_vec = pos_vector/mag_of_vec;
  T d_x = mag_of_vec/no_divisions;
  
  Mat<T> generated_cords(no_divisions+1,3,fill::zeros);
  for(size_t i = 0; i <= no_divisions; ++i){
    generated_cords.row(i) = start + i * d_x  * unit_pos_vec;
  }
  return(generated_cords);
}



// Circle maker Fucntion


template<typename T>
Mat<T> circle_maker(T d_r,T radii, Col<T> origin_circ , Row<T> normal){

  Row<T> n = normal/norm(normal);
  Row<T> y_axis_vector = {0,1,0};
  Row<T> l = cross(normal,y_axis_vector);
  l = l/norm(l);
  Row<T> m = cross(n, l);

  // Calculating the Rotation matrix

  Mat<T> rot;

  rot.insert_rows(0,l);
  rot.insert_rows(1,m);
  rot.insert_rows(2,n);
  
  Mat<T> pos;
  int row_counter_pos = 0;

  Col<T> loop_radians = regspace<Col<T>>(d_r,d_r,radii);

  for(size_t i = 0; i<loop_radians.n_elem;++i){

    T a_n = std::ceil(2*M_PI*loop_radians[i]/d_r);
    Col<T> no_division = regspace<Col<T>>(0,1,a_n-1);

    for(size_t j = 0;j < no_division.n_elem ; ++j){

      Col<T> dpos = {T(std::cos(2*M_PI*no_division[j]/a_n)),T(std::sin(2*M_PI*no_division[j]/a_n)),0} ;
      dpos = dpos * loop_radians[i];
      Col<T> npos = origin_circ +  rot.t() * dpos ;
      Row<T> npos_rowvec =  conv_to<Row<T>>::from(npos);
      pos.insert_rows(row_counter_pos,npos_rowvec);
      row_counter_pos++;
    }
  }
  return pos;
}


/* Input Explanation 
The Input to the rectangle generator function 
o - Corner
x - Horizantal directional vector
y - Vertical Directional Vector
y_div - No of Divisions in Y directions
x_div - No of Divisions in X Directions

ASCII ART Below
==============================
y
|
|
|y_div
|
|
|
|
o---------------------------x
       x_div
=============================== 
 */

template <typename T>
Mat<T> rectangle_generator(Row<T> dir_vec_1, Row<T> dir_vec_2 , Row<T> origin , int x_div , int y_div){

    Row<T> x_dir_rat =  dir_vec_1- origin ;
    Row<T> x_dir_cos = x_dir_rat / (x_div-1);
    Row<T> y_dir_rat = dir_vec_2 - origin ;
    Row<T> y_dir_cos = y_dir_rat / (y_div-1);
    cout << "X Divisions--"<< x_div << "========" <<" Y Divisions--"<< y_div << endl;
    Mat<T> Rec_Point_Mat(x_div  * y_div ,3,fill::zeros);

    for(int i = 0;i < y_div; ++i){
        Rec_Point_Mat(span(i*x_div,(i+1)*x_div-1),span::all) =  line_generator<T>(origin + i * y_dir_cos  , origin + i * y_dir_cos + x_dir_rat , x_div-1);

    }
    return(Rec_Point_Mat);
}



// ------------------------------------------------------------------//
// --------------Source Point_placer--------------------------------//
// ----------------------------------------------------------------//

template<typename T>
Mat<T> source_point_placer(Mat<T> source_location_mat, Row<T> normal_location, T radiii){
  size_t no_of_sources = source_location_mat.n_rows;
  Mat<T> displaced_sources(no_of_sources,3,fill::zeros);
#pragma omp parallel for
  for (size_t i = 0; i < no_of_sources; ++i) { // openMP is Crying if we compare the double and int in a loop so keep that in mind
   displaced_sources.row(i) = source_location_mat.row(i) + radiii * normal_location;
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
	  //outdata[k*3+j] << std::abs(mat_to_be_saved[i](0,0)) << endl;  
    outdata[k*3+j] << mat_to_be_saved[i](k,j)<< endl;  
	}
      }
    }
  }

  for (int i = 0; i < 9; ++i) {
    outdata[i].close();
  }
  return 0;
}



template <typename T>
Mat<T> CUBE_EXTRACT(Cube<T> incoming_tube,int row_no_ex){
    Mat<T> Extracnted(3,3,fill::zeros);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Extracnted(i,j) = incoming_tube(row_no_ex,j,i);
        }
    }
    return Extracnted;
};


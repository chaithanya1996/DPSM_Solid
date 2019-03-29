#include "data_dec.hpp"


// ------------------------------------------------------------------//
// --------------Basic Arthemetic on Defined Fucntion---------------//
// ----------------------------------------------------------------//

template <typename T>
cx_3d<T> add_cx_3d (cx_3d<T> former, cx_3d<T> latter);
  
template <typename T>
Mat<std::complex<T>> conv_cube_2_mat(Cube<std::complex<T>> to_be_conv);


// ------------------------------------------------------------------//
// --------------Source Modelling Functions-------------------------//
// ----------------------------------------------------------------//

// Line Maker Funcion
template<typename T>
Mat<T> line_generator(Row<T> start, Row<T> end , size_t no_divisions = 20);

// Circle maker Fucntion
template<typename T>
Mat<T> circle_maker(T d_r,T radii, Col<T> origin_circ , Row<T> normal);

// rectangle generator
template <typename T>
Mat<T> rectangle_generator(Row<T> dir_vec_1, Row<T> dir_vec_2 , Row<T> origin , int x_div , int y_div);

// ------------------------------------------------------------------//
// -----------------Mesh Genrator Cube function---------------------//
// ----------------------------------------------------------------//

template <typename T>
Mat<T> grid_target_generator(Row<T> x, Row<T> y);

// ------------------------------------------------------------------//
// --------------vector Magnitude generator fucntion----------------//
// ----------------------------------------------------------------//

template <typename T>
T vec_mag(Row<T> pos_vector);

// ------------------------------------------------------------------//
// --------------Source Point_placer--------------------------------//
// ----------------------------------------------------------------//

template<typename T>
Mat<T> source_point_placer(Mat<T> source_location_mat, Row<T> normal_location, T radiii);

// ------------------------------------------------------------------//
// --------------Save Stress To Disk Function-----------------------//
// ----------------------------------------------------------------//

template <typename T>
int save_cx_3d(std::vector<Mat<T>> mat_to_be_saved, const std::string &PATH_TO_SAVE);

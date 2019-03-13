#include "data_dec.hpp"


cx_3d add_cx_3d (cx_3d former, cx_3d latter);
cx_mat conv_cube_2_mat(cx_cube to_be_conv);
mat rectangle_generator(rowvec x, rowvec y , rowvec o , double x_div , double y_div);
std::tuple<mat,mat> circle_maker(double d_r,double radii, vec origin_circ , rowvec normal);
mat line_generator(rowvec start, rowvec end , int no_divisions );
cube grid_cube_generator(rowvec x, rowvec y);
mat grid_target_generator(rowvec x, rowvec y);
double vec_mag(rowvec pos_vector);
mat source_point_placer(mat source_location_mat, rowvec normal_location, double radiii);
template <typename T>
int save_cx_3d(std::vector<Mat<T>> mat_to_be_saved, const std::string &PATH_TO_SAVE);

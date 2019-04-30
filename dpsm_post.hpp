#include "dpsm_core.hpp"
#include<string>

using std::string;
using std::to_string;
using std::vector;

template<typename T>
int DPSM_POST_PROCESS(int x_vec_divs, int y_vec_divs,Mat<T> GEO_MAT, Mat<T> TARGET_MAT,string PATH_TOSAVE,string FILE_NAME_TO_SAVE){

  Mat<T> MESH_GRID_X(y_vec_divs,x_vec_divs ,fill::zeros);
  Mat<T> MESH_GRID_Y(y_vec_divs,x_vec_divs ,fill::zeros);
  Mat<T> MESH_GRID_Z(y_vec_divs,x_vec_divs ,fill::zeros);

  for (int i= 0; i < y_vec_divs; ++i) {
    MESH_GRID_X.row(i) = trans(GEO_MAT(span(i*x_vec_divs,(i+1) * x_vec_divs - 1),0));
    MESH_GRID_Y.row(i) = trans(GEO_MAT(span(i*x_vec_divs,(i+1) * x_vec_divs -1 ),1));
    MESH_GRID_Z.row(i) = trans(GEO_MAT(span(i*x_vec_divs,(i+1) * x_vec_divs -1 ),2));
  }

  string X_NAME = PATH_TOSAVE + FILE_NAME_TO_SAVE + "MESH_GRID_X.csv";
  string Y_NAME = PATH_TOSAVE + FILE_NAME_TO_SAVE + "MESH_GRID_Y.csv";
  string Z_NAME = PATH_TOSAVE + FILE_NAME_TO_SAVE + "MESH_GRID_Z.csv";
  
  MESH_GRID_X.save(X_NAME,csv_ascii);
  MESH_GRID_Y.save(Y_NAME,csv_ascii);
  MESH_GRID_Z.save(Z_NAME,csv_ascii);
  
  int no_of_target_to_save = TARGET_MAT.n_cols;
  
  for (int i = 0 ; i < no_of_target_to_save; ++i) {
    
    // Creating Placeholder For Target Values
    Mat<T> TARGET_VAL(y_vec_divs,x_vec_divs,fill::zeros);
    Col<T> COL_vals = TARGET_MAT.col(i);
    for (int j= 0; j < y_vec_divs; ++j) {
      TARGET_VAL.row(j) = trans(COL_vals(span(j*x_vec_divs,(j+1) * x_vec_divs - 1)));
    }
    string TARGET_NAME = PATH_TOSAVE + FILE_NAME_TO_SAVE + "TARGET-" + to_string(i)+".csv";
    TARGET_VAL.save(TARGET_NAME,csv_ascii);
  }
  return(1);
}


template<typename T>
Mat<T> STRESS_POST_PROCESS(vector<Mat<T>> STRESS_VEC){
  int vec_size = STRESS_VEC.size();
  Mat<T> OUTPUT(vec_size,9,fill::zeros);
  for (int i = 0; i < vec_size; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3 ; ++k) {
	OUTPUT(i,j*3+k) = STRESS_VEC[i](j,k);
      }
    }
  }
  return(OUTPUT);
}

#include<armadillo>
using namespace arma;
int main(){
  cx_mat Displ;
  Displ.load("DBUG_DISP.csv",csv_ascii);
  mat Displ_abs = abs(Displ);
  Displ_abs.save("DBUG_DISP_single_val.csv",csv_ascii);
}

#include<armadillo>
#include<complex>
#include<vector>
#include <string>

#include <cstdlib>
#include <fstream>



using namespace arma;


using std::cerr;
using std::cout;
using std::endl;
using std::size_t;
using std::ofstream;
using std::complex;
typedef std::complex<double> cx_db; // Redundant Stuff

// Implemting Templated Alias Type Def Varibles 

/* ---------------------------------

typedef std::complex<double> cx_db;
typedef std::vector<cx_mat> cx_3d;
typedef std::vector<cx_cube> cx_4d;
typedef std::vector<cx_3d> cx_4d_mat;
typedef std::vector<cx_4d> cx_5d;;

---------------------------------*/

template<typename T>
using cx_3d = std::vector<Mat<std::complex<T>>>;

template<typename T>
using cx_4d = std::vector<Cube<std::complex<T>>>;

template<typename T>
using cx_4d_mat = std::vector<cx_3d<T>>;

template<typename T>
using cx_5d = std::vector<cx_4d<T>>;

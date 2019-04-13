#include "dpsm_core.hpp"

namespace dpsm{
  template<class T>
  class solid_material{
  public:
    T mu, lamda;
    //! Default constructor
    solid_material(){};
    solid_material(T mu, T lamda);
    solid_material(const solid_material &other);
    
  };

  template<class T>
  class transducer{
  public:
    Mat<T> transducer_location;
    Row<T> normal_to_transducer_storage;
    transducer(){};
    transducer(Mat<T> transducer_location_rowvec , Row<T> normal_to_transducer);
    
  };

  /* The Interface Will have Two sides to the Normal created even though one normal is defined 
     The Other normal is taken as the opposite direction by default as this makes sense
  */
  template<class T>
  class solid_interface{
  public:
    Mat<T> interface_source_1;
    Row<T> normal_to_solid_interface;
    //! Default constructor
    solid_interface(){};
    solid_interface(Mat<T> interface_location, Row<T> normal_to_solid_interf);
  };


  template<class T>
  class solid_ply{
  public:
    Mat<T> active_sources,passve_sources;
    //mat mesh_location;
    solid_material<T> mat_prop;
    //! Default constructor
    // This Contructoe is Using member list initialization Check to clear Doubts
    solid_ply(){};
    
    // solid_ply (mat active_sources,mat passve_sources, mat mesh_location,solid_material mat_prop)
    //   : active_sources(active_sources), passve_sources(passve_sources) , mesh_location (mesh_location) , mat_prop(mat_prop) {};
    solid_ply (Mat<T> active_sources,Mat<T> passve_sources,solid_material<T> mat_prop)
      : active_sources(active_sources), passve_sources(passve_sources)  , mat_prop(mat_prop) {};
    // Copy Constructor of Solid Ply

    //solid_ply(const solid_ply &other):solid_ply(other.active_sources,other.passve_sources,other.mesh_location,other.mat_prop){};
    solid_ply(const solid_ply<T> &other):solid_ply(other.active_sources,other.passve_sources,other.mat_prop){};
    solid_ply &operator= (const solid_ply<T> &other){
      active_sources = other.active_sources;
      passve_sources = other.passve_sources;
      mat_prop = other.mat_prop;
      return *this;
    };
    // Need a Copy Constructor here
  };



  /*
    
   */

  template<class T>
  class geometry{
  public:
    
    // Member Varibles 
    T no_of_slices;
    T r_s;
    std::vector<transducer<T>> transducer_array;
    std::vector<solid_material<T>> material_array;
    std::vector<solid_interface<T>> interface_array;

    
    //std::vector<solid_ply> solid_lamina_array;
    //std::vector<rowvec> normal_interface_vectors;
    //mat active_source_locations;
    //mat passive_source_location;
    
    //! Default constructor
    
    /* Assumptions made for the Composite material
       - The Composite is Continuos 
       - The radius id extremely less than the Dimentions supplied and assumed to be correctly while receiving
       - The Normal is same across the material but for sake of extendibility we assume that the 
    */
    geometry(double slice_count,const std::vector<solid_material<T>> &mat_arr, const std::vector<transducer<T>> &trasducer_arr_location,
	     const std::vector<solid_interface<T>> &interface__arr,T radius)
      : no_of_slices(slice_count) , material_array(mat_arr) , transducer_array(trasducer_arr_location) , interface_array(interface__arr) , r_s(radius){}
    
    
    std::vector<solid_ply<T>> get_the_material_grid();  
  };

  
};


template<typename T>
dpsm::solid_interface<T> combine_interface(dpsm::solid_interface<T> interface1,dpsm::solid_interface<T> interface2){
  Mat<T> joined_mats = join_cols(interface1.interface_source_1,interface2.interface_source_1);
  dpsm::solid_interface <T> joint_interface_solid(joined_mats,interface1.normal_to_solid_interface);
  return(joint_interface_solid);
};

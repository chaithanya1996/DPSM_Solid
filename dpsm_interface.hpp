#include "dpsm_core.hpp"

namespace dpsm{

  class solid_material{
  public:
    double mu, lamda;
    //! Default constructor
    solid_material(){};
    solid_material(double mu, double lamda);
    solid_material(const solid_material &other);
    
  };


  class transducer{
  public:
    mat transducer_location;
    rowvec normal_to_transducer_storage;
    transducer(){};
    transducer(mat transducer_location_rowvec , rowvec normal_to_transducer);
    
  };

  /* The Interface Will have Two sides to the Normal created even though one normal is defined 
     The Other normal is taken as the opposite direction by default as this makes sense
  */
  class solid_interface{
  public:
    mat interface_source_1;
    rowvec normal_to_solid_interface;
    //! Default constructor
    solid_interface(){};
    solid_interface(mat interface_location, rowvec normal_to_solid_interf);
  };


  
  class solid_ply{
  public:
    mat active_sources,passve_sources;
    //mat mesh_location;
    solid_material mat_prop;
    //! Default constructor
    // This Contructoe is Using member list initialization Check to clear Doubts
    solid_ply(){};
    
    // solid_ply (mat active_sources,mat passve_sources, mat mesh_location,solid_material mat_prop)
    //   : active_sources(active_sources), passve_sources(passve_sources) , mesh_location (mesh_location) , mat_prop(mat_prop) {};
    solid_ply (mat active_sources,mat passve_sources,solid_material mat_prop)
      : active_sources(active_sources), passve_sources(passve_sources)  , mat_prop(mat_prop) {};
    // Copy Constructor of Solid Ply

    //solid_ply(const solid_ply &other):solid_ply(other.active_sources,other.passve_sources,other.mesh_location,other.mat_prop){};
    solid_ply(const solid_ply &other):solid_ply(other.active_sources,other.passve_sources,other.mat_prop){};
    // Need a Copy Constructor here
  };



  /*
    
   */

  
  class geometry{
  public:
    
    // Member Varibles 
    double no_of_slices;
    double r_s;
    std::vector<transducer> transducer_array;
    std::vector<solid_material> material_array;
    std::vector<solid_interface> interface_array;

    
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
    geometry(double slice_count,const std::vector<solid_material> &mat_arr, const std::vector<transducer> &trasducer_arr_location,
	     const std::vector<solid_interface> &interface__arr,double radius)
      : no_of_slices(slice_count) , material_array(mat_arr) , transducer_array(trasducer_arr_location) , interface_array(interface__arr) , r_s(radius){}
    
    
    std::vector<solid_ply> get_the_material_grid();  
  };

  
};

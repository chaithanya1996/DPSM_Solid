#include "dpsm_interface.hpp"

// Constructor For Solid Material
dpsm::solid_material::solid_material(double mu_val, double lamda_val){
  mu = mu_val;
  lamda = lamda_val;
}

// Solid Material Copy Constructor
dpsm::solid_material::solid_material(const solid_material &other){
  mu = other.mu;
  lamda = other.lamda;
}

// Transducer Constructor
dpsm::transducer::transducer(mat transducer_location_rowvec,rowvec normal_to_transducer){
  transducer_location = transducer_location_rowvec;
  normal_to_transducer_storage = normal_to_transducer;
}

dpsm::solid_interface::solid_interface(mat interface_location,rowvec normal_to_solid_interf){
  interface_source_1 = interface_location;
  normal_to_solid_interface = normal_to_solid_interf;
}
/* The class performs action acting in good faith of the Fixed Co ordinate system not a generic one
   we will Extend this in the Next iteration when things are working perfectly  */

// dpsm::geometry::geometry(double slice_count,const std::vector<solid_material> &mat_arr, const std::vector<transducer> &trasducer_arr_location, const std::vector<solid_interface> &interface__arr,double radius): no_of_slices(slice_count){
//   solid_lamina_array = new std::vector<solid_ply> (slice_count,solid_ply);
  
// }


// Solid Ply Copy Constructor to Enable Working In Geometry
// Solid Ply Deafult constructor is added to facilitate uninitilizatied declaration of solid ply in the geometry

/*----------------------------------------------
Steps To Follow to create a Composite Hetrogeneous Medium
-----------------------------------------------
1. check Transducer is present in the material if present inset into active sources
2. Calculate the PAssive Source 
------------------------------------------------ */

std::vector<dpsm::solid_ply> dpsm::geometry::get_the_material_grid(){
  std::vector<dpsm::solid_ply> resulting_geometry(no_of_slices,dpsm::solid_ply());
  for (int i = 0 ; i < no_of_slices; ++i) {
    
    if (!transducer_array[i].transducer_location.is_empty()) {
      // std::cout << " Transducer -" << i << endl;
      // std::cout << transducer_array[i].transducer_location << std::endl;
      // std::cout << transducer_array[i].normal_to_transducer_storage << std::endl;
      resulting_geometry[i].active_sources =  source_point_placer(transducer_array[i].transducer_location, transducer_array[i].normal_to_transducer_storage,r_s);
    }
    // std::cout << " Interfaec -" << i << endl;
    // std::cout << interface_array[i].interface_source_1 << std::endl;
    // std::cout << interface_array[i].normal_to_solid_interface << std::endl;
      
    resulting_geometry[i].passve_sources = source_point_placer(interface_array[i].interface_source_1, interface_array[i].normal_to_solid_interface,r_s);
    resulting_geometry[i].mat_prop = material_array[i];
  }
  return(resulting_geometry);
}

// This Only Works for the Solid Composite case as the Interface is similar direction
// Ideally we would need the get fucntion to handle this situation

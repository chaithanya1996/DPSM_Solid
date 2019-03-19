#include "dpsm_interface.hpp"

// Constructor For Solid Material
template<typename T>
dpsm::solid_material<T>::solid_material(T mu_val, T lamda_val){
  mu = mu_val;
  lamda = lamda_val;
}


// Solid Material Copy Constructo
template<typename T>
dpsm::solid_material<T>::solid_material(const solid_material<T> &other){
  mu = other.mu;
  lamda = other.lamda;
}

// Transducer Constructor
template<typename T>
dpsm::transducer<T>::transducer(Mat<T> transducer_location_rowvec,Row<T> normal_to_transducer){
  transducer_location = transducer_location_rowvec;
  normal_to_transducer_storage = normal_to_transducer;
}

template<typename T>
dpsm::solid_interface<T>::solid_interface(Mat<T> interface_location,Row<T> normal_to_solid_interf){
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
template<typename T>
std::vector<dpsm::solid_ply<T>> dpsm::geometry<T>::get_the_material_grid(){
  
  std::vector<dpsm::solid_ply<T>> resulting_geometry(no_of_slices,dpsm::solid_ply<T>());
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


// Explicit Instantiations

template class dpsm::solid_material<float>;
template class dpsm::solid_material<double>;

template class dpsm::transducer<float>;
template class dpsm::transducer<double>;

template class dpsm::solid_interface<float>;
template class dpsm::solid_interface<double>;

template class dpsm::solid_ply<float>;
template class dpsm::solid_ply<double>;

template class dpsm::geometry<float>;
template class dpsm::geometry<double>;

#ifndef _MODELSTRUCTURE_H_
#define _MODELSTRUCTURE_H_


#include <iostream>
#include <fstream>

#include <map>
#include <vector>
#include <string>
#include "TypeDefs.h"
#include "ModelOptions.h"
#include "RegionStructure.h"



//!  Class  containing  input-defined model features.
/*!
 * Contains  name, associated physical region IDs, associated boundary conditions 
 * region IDs and  their definition.
 */

class ModelStructure{


 public:

  //!  Constructor. 
  /*!
   * Creates a new structure associated to the model "model_name".
   */
  ModelStructure( const std::string& model_name);

  //!  Destructor 
  /*!
   *
   */
  ~ModelStructure();



  //!  Gets a map with Boundary Conditions definitions associated to each BC ID. 
  /*!
   * Returns a map which associates Boundary Conditions IDs with the a RegionStructure
   * containing description  of
   * each   Boundary Condition  region.
   */
  std::map <ID, RegionStructure>& get_model_BC_map();

  //!  Gets a map  with physical model info 
  /*!
   * Returns a map which associates physical model name  with properties contained in a 
   * ModelOptions object
   */
  std::map<const std::string,ModelOptions>& get_physical_model_map();



  //!  Sets an internal map with Boundary Conditions definitions associated to each BC ID. 
  /*!
   * Creates  a map which associates Boundary Conditions IDs with the RegionStructure 
   * of the BC.
   */
  void set_model_BC_map( std::map <ID, RegionStructure>& mod_BC_map);


  //!  Sets an internal map with physical model info 
  /*!
   * Creates  a map which associates a physical model name  with a ModelOptions 
   * object containing its  properties. 
   */
  void set_physical_model_map(  std::map <const std::string,ModelOptions>& phys_model_map );



  //!  Gets a  vector with the physical regions IDs of this model. 
  /*!
   * Returns a vector with the physical regions IDs   associated to the  current model.
   */
  std::vector<std::string> get_physical_regions();


  //!  Sets a  vector with the physical regions IDs of this model. 
  /*!
   * Creates  a vector with the physical regions IDs   associated to the  current model.
   */
  void set_physical_regions( std::vector<std::string>& list_phys_regions);

 
  //!  Gets a  string with the name of the model.
  std::string  get_model_name();

  //----------------------------------------------------------------------------------


 private:

  //! Model name.
  std::string model_name;

  //!  Physical regions associated to this  model
  std::vector<std::string>  physical_regions;   

  /*! Map which associates Boundary Conditions IDs with the RegionStructure  of
   * the   Boundary Condition  region.
   */
  std::map <ID, RegionStructure> model_BC_map;





  /*! Map of  the  physical models associated with the  current simulation model:
   *  physical model name is  associated with a ModelOptions object.
   */
  std::map <const std::string,ModelOptions>  physical_model_map ; 


  //!  Sets  the name of the model.
  void set_model_name(const std::string& model);


};

#endif // endif define   _INPUTPARSER_H_

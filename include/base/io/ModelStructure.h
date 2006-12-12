#ifndef _MODELSTRUCTURE_H_
#define _MODELSTRUCTURE_H_


#include <iostream>
#include <fstream>

#include <map>
#include <vector>
#include <string>
#include "TypeDefs.h"

using namespace std;

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
  ModelStructure(string& model_name);

  //!  Destructor 
  /*!
   *
   */
  ~ModelStructure();



  //!  Gets a map with Boundary Conditions definitions associated to each BC ID. 
  /*!
   * Returns a map which associates Boundary Conditions IDs with the property map of
   * each   Boundary Condition  region.
   */
  map  <ID,  map <string,string> >&  get_model_BC_map();


  //!  Sets an internal map with Boundary Conditions definitions associated to each BC ID. 
  /*!
   * Creates  a map which associates Boundary Conditions IDs with the property map of
   * each   Boundary Condition  region.
   */
  void set_model_BC_map( map  <ID,  map <string,string> >& id_BC_regions_map   );


  //!  Gets a  vector with the physical regions IDs of this model. 
  /*!
   * Returns a vector with the physical regions IDs   associated to the  current model.
   */
  vector<string> get_physical_regions();


  //!  Sets a  vector with the physical regions IDs of this model. 
  /*!
   * Creates  a vector with the physical regions IDs   associated to the  current model.
   */
  void set_physical_regions( vector<string>& list_phys_regions);

 
  //!  Gets a  string with the name of the model.
  string  get_model_name();

  //----------------------------------------------------------------------------------


 private:

  //! Model name.
  string model_name;

  //!  Physical regions associated to this  model
  vector<string>  physical_regions;   

  /*! Map which associates Boundary Conditions IDs with the property map of
   * each   Boundary Condition  region.
   */
  map  <unsigned int,  map <string,string>   >  model_BC_map;

  //!  Sets  the name of the model.
  void set_model_name(string& model);


};

#endif // endif define   _INPUTPARSER_H_

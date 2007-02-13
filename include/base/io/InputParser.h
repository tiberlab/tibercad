/*=============================================================================
  Copyright (c) 2002-2003 Joel de Guzman
  http://spirit.sourceforge.net/

  Use, modification and distribution is subject to the Boost Software
  License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
  http://www.boost.org/LICENSE_1_0.txt)
  =============================================================================*/


#ifndef _INPUTPARSER_H_
#define _INPUTPARSER_H_

#include <iostream>
#include <fstream>

#include <map>
#include <vector>
#include <string>

#include "TypeDefs.h"
#include "ModelStructure.h"
#include "RegionStructure.h"
#include "ModelOptions.h"



//!  A parser  for  TIBERCAD input  text  file. 
/*!
 * Parses an input text file composed by sections:  each  section begins with a name  
 * preceded by "$" (e.g. $Physics ) and it is composed by a block between "{" and  "}" 
 * parenthesis.
 * This block can be possibly composed by one or more subblocks, each preceded  by a model name.
 * Two special sections are "Device" and "Models"  sections: in  "Device" section, each subblock
 * must be preceded by the keyword "Region", followed by the name of the physical region.
 * In  "Models"  section, one or  more special model-blocks must be  present: each model-block
 * must be preceded by the keyword "Model", followed by the model name. Each model-block must 
 * contain the tagname "phys_regions = ( list of physical regions of the model)" and one 
 * BC_regions-block. The BC_regions-block must be preceded by the  keyword  "BC_Regions" and 
 * it is  composed by one or more subblocks,
 * each  preceded by  the keyword "BC_Region" followed by the name of the boundary condition region.
 * "Models"  section can  contain also  one or more optionals physical model blocks:
 * each physical model-block must be preceded by the keyword "physical_model", followed by 
 * the name of physical model.
 * 
 *  
 * Except the model-blocks and the BC_regions-blocks, each block and  subblock can 
 *  contain zero or  any number of  parameters in the assignement form "tagname = tagvalue", 
 *  where "tagname" is  a  string and "tagvalue" is  a single 
 * numerical or string item or a list of  items between  "(" and ")" parenthesis.
 * Format is free for these assignements, provided they are separated by spaces. 
 * Everything following a '#' is  a  comment and  is  disregarded.
 *
 * Public method read_device()  parses and extracts information 
 * from device section.
 * read_models()  method  parses and extracts information from Models section.
 *  Models data are stored in the class \c ModelStructure.
 * Public methods read_parameters()   parse and extract 
 * information from the other parameters sections.
 */


class InputParser{


 public:


  //!  Constructor 
  /*!
    Defines  name  of input  file ("input_file_name").
  */
  InputParser(const std::string& input_file_name);

  //!  Destructor 
  /*!
  
  */
  ~InputParser(void);




  //!   Parses the  device section of  input  file 
  /*!
   * Method to  read   the  device section of  input  file: returns a map 
   * between region ID and its RegionStructure.
   */
  const std::map <ID, RegionStructure>& read_device(void);


  //!  Parses the "Models"  section  of the  input file. 
  /*!
   * Method to  parse the "Models"  section  of the  input file.
   *  Returns a map between Model name and a pointer to ModelStructure 
   * associated to the model. 
   */
  const std::map <const std::string, ModelStructure*>& read_models(void);

 


  //!   Reads parameters of  a  given  section   for a given  model. 
  /*!
   * Method to   read  parameters of the  section "section_name"  for a given model "model_name".
   * Returns parameters in a ModelOptions object. If "model_name" section is absent, returns 
   *  am empty  ModelOptions.
   */
  const  ModelOptions&  read_parameters(std::string section_name,const std::string& model_name);

 
  //!   Reads parameters of  a  given  section  
  /*!
   * Overloaded method to   read  model-independent parameters of the  section "section_name".
   * Returns parameters in a ModelOptions object. 
   */
  const  ModelOptions& read_parameters(std::string section_name);





 

  //----------------------------------------------------------------------------------

 private:

  //  members data:

  //   control  chars:

  /*!
   *  Char to  start a block ("{").
   */
  std::string  start_symb ;


  /*!
   *  Char to  end a block ("}")
   */
  std::string  end_symb ;
 

 
  /*!
   *  Name of the  input  file.
   */
  std::string  filename;

 

  /*!
   *  Vector of  vectorial-properties labels (presently not used).
   */
  std::vector<std::string> vector_prop_labels;


  /*!
   *  Vector of  properties labels.
   */
  std::vector<std::string> string_prop_labels;





  /*!
   *   Map between label and property in the  parameter section; 
   * property is  read as a string and its correct type is checked elsewhere.
   */
  std::map <const std::string,std::string>  string_prop_labels_map;

  /*!
   *   Map between label and a vector of properties in the  parameter section; 
   * properties  are  read as a string and their  correct type is checked elsewhere.
   */
  std::map <const std::string, std::vector<std::string> >  vector_string_prop_labels_map;



  /*!
   *   Map between physical region number  and  the Region Structure associated to it.
   */
  std::map <ID, RegionStructure> device_map;
  
  /*!
   *   Map between  model name and  \c ModelStructure object associated.
   */
  std::map <const std::string, ModelStructure*>  model_structure_map;

  /*!
   *   Map between  BC region number and the Region Structure  for the associated BC region.
   */

  std::map <ID, RegionStructure> model_BC_map;
  

  /*!
   *   Map between   physical model name and its options (stored in ModelOptions).
   */
  std::multimap <const std::string,ModelOptions> physical_model_map;



  /*!
   *   Pointer to \c ModelStructure  object.
   */
  ModelStructure* current_model_point;

  /*!
   *   ModelOptions  object for  the  options read in  each  section.
   */
  ModelOptions temp_options;



  //  private  methods

 
  /*!
   *  Utility  to  find a  keyword in the input file.
   *    
   */
  void find_keyword(std::ifstream& in_stream, const std::string& keyword);

  /*!
   *  Utility  to  find a  keyword in a section
   *    
   */
//  void find_keyword_in_section(std::ifstream& in_stream, const std::string& keyword);
  bool find_keyword_in_section(std::ifstream& in_stream, const std::string& keyword);

  /*!
   *   Method   to  clear  all  maps 
   *    
   */
  void reset_all_maps(void);

  /*!
   *   Method   to  read free-format parameter  section. 
   *    
   */
  //  void parse_options(ifstream& in_stream );
  void parse_options(std::ifstream& in_stream, ModelOptions& options );



  /*!
   *   Method   to  parse "Models" section. 
   *    
   */
  void parse_model(std::ifstream& in_stream);




  /*!
   *   Utility  to  skip  comments (everything on a line, after "#" ). 
   *   Returns true if  the following part of the line is  a  comment.  
   */
  bool skip_comments(std::ifstream& in_stream, const std::string& item);

 

  //!  Method   to  read a list of strings.
  /*!  Parses an assignement "phys_regions = ( ....)  
   */
  void parse_list_phys_ID(std::ifstream& in_stream,std::vector<std::string>& list_regions   );


  //!   Returns options associated to a  section specified in read_parameters. 
  /*!
   * Returns  a reference to an object ModelOptions, containing information 
   * read from a section of input file, and possibly for a given model, as specified previously by 
   * \c read_parameters.
   */
  const  ModelOptions& get_options(void);


  //!   Gets the map of   device description. 
  /*!
   * Returns   a  map which associates a physical region number 
   * with the associated RegionStructure object.
   */
  std::map <ID, RegionStructure>& get_device_map(void); 



  //!   Gets the map of   model definition. 
  /*!
   * Returns   a  map which associates the name of  a model to a pointer 
   * to the class  \c  ModelStructure, which contains a full
   * description of model parameters:  physical regions associated to it, 
   * list of boundary conditions defined for the model. 
   */
  std::map <const std::string, ModelStructure*>& get_model_structure_map(void);



};



#endif // endif define   _INPUTPARSER_H_

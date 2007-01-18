/*=============================================================================
  Copyright (c) 2002-2003 Joel de Guzman
  http://spirit.sourceforge.net/

  Use, modification and distribution is subject to the Boost Software
  License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
  http://www.boost.org/LICENSE_1_0.txt)
  =============================================================================*/


#ifndef _INPUTPARSER_H_
#define _INPUTPARSER_H_

#include <boost/spirit/core.hpp>
#include <boost/spirit/actor/push_back_actor.hpp>
#include <boost/spirit/dynamic.hpp>
#include <boost/spirit/utility/confix.hpp>
#include <iostream>
#include <fstream>

#include <map>
#include <vector>
#include <string>

#include "TypeDefs.h"
#include "ModelStructure.h"




///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace boost::spirit;



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
 * 
 *  
 * Except the model-blocks and the BC_regions-blocks, each block and  subblock can 
 *  contain zero or  any number of  parameters in the assignement form "tagname = tagvalue", 
 *  where "tagname" is  a  string and "tagvalue" is  a single 
 * numerical or string item or a list of  items between  "(" and ")" parenthesis.
 * Format is free for these assignements, provided they are separated by spaces. 
 * Everything following a '#' is  a  comment and  is  disregarded.
 *
 * Public methods read_device() and get_device_map() parse and extract information 
 * from device section.
 * read_models() and get_model_structure_map() methods  parse and extract information 
 * from Models section. Models data are stored in the class \c ModelStructure.
 * Public methods read_parameters() and get_parameters_map()  parse and extract 
 * information from the other parameters sections.
 */


class InputParser{


 public:


  //!  Constructor 
  /*!
    Defines  name  of input  file ("input_file_name").
  */
  InputParser(string& input_file_name);

  //!  Destructor 
  /*!
  
  */
  ~InputParser();




  //!   Parses the  device section of  input  file 
  /*!
   * Method to  read   the  device section of  input  file in an internal map, 
   * which is  returned by \c  get_device_map()
   */
  void read_device();

  //!   Gets the map of   device description. 
  /*!
   * Returns   a  map which associates a physical region number 
   * with the map (of the kind "label = property" ) of the appropriate physical region.
   */
  map <ID,  map <string,string> >& get_device_map();


  //!  Parses the "Models"  section  of the  input file. 
  /*!
   * Method to  parse the "Models"  section  of the  input file.
   *  Fields are  read in an internal  map and can be  extracted by means 
   * of \c get_model_structure_map() 
   * and the methods of class \c  ModelStructure.
   */
  void read_models();

  //!   Gets the map of   model definition. 
  /*!
   * Returns   a  map which associates the name of  a model to a pointer 
   * to the class  \c  ModelStructure, which contains a full
   * description of model parameters:  physical regions associated to it, 
   * list of boundary conditions defined for the model. 
   */
  map <string, ModelStructure*>& get_model_structure_map();


  //!   Reads parameters of  a  given  section   for a given  model. 
  /*!
   * Method to   read  parameters of the  section "section_name"  for a given model " model_name".
   * Parameters are read in an internal  map and can be  extracted by means of 
   *  \c get_parameters_map().
   */
  void read_parameters( string& section_name, string& model_name);
 
  //!   Reads parameters of  a  given  section  
  /*!
   * Overloaded method to   read  model-independent parameters of the  section "section_name".
   * Parameters are read in an internal  map and can be  extracted by means of 
   *  \c get_parameters_map(). 
   */
  void read_parameters( string& section_name);

  //!   Gets a map of  parameters. 
  /*!
   * Returns   a  map of the kind "label = property" containing information 
   * read in a section of input file, and possibly for a given model, as specified previously by 
   * \c read_parameters.
   */
  map <string,string>& get_parameters_map();

 

  //----------------------------------------------------------------------------------

 private:

  //  members data:

  //   control  chars:

  /*!
   *  Char to  start a block ("{").
   */
  string  start_symb ;


  /*!
   *  Char to  end a block ("}")
   */
  string  end_symb ;
 

 
  /*!
   *  Name of the  input  file.
   */
  string  filename;

 

  /*!
   *  Vector of  vectorial-properties labels (presently not used).
   */
  vector<string> vector_prop_labels;


  /*!
   *  Vector of  properties labels.
   */
  vector<string> string_prop_labels;





  /*!
   *   Map between label and property in the  parameter section; 
   * property is  read as a string and its correct type is checked elsewhere.
   */
  map <string,string>  string_prop_labels_map;

  /*!
   *   Map between label and a vector of properties in the  parameter section; 
   * properties  are  read as a string and their  correct type is checked elsewhere.
   */
  map <string, vector<string> >  vector_string_prop_labels_map;



  /*!
   *   Map between physical region number  and  the map associated to it.
   */
  map <ID,  map <string,string> > device_map;
  
  /*!
   *   Map between  model name and  \c ModelStructure object associated.
   */
  map <string, ModelStructure*>  model_structure_map;

  /*!
   *   Map between  BC region number and the map for the BC region associated.
   */
  map <ID,  map <string,string> > model_BC_map;

  /*!
   *   Pointer to \c ModelStructure  object.
   */
  ModelStructure* current_model_point;


  //  private  methods

  //  void initialize_vectors();


 
 
  /*!
   *  Utility  to  find a  keyword in the input file.
   *    
   */
  void find_keyword(ifstream& in_stream, string& keyword);

  /*!
   *  Utility  to  find a  keyword in a section
   *    
   */
  void find_keyword_in_section(ifstream& in_stream, string& keyword);


  /*!
   *   Method   to  clear  all  maps 
   *    
   */
  void reset_all_maps();

  /*!
   *   Method   to  read free-format parameter  section. 
   *    
   */
  void parse_options(ifstream& in_stream );

  /*!
   *   Method   to  parse "Models" section. 
   *    
   */
  void parse_model(ifstream& in_stream);




  /*!
   *   Utility  to  skip  comments (everything on a line, after "#" ). 
   *   Returns true if  the following part of the line is  a  comment.  
   */
  bool skip_comments(ifstream& in_stream, string& item);

 

  //!  Method   to  read a list of strings.
  /*!  Parses an assignement "phys_regions = ( ....)  
   */
  void parse_list_phys_ID(ifstream& in_stream,vector<string>& list_regions   );




};



#endif // endif define   _INPUTPARSER_H_

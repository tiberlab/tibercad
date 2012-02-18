// $Id$


#ifndef _INPUTPARSER_H_
#define _INPUTPARSER_H_

#include <iostream>
#include <fstream>

#include <string>

#include "TypeDefs.h"
#include "ModelOptions.h"



//!  A parser  for  TIBERCAD input  text  file.
/*!
 * Parses an input text file composed by sections:
 *
 * key [name]
 * {
 *   key = value
 *   ...
 *
 *   key [name]
 *   {
 *     ...
 *   }
 * }
 *
 * special keywords: \c @include
 */
class InputParser
{

 public:

  //! Constructor
  InputParser(void);

  //! Destructor
  ~InputParser(void);


  //! Parse the file
  void parse_file(const std::string file, ModelOptions& data);


  //! Define a macro
  static void add_defined(const std::string& name, const std::string& value = "");

  //! Check if a macro is defined
  static bool defined(const std::string& name);

  //! Get the value of a defined macro
  std::string get_defined(const std::string& name);

  //! In-place expansion of macro
  void expand_macro(std::string& in);


 private:


  //! A line counter
  unsigned int line_counter;

  //! A block_counter
  int block_counter;

  //! The current input file
  std::string _file;

  //! Map of all defined macros
  static std::map<std::string, std::string> _defined;


  void  read_block_no_boost(std::istream& in_stream, ModelOptions& options);

  std::string get_token(std::istream& in_stream, bool expand = true);

  void skip_whitespaces(std::istream& in_stream);

  std::string get_until_closing_brace(std::istream& in_stream);

  std::string get_until_closing_quotes(std::istream& in_stream);

  std::string get_until_eol(std::istream& in_stream);

  std::string skip_until_else_or_endif(std::istream& in_stream);

  bool check_validity(std::string& token);

};



#endif // endif define   _INPUTPARSER_H_

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



 private:


  //! A line counter
  unsigned int line_counter;

  //! A block_counter
  int block_counter;

  //! The current input file
  std::string _file;


  void  read_block_no_boost(std::ifstream& in_stream, ModelOptions& options);

  const std::string get_token(std::ifstream& in_stream);

  void skip_whitespaces(std::ifstream& in_stream);

  const std::string get_until_closing_brace(std::ifstream& in_stream);

  const std::string get_until_closing_quotes(std::ifstream& in_stream);

  bool check_validity(std::string& token);

};



#endif // endif define   _INPUTPARSER_H_

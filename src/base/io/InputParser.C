// $Id$

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <string>

#include "InputParser.h"
#include "Messages.h"
#include "Utils.h"
#include "InitFailedException.h"




using namespace std;




InputParser::InputParser(void)
{
}



void
InputParser::parse_file(const string file, ModelOptions& data)
{
  _file = file;
  ifstream in_stream(file.c_str());
  if (in_stream.fail() || !in_stream.good())
    throw InitFailedException("Cannot open input file for reading: " + file);

  // reset block and line counters
  line_counter = 1;
  block_counter = 0;

  read_block_no_boost(in_stream, data);

}



InputParser::~InputParser(void)
{
}



void InputParser::read_block_no_boost(ifstream& in_stream, ModelOptions& options)
{

  string token_1 =  "";
  string token_2 =  "";
  string token_3 =  "";
  string model_name;

  skip_whitespaces(in_stream);

  token_1 = get_token(in_stream);

  if ( token_1 == "}")
  {
    block_counter--;
  }


  //    (check_validity(token_1)  ==
//   if  (check_validity(token_1)  ==   false)
//   {

//     std::ostringstream stm;

//     stm << " ERROR in  line " <<  line_counter +1 <<   endl;

//     throw InitFailedException(stm.str());
//   }


  while ((token_1 != "}") && (!in_stream.eof()))
  {

    if (token_1 == "{")
    {
      ostringstream stm;
      stm << "in " << _file << " on line " <<  line_counter
          << " : missing keyword before \'{\'" <<  endl;
      throw InitFailedException(stm.str());
    }

    skip_whitespaces(in_stream);

    token_2 = get_token(in_stream);

    if (token_1 == "@include")
    {
      ModelOptions tmp;
      Messages::info("Including file " + token_2);
      InputParser().parse_file(token_2, tmp);
      options += tmp;
    }
    else if (token_2 == "=" )
    {

      skip_whitespaces(in_stream);

      token_3 =  get_token(in_stream);

      if  (( token_3 == "{" ) || ( token_3 == "=" ) || ( token_3 == "}") )
      {
        ostringstream stm;
        stm << "in " << _file << " on line " <<  line_counter <<   endl;
        throw InitFailedException(stm.str());
      }

      options.set_option(token_1,token_3);
    }
    else if  (token_2 == "{" )
    {
      //  1 keyword  header -  block

      model_name = token_1;

      //  create a  ModelOptions object  to  contain the  next subblock
      ModelOptions submodel;
      submodel.set_key(model_name);

      block_counter++;

      // recursively  call  to  read  the  subblock and  put the  data in  'submodel' ModelOptions
      read_block_no_boost(in_stream, submodel);

      // add  the  subblock ModelOptions to  the current-level  options, as  a  submodel
      options.add_submodel(model_name,submodel);

    }
    else
    {
      // 2  keyword-block

      skip_whitespaces(in_stream);

      token_3 =  get_token(in_stream);

      if  (token_3 != "{")
      {
        ostringstream stm;
        stm << "in " << _file << " on line " <<  line_counter <<   endl;
        throw InitFailedException(stm.str());
      }

      //  else   we  read  2 names header !

      model_name = token_1;

      //  create a  ModelOptions object  to  contain the  next subblock
      ModelOptions submodel;
      submodel.set_key(model_name);

      block_counter++;

      // recursively  call  to  read  the  subblock and  put the  data in  'submodel' ModelOptions
      read_block_no_boost(in_stream, submodel);

      // we allow several modifiers as vector
      vector<string> tok2;
      Utils::extract_vector(token_2, tok2);
      for (size_t i = 0; i < tok2.size(); ++i)
      {
        ModelOptions opts(submodel);
        opts.set_name(tok2[i]);
        // add  the  subblock ModelOptions to  the current-level  options, as  a  submodel
        options.add_submodel(model_name, opts);
      }

    }

    //  else  error  ???

    skip_whitespaces(in_stream);

    token_1 = get_token(in_stream);

    if (token_1 == "}") block_counter--;

  } //  end while

  if (block_counter < 0)
  {
    ostringstream stm;
    stm << "in " << _file << " on line " << line_counter
        << " : unexpected \'}\'" <<   endl;
    throw InitFailedException(stm.str());
  }

  if (in_stream.eof())
  {
    if (block_counter != 0)
    {
      ostringstream stm;
      stm << "in " << _file
          << " : unexpected end of file (hint: check curly braces)" <<   endl;
      throw InitFailedException(stm.str());
    }
  }

} // end  method





const string InputParser::get_token(ifstream& in_stream)
{

  string token;
  token = "";

  int temp = 0;

  temp = in_stream.get();

  if (temp == '{' || temp == '}' ||  temp  == '=')
  {
    token += temp;
    return token;
  }
  else if (temp == '(')
  {
    token += '(' + get_until_closing_brace(in_stream);
    return token;
  }
  else if (temp == '"')
  {
    //token += '"' + get_until_closing_quotes(in_stream);
    token += get_until_closing_quotes(in_stream);
    return token;
  }
  else if (temp == EOF)
  {
    return token;
  }
  else token +=  temp;


  while (true)
  {
    temp = in_stream.get();

    if ((temp == EOF) || (temp == '{' || temp == '}' || temp  == '='
        || temp  == ' ' || temp  == '\n') || (temp == '\r') || (temp == '#' ))
    {
      in_stream.unget();
      return token;
    }

    token +=  temp;
  }
}



// find next non-whitespace and  skip  comments
//  "space"  include  \n, \r, and  tab
void InputParser::skip_whitespaces(ifstream& in_stream)
{

  int tmp = in_stream.get();
  do {

    while (isspace(tmp))
    {
      if (tmp == '\n')
        line_counter++;

      tmp = in_stream.get();
      if (!in_stream) return;
    }

    // found a non whitespace
    if (tmp == '#')
    {
      // comment => skip until end of line

      while ((tmp != '\n') && (tmp != '\r') )
      {
	if (!in_stream)
	{
	  in_stream.unget();
	  return;
	}
	tmp = in_stream.get();
      }

      continue;
    }
    else
    {
      in_stream.unget();
      return;
    }
  } while (in_stream);

}



const string InputParser::get_until_closing_brace(ifstream& in_stream)
{

  string str = "";
  int temp = 0;

  while (true)
  {
    temp = in_stream.get();
    if (temp == EOF)
      return str;

    str += temp;
    if (temp == ')')
      return str;

    if (temp == '\n')
      line_counter++;
  }

}



const string InputParser::get_until_closing_quotes(ifstream& in_stream)
{

  string str = "";
  int    temp = 0;

  while (true)
  {

    temp = in_stream.get();
    if (temp == EOF) return str;
    else if(temp == '"')
    {
      return str ;
      //return (str + '"') ;
    }
    else if (temp == '\n')
      line_counter++;
    else
      str += temp;

  }

}



bool InputParser::check_validity(string& token)
{

  bool  check;
  check =  true;
  int i;
  i=0;
  string::iterator it;

  for ( it=token.begin() ; it < token.end(); it++ )
  {
    if ( (!isalnum(*it)  ) && ( *it != '_') && ( *it != '-')   )
      check =  false;

  }

  return  check;
}

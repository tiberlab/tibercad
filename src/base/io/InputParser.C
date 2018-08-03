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



map<string, string>
InputParser::_defined;



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

  read_block(in_stream, data);

}



InputParser::~InputParser(void)
{
}


void
InputParser::add_defined(const string& name, const string& value,
    bool warn_on_redefine)
{
  if (warn_on_redefine && (_defined.find(name) != _defined.end()))
    Messages::warning("Redefining input file macro '" + name + "'");
  _defined[name] = value;
}

bool
InputParser::defined(const string& name)
{
  return (_defined.find(name) != _defined.end());
}

string
InputParser::get_defined(const string& name)
{
  string val("");
  if (_defined.find(name) != _defined.end())
    val = _defined[name];

  return val;
}


void
InputParser::expand_macro(string& in)
{
  map<string, string>::iterator it(_defined.find(in));
  if (it != _defined.end())
    in = it->second;
}

void InputParser::read_block(istream& in_stream, ModelOptions& options)
{

  string token_1 =  "";
  string token_2 =  "";
  string token_3 =  "";
  string model_name;

  token_1 = get_token(in_stream);

  if ( token_1 == "}")
  {
    block_counter--;
  }


  while ((token_1 != "}") && (!in_stream.eof()))
  {

    if (token_1 == "{")
    {
      ostringstream stm;
      stm << "in " << _file << " on line " <<  line_counter
          << " : missing keyword before \'{\'" <<  endl;
      throw InitFailedException(stm.str());
    }


    if (token_1.at(0) == '@')
    {

      if (token_1 == "@include")
      {
        ModelOptions tmp;
        InputParser().parse_file(get_token(in_stream), tmp);
        options += tmp;
      }
      else if (token_1 == "@define")
      {
        // the next token is a new status variable
        string tok1(get_token(in_stream, false));

        istringstream is(get_until_eol(in_stream));
        string value(get_token(is));
        while (is.good())
          value += " " + get_token(is);
        add_defined(tok1, value);
      }
      else if (token_1 == "@ifdef")
      {
        // if it is defined we go on, otherwise we skip
        if (!defined(get_token(in_stream, false)))
          skip_until_else_or_endif(in_stream);
      }
      else if (token_1 == "@ifndef")
      {
        // if it is not defined we go on, otherwise we skip
        if (defined(get_token(in_stream, false)))
          skip_until_else_or_endif(in_stream);
      }
      else if (token_1 == "@else")
      {
        // we can skip until @endif
        skip_until_else_or_endif(in_stream);
      }
      else if (token_1 == "@endif")
      {
      }

      token_1 = get_token(in_stream);

      if (token_1 == "}") block_counter--;
      continue;
    }

    token_2 = get_token(in_stream);

    if (token_2 == "=" )
    {
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
      read_block(in_stream, submodel);

      // add  the  subblock ModelOptions to  the current-level  options, as  a  submodel
      options.add_submodel(model_name,submodel);

    }
    /*
    else if  (token_2 == "}" )
    {
      options.set_option(token_1, true);
      token_1 = token_2;
      in_stream.putback('}');
      //for (int i = 0; i < token_2.size(); i++)
      //  in_stream.unget();
    }
    */
    else
    {
      // 2  keyword-block

      token_3 =  get_token(in_stream);

      if  (token_3 != "{")
      {
        /*
        options.set_option(token_1, true);
        token_1 = token_2;
        in_stream.putback(' ');
        for (int i = 0; i < token_3.size(); i++)
          in_stream.unget();

        continue;
        */
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
      read_block(in_stream, submodel);

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





string InputParser::get_token(istream& in_stream, bool expand)
{
  skip_whitespaces(in_stream);

  string token;

  int temp = 0;
  switch (temp = in_stream.get())
  {
    case '{':
    case '}':
    case '=':
      token += temp;
      return token;
      break;

    case '(':
      token += '(' + get_until_closing_brace(in_stream);
      return token;
      break;

    case '"':
      token += get_until_closing_quotes(in_stream);
      break;

    case EOF:
      return token;
      break;

    default:
      in_stream.unget();
      //token += temp;
      break;
  }


  string tok;
  while (in_stream.good())
  {
    switch (temp = in_stream.get())
    {
      case '"':
        expand_macro(tok);
        token += tok + get_until_closing_quotes(in_stream);
        tok.clear();
        break;

      case '{':
      case '}':
      case '(':
      case '=':
      case ' ':
      case '\n':
      case '\r':
      case '\t':
      case '#':
      case EOF:
        in_stream.unget();
        if (expand)
          expand_macro(tok);
        token += tok;
        return token;
        break;

      default:
        tok += temp;
        break;
    }
  }
}


string
InputParser::skip_until_else_or_endif(std::istream& in_stream)
{
  string tok;
  do
  {
    tok = get_token(in_stream);
  } while ((tok != "@else") && (tok != "@endif"));

  return tok;
}

// find next non-whitespace and  skip  comments
//  "space"  include  \n, \r, and  tab
void InputParser::skip_whitespaces(istream& in_stream)
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



string InputParser::get_until_closing_brace(istream& in_stream)
{

  string str = "";
  int temp = 0;

  skip_whitespaces(in_stream);

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



string InputParser::get_until_closing_quotes(istream& in_stream)
{

  string str = "";
  int temp = 0;
  int last = 0;

  while (true)
  {
    last = temp;
    temp = in_stream.get();
    if (temp == EOF) return str;
    else if ((temp == '"') && (last != '\\'))
      return str ;
    else if (temp == '\n')
      line_counter++;
    else if ((temp == '\\') && (last != '\\'))
      continue;
    else
      str += temp;

  }

}



string InputParser::get_until_eol(istream& in_stream)
{
  string str;
  int temp = 0;

  while (in_stream)
  {
    temp = in_stream.get();

    if ((temp == '\n') || (temp == '#'))
    {
      in_stream.unget();
      break;
    }
    else if (temp == '\\')
      skip_whitespaces(in_stream);
    else
      str += temp;
  }

  return str;
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

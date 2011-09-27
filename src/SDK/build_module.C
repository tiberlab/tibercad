// $Id: main.C 2874 2011-07-15 21:12:36Z maufder $

#include "boost/algorithm/string/trim.hpp"
#include "boost/regex.hpp"

#include "TiberCad.h"
#include "Utils.h"
#include "InputParser.h"
#include "ModelOptions.h"

#include "tiber_config.h"


#include <iostream>
#include <sstream>
#include <fstream>

#ifdef HAVE_LIBREADLINE
# ifdef HAVE_READLINE_READLINE_H
#  include <readline/readline.h>
# else
#  include <readline.h>
# endif
#endif

#if defined(_WIN32)
# include <windows.h>
#endif

#include <cstdio>
#include <getopt.h>



using namespace std;

#if defined(_WIN32)
// Return filename from file open dialog
static string open_file(const char *filter = "All Files (*.*)\0*.*\0", HWND owner = NULL)
{
  char fileName[MAX_PATH] = "";
  OPENFILENAME ofn;
  ZeroMemory(&ofn, sizeof(ofn));

  ofn.lStructSize = sizeof(OPENFILENAME);
  ofn.hwndOwner = owner;
  ofn.lpstrFilter = filter;
  ofn.lpstrFile = fileName;
  ofn.nMaxFile = MAX_PATH;
  ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
  ofn.lpstrDefExt = "";

  string fileNameStr;
  if (GetOpenFileName(&ofn))
    fileNameStr = fileName;

  return fileNameStr;
}
#endif




namespace
{
  //! Interactive or batch mode (only useful for Win)
  bool interactive;

  //! Content of TIBERCADROOT
  string tc_root;

  void usage(void)
  {
#if defined(_WIN32)
    cout << endl << "Usage:" << endl
      << "  from command line: build_module configfile" << endl
      << "  or double click on configfile" << endl << endl;
    cout << "press Enter ...";
    if (interactive) cin.get();
# else
    cout << endl << "Usage: tibercad configfile" << endl << endl;
# endif
  }

  void replace(std::string& s, const std::string& from, const std::string& to)
  {
    const boost::regex e(from);
    string news = regex_replace(s, e, to, boost::match_default | boost::format_sed);
    s = news;
  }
}


class Compiler
{
  public:

    static void setup(const ModelOptions& options);

    static int compile(const std::string& source);

    static int compile(const std::string& compiler, const std::string& source);

    static int link(const std::string& target, const std::string& sources);

    //bool static_link(const std::string& target, const std::string& sources);

  private:

    Compiler(const ModelOptions& options);

    ModelOptions _options;

    std::string _preprocessor_flags;

    std::string _compiler_flags;

    std::string _linker_flags;

    std::string _executable(void);

    int _compile(const std::string& source);

    //! Compiler name to compiler object
    static std::map<std::string, Compiler*> _compilers;

    //! File extension to compiler object
    static std::map<std::string, string> _file_association;

};


void process_module(const string& name, const ModelOptions& options);



// Will be extended with tools for command line argument parsing
// and so on
int main (int argc, char** argv)
{

  interactive = true;

  opterr = 0;
  int c;
  while ((c = getopt(argc, argv, "bv")) != -1)
    switch (c)
    {
      case 'b':
        interactive = false;
        break;

      case '?':
        cout << "Unknown option: -" << (char) optopt << endl;
      default:
        usage();
        return 1;
    }

#ifdef _WIN32
  if (!interactive)
#endif
    if (optind >= argc)
    {
      usage();
      return 1;
    }



  // take input file from command line or ask for it
  string inputfile;
  if (optind < argc)
    inputfile = string(argv[optind]);

//#ifdef HAVE_LIBREADLINE
//    char *line = readline ("Enter input file: ");
//    inputfile = string(line);
//    free(line);
//    boost::algorithm::trim(inputfile);
//    cout << endl;
//#else
//#endif


  // do some preparation
  {
#if defined(_WIN32)
#if defined(__CYGWIN__)
    // we first convert the filename to something more UNIX like
    Utils::convert_win32_path_to_posix(inputfile);
#endif

    if (inputfile.empty())
      inputfile = open_file();


    char* root = getenv("TIBERCADROOT");
    if (root == NULL)
    {
      const size_t bufsize = 1024;
      char buffer[bufsize];
      if (!GetModuleFileName(NULL, buffer, bufsize))
        cerr << "Problems detecting installation path." << endl;
      string program(buffer);
#if defined(__CYGWIN__)
      Utils::convert_win32_path_to_posix(program);
#endif
      string exepath(Utils::dirname(program));
# ifdef HAVE_SETENV
      setenv("TIBERCADROOT", exepath.c_str(), 1);
# else
#  ifdef HAVE_PUTENV
      string newroot("TIBERCADROOT=" + exepath);
      putenv(newroot.c_str());
#  else
#   error "Neither setenv nor putenv available"
#  endif
# endif
    }
#endif

    {
      // we check here if the input file exists
      ifstream infile;
      infile.open(inputfile.c_str());
      if (infile.fail() || !infile.good())
      {
        infile.close();
        cerr << "TiberCAD: Cannot open file " << inputfile <<  " for reading." << endl;
        if (interactive) cin.get();
        return 1;
      }
      infile.close();
    }

#if defined(_WIN32)
    {
      //TCHAR cwd[MAX_PATH] = "";

      //if(!GetCurrentDirectory(sizeof(cwd) - 1, cwd))
      //{
      //  cerr << "TiberCAD: cannot get working directory." << endl;
      //  if (interactive) cin.get();
      //  return 1;
      //}

      //string dirname(cwd);

      string dirname = Utils::dirname(inputfile);
      if (!dirname.empty())
      {
        // try to set the working directory
        if (SetCurrentDirectory(dirname.c_str()) == 0)
        {
          cerr << "TiberCAD: cannot set " << dirname << " as working directory." << endl;
          if (interactive) cin.get();
          return 1;
        }
      }
    }
#endif


  }

  //
  // here begins real task
  //
  int error = 1;

  {
    char* root = getenv("TIBERCADROOT");
    if (root == NULL)
    {
      cerr << "TIBERCADROOT environment variable is not set.\n";
      return 1;
    }
    tc_root = string(root);
  }

  InputParser parser;
  ModelOptions global_config;
  parser.parse_file(tc_root + "/share/modules.conf", global_config);

  ModelOptions::submodel_iterator it(global_config.submodels_begin("Compiler"));
  if (it == global_config.submodels_end("Compiler"))
  {
    cerr << "No compiler defined in configuration file.\n";
    return 1;
  }

  Compiler::setup(it->second);


  ModelOptions config;
  parser.parse_file(inputfile, config);

  it = config.submodels_begin("Module");
  ModelOptions::submodel_iterator end(config.submodels_end("Module"));

  for ( ; it != end; ++it)
  {
    process_module(it->second.get_name(), it->second);
  }
  
  return error;
}






map<string, Compiler*>
Compiler::_compilers;

map<string, string>
Compiler::_file_association;

void
Compiler::setup(const ModelOptions& options)
{

  _file_association["C"] = "C++";
  _file_association["cpp"] = "C++";
  _file_association["cxx"] = "C++";
  _file_association["c"] = "C";
  _file_association["f"] = "Fortran";
  _file_association["F"] = "Fortran";
  _file_association["f90"] = "Fortran";


  ModelOptions::const_submodel_iterator it(options.submodels_begin());
  ModelOptions::const_submodel_iterator end(options.submodels_end());
  for ( ; it != end; ++it)
  {
    //const string& name = it->first;
    //const ModelOptions& opts = it->second;

    ModelOptions opts;
    if (options.find_option("path"))
      opts.set_option("path", options.get_option("path", ""));

    if (options.find_option("includes"))
      opts.set_option("includes", options.get_option("includes", ""));

    opts += it->second;
    opts.set_key(it->second.get_key());

    Compiler* comp = new Compiler(opts);

    _compilers[opts.get_key()] = comp;

  }
}


int
Compiler::compile(const std::string& compiler, const std::string& source)
{
  map<string, Compiler*>::iterator it(_compilers.find(compiler));
  if (it == _compilers.end())
  {
    cerr << "Compiler " << compiler << " not configured." << endl;
    exit(1);
  }

  return it->second->_compile(source);
}


int
Compiler::compile(const std::string& source)
{
  string suffix = Utils::file_extension(source);
  map<string, string>::iterator it(_file_association.find(suffix));
  if (it == _file_association.end())
  {
    cerr << "Unknwon file type: " << suffix << endl;
    exit(1);
  }

  return compile(it->second, source);
}


Compiler::Compiler(const ModelOptions& options) :
  _options(options)
{

  ostringstream pre;
  pre << "-DARCH=" + string(ARCH);

  string path = options.get_option("sdk_path", tc_root);
  replace(path, "@ARCH", ARCH);

  vector<string> includes;
  options.get_option("includes", includes);
  for (size_t i = 0; i < includes.size(); ++i)
  {
    string inc = includes[i];
    replace(inc, "@ARCH", ARCH);
    pre << " -I" << path << "/" << inc;
  }

  _preprocessor_flags = pre.str();

  _compiler_flags = options.get_option("compiler_flags", "");
#if defined(_WIN32)
  replace(_compiler_flags, "-fPIC", "");
#endif

  _linker_flags = options.get_option("linker_flags", "");
}


inline
string
Compiler::_executable(void)
{
  string exe(_options["executable"]);
  replace(exe, "@ARCH", ARCH);
  return tc_root + "/" + exe;
}


int
Compiler::_compile(const std::string& source)
{
  cout << "Compiling " << source << " (" << _options.get_key() << ") ...\n";

  string basename = Utils::basename(source);

  ostringstream cmdline;
  cmdline << _executable() << " " << _preprocessor_flags << " " <<
    _compiler_flags << " " << source;

  return system(cmdline.str().c_str());
}


void process_module(const string& name, const ModelOptions& options)
{
  vector<string> sources;
  options.get_option("sources", sources);

  string creatable = options.get_option("creatable", "");

  char* root = getenv("TIBERCADROOT");

  string binsuffix;
#if defined(_WIN32)
  string libsuffix(".dll");
  binsuffix = ".exe";
#elif defined(_APPLE_)
  string libsuffix(".dylib");
#else
  string libsuffix(".so");
#endif

  string tcroot(root);
  string sdkdir(tcroot + "/SDK");

  string cpp(sdkdir);
  cpp += "/compiler/" + string(ARCH) + "/bin/g++" + binsuffix;

  string cppflags;
  cppflags += "-DARCH=" + string(ARCH) + " ";
  cppflags += "-I" + sdkdir + "/petsc-3.0.0-p12/include ";
  cppflags += "-I" + sdkdir + "/petsc-3.0.0-p12/" + ARCH + "/include ";
  cppflags += "-I" + sdkdir + "/slepc-3.0.0-p7/include ";
  cppflags += "-I" + sdkdir + "/slepc-3.0.0-p7/" + ARCH + "/include ";
  cppflags += "-I" + string(root) + "/include/base/common";

  string cxxflags;
  cxxflags += "-std=gnu++0x -pthread ";
#if !defined(_WIN32)
  cxxflags += "-fPIC ";
#endif

  string ldflags;
  ldflags += "-shared ";
  ldflags += "-L" + tcroot + "/" + ARCH + "/lib ";
  ldflags += "-Wl,-rpath,\'$$ORIGIN\' -Wl,-rpath,\'$$ORIGIN/../lib\' ";

  string libfile = name + libsuffix;

  ostringstream cmdline;
  cmdline << cpp;
  cmdline << " " << cppflags << " " << cxxflags << " " << ldflags
      << " -o " << libfile;
  for (size_t i = 0; i < sources.size(); ++i)
  {
    cmdline << " " << sources[i];
    Compiler::compile("C++", sources[i]);
  }

  //cout << "Compiling " << libfile << " ...";
  //system(cmdline.str().c_str());
 
}

// $Id: main.C 2874 2011-07-15 21:12:36Z maufder $

#include "boost/algorithm/string/trim.hpp"
#include "boost/regex.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/operations.hpp"

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




namespace BuildModule
{
  //! Interactive or batch mode (only useful for Win)
  bool interactive;

  //! Content of TIBERCADROOT
  std::string tc_root;

  //! Verbose mode
  bool verbose = false;

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

    static std::string compile(const std::string& source, const std::string& flags = "",
        const std::string& compiler = "");

    static int link(const std::string& target, const std::vector<std::string>& objects,
        const std::string& flags = "", const std::string& linker = "C++");

    static void add_library(const ModelOptions& options, const std::string& compiler = "");

    //bool static_link(const std::string& target, const std::string& sources);

  private:

    Compiler(const ModelOptions& options);

    ModelOptions _options;

    std::string _outdir;

    std::string _preprocessor_flags;

    std::string _compiler_flags;

    std::string _linker_flags;

    std::string _executable(void);

    std::string _compile(const std::string& source, const std::string& flags);

    int _link(const std::string& target, const std::vector<std::string>& objects,
        const std::string& flags);

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

  BuildModule::interactive = true;

  opterr = 0;
  int c;
  while ((c = getopt(argc, argv, "bv")) != -1)
    switch (c)
    {
      case 'b':
        BuildModule::interactive = false;
        break;

      case '?':
        cout << "Unknown option: -" << (char) optopt << endl;
      default:
        BuildModule::usage();
        return 1;
    }

#ifdef _WIN32
  if (!BuildModule::interactive)
#endif
    if (optind >= argc)
    {
      BuildModule::usage();
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
        if (BuildModule::interactive) cin.get();
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
      //  if (BuildModule::interactive) cin.get();
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
          if (BuildModule::interactive) cin.get();
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
    string rs(root);
    boost::filesystem::path tcroot(rs);
    tcroot = boost::filesystem::system_complete(tcroot);

    BuildModule::tc_root = tcroot.string();
  }

  InputParser parser;
  ModelOptions global_config;
  parser.parse_file(BuildModule::tc_root + "/share/modules.conf", global_config);

  ModelOptions::submodel_iterator it(global_config.submodels_begin("Compiler"));
  if (it == global_config.submodels_end("Compiler"))
  {
    cerr << "No compiler defined in configuration file.\n";
    return 1;
  }

  Compiler::setup(it->second);

  it = global_config.submodels_begin("Library");
  for ( ; it != global_config.submodels_end("Library"); ++it)
  {
    Compiler::add_library(it->second);
  }


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

  _file_association[".C"] = "C++";
  _file_association[".cpp"] = "C++";
  _file_association[".cxx"] = "C++";
  _file_association[".c"] = "C";
  _file_association[".f"] = "Fortran";
  _file_association[".F"] = "Fortran";
  _file_association[".f90"] = "Fortran";


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

    if (options.find_option("linker_flags"))
      opts.set_option("linker_flags", options.get_option("linker_flags", ""));

    opts += it->second;
    opts.set_key(it->second.get_key());

    Compiler* comp = new Compiler(opts);

    _compilers[opts.get_key()] = comp;

  }
}



std::string
Compiler::compile(const std::string& source, const std::string& flags,
    const std::string& compiler)
{
  map<string, Compiler*>::iterator comp(_compilers.end());

  if (compiler.empty())
  {
    string suffix = Utils::file_extension(source);
    map<string, string>::iterator it(_file_association.find(suffix));
    if (it == _file_association.end())
    {
      cerr << "Unknwon file type: " << suffix << endl;
      exit(1);
    }

    comp = _compilers.find(it->second);
  }
  else
    comp = _compilers.find(compiler);

  if (comp == _compilers.end())
  {
    cerr << "Compiler " << compiler << " not configured." << endl;
    exit(1);
  }

  return comp->second->_compile(source, flags);
}



int
Compiler::link(const std::string& target, const std::vector<std::string>& objects,
    const std::string& flags, const std::string& linker)
{
  map<string, Compiler*>::iterator comp(_compilers.find(linker));

  if (comp == _compilers.end())
  {
    cerr << "Compiler " << linker << " not configured." << endl;
    exit(1);
  }

  return comp->second->_link(target, objects, flags);
}




void
Compiler::add_library(const ModelOptions& options, const std::string& compiler)
{
  ostringstream pre;

  string path = options.get_option("path", BuildModule::tc_root);
  BuildModule::replace(path, "@ROOT", BuildModule::tc_root);
  BuildModule::replace(path, "@ARCH", ARCH);
  if (path[0] != '/')
    path = BuildModule::tc_root + "/" + path;

  vector<string> includes;
  options.get_option("includes", includes);
  for (size_t i = 0; i < includes.size(); ++i)
  {
    string inc = includes[i];
    BuildModule::replace(inc, "@ROOT", BuildModule::tc_root);
    BuildModule::replace(inc, "@ARCH", ARCH);
    if (inc[0] != '/')
      pre << " -I" << path << "/" << inc;
    else
      pre << " -I" << inc;
  }

  ostringstream linkerflags;

  string libpath;
  libpath = options.get_option("libpath", libpath);
  if (!libpath.empty())
  {
    BuildModule::replace(libpath, "@ROOT", BuildModule::tc_root);
    BuildModule::replace(libpath, "@ARCH", ARCH);
    if (libpath[0] != '/')
      linkerflags << "-L " << BuildModule::tc_root << "/" << libpath << " ";
    else
      linkerflags << "-L " << libpath << " ";
  }

  vector<string> ldflags;
  options.get_option("linker_flags", ldflags);
  for (size_t i = 0; i < ldflags.size(); ++i)
  {
    string flags = ldflags[i];
    BuildModule::replace(flags, "@ROOT", BuildModule::tc_root);
    BuildModule::replace(flags, "@ARCH", ARCH);
    linkerflags << flags;
  }

  if (compiler.empty())
  {
    map<string, Compiler*>::iterator it(_compilers.begin());
    for ( ; it != _compilers.end(); ++it)
    {
      it->second->_preprocessor_flags += " " + pre.str();
      it->second->_linker_flags += " " + linkerflags.str();
    }
  }

}



Compiler::Compiler(const ModelOptions& options) :
  _options(options)
{

  ostringstream pre;
  pre << "-DARCH=" + string(ARCH);

  string path = options.get_option("sdk_path", BuildModule::tc_root);
  BuildModule::replace(path, "@ROOT", BuildModule::tc_root);
  BuildModule::replace(path, "@ARCH", ARCH);

  vector<string> includes;
  options.get_option("includes", includes);
  for (size_t i = 0; i < includes.size(); ++i)
  {
    string inc = includes[i];
    BuildModule::replace(inc, "@ARCH", ARCH);
    pre << " -I" << path << "/" << inc;
  }

  pre << " -I " << boost::filesystem::current_path();

  _preprocessor_flags = pre.str();

  _compiler_flags = options.get_option("compiler_flags", "");

  _outdir = "." + string(ARCH);

#if defined(_WIN32)
  BuildModule::replace(_compiler_flags, "-fPIC", "");
#endif

#if defined(_LINUX)
  _linker_flags = "-Wl,--as-needed ";
#endif
  _linker_flags += options.get_option("linker_flags", "");
}


inline
string
Compiler::_executable(void)
{
  string exe(_options["executable"]);
  BuildModule::replace(exe, "@ARCH", ARCH);
  if (exe[0] != '/')
    exe = BuildModule::tc_root + "/" + exe;

  return exe;
}


std::string
Compiler::_compile(const std::string& source, const std::string& flags)
{
  cout << "Compiling " << source << " (" << _options.get_key() << ") ...\n";

  boost::filesystem::path outdir(_outdir);
  if (!boost::filesystem::exists(outdir))
    boost::filesystem::create_directory(outdir);

  string basename = Utils::basename(source);
  string target = _outdir + "/" + basename + ".o";

  ostringstream cmdline;
  cmdline << _executable() << " " << _preprocessor_flags << " " <<
    _compiler_flags << " " << flags << " -o " <<
    target << " " << source;


  if (BuildModule::verbose) cout << cmdline.str() << endl;
  if (system(cmdline.str().c_str()) != 0)
    target = "";

  return target;
}



int
Compiler::_link(const std::string& target, const std::vector<std::string>& objects,
    const std::string& flags)
{
  cout << "Linking " << target << " (" << _options.get_key() << ") ...\n";

  ostringstream cmdline;
  cmdline << _executable();
  for (size_t i = 0; i < objects.size(); ++i)
    cmdline << " " << objects[i];

  cmdline << " " << flags << " " <<_linker_flags << " " << " -o " << target;

  if (BuildModule::verbose) cout << cmdline.str() << endl;

  return system(cmdline.str().c_str());
}



void process_module(const string& name, const ModelOptions& options)
{
  vector<string> sources;
  options.get_option("sources", sources);

  string creatable = options.get_option("createable", "");

  if (creatable.empty())
  {
    // guess it from the first source file
    creatable = Utils::basename(sources[0]);
  }

  string module = options.get_option("module", "");

  string modulename = options.get_key();
  string type = options.get_name();
  if (modulename == "Module")
  {
    modulename = options.get_name();
    module = modulename;
    type = options.get_option("type", "");
  }

  if (!type.empty())
    modulename += "_" + type;

  cout << "Processing module " << modulename << " ...\n";

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


  string compileflags = options.get_option("compiler_flags", "");
  compileflags += " -DMODULE_NAME=" + module;
  compileflags += " -DCREATABLE=" + creatable;

  string module_code = options.get_option("creator_code", "");
  if (!module_code.empty())
  {
    ofstream of(module + "_creator.h");
    of << module_code;
    compileflags += " -DCREATORCODE=" + module + "_creator.h";
  }

  vector<string> objects;
  for (size_t i = 0; i < sources.size(); ++i)
  {
    string obj = Compiler::compile(sources[i], compileflags);
    if (obj.empty())
    {
      cerr << "Error in compilation: Could not compile "
          << sources[i] << endl;
      exit(1);
    }
    objects.push_back(obj);
  }

  string linkflags = options.get_option("linker_flags", "");

  string target = modulename + libsuffix;
  Compiler::link(target, objects, linkflags);

  string instpath = options.get_option("installpath", "");
  BuildModule::replace(instpath, "@ARCH", ARCH);
  BuildModule::replace(instpath, "@ROOT", BuildModule::tc_root);
  BuildModule::replace(instpath, "@MODULE", modulename);
  if (instpath[0] != '/')
    instpath = BuildModule::tc_root + "/" + instpath;

  using namespace boost::filesystem;
  path from("./" + target);
  path to(instpath + "/" + target);
  copy_file(from, to, copy_option::overwrite_if_exists);
  remove(from);


  // recursively process submodules
  ModelOptions::const_submodel_iterator it(options.submodels_begin());
  ModelOptions::const_submodel_iterator end(options.submodels_end());
  for ( ; it != end; ++it)
  {
    if (it->first == "Dependency") continue;

    ModelOptions opts(it->second);
    if (!opts.find_option("installpath"))
      opts["installpath"] = instpath;

    //string cpflags = opts.get_option("compiler_flags", "");
    string ldflags = opts.get_option("linker_flags", "");
#if defined(_LINUX)
    ldflags += "-Wl,--as-needed ";
#endif
    ldflags += instpath + "/" + target;
    opts["linker_flags"] = ldflags;

    opts["module"] = module;

    process_module(it->first, opts);
  }


}

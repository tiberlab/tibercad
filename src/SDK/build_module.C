// $Id: main.C 2874 2011-07-15 21:12:36Z maufder $

#include "boost/algorithm/string/trim.hpp"
#include "boost/regex.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/operations.hpp"

#include "tibercad/base/TiberCad.h"
#include "tibercad/utils/Utils.h"
#include "base/io/InputParser.h"
#include "tibercad/base/ModelOptions.h"

#include "tibercad/base/tiber_config.h"


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
  bool interactive = true;

  //! Content of TIBERCADROOT
  string tc_root;

  //! Verbose mode
  bool verbose = false;

  //! Compile method
  string compile_method = "devel";

  void usage(void)
  {
#if defined(_WIN32)
    cout << endl << "Usage:" << endl
      << "  from command line: build_module configfile" << endl
      << "  or double click on configfile" << endl << endl;
    cout << "press Enter ...";
    if (interactive) cin.get();
# else
    cout << endl << "Usage: build_module [-b] [-v] [-m method] [-c globalconfig] moduleconfig" << endl << endl;
# endif
  }

  void replace(string& s, const string& from, const string& to)
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

    static string compile(const string& source, const string& flags = "",
        const string& compiler = "");

    static int link(const string& target, const vector<string>& objects,
        const string& flags = "", const string& linker = "C++");

    static void add_library(const ModelOptions& options, const string& compiler = "");

    //bool static_link(const string& target, const string& sources);

  private:

    Compiler(const ModelOptions& options);

    ModelOptions _options;

    string _outdir;

    string _preprocessor_flags;

    string _compiler_flags;

    string _linker_flags;

    string _executable;

    //! Check dependency timestamps
    /*!
     * If any of the dependencies has a newer (write) timestamp
     * the method returns \c true, \c false otherwise.
     */
    static bool _needs_build(const string& target,
        const vector<string>& dependencies);

    string _compile(const string& source, const string& flags);

    int _link(const string& target, const vector<string>& objects,
        const string& flags);

    //! Compiler name to compiler object
    static map<string, Compiler*> _compilers;

    //! File extension to compiler object
    static map<string, string> _file_association;

};


void process_module(const string& name, const ModelOptions& options);

void check_dependency(const ModelOptions& options, const string& modulename,
    string& compileflags, string& linkflags);


// Will be extended with tools for command line argument parsing
// and so on
int main(int argc, char** argv)
{

  //MPI_Init(&argc, &argv);

  string global_config_file;
  opterr = 0;
  int c;
  while ((c = getopt(argc, argv, "bvc:m:")) != -1)
    switch (c)
    {
      case 'b':
        BuildModule::interactive = false;
        break;

      case 'v':
        BuildModule::verbose = true;
        break;

      case 'c':
        global_config_file = string(optarg);
        break;

      case 'm':
        BuildModule::compile_method = string(optarg);
        break;

      case '?':
        cout << "Unknown option: -" << (char) optopt << endl;
        [[fallthrough]];
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
    if (root == nullptr)
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
  int error = 0;

  {
    char* root = getenv("TIBERCADROOT");
    if (root == nullptr)
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
  if (global_config_file.empty())
    global_config_file = BuildModule::tc_root + "/SDK/etc/modules.conf";
  parser.parse_file(global_config_file, global_config);

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



string
Compiler::compile(const string& source, const string& flags,
    const string& compiler)
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
Compiler::link(const string& target, const vector<string>& objects,
    const string& flags, const string& linker)
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
Compiler::add_library(const ModelOptions& options, const string& compiler)
{
  ostringstream pre;

  string path = options.get_option("path", BuildModule::tc_root + "/SDK");
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
      pre << " -I" << path << "/include/" << inc;
    else
      pre << " -I" << inc;
  }

  ostringstream linkerflags;

  string libpath(path + "@ARCH/lib");
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
    BuildModule::replace(flags, "@METHOD", BuildModule::compile_method);
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

  string path = options.get_option("path", "");
  BuildModule::replace(path, "@ROOT", BuildModule::tc_root);
  BuildModule::replace(path, "@ARCH", ARCH);

  _executable = options.get_option("executable", _executable);
  BuildModule::replace(_executable, "@ARCH", ARCH);
  BuildModule::replace(_executable, "@ROOT", BuildModule::tc_root);
  if (_executable.empty())
  {
    cerr << "You must give executable name for compilers (in "
      << options.get_name() << ")";
    exit(1);
  }

  if ((_executable[0] != '/') && !path.empty())
    _executable = path + "/" + _executable;

  vector<string> includes;
  options.get_option("includes", includes);
  for (size_t i = 0; i < includes.size(); ++i)
  {
    string inc = includes[i];
    BuildModule::replace(inc, "@ROOT", BuildModule::tc_root);
    BuildModule::replace(inc, "@ARCH", ARCH);

    pre << " -I" << inc;
  }

  pre << " -I " << boost::filesystem::current_path();

  _preprocessor_flags = pre.str();

  _compiler_flags = options.get_option("compiler_flags", "");

  string arch_str(ARCH);
  _outdir = "." + arch_str;

  if ((arch_str == "i686-linux") || (arch_str == "x86_64-linux"))
    _linker_flags = "-Wl,--as-needed ";

  _linker_flags += options.get_option("linker_flags", "");

  BuildModule::replace(_linker_flags, "@ROOT", BuildModule::tc_root);
  BuildModule::replace(_linker_flags, "@ARCH", ARCH);

  if ((arch_str == "i686-w64-mingw32") || (arch_str == "x86_64-w64-mingw32"))
  {
    BuildModule::replace(_compiler_flags, "-fPIC", "");
    BuildModule::replace(_compiler_flags, "-pthread", "");
    BuildModule::replace(_linker_flags, "-pthread", "");
  }
}



bool
Compiler::_needs_build(const string& target,
    const vector<string>& dependencies)
{
  bool need_compilation = false;

  boost::filesystem::path target_p(target);
  if (!boost::filesystem::exists(target_p))
  {
    need_compilation = true;
  }
  else
  {
    time_t target_ts = boost::filesystem::last_write_time(boost::filesystem::path(target));
    for (size_t i = 0; i < dependencies.size(); ++i)
    {
      boost::filesystem::path p(dependencies[i]);
      if (boost::filesystem::exists(p))
      {
        time_t dep_ts = boost::filesystem::last_write_time(p);
        if (difftime(target_ts, dep_ts) < 0.0)
        {
          need_compilation = true;
          break;
        }
      }
    }
  }

  return need_compilation;
}




string
Compiler::_compile(const string& source, const string& flags)
{
  cout << "  Compiling " << source << " (" << _options.get_key() << ") ";

  boost::filesystem::path outdir(_outdir);
  if (!boost::filesystem::exists(outdir))
    boost::filesystem::create_directories(outdir);

  string basename = Utils::basename(source);
  string target = _outdir + "/" + basename + ".o";

  // check dependencies
  string depfile = _outdir + "/" + basename + ".d";
  ostringstream cmdline;
  cmdline << _executable << " " << _preprocessor_flags
    << " " << _compiler_flags
      << " -MG -MM -MF " << depfile << " " << source;
  system(cmdline.str().c_str());

  // the source is the main dependency of the target
  vector<string> dependencies(1, source);

  ifstream dep(depfile.c_str());
  while (dep.good() && !dep.eof())
  {
    string buf;
    getline(dep, buf);
    vector<string> deplist;
    Utils::tokenize(buf, deplist, " \\");
    if (deplist.size() == 3)
    {
      dependencies.push_back(deplist[1]);
      dependencies.push_back(deplist[2]);
    }
    else if (deplist.size() == 1)
      dependencies.push_back(deplist[0]);
  }



  if (_needs_build(target, dependencies))
  {
    ostringstream cmdline;
    cmdline << _executable << " " << _preprocessor_flags << " " <<
        _compiler_flags << " " << flags << " -o " <<
        target << " " << source;


    cout << "...\n";
    if (BuildModule::verbose) cout << cmdline.str() << endl;
    if (system(cmdline.str().c_str()) != 0)
      target = "";
  }
  else
    cout << "- nothing to be done\n";

  return target;
}



int
Compiler::_link(const string& target, const vector<string>& objects,
    const string& flags)
{
  cout << "  Linking " << target << " (" << _options.get_key() << ")";

  int returnval = 0;
  if (_needs_build(target, objects))
  {
    cout << " ...\n";
    ostringstream cmdline;
    cmdline << _executable;
    for (size_t i = 0; i < objects.size(); ++i)
      cmdline << " " << objects[i];

    cmdline << " " << " " <<_linker_flags << " " << flags << " -o " << target;

    if (BuildModule::verbose) cout << cmdline.str() << endl;

    returnval = system(cmdline.str().c_str());
  }
  else
    cout << "- nothing to be done\n";

  return returnval;
}



void
check_dependency(const ModelOptions& options, const string& modulename,
    string& compileflags, string& linkflags)
{
  using namespace boost::filesystem;

  ostringstream pre;

  // the path is relative or absolute
  //string path = options.get_option("path", BuildModule::tc_root);
  string path = options.get_option("path", ".");
  BuildModule::replace(path, "@ROOT", BuildModule::tc_root);
  BuildModule::replace(path, "@ARCH", ARCH);
  BuildModule::replace(path, "@MODULE", modulename);
  //if (path[0] != '/')
  //  path = BuildModule::tc_root + "/" + path;

  vector<string> includes;
  options.get_option("includes", includes);
  for (size_t i = 0; i < includes.size(); ++i)
  {
    string inc = includes[i];
    BuildModule::replace(inc, "@ROOT", BuildModule::tc_root);
    BuildModule::replace(inc, "@ARCH", ARCH);
    BuildModule::replace(inc, "@MODULE", modulename);
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
    BuildModule::replace(libpath, "@MODULE", modulename);
    //if (libpath[0] != '/')
    //  linkerflags << "-L " << BuildModule::tc_root << "/" << libpath << " ";
    //else
      linkerflags << "-L " << libpath << " ";
  }

  vector<string> ldflags;
  options.get_option("linker_flags", ldflags);
  for (size_t i = 0; i < ldflags.size(); ++i)
  {
    string flags = ldflags[i];
    BuildModule::replace(flags, "@ROOT", BuildModule::tc_root);
    BuildModule::replace(flags, "@ARCH", ARCH);
    BuildModule::replace(flags, "@MODULE", modulename);
    linkerflags << flags << " ";
  }

  // we do not use this for now, but we assume that all
  // dependencies are already installed to their final
  // destination
  bool install = options.get_option("install", "true");
  options.get_option("libraries", ldflags);
  for (size_t i = 0; i < ldflags.size(); ++i)
  {
    string lib = ldflags[i];
    if (Utils::file_extension(lib).empty())
    {
      lib = " -l" + lib;
      if (install)
      {
        //copy_file(
      }
    }
    else
    {
      BuildModule::replace(lib, "@ROOT", BuildModule::tc_root);
      BuildModule::replace(lib, "@ARCH", ARCH);
      BuildModule::replace(lib, "@MODULE", modulename);
    }
    linkerflags << lib << " ";
  }

  compileflags += " " + pre.str();
  linkflags += " " + linkerflags.str();

}



void process_module(const string& , const ModelOptions& options)
{


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

  // the module library, if present
  string modulelib;

  string instpath = options.get_option("installpath", "@ROOT/@ARCH/modules/@MODULE");
  BuildModule::replace(instpath, "@ARCH", ARCH);
  BuildModule::replace(instpath, "@ROOT", BuildModule::tc_root);
  BuildModule::replace(instpath, "@MODULE", modulename);

  string arch_str(ARCH);

  // the compile and link flags
  string cpflags = options.get_option("compiler_flags", "");
  string ldflags = options.get_option("linker_flags", "");
  if ((arch_str == "i686-linux") || (arch_str == "x86_64-linux"))
    ldflags = " -Wl,--as-needed " + ldflags + " ";


  // process dependencies
  ModelOptions::const_submodel_iterator it(options.submodels_begin("Dependency"));
  ModelOptions::const_submodel_iterator end(options.submodels_end("Dependency"));
  for ( ; it != end; ++it)
  {
    // add installpath to options
    ModelOptions opts(it->second);
    string installpath = opts.get_option("installpath", instpath);
    opts.set_option("installpath", installpath);
    check_dependency(opts, modulename, cpflags, ldflags);
  }


  vector<string> sources;
  options.get_option("sources", sources);

  if (sources.size() > 0)
  {

    string creatable = options.get_option("createable", "");

    if (creatable.empty())
    {
      // guess it from the first source file
      creatable = Utils::basename(sources[0]);
    }

    //char* root = getenv("TIBERCADROOT");


    string compileflags = options.get_option("compiler_flags", "");
    compileflags += cpflags;
    compileflags += " -DMODULE_NAME=" + module;
    compileflags += " -DCREATABLE=" + creatable;

    string module_code = options.get_option("creator_code", "");
    if (!module_code.empty())
    {
      ofstream of(string(module + "_creator.h").c_str());
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

    string binsuffix;
    string libsuffix(".so");
    if ((arch_str == "i686-w64-mingw32") || (arch_str == "x86_64-w64-mingw32"))
    {
      libsuffix = ".dll";
      binsuffix = ".exe";
    }
    else if ((arch_str == "i686-darwin") || (arch_str == "x86_64-darwin"))
      libsuffix = ".dylib";

    modulelib = modulename + libsuffix;

    string linkflags(ldflags);
    if ((arch_str == "i686-linux") || (arch_str == "x86_64-linux"))
      linkflags += " -Wl,-soname," + modulelib + " ";

    using namespace boost::filesystem;

    // build full installation path
    path instpath_p(instpath);
    if (!exists(instpath_p))
      create_directories(instpath_p);
    modulelib = instpath + "/" + modulelib;

    Compiler::link(modulelib, objects, linkflags);
  }

  // recursively process submodules
  it = options.submodels_begin();
  end = options.submodels_end();
  for ( ; it != end; ++it)
  {
    if (it->first == "Dependency") continue;

    ModelOptions opts(it->second);
    if (!opts.find_option("installpath"))
      opts["installpath"] = instpath;

    opts["compiler_flags"] = cpflags;

    opts["linker_flags"] = ldflags + modulelib;

    opts["module"] = module;

    process_module(it->first, opts);
  }


}

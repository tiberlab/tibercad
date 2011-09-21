// $Id: main.C 2874 2011-07-15 21:12:36Z maufder $

#include "boost/algorithm/string/trim.hpp"

#include "TiberCad.h"
#include "Utils.h"
#include "InputParser.h"
#include "ModelOptions.h"

#include "tiber_config.h"


#include <iostream>
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
  bool interactive;

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
}

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
      string tc_root("TIBERCADROOT=" + exepath);
      putenv(tc_root.c_str());
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

  char* root = getenv("TIBERCADROOT");
  if (root == NULL)
  {
    cerr << "TIBERCADROOT environment variable is not set.\n";
    return 1;
  }


  ModelOptions config;
  InputParser parser;
  parser.parse_file(inputfile, config);
  config.print_all();

  string name = "test";

#if defined(_WIN32)
  string binsuffix(".exe");
#else
  string binsuffix("");
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

  string cxxflags;
  cxxflags += "-std=gnu++0x -pthread ";
#if !defined(_WIN32)
  cxxflags += "-fPIC ";
#endif

  string ldflags;
  ldflags += "-shared ";
  ldflags += "-L" + tcroot + "/" + ARCH + "/lib ";
  ldflags += "-Wl,-rpath,\'$$ORIGIN\' -Wl,-rpath,\'$$ORIGIN/../lib\' ";

  cout << cpp << endl << cppflags << endl << ldflags << endl;

  string cmdline(cpp);
  cmdline += " " + cppflags + " " + cxxflags + " " + ldflags + " -o " + name + " " + name + ".C";
  cout << "Compiling " << name << ".C ...";
  //cout << cmdline << endl;
  system(cmdline.c_str());
 
  return error;
}


// $Id$

#include "boost/algorithm/string/trim.hpp"

#include "DLLoader.h"
#include "TiberCad.h"
#include "License.h"
#include "Utils.h"
#include "Messages.h"

#include "tiber_config.h"

#include "auto_ptr.h"

#include <iostream>
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
static string open_file(const char *filter = "tiberCAD input files (*.tib)\0*.tib\0", HWND owner = NULL)
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
  bool stop_on_warning = false;

  // TODO doe sthis work in windows?
  //AutoPtr<ofstream> nullstream(NULL);

  void usage(void)
  {
#if defined(_WIN32)
    cout << endl << "Usage:" << endl
      << "  from command line: tibercad [options] inputfile" << endl
      << "  or double click on inputfile" << endl << endl;
# else
    cout << endl << "Usage: tibercad [options] inputfile" << endl << endl;
# endif

    cout << "Options:" << endl
         << "  -b      batch mode" << endl
         << "  -i      interactive mode" << endl
         << "  -s      stop at warnings (disabled in batch mode)" << endl
         << "  -v      version info" << endl;

#if defined(_WIN32)
    cout << "press Enter ...";
    if (interactive) cin.get();
# endif
  }
}

// Will be extended with tools for command line argument parsing
// and so on
int main (int argc, char** argv)
{

  /*
   * As first thing, we start MPI
   */
  MPI_Init(&argc, &argv);

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  // on rank other than 0 we set stdout to the nullstream
  if (my_rank != 0)
  {
    //nullstream.reset(new ofstream("/dev/null"));
    Messages::set_stdout();
  }


#ifdef _WIN32
  interactive = true;
#else
  interactive = false;
#endif

  opterr = 0;
  int c;
  while ((c = getopt(argc, argv, "bivs")) != -1)
    switch (c)
    {
      case 'v':
        cout << "tiberCAD release " << TiberCad::version_string()
          << "(rev. " << TiberCad::software_revision() << ", "
                  << TiberCad::arch_string() << ")" << endl;
        return 0;
        break;

      case 'b':
        interactive = false;
        stop_on_warning = false;
        break;

      case 'i':
        interactive = true;
        break;

      case 's':
        stop_on_warning = true;
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


#ifdef LICENSE_CHECK
    // check the license
    if (!License::check_license())
    {
      cerr << "Sorry, cannot start TiberCAD as I could not find "
          "a valid license." << endl;
#if defined(_WIN32)
      cout << endl << "press Enter ...";
      if (interactive) cin.get();
# endif
      return 1;
    }
#endif

  }

  //
  // here begins real TiberCAD
  //
  int error = 1;


  Messages::interactive() = interactive;
  Messages::stop_on_warning() = stop_on_warning;

  {
    /*
     * We have this inside a defined scope so that no
     * MPI-using objects will be destroyed after MPI_Finalize
     */

    // Create the entry point object
    // here we pass MPI_COMM_WORLD because we are in the main
    TiberCad tibercad(MPI_COMM_WORLD);

    try {

      tibercad.init(inputfile);

      tibercad.run();

      Messages::info("Simulation finished.");

      error = 0;
    }
    catch (exception& e)
    {
      Messages::error(e.what());
    }
    catch (...)
    {
      Messages::error("TiberCAD crashed for unknown reason.");
    }
#if defined(_WIN32)
    cout << "press Enter ...";
    if (interactive) cin.get();
#endif

    Messages::print_statistics();
    Messages::close_log_file();

    Messages::info("Goodbye");
  }


  /*
   * As last thing, finalize MPI
   */
  MPI_Finalize();

  return error;
}


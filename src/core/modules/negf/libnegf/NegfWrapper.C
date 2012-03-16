#include "NegfWrapper.h"
#include <cstring>

//---------------------------------------------------------------------
NegfWrapper* NegfWrapper::create()
{
  return new NegfWrapper();
}



NegfWrapper::NegfWrapper()
{

    std::cout << "\nCreating libNEGF instance... ";
    f77_negf_init_session(_handler);
    std::cout << "done." << std::endl;

    for  (int ii = 0; ii < NEGF_HSIZE; ++ii) {
      std::cout << _handler[ii] << " ";
    }
    std::cout << std::endl;
    f77_negf_getversion(_handler);

}


NegfWrapper::~NegfWrapper()
{
    f77_negf_destruct_libnegf(_handler);
    f77_negf_destruct_session(_handler);
}

int
NegfWrapper::init()
{
  f77_negf_init(_handler);
  return 0;
}



void NegfWrapper::set_verbose(int verbose_lev)
{

  f77_negf_set_verbosity(_handler,verbose_lev);

}

int
NegfWrapper::current()
{
  f77_negf_current(_handler);
  return 0;
}

//Compute charge density
int
NegfWrapper::density(std::vector<double>& density)
{
   int size = density.size();
   f77_negf_density(_handler, size, &density.front());
}

//Set SC iteration
void
NegfWrapper::set_iteration(int iter)
{
  f77_negf_set_iteration(_handler, iter);
}

void
NegfWrapper::set_scratch_path(std::string path)
{
  char *cpath = new char[NEGF_LC];
  memset(cpath,NEGF_PADCHAR,NEGF_LC);
  path.copy(cpath,path.size());

  f77_negf_set_scratch(_handler, cpath);
}

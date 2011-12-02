#include "NegfWrapper.h"


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



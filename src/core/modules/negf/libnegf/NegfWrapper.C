#include "NegfWrapper.h"

#include <cstring>
#include <iostream>

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
    f77_negf_get_version(_handler);

}


NegfWrapper::~NegfWrapper()
{
    f77_negf_destruct_libnegf(_handler);
    f77_negf_destruct_session(_handler);
}

void
NegfWrapper::init()
{
  f77_negf_init(_handler);
}


int
NegfWrapper::read_HS()
{
  // TODO new api needs file name
  //f77_negf_read_hs(_handler);

  return 0;
}

int
NegfWrapper::read_input()
{
  f77_negf_read_input(_handler);

  return 0;
}

void NegfWrapper::set_verbose(int verbose_lev)
{

  f77_negf_set_verbosity(_handler,verbose_lev);

}

double
NegfWrapper::current(std::string unitsOfH, std::string unitsOfJ)
{
  double current;
  char *unitH = new char[NEGF_SC];
  char *unitJ = new char[NEGF_SC];
  memset(unitH,NEGF_PADCHAR,NEGF_SC);
  unitsOfH.copy(unitH,unitsOfH.size());
  unitsOfJ.copy(unitJ,unitsOfJ.size());

  f77_negf_current(_handler, current, *unitH, *unitJ);

  return current;
}

void
NegfWrapper::density(std::vector<double>& density, std::string particle)
{
   int size = density.size();
   int p = 0;

   if (particle == "el"){ p = +1;}
   if (particle == "hl"){ p = -1;}
   
   f77_negf_density_efa(_handler, size, &density.front(), p);
}


void
NegfWrapper::ldos(std::vector<double>& ldos, int& esteps, int& npoints)
{
  int shape[2];
  double *d = ldos.data();
  f77_negf_associate_ldos(_handler, shape, &d);
  esteps = shape[0];
  npoints = shape[1];
}


void
NegfWrapper::set_iteration(int iter)
{
  // TODO this is missing in current API
  //f77_negf_set_iteration(_handler, iter);
}

void
NegfWrapper::set_kpoint(int kpoint)
{
  f77_negf_set_kpoint(_handler, kpoint);
}

void
NegfWrapper::set_scratch_path(std::string path)
{
  char *cpath = new char[NEGF_LC];
  memset(cpath,NEGF_PADCHAR,NEGF_LC);
  path.copy(cpath,path.size());

  f77_negf_set_scratch(_handler, cpath);
}

void
NegfWrapper::set_output_path(std::string path)
{
  char *cpath = new char[NEGF_LC];
  memset(cpath,NEGF_PADCHAR,NEGF_LC);
  path.copy(cpath,path.size());

  f77_negf_set_output(_handler, cpath);
}

void
NegfWrapper::set_reference(int minmax)
{
  f77_negf_set_reference(_handler,minmax);
}


void
NegfWrapper::clean_libnegf(void)
{
   f77_negf_destruct_libnegf(_handler);
}

void
NegfWrapper::device_contact_dm(int outer)
{
   f77_negf_set_outer(_handler, outer);
}

void
NegfWrapper::set_write_ldos(bool flag)
{
  if (flag) f77_negf_set_writeldos(_handler,1);
  else f77_negf_set_writeldos(_handler,0);
}

void
NegfWrapper::set_write_tunn(bool flag)
{
  if (flag) f77_negf_set_writetunn(_handler,1);
  else f77_negf_set_writetunn(_handler,0);
}


void
NegfWrapper::partition_info(void)
{
   f77_negf_write_partition(_handler);
}

void 
NegfWrapper::set_H_csr(int nrow, std::vector<Complex >& A, 
                           std::vector<int>& JA, std::vector<int>& IA)
{
  f77_negf_set_h(_handler,nrow,&A.front(),&JA.front(),&IA.front());
}

void 
NegfWrapper::set_S_csr(int nrow, std::vector<Complex >& A, 
                           std::vector<int>& JA, std::vector<int>& IA)
{
  f77_negf_set_s(_handler,nrow,&A.front(),&JA.front(),&IA.front());
}


void 
NegfWrapper::set_S_id(int nrow) 
{
  f77_negf_set_s_id(_handler,nrow);
}


void 
NegfWrapper::print_mat(void) 
{
  f77_negf_print_mat(_handler);
}


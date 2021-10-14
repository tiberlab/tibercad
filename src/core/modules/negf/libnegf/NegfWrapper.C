#include "NegfWrapper.h"

#include <cstring>
#include <iostream>

//---------------------------------------------------------------------
NegfWrapper* NegfWrapper::create(void)
{
  return new NegfWrapper();
}



NegfWrapper::NegfWrapper(void)
 : _is_initialized(false)
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


NegfWrapper::~NegfWrapper(void)
{
    f77_negf_destruct_libnegf(_handler);
    f77_negf_destruct_session(_handler);
}

void
NegfWrapper::init(void)
{
  if (!_is_initialized)
  {
    f77_negf_init(_handler);
    _is_initialized = true;
  }
}

void
NegfWrapper::force_reinit(void)
{
  if (_is_initialized)
    f77_negf_destruct_libnegf(_handler);

  init();
}



int
NegfWrapper::read_HS(void)
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

  f77_negf_current(_handler, current, unitsOfH.c_str(), unitsOfJ.c_str());

  return current;
}


void
NegfWrapper::get_energies(std::vector<double>& energies)
{
  int n_erg;
  std::vector<double> dummy;

  f77_negf_get_energies(_handler, n_erg, energies.data(), dummy.data(), 0);

  energies.resize(n_erg);
  dummy.resize(n_erg);

  f77_negf_get_energies(_handler, n_erg, energies.data(), dummy.data(), 1);
}


void
NegfWrapper::get_transmission(std::vector<std::vector<double>>& transmission)
{
  // # erg points, # T
  int shape[2];
  double *data;

  f77_negf_associate_transmission(_handler, shape, &data);

  transmission.resize(shape[1]);

  for (size_t i = 0; i< shape[1]; ++i)
  {
    transmission[i].resize(shape[0]);
    for (size_t j = 0; j < shape[0]; ++j)
      transmission[i][j] = data[i*shape[0] + j];
  }
}

void
NegfWrapper::get_energy_current(std::vector<std::vector<double>>& current)
{
  // # erg points, # J
   int shape[2];
   double *data;

   f77_negf_associate_energy_current(_handler, shape, &data);

   current.resize(shape[1]);

   for (size_t i = 0; i< shape[1]; ++i)
   {
     current[i].resize(shape[0]);
     for (size_t j = 0; j < shape[0]; ++j)
       current[i][j] = data[i*shape[0] + j];
   }
}

void
NegfWrapper::density(std::vector<double>& density, std::string particle)
{
   int size = density.size();
   int p = 0;

   if (particle == "el"){ p = +1;}
   if (particle == "hl"){ p = -1;}
   
   f77_negf_density_efa(_handler, size, density.data(), p);
}


void
NegfWrapper::get_ldos(std::vector<std::vector<double>>& ldos)
{
  int shape[2];
  double *data;
  f77_negf_associate_ldos(_handler, shape, &data);

  ldos.resize(shape[0]);

   for (size_t i = 0; i< shape[0]; ++i)
   {
     ldos[i].resize(shape[1]);
     for (size_t j = 0; j < shape[1]; ++j)
       ldos[i][j] = data[j*shape[0] + i];
   }
}

void
NegfWrapper::init_ldos(unsigned int nldos)
{
  f77_negf_init_ldos(_handler, nldos);
}


void
NegfWrapper::init_ldos(const std::vector<int>& start,
                       const std::vector<int>& end)
{
  int nldos = start.size();
  if (end.size() != nldos)
    throw std::runtime_error("trying to initialize ldos in libnegf "
        "with incompatible index arrays");

  std::vector<int> fstart(start);
  std::vector<int> fend(end);
  for (size_t i = 0; i < nldos; ++i)
  {
    ++fstart[i];
    ++fend[i];
  }

  f77_negf_init_ldos(_handler, nldos);
  f77_negf_set_ldos_intervals(_handler, nldos, start.data(), end.data());
}

void
NegfWrapper::set_ldos_indices(unsigned int dos_index,
                              const std::vector<int>& indices)
{
  // adjust to Fortran indices
  std::vector<int> ids(indices);
  for (auto&& id : ids) ++id;

  f77_negf_set_ldos_indexes(_handler, dos_index+1, ids.size(), ids.data());
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
  f77_negf_set_scratch(_handler, path.c_str());
}

void
NegfWrapper::set_output_path(std::string path)
{
  f77_negf_set_output(_handler, path.c_str());
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
NegfWrapper::init_contacts(int n_cont)
{
  f77_negf_init_contacts(_handler, n_cont);
}

void
NegfWrapper::init_structure(int ncont, const std::vector<int>& surfstart,
    const std::vector<int>& surfend, const std::vector<int>& contend, int npl,
    const std::vector<int>& plend)
{
  f77_negf_init_structure(_handler, ncont, surfstart.data(),
      surfend.data(), contend.data(), npl, plend.data());
}

void
NegfWrapper::get_parameters(NegfWrapper::Parameters& params)
{
  f77_negf_get_params(_handler, params);
}

void
NegfWrapper::set_parameters(const NegfWrapper::Parameters& params)
{
  f77_negf_set_params(_handler, params);
}



void
NegfWrapper::partition_info(void)
{
   f77_negf_write_partition(_handler);
}

void 
NegfWrapper::set_H_csr(int nrow, std::vector<Complex>& A,
                           std::vector<int>& JA, std::vector<int>& IA)
{
  f77_negf_set_h(_handler,nrow,A.data(),JA.data(),IA.data());
}

void 
NegfWrapper::set_S_csr(int nrow, std::vector<Complex >& A, 
                           std::vector<int>& JA, std::vector<int>& IA)
{
  f77_negf_set_s(_handler,nrow,A.data(),JA.data(),IA.data());
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


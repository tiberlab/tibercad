#include "NegfWrapper.h"

#include <cstring>
#include <iostream>

//---------------------------------------------------------------------

namespace
{
  void TwoDArray2Fortran(const std::vector<std::vector<double>> array_2d, std::vector<double>& out_array, 
       int& rows, int& columns, bool transpose=false)
  {
    // Note: Fortran uses column-major order  
    out_array.resize(rows*columns);
    if (!transpose)
    {
      for (int i = 0; i<rows; i++) 
      {
        for (int j = 0; j < columns; j++) 
        {
            out_array[j * rows + i] = array_2d[i][j];    
        }
      }
    }
    else
    {
      for (int i = 0; i<rows; i++) 
      {
        for (int j = 0; j < columns; j++) 
        {
            out_array[i * columns + j] = array_2d[i][j];    
        }
      }
      std::swap(rows, columns);
    }

    return;
  }
}


namespace
{
  void test_mpi_reduce(MPI_Comm comm, std::string name)
  {
    int rank; int size; int err;
    
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    int ARRAY_SIZE = 5;
    double test_array[ARRAY_SIZE];
    for (int i = 0; i < ARRAY_SIZE; ++i) {
        test_array[i] = 1.01;
    }
    std::cout << std::endl << "TEST REDUCE (C++), communicator: " << name << std::endl;
    std::cout << name << ".id = " << comm <<std::endl;
    std::cout << name << ".size = " << size <<std::endl;
    std::cout << name << ".rank = " << rank <<std::endl;
    std::cout << "Array before calling reduce:" << std::endl;
    for (int i = 0; i < ARRAY_SIZE; i++) std::cout << test_array[i] << ", ";
    std::cout << std::endl;
    std::cout << "Calling reduce..." << std::endl;
    err = MPI_Reduce(MPI_IN_PLACE, test_array, ARRAY_SIZE, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm);
    std::cout << "Array after calling reduce:" << std::endl;
    for (int i = 0; i < ARRAY_SIZE; i++) std::cout << test_array[i] << ", ";
    std::cout<<std::endl;
    std::cout << "MPI_IN_PLACE:" << MPI_IN_PLACE <<std::endl;
    std::cout<< "ERROR FLAG: " << err << std::endl;
    std::cout << "END TEST REDUCE" << std::endl << std::endl;

  }
}


NegfWrapper* NegfWrapper::create(void)
{
  return new NegfWrapper();
}



NegfWrapper::NegfWrapper(void)
 : _is_initialized(false)
{

    std::cout << "\nCreating libNEGF instance... ";
    negf_init_session(_handler);
    std::cout << "done." << std::endl;

    for  (int ii = 0; ii < NEGF_HSIZE; ++ii) {
      std::cout << _handler[ii] << " ";
    }
    std::cout << std::endl;
    negf_get_version(_handler);

}


NegfWrapper::~NegfWrapper(void)
{
    negf_destruct_libnegf(_handler);
    negf_destruct_session(_handler);
}

void
NegfWrapper::init(void)
{
  if (!_is_initialized)
  {
    negf_init(_handler);
    _is_initialized = true;
  }
}

int NegfWrapper::set_mpi_comm(MPI_Comm comm)
{
    {
      MPI_Fint global_comm_f = MPI_Comm_c2f(comm);

      negf_set_mpi_fcomm(_handler, global_comm_f);
      return 0;
    }
}

void NegfWrapper :: mpi_cart_init(MPI_Comm inComm, int nGroups, MPI_Comm& cartComm, MPI_Comm& kComm, MPI_Comm& enComm)
{   
    MPI_Fint cartComm_F;
    MPI_Fint kComm_F;
    MPI_Fint enComm_F;
    MPI_Fint global_comm_f = MPI_Comm_c2f(inComm);

    negf_cartesian_init(_handler, global_comm_f, nGroups, cartComm_F, kComm_F, enComm_F);

    cartComm = MPI_Comm_f2c(cartComm_F);
    kComm = MPI_Comm_f2c(kComm_F);
    enComm = MPI_Comm_f2c(enComm_F);

}

void
NegfWrapper::force_reinit(void)
{
  if (_is_initialized)
    negf_destruct_libnegf(_handler);

  init();
}



int
NegfWrapper::read_HS(std::string real_path, std::string imag_path, std::string H_or_S)
{
  if (H_or_S == "H")
  {
    negf_read_hs(_handler, real_path.c_str(), imag_path.c_str(), 0);
  }
  else if (H_or_S == "S")
  {
    negf_read_hs(_handler, real_path.c_str(), imag_path.c_str(), 1);
  }
  else
  {
    throw std::runtime_error("In read_HS: argument H_or_S must be either string 'H' or string 'S'.");
  }

  return 0;
}

int
NegfWrapper::read_input()
{
  negf_read_input(_handler);

  return 0;
}

void NegfWrapper::set_verbose(int verbose_lev)
{

  negf_set_verbosity(_handler,verbose_lev);

}

void NegfWrapper::print_memory_statistics(void)
{
  negf_mem_stats(_handler);
}


double
NegfWrapper::current(std::string unitsOfH, std::string unitsOfJ)
{
  double current;

  negf_current(_handler, current, unitsOfH.c_str(), unitsOfJ.c_str());

  return current;
}


void
NegfWrapper::layer_current(std::vector<double>& layer_current, std::string unitsOfH, std::string unitsOfJ)
{
   int size = layer_current.size();
   
   negf_layer_current(_handler, size, layer_current.data(), unitsOfH.c_str(), unitsOfJ.c_str());
}


void
NegfWrapper::get_energies(std::vector<double>& energies)
{
  int n_erg;
  std::vector<double> dummy;

  negf_get_energies(_handler, n_erg, energies.data(), dummy.data(), 0);

  energies.resize(n_erg);
  dummy.resize(n_erg);

  negf_get_energies(_handler, n_erg, energies.data(), dummy.data(), 1);
}


void
NegfWrapper::get_transmission(std::vector<std::vector<double>>& transmission)
{
  // # erg points, # T
  int shape[2];
  double *data;

  negf_associate_transmission(_handler, shape, &data);

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

   negf_associate_energy_current(_handler, shape, &data);

   current.resize(shape[1]);

   for (size_t i = 0; i< shape[1]; ++i)
   {
     current[i].resize(shape[0]);
     for (size_t j = 0; j < shape[0]; ++j)
       current[i][j] = data[i*shape[0] + j];
   }
}

void
NegfWrapper::density(std::vector<double>& density, std::string particle, bool contact_calc)
{
   int size = density.size();
   int p = 0;

   if (particle == "el"){ p = +1;}
   if (particle == "hl"){ p = -1;}
   
   negf_density_efa(_handler, size, density.data(), p, contact_calc);
}

void
NegfWrapper::quasi_equilibrium_density(std::vector<double>& density, std::string particle,
                                       std::vector<double>& Ec, std::vector<double>& Ev,
                                       std::vector<double>& muN, std::vector<double>& muP)
{

   int size = density.size();
   int p = 0;

   if (particle == "el"){ p = +1;}
   if (particle == "hl"){ p = -1;}

   negf_density_quasi_equilibrium(_handler, size, density.data(), p, Ec.data(),
      Ev.data(), muN.data(), muP.data());
}

void
NegfWrapper::get_ldos(std::vector<std::vector<double>>& ldos)
{
  int shape[2];
  double *data;
  negf_associate_ldos(_handler, shape, &data);

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
  negf_init_ldos(_handler, nldos);
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

  negf_init_ldos(_handler, nldos);
  negf_set_ldos_intervals(_handler, nldos, start.data(), end.data());
}

void
NegfWrapper::set_ldos_indices(unsigned int dos_index,
                              const std::vector<int>& indices)
{
  // adjust to Fortran indices
  std::vector<int> ids(indices);
  for (auto&& id : ids) ++id;

  negf_set_ldos_indexes(_handler, dos_index+1, ids.size(), ids.data());
}


void
NegfWrapper::set_iteration(int iter)
{
  // TODO this is missing in current API
  //negf_set_iteration(_handler, iter);
}

void
NegfWrapper::set_kpoint(int kpoint)
{
  negf_set_kpoint(_handler, kpoint);
}

void
NegfWrapper::set_kpoints(std::vector<std::vector<double>> kpoints, std::vector<double> kweights,  
             std::vector<int> local_k_indices, std::vector<std::vector<double>> equivalent_points, 
             std::vector<int> n_equivalent, bool reduced_BZ)
{
  int set_eq_pts = (reduced_BZ) ? 1 : 0;

  std::vector<double> kpoints_f;
  int rows = kpoints.size(); 
  int columns = kpoints[0].size();

  // We need the number of kpoints to be the second dimension, so we transpose
  TwoDArray2Fortran(kpoints, kpoints_f, rows, columns, true);

  std::vector<double> equivalent_points_f;
  int N = equivalent_points.size();
  int M = kpoints[0].size(); // equivalent_points[0] could be empty. Dimension of a point should be 3 no matter the point.
  TwoDArray2Fortran(equivalent_points, equivalent_points_f, N, M, true);

  //Fortran is 1-indexed: shift local_k_indices of +1
  for (auto&& ind : local_k_indices) ind += 1;

  negf_set_kpoints(_handler, kpoints_f.data(), rows, columns, kweights.data(), local_k_indices.data(),
                      local_k_indices.size(), equivalent_points_f.data(), M, n_equivalent.data(), set_eq_pts);
}

void
NegfWrapper::init_basis(std::vector<std::vector<double>> coordinates, std::vector<int> matrix_indices,
             std::vector<std::vector<double>> lattice_vectors, int transport_direction)
{
  std::vector<double> coords_f;
  int dims = coordinates.size();
  int n_atoms = coordinates[0].size();
  TwoDArray2Fortran(coordinates, coords_f, dims, n_atoms);

  std::vector<double> lattice_vecs_f;
  int n_vecs = lattice_vectors[0].size();
  TwoDArray2Fortran(lattice_vectors, lattice_vecs_f, dims, n_vecs);

  negf_init_basis(_handler, coords_f.data(), n_atoms, dims, matrix_indices.data(), 
                     lattice_vecs_f.data(), n_vecs, transport_direction);
}

void
NegfWrapper::set_scratch_path(std::string path)
{
  negf_set_scratch(_handler, path.c_str());
}

void
NegfWrapper::set_output_path(std::string path)
{
  negf_set_output(_handler, path.c_str());
}

void
NegfWrapper::set_reference(int minmax)
{
  negf_set_reference(_handler,minmax);
}


void
NegfWrapper::clean_libnegf(void)
{
   negf_destruct_libnegf(_handler);
}

void
NegfWrapper::device_contact_dm(int outer)
{
   negf_set_outer(_handler, outer);
}


void
NegfWrapper::init_contacts(int n_cont)
{
  negf_init_contacts(_handler, n_cont);
}

void
NegfWrapper::init_structure(int ncont, const std::vector<int>& surfstart,
    const std::vector<int>& surfend, const std::vector<int>& contend, int npl,
    const std::vector<int>& plend,const std::vector<int>& cblks)
{
  negf_init_structure(_handler, ncont, surfstart.data(),
      surfend.data(), contend.data(), npl, plend.data(), cblks.data());
}

std::vector<int>
NegfWrapper::contact_blocks(int ncont, const std::vector<int>& surfstart,
    const std::vector<int>& surfend, const std::vector<int>& contend, int npl,
    const std::vector<int>& plend)
{
  std::vector<int> cblks(ncont);
  negf_contact_blocks(_handler, ncont, surfstart.data(),
      surfend.data(), contend.data(), npl, plend.data(), cblks.data());

  return cblks;
}

void
NegfWrapper::get_parameters(NegfWrapper::Parameters& params)
{
  negf_get_params(_handler, params);
}

void
NegfWrapper::set_parameters(const NegfWrapper::Parameters& params)
{
  negf_set_params(_handler, params);
}



void
NegfWrapper::partition_info(void)
{
   negf_write_partition(_handler);
}

void 
NegfWrapper::set_H_csr(int nrow, std::vector<Complex>& A,
                           std::vector<int>& JA, std::vector<int>& IA, int iK)
{
  negf_set_h(_handler,nrow,A.data(),JA.data(),IA.data(), iK);
}

void 
NegfWrapper::create_HS_container(int n_Hk)
{
  negf_create_hs(_handler, n_Hk);
}

void 
NegfWrapper::set_S_csr(int nrow, std::vector<Complex >& A, 
                           std::vector<int>& JA, std::vector<int>& IA, int iK)
{
  negf_set_s(_handler,nrow,A.data(),JA.data(),IA.data(),iK);
}


void 
NegfWrapper::set_S_id(int nrow, int iK) 
{
  negf_set_s_id(_handler,nrow,iK);
}


void 
NegfWrapper::print_mat(void) 
{
  negf_print_mat(_handler);
}


void
NegfWrapper::set_elph_dephasing(std::vector<double> coupling, int scba_niter)
{
  negf_set_elph_dephasing(_handler, coupling.data(), coupling.size(), scba_niter);
}


void
NegfWrapper::set_elph_block_dephasing(std::vector<double> coupling, std::vector<int> orbsperatm, int scba_niter)
{
  negf_set_elph_block_dephasing(_handler, coupling.data(), coupling.size(), orbsperatm.data(), orbsperatm.size(), scba_niter);
}


void
NegfWrapper::set_elph_s_dephasing(std::vector<double> coupling, std::vector<int>orbsperatm, int scba_niter)
{
  negf_set_elph_s_dephasing(_handler, coupling.data(), coupling.size(), orbsperatm.data(), orbsperatm.size(), scba_niter);
}


void
NegfWrapper::set_elph_polaroptical(std::vector<double> coupling,  double wq, double kbT, double deltaz, double eps_r,  double eps_inf, 
             double q0, double cell_area, int scba_niter,  bool tTridiagonal)
{
  negf_set_elph_polaroptical(_handler, coupling.data(), coupling.size(), wq, kbT, deltaz, eps_r, eps_inf, q0, cell_area, 
                                scba_niter, tTridiagonal);
}


void
NegfWrapper::set_elph_nonpolaroptical(std::vector<double> coupling, double wq, double kbT, double deltaz, double D0, double cell_area, 
             int scba_niter, bool tTridiagonal)
{
  negf_set_elph_nonpolaroptical(_handler, coupling.data(), coupling.size(), wq, kbT, deltaz, D0, cell_area, scba_niter, tTridiagonal);
}


void
NegfWrapper::set_scba_tolerances(double elastic_tol, double inelastic_tol)
{
  negf_set_scba_tolerances(_handler, elastic_tol, inelastic_tol);
}
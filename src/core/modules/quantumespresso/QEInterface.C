// $Id$
//
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>

#include "QEInterface.h"
#include "AtomisticStructure.h"
#include "TeeStream.h"

#include "TiberModule.h"

using namespace std;

QEInterface::QEInterface(const ModelOptions& options) :
  SimulationInterface(options),
  _qe_pseudo_dir("."),
  _qe_ecutwfc(20.0),
  _qe_conv_thr(1e-6),
  _qe_nbnd(0),
  _qe_k_points({"automatic", "1 1 1 0 0 0"}),
  _qe_outdir("out")
{
}


QEInterface::~QEInterface(void)
{
}


QEInterface*
QEInterface::create(const ModelOptions& options)
{
  return(new QEInterface(options));
}


void
QEInterface::do_init(void)
{
  if (this->get_atomistic_structure() == nullptr)
    throw InitFailedException("QuantumEspresso interface needs an atomistic structure.");

  _qe_pseudo_dir = get_option("pseudo_dir", _qe_pseudo_dir);
  _qe_ecutwfc = get_option("ecutwfc", _qe_ecutwfc);
  _qe_conv_thr = get_option("conv_thr", _qe_conv_thr);
  _qe_nbnd = get_option("nbnd", _qe_nbnd);

  get_option("k_points", _qe_k_points);
  if (_qe_k_points.size() < 2)
    throw InitFailedException("You need to provide the definitions for the kpoints.");

  vector<string> pseudos;
  get_option("pseudopotentials", pseudos);
  for (unsigned int i = 0; i < pseudos.size(); i = i + 2)
  {
    _qe_pseudos[pseudos[i]] = pseudos[i+1];
    cerr << pseudos[i] << "  " << pseudos[i+1] << endl;
  }

  _qe_outdir = get_option("outdir", _qe_outdir);
}


void
QEInterface::do_solve(void)
{
  string prefix = this->get_atomistic_structure()->get_name();
  string outdir = this->get_output_directory() + "/QE";
  string tmpdir = this->get_scratch_directory() + "/QE";

  using namespace boost::filesystem;
  
  // check output and scratch path
  path outpath(outdir);
  if (!exists(outpath))
  {
    // we catch any error here without doing anything yet
    try {
      create_directories(outpath);
    }
    catch (...) {}
  }

  if (!(exists(outpath) && is_directory(outpath)))
  {
    string msg("Cannot create or use '");
    msg += outpath.string() + "' as output directory.";
    throw InitFailedException(msg);
  }


  string qe_file_base = outdir + "/" + prefix;

  ofstream qe_scf(qe_file_base + ".scf.in");
  ofstream qe_nscf(qe_file_base + ".nscf.in");

  TeeStream qe_both(qe_scf, qe_nscf);

  qe_both << "&control" << endl;
  qe_scf  << "   calculation = 'scf'," << endl
          << "   restart_mode='from_scratch'," << endl;
  qe_nscf << "   calculation = 'nscf'," << endl;
  qe_both << "   prefix = '" << prefix << "'," << endl
          << "   pseudo_dir = '" << _qe_pseudo_dir << "'," << endl
          << "   outdir = '" << _qe_outdir << "'," << endl
          << "/" << endl;

  // for now this is fixed to tetragonal cells
  int ibrav = 8;
  libMesh::RealVectorValue a, b, c;
  this->get_atomistic_structure()->get_lattice_vectors(a, b, c);
  int nat = this->get_atomistic_structure()->get_N_atoms();
  int ntyp = this->get_atomistic_structure()->get_N_types(); 
  qe_both << "&system" << endl
          << "   ibrav = " << ibrav << "," << endl
          << "   a = " << a(0) << "," << endl
          << "   b = " << b(1) << "," << endl
          << "   c = " << c(2) << "," << endl
          << "   nat = " << nat << ", ntyp = " << ntyp << "," << endl
          << "   ecutwfc = " << _qe_ecutwfc << "," << endl;
  if (_qe_nbnd > 0)
    qe_nscf << "   nbnd = " << _qe_nbnd << "," << endl;
  qe_both << "/" << endl;

  qe_both << "&electrons" << endl
          << "   conv_thr = " << _qe_conv_thr << "," << endl
          << "/" << endl;

  const vector<string>& atom_types = this->get_atomistic_structure()->get_atom_types();

  qe_both << "ATOMIC_SPECIES" << endl;
  for (unsigned int i = 0; i < ntyp; ++i)
  {
    qe_both << atom_types[i] << " " << Specie(atom_types[i]).get_mass()
      << " " << _qe_pseudos[atom_types[i]] << endl;
  }

  const vector<Atom> atoms = this->get_atomistic_structure()->get_structure_atoms();
  double dx = 0,
         dy = 0,
         dz = 0;

  if (get_option("shift_origin", true))
  {
    for (unsigned int i = 0; i < nat; ++i)
    {
      if (atoms[i].get_position(0) < dx)
        dx = atoms[i].get_position(0);

      if (atoms[i].get_position(1) < dy)
        dy = atoms[i].get_position(1);

      if (atoms[i].get_position(2) < dz)
        dz = atoms[i].get_position(2);
    }

    dx -= 1e-6;
    dy -= 1e-6;
    dz -= 1e-6;
  }

  qe_both << "ATOMIC_POSITIONS angstrom" << endl;
  for (unsigned int i = 0; i < nat; ++i)
  {
    qe_both << atoms[i].get_specie() << "  "
            << atoms[i].get_position(0) - dx << " "
            << atoms[i].get_position(1) - dy << " "
            << atoms[i].get_position(2) - dz << " "
            << endl;
  }
  
  qe_both << "K_POINTS ";
  for (unsigned int i = 0; i < _qe_k_points.size(); ++i)
  {
    qe_scf << _qe_k_points[i] << endl;
  }
  qe_nscf << "crystal_b\n1\n0 0 0 1\n";
}

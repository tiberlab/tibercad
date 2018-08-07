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
  _qe_mixing_beta(0),
  _qe_nbnd(0),
  _qe_k_points_scf({"Gamma"}), //{"automatic", "1 1 1 0 0 0"}),
  _qe_k_points_nscf({"Gamma"}), //{"automatic", "1 1 1 0 0 0"}),
  _qe_outdir("out"),
  _qe_degauss(0.001), //
  _qe_DeltaE(0.01), //
  _qe_Emin(4.0), //
  _qe_Emax(15.0), //
  _qe_filplot("psi"), //
  _qe_fileout("el1.cube"), //
  _qe_plot_num(7), //
  _qe_kpoint(1), //
  _qe_nfile(1), //
  _qe_iflag(3), //
  _qe_nx(20), //
  _qe_ny(20), //
  _qe_nz(20), //
  _qe_pp_outputformat(6), //
  _qe_weight(1.0) //
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
  _qe_mixing_beta = get_option("mixing_beta", _qe_mixing_beta);
  _qe_nbnd = get_option("nbnd", _qe_nbnd);
  _qe_degauss = get_option ("degauss", _qe_degauss); //
  _qe_DeltaE = get_option ("DeltaE", _qe_DeltaE); //
  _qe_Emin = get_option ("Emin", _qe_Emin); //
  _qe_Emax = get_option ("Emax", _qe_Emax); //
  _qe_filplot = get_option ("filplot", _qe_filplot); //
  _qe_fileout = get_option ("fileout", _qe_fileout); //
  _qe_plot_num = get_option ("plot_num", _qe_plot_num); //
  _qe_kpoint = get_option ("kpoint", _qe_kpoint); //
  _qe_nfile = get_option ("nfile", _qe_nfile); //
  _qe_iflag = get_option ("iflag", _qe_iflag); //
  _qe_nx = get_option ("nx", _qe_nx); //
  _qe_ny = get_option ("ny", _qe_ny); //
  _qe_nz = get_option ("nz", _qe_nz); //
  _qe_pp_outputformat = get_option ("pp_outputformat", _qe_pp_outputformat); //
  _qe_weight = get_option ("weight(1)", _qe_weight); //

  get_option("k_points_scf", _qe_k_points_scf);
  get_option("k_points_nscf", _qe_k_points_nscf);
  
  //if (_qe_k_points.size() < 2)
  //  throw InitFailedException("You need to provide the definitions for the kpoints.");

  vector<string> pseudos;
  get_option("pseudopotentials", pseudos);
  for (unsigned int i = 0; i < pseudos.size(); i = i + 2)
  {
    _qe_pseudos[pseudos[i]] = pseudos[i+1];
    //cerr << pseudos[i] << "  " << pseudos[i+1] << endl;
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
  ofstream qe_pdos(qe_file_base + ".pdos.in"); //
  ofstream qe_ppel(qe_file_base + ".pp.el1.in"); //

  TeeStream qe_both(qe_scf, qe_nscf);

  qe_both << "&control" << endl;
  qe_pdos << "&projwfc" << endl; //
  qe_ppel << "&INPUTPP" << endl; //
  if (get_option("relax", false))
    qe_scf  << "   calculation = 'relax'," << endl;
  else
    qe_scf  << "   calculation = 'scf'," << endl;

  qe_scf  << "   restart_mode='from_scratch'," << endl;
  qe_nscf << "   calculation = '"
          << get_option("nscf_calculation", "nscf")
          << "'," << endl;
  qe_both << "   prefix = '" << prefix << "'," << endl
          << "   pseudo_dir = '" << _qe_pseudo_dir << "'," << endl
          << "   outdir = '" << _qe_outdir << "'," << endl
          << "   verbosity = 'high'," << endl
          << "   tstress = .true.," << endl
          << "   tprnfor = .true.," << endl;
  if (get_option("wf_collect", false))
    qe_both << "   wf_collect = .true.," << endl;
  qe_both << "/" << endl;
  qe_pdos << "   prefix = '" << prefix << "'," << endl //
          << "   outdir = '" << _qe_outdir << "'," << endl //
          << "   degauss = " << _qe_degauss << "," << endl //
          << "   DeltaE = " << _qe_DeltaE << "," << endl //
          << "   Emin = " << _qe_Emin << "," << endl //
          << "   Emax = " << _qe_Emax << "," << endl //
          << "/" << endl;
  qe_ppel << "   prefix = '" << prefix << "'," << endl //
          << "   outdir = '" << _qe_outdir << "'," << endl //
          << "   filplot = '" << _qe_filplot << "'," << endl //
          << "   plot_num = " << _qe_plot_num << "," << endl //
          << "   kpoint = " << _qe_kpoint << "," << endl //
          << "   kband = " << _qe_nbnd << "," << endl //
          << "/" << endl //
          << "&PLOT" << endl //
          << "   nfile = " << _qe_nfile << "," << endl //
          << "   weight(1) = " << _qe_weight << "," << endl //
          << "   iflag = " << _qe_iflag << "," << endl //
          << "   nx = " << _qe_nx << "," << endl //
          << "   ny = " << _qe_ny << "," << endl //
          << "   nz = " << _qe_nz << "," << endl //
          << "   output_format = " << _qe_pp_outputformat << "," << endl //
          << "   fileout = '" << _qe_fileout << "'," << endl //
          << "/" << endl; //

  // for now this is fixed to tetragonal cells or hexagonal
  int ibrav = 8;

  libMesh::RealVectorValue a, b, c;
  this->get_atomistic_structure()->get_lattice_vectors(a, b, c);
  Messages m;
  m.info("Periodicity vectors:");
  m.indent();
  ostringstream os;
  os << a << endl
    << b << endl
    << c << endl;
  m.info(os.str());
  
  double ab = a*b;
  double ac = a*c;
  double bc = b*c;
  //cerr << ab << " " << ac << " " << bc << endl;

  int nat = this->get_atomistic_structure()->get_N_atoms();
  int ntyp = this->get_atomistic_structure()->get_N_types(); 
  qe_both << "&system" << endl;

  if ((fabs(ab) < 1e-12) && (fabs(ac) < 1e-12) && (fabs(bc) < 1e-12))
  {
    m.info("Using orthorhombic Bravais lattice (ibrav = 8)");
  }
  else
  {
    //m.info("Using hexagonal Bravais lattice (ibrav = 4)");
    //ibrav = 4;
    m.info("Using generic lattice definition (ibrav = 0)");
    ibrav = 0;
  }

  qe_both << "   ibrav = " << ibrav << "," << endl;

  if (ibrav == 8)
  {
    qe_both << "   a = " << a.norm() << "," << endl
            << "   b = " << b.norm() << "," << endl
            << "   c = " << c.norm() << "," << endl;
  }
  else if (ibrav == 4)
  { 
    double alat = 0, c_o_a = 0;
    // detect the c axis
    if (fabs(ab) < 1e-12)
    {
      // a or b
      if (fabs(ac) < 1e-12)
      {
        alat = b.norm();
        c_o_a = a.norm() / alat;
      }
      else
      {
        alat = a.norm();
        c_o_a = b.norm() / alat;
      }
    }
    else
    {
      alat = a.norm();
      c_o_a = c.norm() / alat;
    }

    alat /= Constants::bohr_radius / 1e-10;

    qe_both << "   celldm(1) = " << alat << "," << endl
            << "   celldm(3) = " << c_o_a << "," << endl;

    /*
    qe_both << "   a = " << a.norm() << "," << endl
            << "   b = " << b.norm() << "," << endl
            << "   c = " << c.norm() << "," << endl
    */

  }

  qe_both << "   nat = " << nat << ", ntyp = " << ntyp << "," << endl
          << "   ecutwfc = " << _qe_ecutwfc << "," << endl;
  if (_qe_nbnd > 0)
    qe_nscf << "   nbnd = " << _qe_nbnd << "," << endl;
  qe_both << "/" << endl;

  qe_both << "&electrons" << endl
          << "   conv_thr = " << _qe_conv_thr << "," << endl;
  if (_qe_mixing_beta > 0)
    qe_both << "   mixing_beta = " << _qe_mixing_beta << "," << endl;
  qe_both << "/" << endl;

  if (get_option("relax", false))
  {
    qe_scf << "&ions" << endl
      << "/" << endl;
  }

  if (ibrav == 0)
  {
    qe_both << "CELL_PARAMETERS angstrom" << endl;
            //<< a.print_unformatted() << endl
    qe_both << "  "
      << a(0) << " " << a(1) << " " << a(2) << endl;
    qe_both << "  "
      << b(0) << " " << b(1) << " " << b(2) << endl;
    qe_both << "  "
      << c(0) << " " << c(1) << " " << c(2) << endl;
  }



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
  for (unsigned int i = 0; i < _qe_k_points_scf.size(); ++i)
  {
    qe_scf << _qe_k_points_scf[i] << endl;
  }
  //qe_nscf << "crystal_b\n1\n0 0 0 1\n";
  for (unsigned int i = 0; i < _qe_k_points_nscf.size(); ++i)
  {
    qe_nscf << _qe_k_points_nscf[i] << endl;
  }
}

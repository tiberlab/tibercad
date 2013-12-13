#include "AtomisticBasis.h"
#include "Device.h"


//STD library includes
#include<iostream>
#include<fstream>
#include<sstream>
//-------------------

AtomisticBasis::~AtomisticBasis(void)
{
  if (_bondmap != NULL) delete _bondmap;
}

AtomisticBasis::AtomisticBasis(void)
:_bondmap(NULL),
_lattice_vectors(9,0.0)
{

}


void
AtomisticBasis::set_ttype_lattice_vectors(const Tensor2Gen& T)
{
  unsigned int count = 0;
  for (int i = 0; i < 3 ; i++)
    {
      for (int j = 0; j < 3 ; j++)
        {
          _lattice_vectors[count] = T(j+1,i+1);
          count++;
        }
    }
  _is_periodic = true;
}

Tensor2Gen 
AtomisticBasis::get_ttype_lattice_vectors(void)
{
  Tensor2Gen T;
  unsigned int count = 0;
  for (int i = 0; i < 3 ; i++)
  {
      for (int j = 0; j < 3 ; j++)
      {
          T(j+1,i+1) = _lattice_vectors[count];
          count++;
      }
  }
  return T;

}

int
AtomisticBasis::get_type_index(const std::string& type)
{
  int result = 0;
  for (int i = 0; i < N_types; i++){
      if ( (type.compare( _atom_types[i] ) == 0) ) result = i + 1;
  }

  return result;

}


void
AtomisticBasis::set_atom_types(const std::set<std::string>& atom_types)
{

  for (std::set<std::string>::iterator types = atom_types.begin(); types != atom_types.end(); types++)
    {
      _atom_types.push_back( *types );
    }

}



void
AtomisticBasis::build_bond_map(void)
{
  if (_bondmap != NULL)
  {
    delete _bondmap;
  }
  _bondmap = new BondMap(_atoms.size());

  Tensor2Gen period;
  for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < 3; j++)
        {
          period(j + 1, i + 1) = _lattice_vectors[i*3 + j];
        }
    }
  _bondmap->do_solve(_atoms, period);

}


void
AtomisticBasis::refresh(void)
{
//This should contain all important operation which build additional infos
//(bondmap, atom types, N atoms) starting from lattice vectors and
//atom vector only

  build_bond_map();
  N_atoms = _atoms.size();
  std::set<std::string> types;
  for (std::vector<Atom>::iterator it = _atoms.begin(); it != _atoms.end(); ++it)
  {
    types.insert(it->get_specie().get_string());
  } 
  set_atom_types(types);
}

//-------------------------------------------------------------------
//Print utilities
//-------------------------------------------------------------------

void
AtomisticBasis::print_xyz(const std::string& path) const
{
  std::ofstream file;
  // -------------------------------------------

  std::string outdir = TiberCad::get_output_dir();
  std::string file_name = outdir + "/" + path;

  // --------------------------------------------
  file.open(file_name.c_str());
  file << _atoms.size() << std::endl << std::endl;
  for (unsigned int i = 0; i < _atoms.size(); i++)
  {
    file << std::setw(2) << _atoms[i].get_specie()
    << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(0))
    << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(1))
    << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(2)) << "\n";
  }
  file.close();
}

void
AtomisticBasis::print_xyb(const std::string& path) const
{

  std::ofstream file;
  // -------------------------------------------

  std::string outdir = TiberCad::get_output_dir();
  std::string file_name = outdir + "/" + path;

  // --------------------------------------------
  file.open(file_name.c_str());
  file << _atoms.size() << std::endl << std::endl;

      for (unsigned int i = 0; i < _atoms.size(); i++)
        {
          file << std::setw(2) << _atoms[i].get_specie()
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(0))
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(1))
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(2));

          if (_bondmap != NULL)
            {

              file << std::setw(5) << (*_bondmap)[i].size();

              // N.B. Indexing is in Fortran notation (first atom is labelled as 1) !!!!!!!!!!!!!!!!
              for (unsigned int j = 0; j < (*_bondmap)[i].size(); j++)
                {
                  file << std::setw(10) << (*_bondmap)[i][j] + 1;
                }
              ///////////////////////////////////////////

            }

          file << std::endl;
        }
  file.close();
}


void
AtomisticBasis::print_gen(const std::string& path) const
{

  std::ofstream file;
  // -------------------------------------------

  std::string outdir = TiberCad::get_output_dir();
  std::string file_name = outdir + "/" + path;

  // --------------------------------------------
  file.open(file_name.c_str());
     file << _atoms.size();

      if (_is_periodic) file << std::setw(10) << "S \n";
      else file << std::setw(10) << "C \n";

      for (unsigned int i = 0; i < _atom_types.size(); i++)
        {
          file << std::setw(6) << _atom_types[i];
        }
      file << std::endl;

      for (unsigned int i = 0; i < _atoms.size(); i++)
        {
          unsigned int n_specie;
          for (n_specie = 0; n_specie < _atom_types.size(); n_specie++)
            {
              if (_atom_types[n_specie] == _atoms[i].get_specie() ) break;
            }
          file << std::setw(10) << i + 1 << std::setw(5) << n_specie + 1
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(0))
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(1))
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(2)) << "\n";
        }

      //A line of zeros is put here (coordinates origin)
      file <<  std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0)
        << std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0)
        << std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0) << "\n";

      // Periodicity vectors at the bottom
      unsigned int count = 0;
      for (unsigned int i = 0; i < 3; i++)
        {
          for (unsigned int j = 0; j < 3; j++)
            {
              file << std::setw(20) << std::setprecision(10) << std::fixed <<
                  _lattice_vectors[count];
              count++;
            }
          file << "\n";
        }
  file.close();

}

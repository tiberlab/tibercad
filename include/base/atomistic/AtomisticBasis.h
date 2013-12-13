// $Id: AtomisticStructure.h 3579 2013-04-22 12:25:47Z gpenazzi $

#ifndef _ATOMISTICBASIS_H_
#define _ATOMISTICBASIS_H_

#include "Atom.h"
#include "BondMap.h"
#include "Messages.h"
#include "Specie.h"
#include "InitFailedException.h"
#include "BondMap.h"

//STD library includes
#include <set>
//-------------------

class AtomisticBasis
{
  public:

  //! Virtual destructor, in case of cast of derived class 
  virtual ~AtomisticBasis();

  //! Get bond map
  const BondMap& get_bond_map(void) const;
    
  //! Get neighbor periodic image translation
  //! (same indexing as bondmap)
  const BondMap::Translation& get_neighbor_translation(void) const;

  //! Return a const reference to structure atoms
  const std::vector<Atom>& get_structure_atoms(void) const;
  
  //! Return a reference to structure atoms
  std::vector<Atom>& get_structure_atoms(void);

  //! Return a cons reference to an atom
  const Atom& get_structure_atom(unsigned int i) const;

  //! Assign structure atoms vector
  void set_structure_atoms(const std::vector<Atom>& atoms);

  //! Get periodicity vectors for the structure:
  //[x1, y1, z1, x2, y2, z2, x3, y3, z3]
  const std::vector<double>& get_lattice_vectors(void) const;

  //! Set the periodicity
  void set_ttype_lattice_vectors(const Tensor2Gen& T);
  
  //! Set the periodicity
  Tensor2Gen get_ttype_lattice_vectors(void);

  //! get the number of atoms in the structure
  int get_N_atoms() const {return N_atoms;}

  //! set number of atoms
  void set_N_atoms(long N) {N_atoms=N;}

  //! get the number of types in the structure
  int get_N_types() const {return N_types;}

  //! set number of atoms
  void set_N_types(int N) {N_types=N;}

  //! Print structure in xyz format
  void print_xyz(const std::string& path) const;

  //! Print structure in xyb format (xyz + bondmap)
  void print_xyb(const std::string& path) const;

  //! Print structure in gen format
  void print_gen(const std::string& path) const;

  //! Return a reference to atom types
  const std::vector<std::string>& get_atom_types (void);

  //! Get atom type index (types are stored in _atom_types)
  int get_type_index   (  const std::string &   type    );

  //! Empty list of atom types
  void clear_atom_types(){_atom_types.clear();};

  //! Set list of atom types
  void set_atom_types(const std::set<std::string>& types);

  //! Set periodicity information
  void set_periodic(bool periodic);

  //!Tells if the structure is periodical
  bool is_periodic() const;

  const Atom& operator[](unsigned int i) const;

  protected:
 
  //! Now we don't allow direct declaration of AtomisticBasis,
  //! maybe this will change in the future
  AtomisticBasis(void);

  //! Build bond map
  void build_bond_map(void);

  //! Refresh (or build) information as bondmap, number of atoms etc.
  //! This can be needed if we move atoms and the bondmap change, or
  //! we change total number of atoms
  void refresh(void);
    
  //! Bond Map object pointer
  BondMap* _bondmap;

  //! Tells if the structure is periodical
  bool _is_periodic; 

  //! Vector containing structure atoms 
  std::vector<Atom>  _atoms;
  
  //! Periodicity vectors in canonical basis
  //TODO: this should be changed back in Tensor2Gen 
  //or all the other Tensor2Gen should be changed.
  //There's no reason why we should use two different 
  //representations
  std::vector<double> _lattice_vectors;
  
  //! Number of atoms in structure
  unsigned int N_atoms;

  //! Number of species
  int N_types;

  //! Set of all atom types in structure
  std::vector<std::string> _atom_types;
  
  private:

};


inline
bool
AtomisticBasis::is_periodic() const 
{
  return _is_periodic;
}

inline
void
AtomisticBasis::set_periodic(bool periodic)
{
  _is_periodic = periodic;
}

inline
const BondMap&
AtomisticBasis::get_bond_map() const
{
  return *_bondmap;
}


inline
const BondMap::Translation&
AtomisticBasis::get_neighbor_translation(void) const
{
  return _bondmap->get_translation();
}


inline
const std::vector<Atom>& 
AtomisticBasis::get_structure_atoms(void) const
{
  return _atoms;
}

inline
std::vector<Atom>&
AtomisticBasis::get_structure_atoms(void)
{
  return _atoms;
}


inline
const Atom& 
AtomisticBasis::get_structure_atom(unsigned int i) const
{
  return _atoms[i];
}

inline
const Atom& 
AtomisticBasis::operator[](unsigned int i) const
{
  return _atoms[i];
}



inline
void 
AtomisticBasis::set_structure_atoms(const std::vector<Atom>& atoms)
{
  _atoms = atoms;
}


inline
const std::vector<double>& 
AtomisticBasis::get_lattice_vectors(void) const
{
  return _lattice_vectors;
}



inline
const std::vector<std::string>& AtomisticBasis::get_atom_types(void)
{
  return _atom_types;
}
#endif // _ATOMISTICBASIS_H_

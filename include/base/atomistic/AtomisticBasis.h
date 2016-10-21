// $Id: AtomisticStructure.h 3579 2013-04-22 12:25:47Z gpenazzi $

#ifndef _ATOMISTICBASIS_H_
#define _ATOMISTICBASIS_H_

#include "Atom.h"
#include "BondMap.h"
#include "Messages.h"
#include "Specie.h"
#include "InitFailedException.h"
#include "BondMap.h"

#include "HashSet.h"
#include "HashMap.h"

//STD library includes
#include <set>
//-------------------


namespace libMesh
{
  template <typename T> class VectorValue;
  typedef VectorValue<double> RealVectorValue;
}

class ModelOptions;



class AtomisticBasis
{
  public:

  //! An iterator to iterate over atoms in neighbor ordered way
  /*!
   * This iterator honours periodicity of the structure
   */
  class neighbor_iterator;

  //class const_neighbor_iterator;

  //! Virtual destructor, in case of cast of derived class 
  virtual ~AtomisticBasis();

  //! Get bond map
  const BondMap& get_bond_map(void) const;
    
  //! Get neighbor periodic image translation
  /*!
   * (same indexing as bondmap)
   */
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
  /*!
   * [x1, y1, z1, x2, y2, z2, x3, y3, z3]
   */
  const std::vector<double>& get_lattice_vectors(void) const;

  void get_lattice_vectors(libMesh::RealVectorValue& a,
                           libMesh::RealVectorValue& b,
                           libMesh::RealVectorValue& c) const;

  //! Set the periodicity
  void set_ttype_lattice_vectors(const Tensor2Gen& T);
  
  //! Set the periodicity
  Tensor2Gen get_ttype_lattice_vectors(void);

  //! Build and return a BondMap object
  /*!
   * \param periodicity indicates periodicity along the three
   *        lattice vectors (1 = periodic)
   */
  BondMap* build_bond_map(bool periodicity[3]) const;

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
  int get_type_index(  const std::string &   type ) const;

  //! Empty list of atom types
  void clear_atom_types(){_atom_types.clear();};
  
  //! Set list of atom types
  void set_atom_types(const std::set<std::string>& types);

  //! Set periodicity information
  void set_periodicity(std::vector<bool> periodicity);

  //! Set periodicity information
  void set_periodicity(bool, bool, bool);

  //! Ask for periodicity along coordinate axes
  /*!
   * 0 = x, 1 = y, 2 = z
   */
  bool is_periodic(unsigned int direction) const;

  //void set_periodic(bool periodic) {};

  //! Tells if the structure is periodic in any direction
  bool is_periodic(void) const;

  const Atom& operator[](unsigned int i) const;


  //! Get a subset of atoms
  /*!
   * \param subset vector containing indices of atoms in the subvolume
   * \param opt ModelOptions object describing the subvolume
   */
  void get_subset(std::vector<unsigned int>& subset,
      const ModelOptions& opt);


  //! Get the neighbour iterator for an atom
  /*!
   * \param index index of the atom
   * \param cutoff cutoff distance in A
   */
  neighbor_iterator neighbors_begin(unsigned int index,
      double length = 5.0, double height = 0.0, double width = 0.0) const;

  //! Get the past-the-end iterator for an atom
  /*!
   * \param index index of the atom
   * \param cutoff cutoff distance in A
   */
  neighbor_iterator neighbors_end(unsigned int index,
      double length = 5.0, double height = 0.0, double width = 0.0) const;


  protected:
 
  /*! Now we don't allow direct declaration of AtomisticBasis,
   * maybe this will change in the future
   */
  AtomisticBasis(void);

  //! Build bond map
  void build_bond_map(void);

  //! Refresh (or build) information as bondmap, number of atoms etc.
  /*! This can be needed if we move atoms and the bondmap change, or
   * we change total number of atoms
   */
  void refresh(void);
    
  //! Bond Map object pointer
  BondMap* _bondmap;

  //! Vector containing structure atoms 
  std::vector<Atom>  _atoms;
  
  //TODO: this should be changed back in Tensor2Gen 
  //or all the other Tensor2Gen should be changed.
  //There's no reason why we should use two different 
  //representations
  //! Periodicity vectors in canonical basis
  std::vector<double> _lattice_vectors;
  
  //! Number of atoms in structure
  unsigned int N_atoms;

  //! Number of species
  int N_types;

  //! Set of all atom types in structure
  std::vector<std::string> _atom_types;

  private:

  //! Periodicity for each direction
  std::vector<bool> _periodicity;


  public:

  class neighbor_iterator
  {
    public:
    neighbor_iterator(const AtomisticBasis& structure,
        unsigned int start,
        double length = 10, double height = 0, double width = 0,
        bool begin = true);

    neighbor_iterator(const neighbor_iterator& rhs);

    neighbor_iterator& operator++(void);


    bool operator==(const neighbor_iterator& rhs)
    {
      bool equal = &_structure == &rhs._structure;
      equal &= _start == rhs._start;
      equal &= _current == rhs._current;
      return(equal);
    }

    bool operator!=(const neighbor_iterator& rhs)
    {
      return(!(*this == rhs));
    }

    const Atom* operator*(void)
    {
      const Atom* at =
          (_current < _structure.get_N_atoms()) ?
              &_structure.get_structure_atom(_current) : NULL;

      return(at);
    }

    unsigned int atom_index(void)
    {
      return _current;
    }


    libMesh::Point atom_translation(void)
    {
      return _image;
    }


    private:

    typedef HashMultiMap<unsigned int, libMesh::Point>::Type HMMap;
    typedef HashMap<unsigned int, libMesh::Point>::Type HMap;

    neighbor_iterator& operator=(const neighbor_iterator& rhs);

    const AtomisticBasis& _structure;
    unsigned int _start;
    unsigned int _current;


    double _length;
    double _height;
    double _width;

    libMesh::Point _image;

    double _min_bond = 4;


    HMMap _visited;
    HMap _setA;
    HMap _setB;
    HMap::iterator _itA;

  };
};



inline
void
AtomisticBasis::set_periodicity(std::vector<bool> periodicity)
{
  _periodicity = periodicity;
}

inline
void
AtomisticBasis::set_periodicity(bool px, bool py, bool pz)
{
  _periodicity.resize(3);
  _periodicity[0] = px;
  _periodicity[1] = py;
  _periodicity[2] = pz;
}


inline
bool
AtomisticBasis::is_periodic(unsigned int direction) const
{
  assert(direction < 3);
  return(_periodicity[direction]);
}




inline
bool
AtomisticBasis::is_periodic(void) const
{
  return(_periodicity[0] || _periodicity[1] || _periodicity[2]);
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

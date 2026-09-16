/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file AtomisticBasis.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_ATOMISTICBASIS_H
#define TC_ATOMISTICBASIS_H

#include "tibercad/atomistic/Atom.h"
#include "tibercad/math/Tensor2.h"
#include "tibercad/base/HashMap.h"

//STD library includes
#include <set>
//-------------------


namespace libMesh
{
  template <typename T> class VectorValue;
  typedef VectorValue<double> RealVectorValue;
}

class ModelOptions;
class BondMap;



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

  //! Get the bounding box
  std::pair<libMesh::Point, libMesh::Point> get_bounding_box(void) const;

    
  //! Return a const reference to structure atoms
  const std::vector<Atom>& get_structure_atoms(void) const;
  
  //! Return a reference to structure atoms
  std::vector<Atom>& get_structure_atoms(void);

  //! Return a cons reference to an atom
  const Atom& get_structure_atom(unsigned int i) const;

  //! Assign structure atoms vector
  void set_structure_atoms(const std::vector<Atom>& atoms);

  //! Get periodicity vectors for the structure:
  const std::vector<libMesh::RealVectorValue>&
    get_lattice_vectors(void) const;

  //! Get the three lattice vectors
  void get_lattice_vectors(libMesh::RealVectorValue& a,
                           libMesh::RealVectorValue& b,
                           libMesh::RealVectorValue& c) const;

  //! Set the lattice vectors
  void set_lattice_vectors(const libMesh::RealVectorValue& a,
                           const libMesh::RealVectorValue& b,
                           const libMesh::RealVectorValue& c);

  //! Set the periodicity
  void set_ttype_lattice_vectors(const Tensor2& T);
  
  //! Set the periodicity
  Tensor2 get_ttype_lattice_vectors(void);

  //! Refresh (or build) information as bondmap, number of atoms etc.
  /*! This can be needed if we move atoms and the bondmap change, or
   * we change total number of atoms
   */
  void refresh(void);

  //! Build or re-build bond map
  void build_bond_map(void);

  //! Set the bond map
  void set_bond_map(const BondMap& bondmap);

  //! Build and return a BondMap object
  /*!
   * \param periodicity indicates periodicity along the three
   *        lattice vectors (1 = periodic)
   */
  BondMap* build_bond_map(bool periodicity[3]) const;

  //! get the number of atoms in the structure
  unsigned int get_N_atoms() const {return _atoms.size();}

  //! get the number of types in the structure
  int get_N_types() const {return _atom_types.size();}

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

  //! Get periodicity vectors for the structure:
  std::vector<libMesh::RealVectorValue>&
    get_lattice_vectors(void);

  //! Ask for periodicity along coordinate axes
  /*!
   * 0 = x, 1 = y, 2 = z
   */
  bool is_periodic(unsigned int direction) const;

  // ! Returns the vector of periodicity along all axes
  std::vector<bool> get_periodicity_vector(void) const; 

  //! Tells if the structure is periodic in any direction
  bool is_periodic(void) const;

  //! Set the origin
  void set_origin(const libMesh::Point& origin);

  //! Get the origin of the structure
  const libMesh::Point& get_origin(void) const;

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

  //! Copy constructor
  AtomisticBasis(const AtomisticBasis& other);

    
  //! Bond Map object pointer
  BondMap* _bondmap;

  //! Vector containing structure atoms 
  std::vector<Atom>  _atoms;
  
  //! Set of all atom types in structure
  std::vector<std::string> _atom_types;

  private:

  //! Periodicity for each direction
  std::vector<bool> _periodicity;

  //! Periodicity vectors in canonical basis
  std::vector<libMesh::RealVectorValue> _lattice_vectors;

  //! The origin for the structure
  libMesh::Point _origin;


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
void
AtomisticBasis::set_origin(const libMesh::Point& origin)
{
  _origin = origin;
}


inline
const libMesh::Point&
AtomisticBasis::get_origin(void) const
{
  return(_origin);
}



inline
bool
AtomisticBasis::is_periodic(unsigned int direction) const
{
  assert(direction < 3);
  return(_periodicity[direction]);
}


inline
std::vector<bool>
AtomisticBasis::get_periodicity_vector(void) const
{
  return(_periodicity);
}



inline
bool
AtomisticBasis::is_periodic(void) const
{
  return(_periodicity[0] || _periodicity[1] || _periodicity[2]);
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
const std::vector<libMesh::RealVectorValue>&
AtomisticBasis::get_lattice_vectors(void) const
{
  return _lattice_vectors;
}


inline
std::vector<libMesh::RealVectorValue>&
AtomisticBasis::get_lattice_vectors(void)
{
  return _lattice_vectors;
}




inline
const std::vector<std::string>& AtomisticBasis::get_atom_types(void)
{
  return _atom_types;
}

#endif // TC_ATOMISTICBASIS_H

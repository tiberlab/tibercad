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
 * \file BondMap.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_BONDMAP_H
#define TC_BONDMAP_H

//-------------------------------------------

#include "tibercad/atomistic/Atom.h"
#include "tibercad/atomistic/Specie.h"
#include "tibercad/math/Tensor2.h"

#include "libmesh/vector_value.h"

#include <vector>
#include <map>
#include <set>



//! A class for managing bond maps
/*!
 * When BondMap is created, it can calculate
 * the bond map of the system with a O(N) algorithm.
 * single atom cutoff distancies (harcoded) are used.
 * Two atoms are neighbours if condition
 * $|p_{1} - p_{2}| \leq  |c_{1} + c{2}|$
 * where p_{i} and c_{i} are position and cutoff parameter
 * of atom i, but only if 2 cannot be reached via another
 * more close first neighbor.
 * Periodical structure bond maps are also calculated.
 * Bond Map are stored in a vector.
 */
class BondMap : public std::vector<std::vector<unsigned int>> 
{

public:
  
  typedef std::vector<std::vector<libMesh::Point>> Translation;


  //! BondMap constructor with size as number of atoms
  BondMap(unsigned int structure_size, unsigned int valence = 4);

  //!BondMap destructor
  ~BondMap();


  //! Calculates bond map
  void solve(const std::vector<Atom>& basis,
      const Tensor2& period, const libMesh::Point& origin);

  //! Set the periodicity vectors
  /*!
   * This method is used e.g. after deforming the structure, assuming
   * that the bond map (i.e. connectivity remains unaltered).
   */
  void set_periodicity(const libMesh::RealVectorValue& a,
                       const libMesh::RealVectorValue& b,
                       const libMesh::RealVectorValue& c);

  //! Define edges of atomic basis
  static
  void define_edges(const std::vector<Atom>& basis, Tensor1& edge_min, Tensor1& edge_max);

  //! Gives translation vector for periodical images
  const Translation& get_translation(void) const;

  //! Gives translation vector for periodical images
  Translation& get_translation(void);

  //! Get the spatial translation vector for a neighbor
  /*
   * Translation vector \c t is defined such that the real position of
   * \c neighbor is \c r \c + \c t.
   */
  libMesh::Point get_translation(unsigned int atom, unsigned int neighbor) const;

  //! print
  void print(void) const;
  void print(const std::vector<Atom>& basis) const;
  
  //! Remove a set of atoms from the bond map
  void remove_atoms(const std::set<unsigned int> ids);


  //! Reserve size
  void reserve(unsigned int size);

 
private:


  //! Process two atoms
  void process_atoms(const std::vector<Atom>& basis, 
                     const unsigned int i,
                     const unsigned int j, 
                     const Tensor1& period);
  
  //! Fix Bondmap by eliminating 2nd NN
  void fix_bondmap(const std::vector<Atom>& basis);


  //! Build cutoff distancies map
  void set_cutoff();

  //! Map for cutoff parameters
  std::map<Specie, double> _cutoff;
  
  
  

  //! Clean informations no more useful after bond map calculation
  void clean();


  //! Translation vector for periodic images
  /*!
   * tells for each neighbour the translation
   * vector for which it's a neighbour
   */
  Translation _translation;

  //! Structure periodicity
  Tensor2 _period;

};


//----------------------------------------------------
// Inline member functions
//----------------------------------------------------


#endif // TC_BONDMAP_H

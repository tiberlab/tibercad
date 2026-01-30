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


#ifndef _BONDMAP_H_
#define _BONDMAP_H_


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
 * In a first pass, possible neighbors are found by comparing
 * inter-atomic distances with a cutoff proportional to the
 * sums of the covalent radii, i.e. two atoms are potentially
 * neighbours if 
 * \f$|r_{1} - r_{2}| \leq \alpha(c_{1} + c_{2})\f$
 * where \f$r_{i}\f$ and \f$c_{i}\f$ are position and covalent radius
 * of atom i, and \f$\alpha\f \sim 1.5$.
 * Then, all neighbors not contributing the Voronoi cell of each atom
 * are eliminated, leading to the geometric bond map. Last, a cone
 * criterion is used to filter out spurious bonds with small but existing
 * Voronoi faces, using a variable angle cutoff based on the coordination
 * number guessed from the Voronoi neighbors.
 */
class BondMap : public std::vector<std::vector<unsigned int>> 
{

public:
  
  //! A typedef for the map of translation vectors
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
  /*!
   * This method checks the distance of a pair of atoms,
   * considering also possible periodic copies.
   */
  void process_atoms(const std::vector<Atom>& basis, 
                     const unsigned int i,
                     const unsigned int j, 
                     const Tensor1& period);
  
  //! Fix Bondmap
  /*!
   * This method eliminates neighbors not contributing to
   * the Voronoi cell, or having too small cone angles (i.e.
   * small Voronoi faces). This should provide a good guess of
   * the actual chemical bonding in most cases.
   */
  void fix_bondmap(const std::vector<Atom>& basis);


  //! Build cutoff distancies map
  void set_cutoff();

  //! Map for cutoff parameters
  std::map<Specie, double> _cutoff;
  
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


#endif // _BONDMAP_H_

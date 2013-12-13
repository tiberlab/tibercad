// $Id$

#ifndef _BONDMAP_H_
#define _BONDMAP_H_

//-------------------------------------------

#include "BondMap.h"
#include "Atom.h"
#include "tensor.h"
#include "Specie.h"
#include "GridCells.h"
#include "tiber_dll.h"

#include <vector>
#include <map>


//! A class for managing bond maps
/*!
 *When BondMap is created, it can calculate
 *the bond map of the system with a O(N) algorithm.
 *single atom cutoff distancies (harcoded) are used.
 *Two atoms are neighbours if condition
 * $|p_{1} - p_{2}| \leq  |c_{1} + c{2}|$
 * where p_{i} and c_{i} are position and cutoff parameter
 * of atom i.
 * Periodical structure bond maps are also calculated.
 * Bond Map are stored in a vector. Last elements of each
 * row is the number of neighbours.
 */


class BondMap : public std::vector<std::vector<unsigned int>> 
{

public:
  
  typedef std::vector<std::vector<Tensor1>> Translation;


  //!BondMap constructor with size as number of atoms
  BondMap(unsigned int structure_size, unsigned int valence = 4);

  //!BondMap destructor
  ~BondMap();


  //!Calculates bond map
  void do_solve(const std::vector<Atom>& basis, const Tensor2Gen& period);

  //! Define edges of atomic basis
  static
  void define_edges(const std::vector<Atom>& basis, Tensor1& edge_min, Tensor1& edge_max);

  //! Gives translation vector for periodical images
  const Translation& get_translation(void) const;

  //!print
  void print(const std::vector<Atom>& basis);


private:


  //! Process two atoms
  void process_atoms(const std::vector<Atom>& basis, 
                     const unsigned int i,
                     const unsigned int j, 
                     const Tensor1& period) TBDLLOCAL;

  //! Build cutoff distancies map
  void set_cutoff() TBDLLOCAL;

  //! Map for cutoff parameters
  std::map<Specie, double> _cutoff;

  //! Clean informations no more useful after bond map calculation
  void clean() TBDLLOCAL;


  //! Translation vector for periodic images (tells for each neighbour the translation
  //! vector for which it's a neighbour)
  Translation _translation;
  //--------------------------------------------------------------------

  //! Structure periodicity
  Tensor2Gen _period;

};


//----------------------------------------------------
// Inline member functions
//----------------------------------------------------

inline
const BondMap::Translation&
BondMap::get_translation(void) const
{
  return _translation;
}



#endif // _BONDMAP_H_

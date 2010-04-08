#ifndef _BONDMAP_H_
#define _BONDMAP_H_

//-------------------------------------------

#include "BondMap.h"
#include "Atom.h"
#include "tensor.h"

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

typedef
std::vector < std::vector < unsigned int > > Bondmap;

class BondMap
{

public:

  //!BondMap constructor
  BondMap();

  //!BondMap destructor
  ~BondMap();

  //!Initialize internal grid and allocate arrays
  void do_init(const unsigned int structure_size);

  //!Calculates bond map
  void do_solve(const std::vector<Atom>& basis, const Tensor2Gen& period);

  //! Define edges of atomic basis
  static
  void define_edges(const std::vector<Atom>& basis, Tensor1& edge_min, Tensor1& edge_max);

  //!Gives pointer to bond map
  Bondmap& get_bond_map();

  //! Gives translation vector for periodical images
  std::vector < std::vector < Tensor1 > >& get_translation(void);

private:


  //! Internally defined parallepipedal grid
  std::vector<std::vector<std::vector<std::vector<unsigned int> > > >  _grid_cell;

  //! Set private members related to grid definition
  void define_grid(const double minimum_spacing, const Tensor1& edge_min, const Tensor1& edge_max);

  //! Include atom indexes in proper cells
  void include_atoms(const std::vector<Atom>& basis);

  //! Process atoms in defined cells. Periodicity informations are used internally
  //! for border cells
  void process_cell(const std::vector<Atom>& basis, const unsigned int x1,
      const unsigned int y1, const unsigned int z1,
      int x2, int y2, int z2);


  //!Process a cell with all surrounding cells (only for existing ones, doen't take in account periodicity)
  void process_surrounding_cell(const std::vector<Atom>& basis, const unsigned int x, const unsigned int y,
      const unsigned int z);


  //! Process a single atom with all atoms of a cell
  void process_atom_with_cell(const std::vector<Atom>& basis, const unsigned int i,
      const unsigned int x, const unsigned int y,
      const unsigned int z, const Tensor1& period);

  ////! Same as above, but without periodicity (period = 0)
  //void process_atom_with_cell(const std::vector<Atom>& basis, const unsigned int i,
  //		const unsigned int x, const unsigned int y, const unsigned int z);

  //! Process two atoms
  void process_atoms(const std::vector<Atom>& basis, const unsigned int i,
      const unsigned int j, const Tensor1& period);

  //! Build cutoff distancies map
  void set_cutoff();

  //! Map for cutoff parameters
  std::map<std::string, double> _cutoff;

  //! Clean informations no more useful after bond map calculation
  void clean();

  //------------------------------------------------------------------
  //TODO: translation vector and double pointer must be substituted
  //with a vector of structures, or something similar
  //! Bond map
  Bondmap _bond_map;

  //! Translation vector for periodic images (tells for each neighbour the translation
  //! vector for which it's a neighbour)
  std::vector<std::vector<Tensor1> > _translation;
  //--------------------------------------------------------------------

  //! Spacing of the grid along x axis
  double _x_spacing;

  //! Spacing of the grid along y axis
  double _y_spacing;

  //! Spacing of the grid along z axis
  double _z_spacing;

  //! Number of grid section in x direction
  unsigned int _n_x;

  //! Number of grid section in y direction
  unsigned int _n_y;

  //! Number of grid section in z direction
  unsigned int _n_z;

  //! Local axis origin
  Tensor1 _origin;

  //! Structure periodicity
  Tensor2Gen _period;

};


//----------------------------------------------------
// Inline member functions
//----------------------------------------------------

inline
Bondmap&
BondMap::get_bond_map(void)
{
  return _bond_map;
}

inline
std::vector < std::vector < Tensor1 > >&
BondMap::get_translation(void)
{
  return _translation;
}




#endif // _BONDMAP_H_

// $Id: GridCells.h 2059 2010-08-31 13:11:02Z pecchia $
#ifndef _GRIDCELLS_H_
#define _GRIDCELLS_H_


#include "Atom.h"
#include "tensor.h"
#include "tiber_dll.h"

#include "libmesh/tensor_value.h"

#include <vector>
#include <map>

//! Class used to define a regular structured grid 

class GridCells
{

  public:
  
  class NeighborIterator;

  GridCells(const std::vector<Atom>& basis, const Tensor2Gen& period,
      const libMesh::Point& origin,
      const double minimum_spacing, unsigned int projected_dim = 3);

  ~GridCells(void);

  //! Number of grid section in x direction
  unsigned int n_x;

  //! Number of grid section in y direction
  unsigned int n_y;

  //! Number of grid section in z direction
  unsigned int n_z;
  
  //! Supercell vectors
  Tensor2Gen _period;
 

  void init(void);  
 
  //! Redefines mygrid(i,j,k)
  const std::vector<unsigned int>& operator()(const unsigned int i,  
                                              const unsigned int j,  
                                              const unsigned int k) const; 

  const std::vector<unsigned int>& operator[](const unsigned int i) const;  

  unsigned int size(void) const;

  //! Returns the cell were a point falls  
  unsigned int index(const unsigned int i,const unsigned int j,const unsigned int k) const;

  //! Returns index conversions 
  void index(const unsigned int c, unsigned int& x, unsigned int& y, unsigned int& z) const;
  
  //! Returns the cell were a point falls  
  void get_cell(const libMesh::Point& p, unsigned int& i,  
                                unsigned int& j, 
                                unsigned int& k) const;


  void print_statistics(void);
  

  NeighborIterator begin(unsigned int c);
  NeighborIterator begin(unsigned int x, unsigned int y, unsigned int z);
  NeighborIterator end(unsigned int c);
  NeighborIterator end(unsigned int x, unsigned int y, unsigned int z);

  private:

  std::vector<std::vector<unsigned int> > _grid_cell; 

  //! Define edges of atomic basis
  void define_edges(const std::vector<Atom>& basis) TBDLLOCAL;

  
  //! Set private members related to grid definition
  void define_grid(void) TBDLLOCAL;


  //! Include atom indexes in proper cells
  void include_atoms(const std::vector<Atom>& basis) TBDLLOCAL;

  //! The transformation matrix to get point coordinates in given basis
  libMesh::RealTensor _inv_transform;

  Tensor1 _edge_min;
  Tensor1 _edge_max;


  //! minimal spacing along x used for the grid
  double _minimum_spacing_x;

  //! minimal spacing along y used for the grid
  double _minimum_spacing_y;

  //! minimal spacing along z used for the grid
  double _minimum_spacing_z;

  //! Spacing of the grid along x axis
  double _x_spacing;

  //! Spacing of the grid along y axis
  double _y_spacing;

  //! Spacing of the grid along z axis
  double _z_spacing;

public:

  class NeighborIterator
  {
     public:

       //! Default constructor using 3 indeces
       NeighborIterator(const GridCells* container, unsigned int i, unsigned int j, unsigned int k,
                       int dx=-1, int dy=-1, int dz=-1) : 
        _container(container), _x(i), _y(j), _z(k), _dx(dx), _dy(dy), _dz(dz), _shift(0.0){};

       //! constructor using 1 index
       NeighborIterator(const GridCells* container, unsigned int c,
                        int dx, int dy, int dz) : 
         _container(container), _dx(dx), _dy(dy), _dz(dz), _shift(0.0)
       {
          _z = (unsigned int) c / (_container->n_x * _container->n_y);
          unsigned int r = c % (_container->n_x * _container->n_y);
          _y = (unsigned int) r / _container->n_x;
          _x = r % _container->n_x;
       }

       
       //! Copy constructor
       NeighborIterator(const NeighborIterator&  it) :
        _container(it._container), _x(it._x), _y(it._y), _z(it._z),
        _shift(it._shift), _dx(it._dx), _dy(it._dy), _dz(it._dz), _pair(it._pair)
       {}; 
       

       ~NeighborIterator(){};

       //! Get the shift in real coordinates
       Tensor1 get_shift(void);

       //! Iteration operator 
       NeighborIterator& operator++(void);
       
       //! = operator 
       NeighborIterator& operator=(const NeighborIterator& it)
       {
         _container = it._container;
         _x = it._x; _y = it._y, _z=it._z;
         _dx = it._dx; _dy = it._dy, _dz=it._dz;
         _shift = it._shift;
         _pair = it._pair;
         return *this;
       } 

       //! == operator 
       bool operator ==(const NeighborIterator& it)
       {
         return ( _x == it._x && _y == it._y && _z == it._z && 
              _dx == it._dx && _dy == it._dy && _dz == it._dz );
       }
       
       bool operator !=(const NeighborIterator& it)
       {
         return(!(*this == it));
       }

       //! dereferencing operator, return cell and shift vector 
       const std::pair<unsigned int , const Tensor1* >& operator *(void)
       {
          return _pair;   
       }

     private: 

       const GridCells* _container;
       unsigned int _x, _y, _z; 
       int _dx, _dy, _dz; 
       Tensor1 _shift; //(0.0);
       std::pair<unsigned int, const Tensor1* > _pair; //(0,_shift);  

  };


};

inline
const std::vector<unsigned int>&
GridCells::operator()(const unsigned int i, const unsigned int j, const unsigned int k) const
{
  return _grid_cell[k * n_x*n_y + j * n_x + i];
}
  
inline
const std::vector<unsigned int>& 
GridCells::operator[](const unsigned int i) const
{
   return _grid_cell[i];
}

inline 
unsigned int 
GridCells::size(void) const
{
  return n_x*n_y*n_z;
}

inline
unsigned int
GridCells::index(const unsigned int i, const unsigned int j, const unsigned int k) const
{
  return k*n_x*n_y + j*n_x + i;
}

inline
void
GridCells::index(const unsigned int c, unsigned int& x, unsigned int& y, unsigned int& z) const
{
  z = (unsigned int) c / (n_x * n_y);
  unsigned int r = c % (n_x * n_y);
  y = (unsigned int) r / n_x;
  x = r % n_x;
}


inline  
GridCells::NeighborIterator 
GridCells::begin(unsigned int c)
{
  GridCells::NeighborIterator b(this,c,-1,-1,-2);
  // dirty trick in order to initialize _shift properly using the code implemented in ++
  ++b;
  return b;
}

inline  
GridCells::NeighborIterator 
GridCells::begin(unsigned int x, unsigned int y, unsigned int z)
{
  return begin(index(x,y,z));
}

inline
GridCells::NeighborIterator
GridCells::end(unsigned int c)
{
  GridCells::NeighborIterator e(this,c,+1,+1,+2);
  // dirty trick in order to set the end just after the last 
  return e;
}

inline  
GridCells::NeighborIterator 
GridCells::end(unsigned int x, unsigned int y, unsigned int z)
{
  return end(index(x,y,z));
}


#endif // _GRIDCELLS_H_

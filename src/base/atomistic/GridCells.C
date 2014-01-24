
#include "GridCells.h"
#include "vector_value.h"




GridCells::GridCells(const std::vector<Atom>& basis, const Tensor2Gen& period,
    const double minimum_spacing, unsigned int projected_dim)
{
  _period = period;

  _minimum_spacing_x = minimum_spacing;
  _minimum_spacing_y = minimum_spacing;
  _minimum_spacing_z = minimum_spacing;
  
  if (projected_dim < 3)
  {
    RealVectorValue period_x(0.0), period_y(0.0), period_z(0.0);
    period_x(0) = _period(1,1);
    period_x(1) = _period(2,1);
    period_x(2) = _period(3,1);
    period_y(0) = _period(1,2);
    period_y(1) = _period(2,2);
    period_y(2) = _period(3,2);
    period_z(0) = _period(1,3);
    period_z(1) = _period(2,3);
    period_z(2) = _period(3,3);
    switch (projected_dim)
    {
      case 0:
        _minimum_spacing_x = period_x.size();

      case 1:
        _minimum_spacing_y = period_y.size();

      case 2:
        _minimum_spacing_z = period_z.size();
    }
  }

  define_edges(basis);

  define_grid();

  //resize grid_cell according to grid dimensions
  _grid_cell.resize(n_x * n_y * n_z);

  //a guess of atom density for cell is used for a first guess memory reservation
  //it could waste some memory but make calculation faster
  double density = basis.size() / (n_x * n_y * n_z);
  int int_density = static_cast<int>(floor(density));

  for (unsigned int i = 0; i < n_x*n_y*n_z; i++)
  {  
    _grid_cell[i].reserve(2 * int_density);
  }
  
  include_atoms(basis);

  //print_statistics();
  //std::cout<<"Number of atoms: "<<basis.size()<<std::endl;

};


GridCells::~GridCells(void)
{
    
};


void
GridCells::define_edges(const std::vector<Atom>& basis)
{
  _edge_min(1) = 1e9; _edge_min(2) = 1e9; _edge_min(3) = 1e9;
  _edge_max(1) =-1e9; _edge_max(2) =-1e9; _edge_max(3) =-1e9;

  for (unsigned int i = 0; i < basis.size(); i++)
    {
      if (basis[i].get_position(0) < _edge_min(1)) _edge_min(1) = basis[i].get_position(0);
      if (basis[i].get_position(0) > _edge_max(1)) _edge_max(1) = basis[i].get_position(0);

      if (basis[i].get_position(1) < _edge_min(2)) _edge_min(2) = basis[i].get_position(1);
      if (basis[i].get_position(1) > _edge_max(2)) _edge_max(2) = basis[i].get_position(1);

      if (basis[i].get_position(2) < _edge_min(3)) _edge_min(3) = basis[i].get_position(2);
      if (basis[i].get_position(2) > _edge_max(3)) _edge_max(3) = basis[i].get_position(2);
    }

  // min max edges must be slightly enlarged to avoid strict inclusion errors
  _edge_min(1) = _edge_min(1) - 0.1;
  _edge_min(2) = _edge_min(2) - 0.1;
  _edge_min(3) = _edge_min(3) - 0.1;
  _edge_max(1) = _edge_max(1) + 0.1;
  _edge_max(2) = _edge_max(2) + 0.1;
  _edge_max(3) = _edge_max(3) + 0.1;

}





void
GridCells::define_grid(void)
{
  double lenght_x, lenght_y, lenght_z;

  lenght_x = _edge_max(1) - _edge_min(1);
  lenght_y = _edge_max(2) - _edge_min(2);
  lenght_z = _edge_max(3) - _edge_min(3);

  n_x = static_cast<unsigned int>((floor(lenght_x / _minimum_spacing_x)) + 1);
  _x_spacing = lenght_x / n_x;

  n_y = static_cast<unsigned int>((floor(lenght_y / _minimum_spacing_y)) + 1);
  _y_spacing = lenght_y / n_y;

  n_z = static_cast<unsigned int>((floor(lenght_z / _minimum_spacing_z)) + 1);
  _z_spacing = lenght_z / n_z;

}




void
GridCells::include_atoms(const std::vector<Atom>& basis)
{
  unsigned int x, y, z;

  for (unsigned int i = 0; i < basis.size(); i++)
  {
    get_cell(basis[i].get_position(),x,y,z);
    _grid_cell[z * n_x*n_y + y * n_x + x].push_back(i);
  }

}

void 
GridCells::get_cell(const Point& p, unsigned int& x,  unsigned int& y, unsigned int& z) const
{
  
  double dx = p(0) - _edge_min(1);
  double dy = p(1) - _edge_min(2); 
  double dz = p(2) - _edge_min(3);
 
  if (dx < 0){ x = 0;}
  else
  {
    x = static_cast<unsigned int>(( floor( dx / _x_spacing ) ));
    if (x>n_x-1){x = n_x-1;}
  }

  if (dy < 0){ y = 0;}
  else
  {
    y = static_cast<unsigned int>(( floor( dy / _y_spacing ) ));
    if (y>n_y-1){y = n_y-1;}
  }

  if (dz < 0){ z = 0;}
  else
  {
    z = static_cast<unsigned int>(( floor( dz / _z_spacing ) ));
    if (z>n_z-1){z = n_z-1;}
  }

}

void
GridCells::print_statistics(void)
{
  std::cout<<"edge min: "<<_edge_min(1)<<" "<< _edge_min(2)<<" "<<_edge_min(3)<<std::endl;
  std::cout<<"edge max: "<<_edge_max(1)<<" "<< _edge_max(2)<<" "<<_edge_max(3)<<std::endl;
  std::cout<<"grid size ("<<size()<<") : "<<n_x<<" "<<n_y<<" "<<n_z<<std::endl;
  std::cout<<"Sides: "<<n_x*_x_spacing<<" "<<n_y*_y_spacing<<" "<<n_z*_z_spacing<<std::endl;
  std::cout<<"cell sizes: "<<_x_spacing<<" "<<_y_spacing<<" "<<_z_spacing<<std::endl;

  unsigned int minatoms=1e9, maxatoms=0, totatoms=0;
  for(unsigned int c=0; c < size(); c++)
  {
     unsigned int natoms= _grid_cell[c].size();
     if (natoms<minatoms){ minatoms = natoms;} 
     if (natoms>maxatoms){ maxatoms = natoms;} 
     totatoms += natoms;
  }

  std::cout<<"Number of atoms: "<<totatoms<<std::endl;
  std::cout<<"Average per cell: "<<totatoms/size()<<std::endl;
  std::cout<<"Minimum per cell: "<<minatoms<<std::endl;
  std::cout<<"Maximum per cell: "<<maxatoms<<std::endl;
  
}

GridCells::NeighborIterator&
GridCells::NeighborIterator::operator ++(void)
{
  //-1 +1 +1
  if (_dz < 1 ){ _dz += 1;}
  else if (_dy < 1){_dy += 1; _dz=-1;}
  else if (_dx < 1){_dx += 1; _dy=-1; _dz=-1;}
  else {_dz+=1; return (*this);}
 
  //std::cout<<"dx="<<_dx<<" dy="<<_dy<<" dz="<<_dz<<std::endl; 
  
  int tmp_x = static_cast<int>(_x);
  int tmp_y = static_cast<int>(_y);
  int tmp_z = static_cast<int>(_z);
  
  int x1 = tmp_x + _dx;
  int y1 = tmp_y + _dy;
  int z1 = tmp_z + _dz;

  unsigned int x2 = static_cast<unsigned int>(x1);
  unsigned int y2 = static_cast<unsigned int>(y1);
  unsigned int z2 = static_cast<unsigned int>(z1);


  Tensor1 period_x(0.0), period_y(0.0), period_z(0.0);
  // Initialize _shift to 0.
  _shift=period_x;

  period_x(1) = _container->_period(1,1); 
  period_x(2) = _container->_period(2,1); 
  period_x(3) = _container->_period(3,1);
  period_y(1) = _container->_period(1,2); 
  period_y(2) = _container->_period(2,2); 
  period_y(3) = _container->_period(3,2);
  period_z(1) = _container->_period(1,3); 
  period_z(2) = _container->_period(2,3); 
  period_z(3) = _container->_period(3,3);
 
  if (x1 == -1)
  {
    x2 = (_container->n_x) - 1; _shift = _shift - period_x;
  }
  else if (x1 == static_cast<int>(_container->n_x))
  {
    x2 = 0; _shift = _shift + period_x;
  }
  
  if (y1 == -1)
  {
    y2 = (_container->n_y) - 1; _shift = _shift - period_y;
  }
  else if (y1 == static_cast<int>(_container->n_y))
  {
    y2 = 0; _shift = _shift + period_y;
  }

  if (z1 == -1)
  {
    z2 = (_container->n_z) - 1; _shift = _shift - period_z;
  }
  else if (z1 == static_cast<int>(_container->n_z))
  {
    z2 = 0; _shift = _shift + period_z;
  }
  //------------------------------------------

  //std::cout<<"shift="<<_shift(1)<<" "<<_shift(2)<<" "<<_shift(3)<<std::endl;

  _pair = std::make_pair(_container->index(x2,y2,z2), &_shift);

  return (*this);

}

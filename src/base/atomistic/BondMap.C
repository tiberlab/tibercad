// $Id$

#include "BondMap.h"
#include "GridCells.h"
#include "Messages.h"
#include "Database.h"
#include "Messages.h"
#include "Utils.h"

//#include "libmesh/type_vector.h"
#include <libmesh/plane.h>

#include <fstream>


BondMap::BondMap(unsigned int structure_size, unsigned int valence)
:_period(0)
{

  (*this).resize(structure_size);
  for (unsigned int i = 0; i < structure_size; i++)
  {
    (*this)[i].reserve(valence);
  }

  _translation.resize(structure_size);
  for (unsigned int i = 0; i < structure_size; i++)
  {
    _translation[i].reserve(valence);
  }

}

void
BondMap::reserve(unsigned int size)
{
  vector<vector<unsigned int>>::reserve(size);
  _translation.reserve(size);
}


BondMap::~BondMap(void)
{
}


void
BondMap::clean()
{
  //_grid_cell.clear();
}

void
BondMap::set_cutoff()
{
  std::ifstream file;
  std::string database_path = Database::get_default_search_path();
  std::string filename = "cutoff.dat";
  std::string line, record;
  std::stringstream line_stream;
  Specie s;

  filename = database_path + "/" + filename;

    
  file.open(filename.c_str());
  if (!file)
  {
    Messages::error("Exception opening/reading file cutoff.dat");
    exit(1);
  }
  
  while (!file.eof())
  {
    getline(file, line);
    line_stream.str(std::string());
    line_stream.clear(std::stringstream::goodbit);
    line_stream << line;
    line_stream >> record;
    s = record;
    line_stream >> record;
    
    _cutoff[s] = atof(record.c_str());
    
  }
  
  file.close();
}


const BondMap::Translation&
BondMap::get_translation(void) const
{
  return _translation;
}


BondMap::Translation&
BondMap::get_translation(void)
{
  return _translation;
}




void
BondMap::set_periodicity(const libMesh::RealVectorValue& a,
                         const libMesh::RealVectorValue& b,
                         const libMesh::RealVectorValue& c)
{
  for (int j = 0; j < 3 ; j++)
  {
    _period(j+1, 1) = a(j);
    _period(j+1, 2) = b(j);
    _period(j+1, 3) = c(j);
  }
}




libMesh::Point
BondMap::get_translation(unsigned int atom, unsigned int neighbor) const
{
  const libMesh::Point& c = _translation[atom][neighbor];

  libMesh::Point shift;
  shift(0) = c(0) * _period(1, 1) + c(1) * _period(1, 2) + c(2) * _period(1, 3);
  shift(1) = c(0) * _period(2, 1) + c(1) * _period(2, 2) + c(2) * _period(2, 3);
  shift(2) = c(0) * _period(3, 1) + c(1) * _period(3, 2) + c(2) * _period(3, 3);

  return(shift);
}





void
BondMap::do_solve(const std::vector<Atom>& basis,
    const Tensor2Gen& period, const libMesh::Point& origin)
{

  Messages::debug("BondMap::do_solve");

  Tensor1 edge_min, edge_max;

  //set cutoff distancies
  set_cutoff();
  
  for (unsigned int i=0; i< basis.size(); i++)
  {
    if (_cutoff[basis[i].get_specie()] == 0.0)
    {
      Messages::error("Cutoff distance for specie "
                      + basis[i].get_specie().get_string()
                      +" is not defined. Add to materials/cutoff.dat");
    }
  }
 
  _period = period;
  
  //define the minimum spacing of the grid. the smaller it is, the faster is bonds calculations
  //cannot be smaller than the largest bond length (in Angstrom)
  GridCells cells(basis, period, origin, 8.0);

  // Loop on all cells
  Utils::Progress prog("BondMap", cells.size());

  for (unsigned int c1 = 0; c1 < cells.size(); c1++)
  {
    GridCells::NeighborIterator it = cells.begin(c1);
    GridCells::NeighborIterator end = cells.end(c1);
 
    prog.progress_message(c1+1);

    // Loop on all 27 neighboring cells (periodicity is taken care by the iterator)
    for ( ; it != end; ++it)
    {
      unsigned int c2 = (*it).first;
      const Tensor1& shift = *((*it).second);

      // Loop over all atoms in each cell
      for (unsigned int i = 0; i < cells[c1].size(); i++)
      {     
        for (unsigned int j = 0; j < cells[c2].size(); j++)
        {
          process_atoms(basis, cells[c1][i], cells[c2][j], shift);
        }  
      }
    }
  }

  fix_bondmap(basis);

}


void
BondMap::process_atoms(const std::vector<Atom>& basis, 
                       const unsigned int i,
                       const unsigned int j, 
                       const Tensor1& period)
{

  unsigned int put_here = 0;
  bool not_already_signed;
  double cutofftmp;
  const libMesh::Point per(period(1), period(2), period(3));
  libMesh::Point shift;
  shift(0) = per(0) * _period(1, 1) + per(1) * _period(1, 2) + per(2) * _period(1, 3);
  shift(1) = per(0) * _period(2, 1) + per(1) * _period(2, 2) + per(2) * _period(2, 3);
  shift(2) = per(0) * _period(3, 1) + per(1) * _period(3, 2) + per(2) * _period(3, 3);
  

  if ((i != j) || (per.norm() > 1e-5))
  {
   
    cutofftmp = _cutoff[basis[i].get_specie()] + _cutoff[basis[j].get_specie()];
    const libMesh::Point& position1 = basis[i].get_position();
    libMesh::Point position2 = basis[j].get_position() + shift;
    libMesh::Point vector = position2 - position1;
    

    if (vector.norm() < cutofftmp)
    {

      bool add_bond = true;
      put_here = (*this)[i].size();
      
      // check if it is already there
      for (unsigned int n = 0; n < (*this)[i].size(); n++){

        if ((*this)[i][n] == j && libMesh::Point(get_translation(i, n) - shift).norm() < 1e-6)
        {
          add_bond = false;
          break;
        }

      }

      if (add_bond)
      {
        (*this)[i].push_back(j);
        _translation[i].push_back(per);
      }
      
    }

  }

}


void
BondMap::fix_bondmap(const std::vector<Atom>& basis)
{
  // we decide whether a bond needs to be kept by checking
  // if the midpoint lies inside the Voronoi box around the
  // center atom. But this seems to work only if we check for
  // something longer than midpoint. Otherwise, some 2nd NN
  // bondings appear as 1st NN.
  double scale = 1.0;

  // loop over all atoms
  for (size_t i = 0; i < this->size(); ++i)
  {

    // first we reorder the bonds in increasing length
    // (using insertion sort algorithm)
    for (int j = 1; j < (*this)[i].size(); ++j)
    {
      unsigned int tmp = (*this)[i][j];
      libMesh::Point tmp_p = _translation[i][j];

      double dist = libMesh::Point(basis[(*this)[i][j]].get_position() +
          get_translation(i, j) - basis[i].get_position()).norm();

      int k;

      for (k = j; k > 0; --k)
      {
        double dist2 = libMesh::Point(basis[(*this)[i][k-1]].get_position() +
          get_translation(i, k-1) - basis[i].get_position()).norm();

        if (dist2 <= dist)
          break;

        (*this)[i][k] = (*this)[i][k-1];
        _translation[i][k] = _translation[i][k-1];
      }
      (*this)[i][k] = tmp;
      _translation[i][k] = tmp_p;

    }
    /*
    std::cerr << i << "  : ";
    for(int j = 0; j < (*this)[i].size(); ++j)
    {
        libMesh::Point normal = basis[(*this)[i][j]].get_position() +
            _translation[i][j] - basis[i].get_position();
      std::cerr <<  normal.norm() << " ";
    }
    std::cerr << "\n";
    */


    // use greater so that erasing happens starting from highest index
    std::set<size_t, std::greater<size_t>> flag;

    // loop over all neighbors
    for (int j = 0; j < (*this)[i].size(); ++j)
    {
      if (flag.count(j) == 0)
      {
        libMesh::Point normal = basis[(*this)[i][j]].get_position() +
            get_translation(i, j) - basis[i].get_position();

        //std::cerr << "norm = " << normal.norm() << "\n";
        libMesh::Plane p(basis[i].get_position() + 0.5*normal, normal);

        for(int k = 0; k < (*this)[i].size(); ++k)
        {
          if ((k != j) && (flag.count(k) == 0))
          {
            libMesh::Point vector = basis[(*this)[i][k]].get_position() +
                get_translation(i, k) - basis[i].get_position();

            if (vector.norm() > 1.1 * normal.norm())
            {
              libMesh::Point testpoint = basis[i].get_position() + scale * vector;

              if (p.above_surface(testpoint))
                flag.insert(k);
            }
          }
        }
      }
    }


    /*
    std::cerr << "\n" << i << ": ";
    for (auto&& n : (*this)[i])
      std::cerr << n << " ";
    std::cerr << std::endl << "delete: ";
    for(std::set<size_t>::iterator it = flag.begin(); it != flag.end(); it++)
        std::cerr << *it << " ";
    std::cerr << std::endl;
    */

    if (!flag.empty())
    {
      for (auto it = flag.begin(); it != flag.end(); ++it)
      {
        unsigned int id = *it;
        unsigned int neighbour = (*this)[i][id];

        (*this)[i].erase((*this)[i].begin() + id);
        _translation[i].erase(_translation[i].begin() + id);
/*
        if (neighbour != i)
        {
          for (int j = (*this)[neighbour].size() - 1; j >= 0; --j)
          {
            if ((*this)[neighbour][j] == i)
            {

              (*this)[neighbour].erase((*this)[neighbour].begin() + j);
              _translation[neighbour].erase(_translation[neighbour].begin() + j);
            }
          }
        }*/
      }
    }
  }

}   

void
BondMap::remove_atoms(const std::set<unsigned int> ids)
{
  if (ids.empty())
    return;

  vector<vector<unsigned int>> newbm;
  vector<vector<libMesh::Point>> newtr;

  unsigned int n = ids.size();

  newbm.resize(this->size() - n);
  newtr.resize(newbm.size());

  vector<unsigned int> idv;
  idv.assign(ids.begin(), ids.end());

  // a lambda function for adjusting atom id
  auto newid = [&](unsigned int id) {
    unsigned int _newid = id;
    if (id > idv[n-1])
    {
      _newid = id - n;
    }
    else if (id > idv[0])
    {
      unsigned int k = 0;
      while (id > idv[k]) { ++k; };
      _newid -= k;
    }
    return(_newid);
  };

  size_t ctr = 0;
  for (size_t i = 0; i < this->size(); ++i)
  {
    if (!ids.count(i))
    {
      unsigned int bonds = (*this)[i].size();
      newbm[ctr].reserve(bonds);
      newtr[ctr].reserve(bonds);
      for (unsigned int j = 0; j < bonds; ++j)
      {
        unsigned int id_b = (*this)[i][j];
        if (!ids.count(id_b))
        {
          unsigned int id_new = newid(id_b);

          newbm[ctr].push_back(id_new);
          newtr[ctr].push_back(_translation[i][j]);
        }
      }
      newbm[ctr].shrink_to_fit();
      newtr[ctr].shrink_to_fit();

      ++ctr;
    }
  }


  this->swap(newbm);
  _translation.swap(newtr);

  newbm.clear();
  newtr.clear();

}



void
BondMap::print(void) const
{
  std::cout << std::endl;
  for (unsigned int i = 0; i < this->size(); i++)
  {
    std::cout<<"BondMap["<<i<<"]= ";
    for (unsigned int j = 0; j < (*this)[i].size(); j++)
    {
      std::cout<<(*this)[i][j];
      if (_translation[i][j].norm() > 1e-6)
        std::cout << "'";
      std::cout<<" ";
    }

    std::cout<<std::endl;
  }
}





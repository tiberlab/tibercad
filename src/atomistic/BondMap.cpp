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
 * \file BondMap.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/atomistic/BondMap.h"
#include "tibercad/atomistic/GridCells.h"
#include "tibercad/io/Messages.h"
#include "tibercad/io/Database.h"
#include "tibercad/io/Messages.h"
#include "tibercad/utils/Utils.h"

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
BondMap::solve(const std::vector<Atom>& basis,
    const Tensor2& period, const libMesh::Point& origin)
{

  Messages::debug("BondMap::solve");

  // nothing to be done if there are no atoms
  if (basis.empty())
    return;

  _basis = &basis;

  clear();
  _translation.clear();

  size_t na = basis.size();
  resize(na);
  _translation.resize(na);

  Tensor1 edge_min, edge_max;

  _period = period;
  
  //define the minimum spacing of the grid. the smaller it is, the faster is bonds calculations
  //cannot be smaller than the largest bond length (in Angstrom)
  GridCells cells(basis, period, origin, 8.0);
  //cells.print_statistics();

  // Loop on all cells
  Utils::Progress* prog = nullptr;
  if (basis.size() > 10000)
    prog = new Utils::Progress("BondMap", cells.size());

  for (unsigned int c1 = 0; c1 < cells.size(); c1++)
  {
    GridCells::NeighborIterator it = cells.begin(c1);
    GridCells::NeighborIterator end = cells.end(c1);
 
    if (prog != nullptr)
      prog->progress_message(c1+1);

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
          process_atoms(cells[c1][i], cells[c2][j], shift);
        }  
      }
    }
  }

  fix_bondmap();

  delete prog;
}


void
BondMap::process_atoms(const unsigned int i,
                       const unsigned int j, 
                       const Tensor1& period)
{
  if (_basis == nullptr)
    return;
  
  const std::vector<Atom>& basis = *_basis;

  const libMesh::Point per(period(1), period(2), period(3));
  libMesh::Point shift;
  shift(0) = per(0) * _period(1, 1) + per(1) * _period(1, 2) + per(2) * _period(1, 3);
  shift(1) = per(0) * _period(2, 1) + per(1) * _period(2, 2) + per(2) * _period(2, 3);
  shift(2) = per(0) * _period(3, 1) + per(1) * _period(3, 2) + per(2) * _period(3, 3);
  

  // 2026-01-26: We now use the covalent radii to decide possible bonding
  if ((i != j) || (per.norm() > 0.1))
  {

    double cutofftmp = basis[i].get_specie().get_covalent_radius() +
                       basis[j].get_specie().get_covalent_radius();
    // we are generous to not miss real bonds
    double f = 2.0;
    
    // tighten cutoff if the two atoms are the same
    // TODO check whether to look also at the label
    //if (basis[i].get_specie() == basis[j].get_specie())
    //  f = 1.2; // or even 1.15 for heavy cations

    cutofftmp *= f;

    const libMesh::Point& position1 = basis[i].get_position();
    libMesh::Point position2 = basis[j].get_position() + shift;
    libMesh::Point vector = position2 - position1;
    
    if (vector.norm() < cutofftmp)
    {

      bool add_bond = true;
      
      // check if it is already there
      for (unsigned int n = 0; n < (*this)[i].size(); n++)
      {
        if (((*this)[i][n] == j) && (libMesh::Point(get_translation(i, n) - shift).norm() < 1e-6))
        {
          add_bond = false;
          break;
        }

      }

      if (add_bond)
      {
        (*this)[i].push_back(j);
        _translation[i].push_back(per);

        (*this)[j].push_back(i);
        _translation[j].push_back(-per);
      }
      
    }
  }
}


void
BondMap::fix_bondmap(void)
{
  // we decide whether a bond should be kept by checking
  // - if the midpoint lies inside the Voronoi box around the
  //   center atom defined by the other atoms with shorter bonds
  // - check bond weight based on ECN

  if (_basis == nullptr)
    return;
  
  const std::vector<Atom>& basis = *_basis;
  
  // loop over all atoms
  for (size_t i = 0; i < this->size(); ++i)
  {
 
    // first we reorder the bonds in increasing length
    // (using insertion sort algorithm)
    for (unsigned int j = 1; j < (*this)[i].size(); ++j)
    {
      unsigned int tmp = (*this)[i][j];
      libMesh::Point tmp_p = _translation[i][j];

      double dist = libMesh::Point(basis[(*this)[i][j]].get_position() +
          get_translation(i, j) - basis[i].get_position()).norm();

      unsigned int k;

      for (k = j; k > 0; --k)
      {
        double dist2 = libMesh::Point(basis[(*this)[i][k-1]].get_position() +
          get_translation(i, k-1) - basis[i].get_position()).norm();

        // strict inequality is necessary for stable sorting
        if (dist2 < dist)
          break;

        (*this)[i][k] = (*this)[i][k-1];
        _translation[i][k] = _translation[i][k-1];
      }
      (*this)[i][k] = tmp;
      _translation[i][k] = tmp_p;

    }
    /*
    std::cerr << i << "  : ";
    for(int j = 0; j < std::min((*this)[i].size(), (size_t)7); ++j)
    {
        libMesh::Point normal = basis[(*this)[i][j]].get_position() +
            get_translation(i,j) - basis[i].get_position();
      std::cerr <<  "(" << (*this)[i][j] << ", " << normal.norm() << ") ";
    }
    std::cerr << "\n";
    */


    // use greater so that erasing happens starting from highest index
    std::set<size_t, std::greater<size_t>> flag;

    // factor for the testpoint for Voronoi box testing
    // (slightly outside of midpoint to prevent false neighbors).
    double alpha = 0.5 + 1e-3;

    // loop over all neighbors
    for (unsigned int j = 0; j < (*this)[i].size(); ++j)
    {
      if (flag.count(j) == 0)
      {
        libMesh::Point normal = basis[(*this)[i][j]].get_position() +
            get_translation(i, j) - basis[i].get_position();

        // the mid-perpendicular plane
        libMesh::Plane p(basis[i].get_position() + 0.5*normal, normal);

        // the shorter bonds cannot be masked by longer ones
        for(unsigned int k = j + 1; k < (*this)[i].size(); ++k)
        {
          if (flag.count(k))
            continue;

          libMesh::Point vector = basis[(*this)[i][k]].get_position() +
                                  get_translation(i, k) - basis[i].get_position();

          libMesh::Point testpoint = basis[i].get_position() + alpha * vector;

          if (p.above_surface(testpoint))
            flag.insert(k);
        }
      }
    }

    // now we delete the flagged bonds
    for (auto it = flag.begin(); it != flag.end(); ++it)
    {
      unsigned int id = *it;

      (*this)[i].erase((*this)[i].begin() + id);
      _translation[i].erase(_translation[i].begin() + id);
    }

    flag.clear();

    // Next step: calculate effective coordination number (ECN) and get a weight
    // to eliminate chemically improbable neighbors
    // see: Hoppe, Z. Kristallogr. 150, 23–52 (1979)

    std::vector<double> weights;
    double ecn = get_effective_coordination_number(i, weights);

    double wmin = 0.4;
    for (unsigned int j = 0; j < (*this)[i].size(); ++j)
    {
      if (weights[j] < wmin)
        flag.insert(j);
    }

  
    // now we delete once again the flagged bonds
    // this leads to an asymmetric bond map we have to adjust afterwards
    for (auto it = flag.begin(); it != flag.end(); ++it)
    {
      unsigned int id = *it;

      (*this)[i].erase((*this)[i].begin() + id);
      _translation[i].erase(_translation[i].begin() + id);
    }
  }

  // Now we have a potentially asymmetric bond map
  // we have to check that a bond seen by one atom also exists for the other

  for (size_t i = 0; i < this->size(); ++i)
  {
    for (unsigned int j = 0; j < (*this)[i].size(); ++j)
    {
      size_t neigh = (*this)[i][j];

      if (neigh == i)
        continue;

      bool found = false;

      for (unsigned k = 0; k < (*this)[neigh].size(); ++k)
      {
        if (((*this)[neigh][k] == i) &&
            ((get_translation(i, j) + get_translation(neigh, k)).norm() < 1e-6))
        {
          found = true;
          break;
        }
      }

      if (!found)
      {
        // we need to add it
        (*this)[neigh].push_back(i);
        _translation[neigh].push_back(-_translation[i][j]);
      }
    }
  }
}   


double
BondMap::get_effective_coordination_number(size_t atom, std::vector<double>& weights) const
{
  if (_basis == nullptr)
    return(0);

  const std::vector<Atom>& basis = *_basis;

  weights.resize((*this)[atom].size());

  double rmin = (basis[(*this)[atom][0]].get_position() +
            get_translation(atom, 0) - basis[atom].get_position()).norm();
  // reference distance, start with smallest bond
  double r0 = rmin;
  int n = 6;
  for (unsigned int it = 0; it < 2; ++it)
  {
    rmin = 0; // we use it in the following loop
    double den = 0;
    for (unsigned int j = 0; j < (*this)[atom].size(); ++j)
    {
      libMesh::Point vec = basis[(*this)[atom][j]].get_position() +
                           get_translation(atom, j) - basis[atom].get_position();
      double rij = vec.norm();

      double weight = std::exp(-std::pow(rij / r0, n));
      rmin += rij * weight;
      den += weight;
    }

    r0 = rmin / den;
  }

  // the effective coordination number
  double ECN = 0.0;
  double wmin = 0.4;

  for (unsigned int j = 0; j < (*this)[atom].size(); ++j)
  {
    libMesh::Point normal = basis[(*this)[atom][j]].get_position() +
                            get_translation(atom, j) - basis[atom].get_position();
    double rij = normal.norm();

    double wij = std::exp(1 - std::pow(rij / r0, n));
    weights[j] = wij;

    ECN += wij;
  }
   
  return(ECN);
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
BondMap::print() const
{
  if (_basis == nullptr)
    return;

  const std::vector<Atom>& basis = *_basis;

  std::cout << std::endl;
  for (unsigned int i = 0; i < this->size(); i++)
  {
    std::cout << i << " (" << basis[i].get_specie() << ") ";
    for (unsigned int j = 0; j < (*this)[i].size(); j++)
    {
      std::cout<<(*this)[i][j];
      if (_translation[i][j].norm() > 1e-6)
        std::cout << "'";
      std::cout<<" ";
    }
    std::cout << "  :  ";
    for (unsigned int j = 0; j < (*this)[i].size(); j++)
    {
      libMesh::Point dist(basis[(*this)[i][j]].get_position());
      dist = dist + get_translation(i, j) - basis[i].get_position();
      std::cout << dist.norm() << " ";
    }

    std::cout<<std::endl;
  }
}






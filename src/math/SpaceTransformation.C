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
 * \file SpaceTransformation.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/math/SpaceTransformation.h"
#include "tibercad/atomistic/CrystalDefs.h"

#include "libmesh/point.h"
#include "libmesh/tensor_value.h"

#include <map>
#include <algorithm>

using namespace std;
using namespace libMesh;

// here we keep the generators
namespace {

  // The basic generator matrices
  // sigma is written as s
  // a trailing p indicates a prime (')
  map<string, RealTensor> generators = {

      {"E",   RealTensor(1, 0, 0,
                         0, 1, 0,
                         0, 0, 1)},

      {"C2",  RealTensor(-1, 0, 0,
                         0, -1, 0,
                         0, 0, 1)},

      {"i",   RealTensor(-1, 0, 0,
                         0, -1, 0,
                         0, 0, -1)},

      {"C2p", RealTensor(1, 0, 0,
                         0, -1, 0,
                         0, 0, -1)},

      {"sh",  RealTensor(1, 0, 0,
                         0, 1, 0,
                         0, 0, -1)},

      {"sn",  RealTensor(-1, 0, 0,
                         0, 1, 0,
                         0, 0, 1)},

      {"snp", RealTensor(1, 0, 0,
                         0, -1, 0,
                         0, 0, 1)},

      {"C4",  RealTensor(0, -1, 0,
                         1, 0, 0,
                         0, 0, 1)},

      {"C4i", RealTensor(0, -1, 0,
                         1, 0, 0,
                         0, 0, -1)},

      {"P",   RealTensor(0, 0, 1,
                         1, 0, 0,
                         0, 1, 0)},

      {"Pi",  RealTensor(0, -1, 0,
                         0, 0, -1,
                         -1, 0, 0)},

      {"C6",  RealTensor(0.5, -sqrt(3)/2.0, 0,
                         sqrt(3)/2.0, 0.5, 0,
                         0, 0, 1)},

      {"C6i", RealTensor(-0.5, -sqrt(3)/2.0, 0,
                         sqrt(3)/2.0, -0.5, 0,
                         0, 0, -1)},

      {"C3",  RealTensor(-0.5, -sqrt(3)/2.0, 0,
                         sqrt(3)/2.0, -0.5, 0,
                         0, 0, 1)},

      {"C3i", RealTensor(0.5, -sqrt(3)/2.0, 0,
                         sqrt(3)/2.0, 0.5, 0,
                         0, 0, -1)},
  };




  // list of generators for known symmetries
  map<string, vector<string>> symmetries = {
      {"1", { } },
      {"-1", { "i" } },
      {"2", { "C2" } },
      {"m", { "sh" } },
      {"2/m", { "C2", "sh" } },

      {"222", { "C2", "C2p" } },
      {"mm2", { "C2", "sn" } },
      {"mmm", { "sh", "sn", "snp" } },

      {"4", { "C4" } },
      {"-4", { "C4i" } },
      {"4/m", { "C4", "sh" } },
      {"422", { "C4", "snp" } },
      {"4mm", { "C4", "sn" } },
      {"-42m", { "C4i", "snp" } },
      {"4/mmm", { "C4", "sn", "sh" } },

      {"3", { "C3" } },
      {"-3", { "C3i" } },
      {"32", { "C3", "snp" } },
      {"3m", { "C3", "sn" } },
      {"-3m", { "C3i", "sn" } },

      {"6", { "C6" } },
      {"-6", { "C6i" } },
      {"6/m", { "C6", "sh" } },
      {"622", { "C6", "C2p" } },
      {"6mm", { "C6", "sn" } },
      {"-6m2", { "C6i", "sn" } },
      {"6/mmm", { "C6", "sn", "sh" } },

      {"23",  { "C2", "P" } },
      {"m3", { "C2", "Pi" } },
      {"432",  { "C4", "P" } },
      {"-43m", { "C4i", "P" } },
      {"m3m", { "C4", "Pi" } },
  };

  // list of transformation matrices for used symmetries
  map<string, vector<RealTensor>> transformations;

  // fuzzy comparison of two points
  bool compare_points(const Point& a, const Point& b)
  {
    return(a.absolute_fuzzy_equals(b));
  }

  // fuzzy compare to RealTensors
  bool compare_tensors(const RealTensor& a, const RealTensor& b)
  {
    RealTensor diff(a - b);
    return(diff.norm() < 1e-9);
  }
}

void
SpaceTransformation::rotate(const RealVectorValue& axis,
                            const double angle,
                            Point& point)
{
  RealVectorValue u(axis.unit());
  double ux = u(0);
  double uy = u(1);
  double uz = u(2);

  double cost = std::cos(angle);
  double sint = std::sin(angle);

  RealTensor R;

  R(0,0) = cost + ux*ux*(1-cost);
  R(0,1) = ux*uy*(1-cost) - uz*sint;
  R(0,2) = ux*uz*(1-cost) + uy*sint;
  R(1,0) = ux*uy*(1-cost) + uz*sint;
  R(1,1) = cost + uy*uy*(1-cost);
  R(1,2) = uy*uz*(1-cost) - ux*sint;
  R(2,0) = ux*uz*(1-cost) - uy*sint;
  R(2,1) = uz*uy*(1-cost) + ux*sint;
  R(2,2) = cost + uz*uz*(1-cost);

  point = R*point;
}


void
SpaceTransformation::create_star(const string& symmetry,
                                 const Point& point,
                                 vector<Point>& star)
{
  // convert given symmetry to International Symbol
  string isym = CrystalDefs::convert_to_international_symbol(symmetry);

  auto& trafos = transformations[isym];


  if (trafos.empty())
    generate_transformations(isym);

  star.clear();

  for (auto&& R : trafos)
    star.push_back(R * point);

  sort(star.begin(), star.end());

  auto ip = unique(star.begin(),
                   star.end(),
                   compare_points);
  star.resize(distance(star.begin(), ip));
}


void
SpaceTransformation::generate_transformations(const std::string& symmetry)
{

  auto& trafos = transformations[symmetry];
  trafos.push_back(generators["E"]);

  vector<RealTensor> gen;

  for (auto&& g : symmetries[symmetry])
  {
    const RealTensor& gt = generators[g];
    gen.push_back(gt);
    trafos.push_back(gt);
  }


  // we implement here a brute force approach, but I think
  // this is ok due to the small number of <= 3 of generators
  int n_trafo = 0;
  do
  {
    n_trafo = trafos.size();
    for (int i = 1; i < n_trafo; ++i)
    {
      for (auto&& G : gen)
      {
        RealTensor R = G * trafos[i];
        bool add = true;
        for (int j = 0; j < trafos.size(); ++j)
        {
          if (compare_tensors(R, trafos[j]))
          {
            add = false;
            break;
          }
        }

        if (add)
          trafos.push_back(R);

        RealTensor R2 = trafos[i] * G;
        if (!compare_tensors(R, R2))
        {
          add = true;
          for (int j = 0; j < trafos.size(); ++j)
          {
            if (compare_tensors(R2, trafos[j]))
            {
              add = false;
              break;
            }
          }

          if (add)
            trafos.push_back(R2);
        }
      }
    }
  }
  while (n_trafo < trafos.size());

  //cerr << symmetry << ": " << trafos.size() << "\n";
  //for (int i = 0; i < trafos.size(); ++i)
  //  cerr << trafos[i] << "\n";
  //cerr << "\n";
}

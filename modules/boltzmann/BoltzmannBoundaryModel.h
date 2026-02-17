/*  
 * This file is part of the tiberCAD module boltzmann.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file BoltzmannBoundaryModel.h
 * \brief tiberCAD boltzmann module header.
 *
 * \note This file is part of module boltzmann.
 */


#ifndef TC_BOLTZMANNBOUNDARYMODEL_H
#define TC_BOLTZMANNBOUNDARYMODEL_H

#include "tibercad/physics/PhysicalModel.h"
#include "vector_value.h"

namespace libMesh
{
  class Point;
  class Elem;
}


//! The base class for thermal balance boundary conditions
class BoltzmannBoundaryModel : public PhysicalModel
{

  public:

    //! Destructor
    ~BoltzmannBoundaryModel(void) {};

    //! Creator function
    static BoltzmannBoundaryModel* create(const MaterialBoundary* boundary,
        const ModelOptions& options);


    //! Calculate for a point on the given side
    virtual void calculate(const Elem* elem, unsigned int side,
        const Point& point) = 0;

    //! Calculate the periodic vector
    libMesh::RealGradient get_periodicity(void);

    //! Calculate the periodic vector
    double get_deltaT(void);

    //! Get coefficients
    void get_coefficients(double& a, double& b, double& c);

   ///!Set the current element
   void set_element(const Elem* elem);

  protected:

  //! Constructor
  BoltzmannBoundaryModel(const ModelOptions& options);

  void set_coefficients(double a, double b, double c);

  //! Set periodicity
  void set_periodicity(const libMesh::RealGradient& periodicity);

  //! Set deltaT
  void set_deltaT(double deltaT);


  private:

  double _alpha;
  double _beta;
  double _gamma;
  double _deltaT;

  libMesh::RealGradient _periodicity;


};


inline
BoltzmannBoundaryModel::BoltzmannBoundaryModel(const ModelOptions& options) :
  PhysicalModel(options)
{
}


inline
void
BoltzmannBoundaryModel::get_coefficients(double& a, double& b, double& c)
{
 a = _alpha;
 b = _beta;
 c = _gamma;
}


inline
double
BoltzmannBoundaryModel::get_deltaT(void)
{
 return _deltaT;
}

inline
libMesh::RealGradient
BoltzmannBoundaryModel::get_periodicity(void)
{
 return _periodicity;
}

inline
void
BoltzmannBoundaryModel::set_coefficients(double a, double b, double c)
{
 _alpha = a;
 _beta  = b;
 _gamma = c;
}


inline
void
BoltzmannBoundaryModel::set_periodicity(const libMesh::RealGradient& periodicity)
{
_periodicity = periodicity;
}

inline
void
BoltzmannBoundaryModel::set_deltaT(double deltaT)
{
_deltaT = deltaT;
}



#endif // TC_THERMALBOUNDARYMODEL_H

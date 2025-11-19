// $Id$

#ifndef _TMMBOUNDARYMODEL_H_
#define _TMMBOUNDARYMODEL_H_

#include "tibercad/physics/PhysicalModel.h"
#include "Tmm.h"


namespace libMesh
{
  class Elem;
  class Point;
}

//! The base class for thermal balance boundary conditions
class TmmBoundaryModel : public  PhysicalModel
{

public:

  //! Destructor
  ~TmmBoundaryModel(void) {};

  //! Creator function
  static TmmBoundaryModel* create(const MaterialBoundary* boundary,
    const ModelOptions& options);

  virtual void Calculate_M_Matrix(void) = 0;
  void set_elements(double, double, double, double);
  virtual double get_element(int);
  virtual double get_kr(void);
  virtual double get_steps(void);
  void set_dipole_elements(double, double);
  std::string read_type(void);








protected:

  //! Constructor
  TmmBoundaryModel(const ModelOptions& options);
  void write_type(std::string);





private:
  std::string typer;
  double _mmm00;
  double _mmm01;
  double _mmm10;
  double _mmm11;
  double __kr;
  double __steps;



};


inline
TmmBoundaryModel::TmmBoundaryModel(const ModelOptions& options) :
  PhysicalModel(options)
{
}



#endif // _TMMBOUNDARYMODEL_H_

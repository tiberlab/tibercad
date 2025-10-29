// $Id: TmmBulkModel.h 4391 2017-04-07 11:16:58Z pecchia $

#ifndef _TmmDipoleSource_H_
#define _TmmDipoleSource_H_

#include "PhysicalModel.h"
#include "elem.h"
// #include "Tmm.h"

using namespace std;

//! This is the base class for the TMM bulk physical model
class TmmDipoleSource : public PhysicalModel
{

public:


  //! Destructor
  virtual ~TmmDipoleSource(void) {};

  //! Creator function
  static TmmDipoleSource* create(const ModelOptions& options);


  const double& get_emission_power(void) const;


  virtual void calculate(const libMesh::Elem* elem, const libMesh::Point& point, double lambda) {};






protected:

  //! Constructor
  TmmDipoleSource(const ModelOptions& options);


  void set_emission_power(const double& emission_power);


private:


  double _emission_power;



};




#endif // _TMMBULKMODEL_H_

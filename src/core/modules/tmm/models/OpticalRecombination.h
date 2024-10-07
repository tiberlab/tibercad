// $Id: PermittivityModel.h 4729 2018-12-05 07:58:16Z maufder $

#ifndef _OpticalRecombination_H_
#define _OpticalRecombination_H_

#include "Tmm.h"
#include "TmmDipoleSource.h"

namespace libMesh
{
  class Elem;
}

// Base class for InCoherence model
class TBDLEXPORT OpticalRecombination : public TmmDipoleSource
{

  public:

  virtual ~OpticalRecombination(void) {};
  static OpticalRecombination* create(const ModelOptions& options);
  virtual void calculate(const libMesh::Elem* elem, const libMesh::Point& point, double lambda);

  
protected:
  
    OpticalRecombination(const ModelOptions& options);
	
	virtual void do_init(void){};





  private:

   double _recombination_rate;
   
       //! The generation model
    std::vector<SimulationInterface*> _recombination_model;

    //! The solution ID of the generation models variable
    std::vector<ID> _recombination_id;

};

inline
OpticalRecombination::OpticalRecombination(const ModelOptions& options) :
  TmmDipoleSource(options),
  _recombination_rate(0.0)
{
}

inline
OpticalRecombination*
OpticalRecombination::create(const ModelOptions& options)
{
  return new OpticalRecombination(options);
}



#endif // _POLARIZATIONMODEL_H_

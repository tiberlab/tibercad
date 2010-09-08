// $Id$

#ifndef _DDHEATSOURCE_H_
#define _DDHEATSOURCE_H_

#include "HeatSourceModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "PhysicalModelInterface.h"
#include "SimulationInterface.h"
#include "tiber_dll.h"

class Elem;


//! The base class for Poisson boundary conditions
class TBDLLOCAL DDHeatSource : public HeatSourceModel
{

  public:
 
     //! Destructor
  ~DDHeatSource(void) {};
  
  //! Creator function
  static DDHeatSource* create(const ModelOptions& options);

  virtual void calculate(const Elem* elem, const Point& point);
    

  protected:

    //! Initialize
    virtual void do_init(void);

    /* In some cases it might be useful to reimplement this: */
    // virtual void do_init_interface(const PhysicalModelInterface* comp_A,
    //         const PhysicalModelInterface* comp_B);


    /* This is not used here: */
    // virtual void read_database(void);


    /* We do not use this here: */
    // virtual void read_interface_database(void);


    //! Create a new object of the same type
    virtual PhysicalModelInterface* create_new(void) const;


  private:
  
    double _heat_source;

  enum heat_variables
    {
      EJOULE = 0,
      HJOULE,
      RECHEAT,
      EPELTH,
      HPELTH,
      WNX,
      WPX
    };
  
  //!Heat source variables for drift diffusion
  std::set<ID> ID_set;
  
  //!Variable map
  std::map<ID,ID> var_map;

  //!Pointer to drift diffusion simulation
  SimulationInterface* _simul;

  //! Constructor
    DDHeatSource(const ModelOptions& options);
  
};


inline
PhysicalModelInterface*
DDHeatSource::create_new(void) const
{
  return new   DDHeatSource(get_options());
}

inline
DDHeatSource*
DDHeatSource::create(const ModelOptions& options)
{
  return new  DDHeatSource(options);
}




#endif // _GRAYMODEL_H_

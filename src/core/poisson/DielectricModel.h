#ifndef _DIELECTRICMODEL_H_
#define _DIELECTRICMODEL_H_

#include "PhysicalModelInterface.h"

#include "Constants.h"

//! This class read the dielectric constant from database
/*!


*/


class DielectricModel : public PhysicalModelInterface
{

public:
  
  //!Constructor 
   DielectricModel(const ModelOptions& options);
 
   //!Destructor
  ~DielectricModel(){}

   //!provides electrons thermoelectric power [V/K]
  Tensor2Sym  get_dielectric_constant(void) const;
  
  static  DielectricModel* create(const ModelOptions& options);


private:
  
   //!dielectric constant scalar
  double _ep_x; 
 
  //!dielectric constant tensor
  Tensor2Sym _dielectric_constant; 
  
protected:

  
  virtual void do_init (void);

  virtual void read_database(void);

  virtual void read_database_alloy(void){};

  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa); 


  virtual PhysicalModelInterface* create_new (void) const;
 
 
 
};



inline
Tensor2Sym
DielectricModel::get_dielectric_constant(void) const
{
  return _dielectric_constant * Constants::epsilon;
}


inline
DielectricModel*  DielectricModel::create(const ModelOptions& options)
{
  return (new  DielectricModel(options));
}

inline
PhysicalModelInterface*  DielectricModel::create_new () const
{
  return (new  DielectricModel(get_options()) );
}




#endif

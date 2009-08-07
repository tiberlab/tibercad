#ifndef _PIEZOELECTRICMODEL_H_
#define _PIEZOELECTRICMODEL_H_


#include "PhysicalModelInterface.h"
#include "SimulationInterface.h"

//! A class that containes Young modules for the elasticity problem

class PiezoelectricModel : public PhysicalModelInterface
{
 public:

  PiezoelectricModel();

  ~PiezoelectricModel(){};

  //! creates new object
  static PiezoelectricModel* create(const std::string& name,  const ModelOptions& options);
 
  void get_piezopolarization(Tensor1& P);

  virtual void calculate_piezopolarization(const Elem* elem) = 0;

  SimulationInterface* _simul;

 protected:

  Tensor1 _P; 

  void rotate_to_calc_system(const Tensor2Gen& RotMatrix);

  virtual void read_database(void)=0;

  virtual void do_init(void) = 0;

  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);

  virtual PhysicalModelInterface* create_new(void) const = 0;


};

//------------------------------------------------------------------------------------//
inline PiezoelectricModel* PiezoelectricModel::create( const std::string& name,  const ModelOptions& options )
{
  return dynamic_cast<PiezoelectricModel*>(PhysicalModelInterface::create("piezoelectric_model_" + name, options));
}

inline
void 
PiezoelectricModel::get_piezopolarization(Tensor1& P)
{

  P = _P;

}

#endif

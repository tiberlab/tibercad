#include "tensor.h"
#include "xtensor.h"
#include <cmath>

#include <vector>

class Piezoelectricity
{
 public:
  Piezoelectricity();
 

  void set_moduli(double e14);

  void set_moduli(double  e31, double e15, double e33);


  void set_pyro_module(double p);

  Tensor1 get_polariz_cryst(Tensor2Sym& strain_cryst);


   
 private:

  std::string             type;

  std::vector<double>   moduli;

  double            pyro_const;            


};

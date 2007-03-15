#ifndef _PIEZOELECTRICITY_H_
#define _PIEZOELECTRICITY_H_

#include "tensor.h"
#include "xtensor.h"
#include <cmath>
#include "PhysicalModel.h"
#include "PhysicalModelInterface.h"
class Piezoelectricity : public PhysicalModelInterface
{
 public:

  //!Empty constructor
  Piezoelectricity();

  //! returns polarization (piezo + pyro) in crystal system
  virtual Tensor1  get_polariz_cryst(Tensor2Sym& strain_cryst) = 0;



  //!calculate in crystal system \f$ r_{jk} = d_{i,jk}f_i  \f$
  /*!
    \param f input vector in crystal system
    \param r result in a crystal system
   */
  virtual void calculate_product_by_vector(const Tensor1& f, Tensor2Sym& r) const {} ;


  //!returns \f$ q = d_{i,jk} f_i e_{jk} \f$
  virtual double calculate_product_by_vector_and_tensor(const Tensor1& f, const Tensor2Sym& e) const {};

  
  static Piezoelectricity* create(const std::string& name,  const ModelOptions& options);

 protected:


  virtual void read_database ( ) = 0;


  virtual void do_init(void) = 0;


  virtual void copy_from (const PhysicalModelInterface *rhs) = 0;


  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa) = 0;

  
  virtual PhysicalModelInterface* create_new(void) const = 0;



};


inline Piezoelectricity* Piezoelectricity::create(const std::string& name,  const ModelOptions& options)
{
  
  return dynamic_cast<Piezoelectricity*>(PhysicalModelInterface::create("piezo_" + name, options));

}


#endif

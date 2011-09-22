// $Id$

#ifndef _PIEZOELECTRICITY_H_
#define _PIEZOELECTRICITY_H_

#include "tensor.h"
#include "xtensor.h"
#include "PhysicalModelInterface.h"


class Piezoelectricity : public PhysicalModelInterface
{
 public:

  //! returns polarization (only piezo) in crystal system
  virtual Tensor1  get_polariz_cryst(Tensor2Sym& strain_cryst) = 0;



  //!calculate in crystal system \f$ r_{jk} = d_{i,jk}f_i  \f$
  /*!
    \param f input vector in crystal system
    \param r result in a crystal system
   */
  virtual void calculate_product_by_vector(const Tensor1& f, Tensor2Sym& r) const {} ;


  //!returns \f$ q = d_{i,jk} f_i e_{jk} \f$
  inline double calculate_product_by_vector_and_tensor(const Tensor1& f, const Tensor2Sym& e) const;

  
  static Piezoelectricity* create(const Material* mat, const ModelOptions& options);


 protected:

  //!Empty constructor
  Piezoelectricity(const ModelOptions& options);


  virtual void read_database ( ) = 0;


  virtual void do_init(void) = 0;


  virtual PhysicalModelInterface* create_new(void) const = 0;



};


inline Piezoelectricity* Piezoelectricity::create(const Material* mat, const ModelOptions& options)
{
  
  std::string structure = mat->get_structure();
  return dynamic_cast<Piezoelectricity*>(PhysicalModelInterface::create("piezo_" + structure,
      mat, options));

}


inline double Piezoelectricity::calculate_product_by_vector_and_tensor(const Tensor1& f, const Tensor2Sym& e) const
{
  Tensor2Sym temp;
  calculate_product_by_vector(f, temp);
  return (doubleContraction(temp, e) );
}

#endif

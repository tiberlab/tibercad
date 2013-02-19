#ifndef _EFABULKMODEL_H_
#define _EFABULKMODEL_H_

#include "PhysicalModel.h"
#include "EFAbulkHamiltonian.h"

//! A container class that contains a dynamical object of bulk Hamiltonian
class EFAbulkModel: public PhysicalModel
{

 public: 

  //!Destructor
  ~EFAbulkModel();

  //!creates a new object
  static EFAbulkModel* create(const ModelOptions& options);

  EFAbulkHamiltonian* get_Hamiltonian_model(void) const;

 protected:

  //!Constructor
  EFAbulkModel(const ModelOptions& options);

  virtual PhysicalModelInterface* create_new (void) const;


  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);

  virtual void do_init();
  
  virtual void do_print_info();



 private:

  //! a pointer to the bulk Hamiltonian
  EFAbulkHamiltonian*  _bulkHamiltonian;


};


inline EFAbulkModel*  EFAbulkModel::create(const ModelOptions& options)
{
  return new EFAbulkModel(options);
}

inline  PhysicalModelInterface* EFAbulkModel::create_new (void) const
{
  return new EFAbulkModel(get_options());
}

inline  EFAbulkHamiltonian* EFAbulkModel::get_Hamiltonian_model(void) const
{
  return  _bulkHamiltonian;
}

#endif

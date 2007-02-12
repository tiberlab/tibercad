#ifndef _EFABULKMODEL_H_
#define _EFABULKMODEL_H_

#include "PhysicalModel.h"
#include "EFAbulkHamiltonian.h"

//! A container class that contains a dynamical object of bulk Hamiltonian
class EFAbulkModel: public PhysicalModel
{

 public: 
  //!Constructor
  EFAbulkModel();

  //!Destructor
  ~EFAbulkModel();

  //!creates a new object
  static EFAbulkModel* create();

  EFAbulkHamiltonian* get_Hamiltonian_model(void) const;

 protected:

  virtual PhysicalModelInterface* create_new (void) const;


  virtual void copy_from(const PhysicalModelInterface *rhs) {};


  virtual void read_database (void) {};
 

  virtual void read_bowing_parameters (void) {};


  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);

  virtual void do_init();
  


 private:

  //! a pointer to the bulk Hamiltonian
  EFAbulkHamiltonian*  _bulkHamiltonian;

  //!copy constructor should not be used
  EFAbulkModel(const EFAbulkModel & t ) {};

};


inline EFAbulkModel*  EFAbulkModel::create()
{
  return new EFAbulkModel();
}

inline  PhysicalModelInterface* EFAbulkModel::create_new (void) const
{
  return new EFAbulkModel();
}

inline  EFAbulkHamiltonian* EFAbulkModel::get_Hamiltonian_model(void) const
{
  return  _bulkHamiltonian;
}

#endif

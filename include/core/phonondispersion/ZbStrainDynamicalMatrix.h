// $Id$

#ifndef _ZB_FREESTRAINDYNAMICALMATRIX_H_
#define _ZB_FREESTRAINDYNAMICALMATRIX_H_


#include "DynamicalMatrix.h"

#include <set>

class PhononModel;
class SimulationInterface;


class ZbStrainDynamicalMatrix: public DynamicalMatrix
{
 public:
  //!constructor
  ZbStrainDynamicalMatrix(const ModelOptions& options) ;

  //!destructor
  ~ZbStrainDynamicalMatrix() {};

  //! Create a ZbLatticeThermalConductivity object
  static  ZbStrainDynamicalMatrix* create(const ModelOptions& options);


  //! Update the lattice thermal conductivity given the Temperature
  virtual void re_init(); 


  virtual void set_phonon_model(PhononModel* phonon_model);
 
 private:

  PhononModel* _phonon_model;
  double p_norm;
  double q_norm;
  double r_norm;
  double w0;
  Tensor4DSym deformation_potential;
  SimulationInterface* _simul;


 protected:

  virtual void read_database(void);

  virtual void do_init(void);

  virtual PhysicalModelInterface*  create_new (void) const;

 private:

  enum strain_variables
    {
      E_XX = 0,
      E_XY,  
      E_XZ,
      E_YY,
      E_YZ,
      E_ZZ
    };

   //!Strain variables 
   std::set<ID> ID_set;

   //!Variable map
   std::map<ID,ID> var_map;



};

inline
ZbStrainDynamicalMatrix* ZbStrainDynamicalMatrix::create(const ModelOptions& options)
{
  return (new ZbStrainDynamicalMatrix(options));
}

inline
PhysicalModelInterface*  ZbStrainDynamicalMatrix::create_new (void) const
{
  return (new  ZbStrainDynamicalMatrix(get_options()) );
}

inline
void ZbStrainDynamicalMatrix::set_phonon_model(PhononModel* phonon_model)
{
  _phonon_model = phonon_model;
}



#endif

#ifndef _KPBULKHAMILTONIAN_H_
#define _KPBULKHAMILTONIAN_H_

#include <complex>
#include <vector>
#include "tensor.h"
#include "EFAbulkHamiltonian.h"
#include "PhysicalModel.h"
#include "Semiconductor.h"
#include "KPparameters.h"

/*!
 * \brief A class that builds kp bulk Hamiltonian
 *
 * This class handles in a unified way different kp models
 * with different number of bands. The Bloch basis functions
 * and there ordering are as follows:
 *
 * 2x2:
 * |CB> |VB>
 *
 * 6x6:
 * |X_u> |Y_u> |Z_u>  |X_d> |Y_d> |Z_d>
 *
 * 8x8:
 * |S_u> |S_d> |X_u> |Y_u> |Z_u>  |X_d> |Y_d> |Z_d>
 *
 */
class KPbulkHamiltonian : public EFAbulkHamiltonian
{
 public:

 
  //! destructor. 
  ~KPbulkHamiltonian(void);
  

  virtual void calculate_Hamiltonian_k_par(void);

 
 
  virtual void calculate_Hamiltonian_gen(void); 



  virtual void apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential);


  //!calculates momentum operator P without k|| application
  void calculate_optical_operator(void);

  //!calculates momentum operator P with k|| application
  void calculate_optical_operator_k_par(void);

 
  void set_parameters(const KPparams&  par1 );


  const std::vector< std::vector <std::vector<MatrixElement> > >& get_optical_operator(void) const;

  //! Get Hamiltonian as coefficients for k-expansion
  void get_hamiltonian_without_k(std::vector<std::vector<MatrixElement> >& ham) const;


  static KPbulkHamiltonian* create(const ModelOptions& options);


  //!set semiconductor
  void set_semiconductor(Semiconductor* semicond);


  //! sets temperature
  virtual void set_temperature(double Temperature);
 

 private:
  //!simmetrize valence-valence term 
  bool kpVVtermSymmetric;
  
  //!simmetrize conduction-valence term 
  bool kpCVtermSymmetric;

  //model_name  name of the model "8x8" or "6x6"
  std::string model_name;
  

  //! a pointer to a semiconductor that contains parameters
  Semiconductor*  semiconductor;

  //! minimal used band in 8x8 Hamiltonian 
  short band_min;

  //! maximal used band in 8x8 Hamiltonian 
  short band_max;

  //! number of bands (8 or 14)
  short num_bands;
 
  //! Hamiltonian without k application and without strain
  std::vector<std::vector<MatrixElement> > Ham; 


  //! P-operator matrixes with k|| applied
  /*
    P[i1][i2][i3]:
    i1 - P-vector component number in crystal system: 0 - "x", 1 - "y", 2 - "z"
    i2, i3 - band indexes like in the Hamiltonian matrix.
  */
  std::vector< std::vector <std::vector<MatrixElement> > > P; 


  //! P-operator matrixes without k|| applied
  std::vector< std::vector <std::vector<MatrixElement> > > P_gen; 

  //! nullify k.p  parameters
  void   nullify_parameters(void); 

  //! k.p wurztzite parameters
  KPparams par;



 protected:

  //! default constructor.
  KPbulkHamiltonian(const ModelOptions& options);

  virtual PhysicalModel* create_new(void) const;

  virtual void do_init(void);

  virtual void do_print_info(void);
  
  virtual void prepare_submodels(void);
 
};


inline   KPbulkHamiltonian* KPbulkHamiltonian::create(const ModelOptions& options)
{
  return new KPbulkHamiltonian(options);
}

inline PhysicalModel* KPbulkHamiltonian::create_new() const
{
  return new KPbulkHamiltonian(get_options());
}

inline
void  KPbulkHamiltonian::set_temperature(double Temperature)
{
  semiconductor->set_temperature(Temperature);
  semiconductor->apply_temperature();
  calculate_Hamiltonian_gen();
}




#endif

#ifndef _KPBULKHAMILTONIAN_H_
#define _KPBULKHAMILTONIAN_H_
//! A class that builds kp bulk Hamiltonian

#include <complex>
#include <vector>
#include "tensor.h"
#include "EFAbulkHamiltonian.h"
#include "PhysicalModelInterface.h"
#include "Semiconductor.h"
#include "KPparameters.h"
typedef std::complex<double> Complex;

class KPbulkHamiltonian : public EFAbulkHamiltonian
{
 public:

  




 


  //! default constructor. 
  KPbulkHamiltonian(void);



 
  //! destructor. 
  ~KPbulkHamiltonian(void);
  

  virtual void calculate_Hamiltonian_k_par(void);

 
 
  virtual void calculate_Hamiltonian_gen(void); 



  virtual void apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential);


  //!calculates momentum operator P without k|| application
  void calculate_optical_operator(void);

  //!calculates momentum operator P with k|| application
  void calculate_optical_operator_k_par(void);

 
  Tensor2Sym strainM;

  void set_parameters(const KPparams&  par1 );


  const std::vector< std::vector <std::vector<MatrixElement> > > & get_optical_operator(void) const;


  static KPbulkHamiltonian* create( );
 

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

 
  //! Hamiltonian without k|| application
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
  KPparams par ;

 
 
 

 protected:

  virtual PhysicalModelInterface* create_new(void) const;

  virtual void copy_from (const PhysicalModelInterface *rhs){};

  virtual void do_init(void);

  virtual void read_database(void){};

  virtual void read_bowing_parameters(void){};
 
  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa) ;
 

 
  
};


inline   KPbulkHamiltonian* KPbulkHamiltonian::create( )
{
  return new KPbulkHamiltonian();
}

inline PhysicalModelInterface* KPbulkHamiltonian::create_new() const
{
  return new KPbulkHamiltonian();
}


#endif

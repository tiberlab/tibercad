#ifndef _KPBULKHAMILTONIAN_H_
#define _KPBULKHAMILTONIAN_H_
//! A class that builds kp bulk Hamiltonian
#include <complex>
#include <vector>
#include "tensor.h"
#include "EFAbulkHamiltonian.h"

typedef std::complex<double> Complex;

class KPbulkHamiltonian : public EFAbulkHamiltonian
{
 public:

  

 struct  KPparams
  {

   

    double L1;
    double L2;
    double M1;
    double M2;
    double M3;
    double N1;
    double N2;    
    double P1;
    double P2;
    double s1;
    double s2;
    double E_c;
    double E_v;
    double d1;
    double d2;
    double d3;
    double N1_xy;
    double N1_yx; 
    double N2_xy; 
    double N2_yx;
    double l1s;
    double l2s;
    double n1s;  
    double n2s;
    double m1s;
    double m2s;
    double m3s;
    double axs;
    double azs;
  };


 


  //! default constructor. 
  KPbulkHamiltonian(void);


  //! constructor
  /*!
    \param model_name  name of the model "8x8" or "6x6"
  */
  KPbulkHamiltonian(const std::string model_name );
  

 

  


  virtual void calculate_Hamiltonian_k_par(void);

 

  virtual void calculate_Hamiltonian_gen(void); 



  virtual void apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential);


  bool kpVVtermSymmetric;
  Tensor2Sym strainM;

  bool kpCVtermSymmetric;

  void set_parameters(const KPparams&  par1 );




 private:

  //! minimal used band in 8x8 Hamiltonian 
  short band_min;


  //! maximal used band in 8x8 Hamiltonian 
  short band_max;

 
  //! Hamiltonian without k|| application
  std::vector<std::vector<MatrixElement> > Ham; 



  void   nullify_parameters(void); 

  KPparams par ;

 
  
};
#endif

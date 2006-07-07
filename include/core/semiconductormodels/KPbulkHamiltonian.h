#ifndef _KPBULKHAMILTONIAN_H_
#define _KPBULKHAMILTONIAN_H_
//! A class that builds kp bulk Hamiltonian

#include <complex>
#include <vector>
#include "tensor.h"
#include "EFAbulkHamiltonian.h"


class DDsemiconductor;

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
  

  //!Copy-constructor
  KPbulkHamiltonian(const KPbulkHamiltonian& kp_ham);
 
  //! destructor. 
  ~KPbulkHamiltonian(void);
  
  virtual void read_database(const Dummy&);

  //! \deprecated { Create parameters for an alloy }
  /*!
   * \deprecated { This method will live as long as the database is
   * not used yet.}
   */

  virtual void build_alloy(const std::string& component2,
			   const std::string& bowing_params, double content);



  void set_data_file(const std::string& filename)
      { _filename = filename; };


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


  //! nullify k.p  parameters
  void   nullify_parameters(void); 

  //! k.p wurztzite parameters
  KPparams par ;

 
  std::string _filename;

  //model_name  name of the model "8x8" or "6x6"
  std::string model_name;

  //! a pointer to a semiconductor that contains parameters
  DDsemiconductor* semiconductor;
 

};
#endif

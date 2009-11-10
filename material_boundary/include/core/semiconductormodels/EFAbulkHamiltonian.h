#ifndef _EFAbulkHamiltonian_H_
#define _EFAbulkHamiltonian_H_
//! A general class for envelope function bulk Hamiltonian

#include <complex>
#include <vector>
#include <map>
#include "tensor.h"

#include "PhysicalModelInterface.h"
class EFAbulkHamiltonian: public PhysicalModelInterface
{

 public:

  typedef std::complex<double> Complex;

  struct MatrixElement
  {
    Complex constant;
    Complex linear_left[3];
    Complex linear_right[3];
    Complex quad[3][3];
    
  };

  //!constructor
  EFAbulkHamiltonian();


  //!destructor
  virtual ~EFAbulkHamiltonian() {};
  

  //!set Bloch k-vector
  /*!
    \param kvector - k vector [a.u.]
  */
  void set_k_vector(const double kvector[3]);
 
  //!set Bloch k-vector
  /*!
    \param k_in - k vector [a.u.]
  */  
  void set_k_vector(Tensor1 k_in); 

 

  //! calculate model Hamiltonian without application of k||
  virtual void calculate_Hamiltonian_gen(void) = 0; 
 
  //! apply k|| to the Hamiltonian
  virtual void calculate_Hamiltonian_k_par(void) = 0;


  //! apply strain and potential to the EFA Hamiltonian
  /*!
    \param strain_crystal strain tensor in crystal system
    \param el_potentia electric potential [V]
  */
  virtual void apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential) = 0;


  //!returns Hamiltonian
  std::vector< std::vector<MatrixElement > >& get_Hamiltonian(void); 

 
  const std::map<short, short>&  get_kp_bands_map(void) const; 

  //!creates a new object
  /*!
    \param structure crystal structure, e.g. "zb", "wz"
    \param options options for the model
  */
  static EFAbulkHamiltonian* create (const std::string& structure,  const ModelOptions& options = ModelOptions());


  //! sets temperature
  virtual void set_temperature(double Temperature) {};


 protected:
  

  virtual PhysicalModelInterface* create_new(void) const = 0;

  virtual void do_init (void);

  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa) = 0;
  

  //!result Hamiltonian in k representation 
  std::vector< std::vector<MatrixElement > > Hamiltonian;


  //! Hamiltonian in k representation that is used by apply_strain_and_potential
  std::vector< std::vector<MatrixElement > > Hamiltonian_without_strain_pot;

  //! k-vector in simualtion system
  double k_vector[3];


  //! Hartree energy in eV
  static const double Hartree;


  //!rotation matrix
  double  rot_matrix[3][3];
  
  //!rotates tensor of rank 1 from crystal to calculation system
  void    rotate_linear(std::complex<double> *vector);

  //!rotates tensor of rank 2 from crystal to calculation system
  void    rotate_quad(std::complex<double> matrix[][3]);

  //!numbers of the bands in a 8x8 k.p basis
  /*
    0,1 - conduction bands;
    2-7 - valence bands
  */
  std::vector<short> kp_bands;


  //!map between band numbers  
  std::map <short, short> kp_bands_map;
  
   void set_rotation_matrix(void);
 

 private:
};


inline  EFAbulkHamiltonian* EFAbulkHamiltonian::create (const std::string& name,  const ModelOptions& options)
{

  if (! (options.find_option("model")) )
  {
    std::cerr << "EFAbulkHamiltonian* EFAbulkHamiltonian::create   model must be specified \n"; 
    options.print_all();
    exit(1);
  }

  const std::string&  model_name = options.get_option("model", ""); 
  
  std::string model;

  if ( model_name == "kp")
    //model = "quantum_kp_" + name;
    model = "quantum_kp";
  else if ( model_name == "sb_user_defined")
    model = "quantum_user";
  else if ( model_name == "conduction_band")
    model ="quantum_cond_band_" + name; 

  return dynamic_cast<EFAbulkHamiltonian*> ( PhysicalModelInterface::create(model, options) );

}

#endif

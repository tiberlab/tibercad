#ifndef _STIFFNESS_H_
#define _STIFFNESS_H_

#include "tensor.h"
#include "xtensor.h"
#include <cmath>
#include "PhysicalModelInterface.h"
#include "PhysicalModel.h"
#include "SimulationInterface.h"


//! A class that containes Young modules for the elasticity problem

class Stiffness : public PhysicalModelInterface
{
 public:

  Stiffness() ;

  void rotate_to_calc_system(const Tensor2Gen& RotMatrix);
  
  //! stiffness tensor in calculation system (rank 4, double symmetric, 21 independent components)
  Tensor4DSym     C_calc  ;   
  
  //!  Returns subtensor \f$ B_{il} = C_{ijkl}, j = j_0, k = k_0 \f$
  Tensor2Gen get_subtensor(int j, int k);


  //!  Returns subtensor \f$ B_{kl} = C_{ijkl}, i = i_0, j = j_0 \f$ 
  Tensor2Sym get_another_subtensor(int i, int j); 

  //! sets stiffness tensor in crystal system [Gpa]  
  void set_C_tensor_crystal(const Tensor4DSym&     C);

  //! creates new object
  static Stiffness* create(const std::string& name,  const ModelOptions& options);
 

 private:

 
 
 protected:

  //! stiffness tensor in crystal system     (rank 4, double symmetric, 21 independent components) [GPa]
  Tensor4DSym     C_cr    ;    


  virtual void read_database ( ) {};


  virtual void do_init(void)=0;


  virtual void copy_from (const PhysicalModelInterface *rhs);


  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);


  virtual PhysicalModelInterface* create_new(void) const = 0;

 


};

//------------------------------------------------------------------------------------//
inline Stiffness* Stiffness::create( const std::string& name,  const ModelOptions& options )
{
  return dynamic_cast<Stiffness*>(PhysicalModelInterface::create("stiffness_" + name, options));
}


//------------------------------------------------------------------------------------//
inline  Tensor2Gen Stiffness::get_subtensor(int j, int k) 
{
  //from a tensor C_ijkl gives a tensor with j and k fixed
  Tensor2Gen temp ;
  const Tensor4DSym&  C1 = C_calc;

  for ( int i = 1; i <=3; i++ ) 
    for ( int l = 1; l <=3; l++ )
      {
	const double  p = C1(i,j,k,l) ; //unfortunatly, I have to work with constant object C1 due to tensor class restriction.
	                                //may be, I'll modify class in future. 	
	temp(i,l) = p; 
      }
  return(temp);
}

//---------------------------------------------------------------------------------------//

inline Tensor2Sym Stiffness::get_another_subtensor(const int i, const int j)
{
  //from a tensor C_ijkl gives a tensor with i and j fixed
  Tensor2Sym temp ;
  const Tensor4DSym&  C1 = C_calc;

  for ( int k = 1; k <=3; k++ ) 
    for ( int l = 1; l <=k; l++ )
      {

	const double  p = C1(i,j,k,l) ; //unfortunatly, I have to work with constant object C1 due to tensor class restriction.
	                                //may be, I'll modify class in future. 	
	temp(k,l) = p; 
      }
  return(temp);
}


#endif

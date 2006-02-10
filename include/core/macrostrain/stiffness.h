#include "tensor.h"
#include "xtensor.h"
#include <cmath>

using std::cout;

class stiffness
{
 public:
  stiffness();
  
  stiffness(double c11, double c12, double c44);
  stiffness(double c11, double c12, double c13, double c33, double c44 );
  
  void set_moduli(double c11, double c12, double c44);
  void set_moduli(double c11, double c12, double c13, double c33, double c44 );
  
  void print_crystall_system();
  

  void rotate_to_calc_system(const Tensor2Gen& RotMatrix);
 
  Tensor4DSym     C_calc  ;    // stiffness tensor in calculation system (rank 4, double symmetric, 21 independent components)
   
  Tensor2Gen get_subtensor(int j, int k);// Returns subtensor B(:,:) = C(:,j,k,:)

  Tensor2Sym get_another_subtensor(int i, int j); //Returns subtensor B(:,:) = C(i,j,:,:)


 private:

  Tensor4DSym     C_cr    ;    // stiffness tensor in crystal system     (rank 4, double symmetric, 21 independent components)
  
  

};

//---------------------------------------------------------------------
inline  Tensor2Gen stiffness::get_subtensor(int j, int k) 
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

//---------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------




inline Tensor2Sym stiffness::get_another_subtensor(const int i, const int j)
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



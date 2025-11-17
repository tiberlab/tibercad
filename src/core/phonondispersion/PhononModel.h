#ifndef _PHONONMODEL_H_
#define _PHONONMODEL_H_

#include "DynamicalMatrix.h"
#include "RamanTensor.h"
#include "PhysicalModel.h"
#include "libMeshDefs.h"

class RamanTensor;
       
class DynamicalMatrix;

//!Class that contains all the object, necessary for phonon dispersion solver
class PhononModel: public PhysicalModel
{
 public:

  //!Destructor
  ~PhononModel();
   
   //! creates a new object
  static  PhononModel* create(const ModelOptions& options);

  //!Get the thermal lattice conductivity
  // void get_dynamical_tensor(Tensor2& dynamical_tensor){};
  void fake(void);

   //! Init all fields
  void re_init(void){};

  //!Set the current element
  void set_element(const Elem* elem);
  
  void clear_dynamical_matrix_models();

  //!get the current element
   const Elem* get_element(void);
  //  void get_free_dynamical_matrix(std::vector<std::vector< double > >& D);
   void get_free_dynamical_matrix(Tensor2& dynamical_matrix);
   void get_full_dynamical_matrix(Tensor2& dynamical_matrix);

   void get_raman_tensor(std::vector<Tensor2> & raman_tensor);

   void get_light_polarization(std::vector<Tensor1>& light_polarization);

  //void get_full_dynamical_matrix(std::vector<std::vector< double > >& D);
    //! \copydoc PhysicalModel::do_print_info(void)
   virtual  void do_print_info(void){};

  
 RamanTensor* raman_tensor_model;


 private:
  
  //LightDirections
  std::vector<Tensor1> _light_polarization;


  ID free_ID;

  ID add_dynamical_matrix_model(const std::string& model_name, 
				  const ModelOptions& options);

  std::map<ID,DynamicalMatrix* > _dynamical_matrix_models;

         //!Iterator for heat source model
  typedef std::map<ID, DynamicalMatrix*>::iterator dyn_mat_iterator; 



    
   //! Current element 
   const Elem* _elem; 
  
  //   model_options model_opt;

   //!Update lattice thermal conductivity
  void  update_dynamical_matrix(void){}; 
  
  
   //! Lattice thermal conductivity
  // Tensor4DSym _dynamical_tensor;  
 
   //! Lattice thermal conductivity model
   //DynamicalTensor* Dt;


   
 protected:

  //!Constructor
  PhononModel(const ModelOptions& options);

  virtual PhysicalModel* create_new (void) const;

  virtual void do_init();

};

inline
PhononModel* PhononModel::create(const ModelOptions& options)
{
  std::cout<<"OK"<<std::endl;
  return new  PhononModel(options);
}


inline
PhysicalModel*
PhononModel::create_new(void) const
{
  return new  PhononModel(get_options());
}



inline
void 
PhononModel::set_element(const Elem* elem)
{

 _elem = elem;

}

inline
const Elem*
PhononModel::get_element(void)
{

  return _elem;

}




#endif

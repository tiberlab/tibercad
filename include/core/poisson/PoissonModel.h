#ifndef _POISSONMODEL_H_
#define _POISSONMODEL_H_


#include "PhysicalModel.h"
#include "elem.h"
#include "SimulationInterface.h"
#include "ChargeDensityModel.h"
//#include "DielectricModel.h"
#include "Macrostrain.h"
#include "SimulationEnvironment.h"
#include "StrainInterface.h"
#include "OptDielectricConstant.h"
#include "PiezoelectricModel.h"
#include "tensor_value.h"

//class PiezoelectricModel;

//!Class that contains all the physical quantities necessary for the POISSON solver
class PoissonModel: public PhysicalModel
{
 public:

  //!Constructor
  PoissonModel(const ModelOptions& options);

  //!Destructor
  ~PoissonModel(void);

  //! creates a new object
  static  PoissonModel* create(const ModelOptions& options);

  void 	re_init(void);

  //!Set the current element
  void set_element(const Elem* elem);

  // //!Return the charge density for the current element
  // double get_charge_density();

   //!Return the charge density for the current element
  std::vector<double> get_charge_density();

  //!Return the charge density for the current element
  void get_dielectric_constant(RealTensor& epsilon);


  //Tensor1 get_built_in_polarization(const std::vector<Point> q_point, std::vector<Tensor1>& built_in_polarization);
  //!charge density (electron/cm^3)
  void get_charge_density(const std::vector<Point> q_point, std::vector<double>& charge_density);



  void  get_total_polarization(std::vector<RealGradient>& pol,const std::vector< Point >& points);
  void  get_pyro_polarization(std::vector<RealGradient>& pol,const std::vector< Point >& points);
  void  get_piezo_polarization(std::vector<RealGradient>& pol,const std::vector< Point >& points);

  //! Get the strain
    const Tensor2Sym& get_strain(void) const;

 PiezoelectricModel* _piezo_model;


 private:

 //----------------


  //!Pointer to drift diffusion simulation
  SimulationInterface* _simul;

  Point _coord;

  OptDielectricConstant* _epsilon_model;

  Tensor2Sym _strain;

    //! The interface to a strain simulation
    StrainInterface _strain_if;

  //  enum strain_variables
  //{
  //   E_XX = 0,
  //   E_XY,
  //   E_XZ,
  //   E_YY,
  //   E_YZ,
  //   E_ZZ
  // };

   //!Strain variables
   //std::set<ID> pol_ID;

   //!Variable map
   //std::map<ID,ID> var_map;

   //!A pointer to cherge density simulation
  SimulationInterface* _chd_sim;

  //!A pointer to cherge piezoelectricity simulation
  Macrostrain* _strain_sim;

  //! ID for charge density simulation
  ID charge_id;

  //Model options structure
  struct model_options
  {
     bool pyro_pol;

     bool piezo_pol;

     std::string  chd_sim_name;

    bool add_doping;

   };

  //!Options for Poisson model
  model_options model_opt;

   //! update the charge density
   void  update_charge_density(void);

  //! Update the built in polarization
  void  update_built_in_polarization(void);

   //! Current element
   const Elem* _elem;

   //!relative dielectric constant
  Tensor2Sym _epsilon;

  // //!charge density (electron/cm^3)
   //std::vector<double> _charge_density;


  // //!A pointer to cherge density object
  // ChargeDensityModel* chd_model;

   //!A pointer to dielectric constant model
  //   DielectricModel* dielectric_model;

  //!Pyropolarization
  //  Tensor1  _pyropolarization;

  //!Piezopolarization
  //Tensor1  _piezopolarization;

 protected:

  virtual PhysicalModelInterface* create_new (void) const;


  virtual void read_database (void){};


  virtual void read_database_alloy (void) {};


  virtual void do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);

  virtual void do_init();

  virtual void create_submodels(void);

};


inline
PoissonModel*
PoissonModel::create(const ModelOptions& options)
{
  return new PoissonModel(options);
}

inline
void
PoissonModel::set_element(const Elem* elem)
{

  _elem = elem;

}


// inline
// double
// PoissonModel::get_charge_density()
// {

//  return _charge_density;

// }


#endif

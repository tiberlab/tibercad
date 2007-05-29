#ifndef _HEATMODEL_H_
#define _HEATMODEL_H_


#include "PhysicalModel.h"
#include "LatticeThermalConductivity.h"
#include "ThermoelectricPower.h"
#include "DriftDiffusion.h"
#include "elem.h"




//!Class that contains all the object, necessary for Heat Transport solver
class HeatModel: public PhysicalModel
{
 public:
  //!Constructor
  HeatModel();

  //!Destructor
  ~HeatModel();
   
   //! creates a new object
  static HeatModel* create();
  
    //! Pointer to a DriftDiffusion simulation
 
  bool get_peltier_thomson_opt();

  //!Get the thermal lattice conductivity
  void get_lattice_thermal_conductivity(Tensor2Sym& lattice_thermal_conductivity);
  
  //! Init all fields
  void re_init();

  void get_dd_solution(std::vector<Point> g_point,
                                std::vector<DriftDiffusion::Solution>& potentials,
				std::vector<DriftDiffusion::Currents>& currents);

   //!Set the current element
   void set_element(const Elem* elem);

   //!Return the electrons thermoelectric_power
   double get_electrons_thermoelectric_power();

   //!Return the holes thermoelectric_power
   double get_holes_thermoelectric_power();

   //!Get 
    bool get_joule_opt();
   
    DriftDiffusion* get_dd_simulation();

   //! Set the temperature 
   void set_temperature(double temperature); 

   //! Get the Simulation Environment of the DD simulation
   SimulationEnvironment& get_dd_environment();


 private:

     struct model_options
   {
       bool joule_effect;

       bool peltier_thomson_effect;

   };
    
   //! Current element 
   const Elem* _elem; 

   DriftDiffusion* _dd_simul;

   model_options model_opt;

   void  update_electrons_thermoelectric_powers();
   
   void  update_holes_thermoelectric_powers();

   double _eTEpower;

   double _hTEpower;
 
   void  update_lattice_thermal_conductivity(void); 
  
   double _temperature;

   Tensor2Sym _lattice_thermal_conductivity; 
  
   // options model_opt;

   LatticeThermalConductivity* kappa;




   
   //!copy constructor should not be used
   HeatModel (const HeatModel &  t) {};
  
 protected:

  virtual PhysicalModelInterface* create_new (void) const;


  virtual void copy_from(const PhysicalModelInterface *rhs){};


  virtual void read_database (void){};
 

  virtual void read_bowing_parameters (void) {};


  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);

  virtual void do_init();

};



inline 
HeatModel* HeatModel::create()
{
  return new  HeatModel();
}



inline 
void
HeatModel::get_lattice_thermal_conductivity(Tensor2Sym& lattice_thermal_conductivity)
{
  lattice_thermal_conductivity = _lattice_thermal_conductivity;
}



inline
void 
HeatModel::set_temperature(double temperature)
{
  _temperature = temperature;
}

inline
double
HeatModel::get_electrons_thermoelectric_power()
{

  return _eTEpower;

}

inline   
double
HeatModel::get_holes_thermoelectric_power()
{

   return _hTEpower; 

}

inline
void 
HeatModel::set_element(const Elem* elem)
{

  _elem = elem;

}

inline
bool
HeatModel::get_joule_opt()
{

  return model_opt.joule_effect;

}

inline
SimulationEnvironment&
HeatModel::get_dd_environment()
{

  return _dd_simul->get_environment();

}

inline
bool
HeatModel::get_peltier_thomson_opt()
{
return model_opt.peltier_thomson_effect;

}

inline
DriftDiffusion* HeatModel::get_dd_simulation()
{
  return _dd_simul;

}

#endif

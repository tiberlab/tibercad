#ifndef _HEATMODEL_H_
#define _HEATMODEL_H_


#include "PhysicalModel.h"
#include "LatticeThermalConductivity.h"
#include "ThermoelectricPower.h"
#include "DriftDiffusion.h"
#include "elem.h"
#include "ParticleThermalConductivity.h"
#include <point.h>
#include "Constants.h"

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
  void get_thermal_conductivity(Tensor2Sym& thermal_conductivity);
  
  //! Init all fields
  void re_init();


  void get_dd_solution_secure( std::vector<Point> g_point,
			       std::vector<double>& QfermiE,
			       std::vector<double>& QfermiH,
			       std::vector<Point>& JE,
			       std::vector<Point>& JH);
  
  void get_ex_solution_secure( std::vector<Point> g_point,
			       std::vector<double>& ex_potential,
			       std::vector<Point>& J_ex,
                               std::vector<double>& radiative_power);
  


   //!Set the current element
   void set_element(const Elem* elem);

   //!Return the electrons thermoelectric_power
   double get_electrons_thermoelectric_power();

   //!Return the holes thermoelectric_power
   double get_holes_thermoelectric_power();

   //!get the information about the drift diffusion simulation
    bool get_joule_opt();

   //!get the information about the excitons simulation
    bool get_excitons_opt();
 
   //! Set the temperature 
   void set_temperature(double temperature); 

   //! Get the Simulation Environment of the DD simulation
   SimulationEnvironment& get_dd_environment();

 //! Get the Simulation Environment of the Excitons simulation
   SimulationEnvironment& get_ex_environment();

 private:
 //! The variables that can be provided
    enum dd_var
    {
      QFERMIE,
      QFERMIH,
      JNX,
      JNY,
      JNZ,
      JPX,
      JPY,
      JPZ
     
    };

    enum dd_var_kpart
    {
      CONDE,
      CONDH
    };

  enum dd_var_TE
    {
      TEPOWERE,
      TEPOWERH
      
    };

    enum ex_var
    {
      JEX_X,
      JEX_Y,
      JEX_Z,
      EX_POTENTIAL,
      RADPOWER
      
    };
    struct model_options
   {
     bool joule_effect;

     bool peltier_thomson_effect;
  
     bool particle_thermal_conductivity;
 
     bool excitons;

   };
    
   //! Current element 
   const Elem* _elem; 


   //! drift diffusion simulation
   SimulationInterface* _dd_simul;

   //! excitons simulation
   SimulationInterface* _ex_simul;

   // For general solution

   std::set< ID > dd_ID_je;

  std::vector<ID> ID_je;

  //!For excitons solution
   std::set< ID > ex_set_ID;

   std::vector<ID> ID_ex;



   //!For particle solution
   std::set< ID >  dd_ID_kpart;

  std::vector<ID> ID_kpart;
  
  //For thermoelectric power

  std::set< ID > dd_ID_TEpower;
  
  std::vector< ID > ID_TEpower;
 


   model_options model_opt;

   void  update_thermoelectric_powers();
   
   void  update_particle_thermal_conductivity();

   

   double _eTEpower;

   double _hTEpower;
 
   void  update_lattice_thermal_conductivity(void); 
  
   double _temperature;

   Tensor2Sym _lattice_thermal_conductivity; 

   Tensor2Sym _electrons_thermal_conductivity; 
   
   Tensor2Sym _holes_thermal_conductivity; 
  
   // options model_opt;

   LatticeThermalConductivity* kappa;

   ParticleThermalConductivity* kappa_carrier;




   
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
HeatModel::get_thermal_conductivity(Tensor2Sym& thermal_conductivity)
{
  thermal_conductivity = _lattice_thermal_conductivity +
    +  _electrons_thermal_conductivity
    +  _holes_thermal_conductivity;
  
    
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
bool
HeatModel::get_excitons_opt()
{

  return model_opt.excitons;

}

inline
SimulationEnvironment&
HeatModel::get_ex_environment()
{

  return _ex_simul->get_environment();

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


#endif

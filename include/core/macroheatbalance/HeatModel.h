#ifndef _HEATMODEL_H_
#define _HEATMODEL_H_


#include "PhysicalModel.h"
#include "LatticeThermalConductivity.h"
#include "SimulationInterface.h"
#include "elem.h"
#include "ParticleThermalConductivity.h"
#include "point.h"
#include "Constants.h"
#include "HeatSourceInterface.h"
#include "ThermalConductivityInterface.h"       

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

  //!Get the thermal lattice conductivity
  void get_thermal_conductivity(Tensor2Sym& thermal_conductivity);
  
  //! Init all fields
  void re_init();

   //!Set the current element
   void set_element(const Elem* elem);

   //! Set the temperature 
   void set_temperature(double temperature); 

  
  //!Get total the heat source model given an ID
  HeatSourceInterface* get_heat_source_model(ID id);
  
  
  //! Get the ids of the heat source models 
  int get_heat_source_IDs(std::vector<ID>& ids);
  

  //!Get total heat source
  void get_total_heat_source(std::vector<Point> h_point,
			     std::vector<double>& total_heat_source);

  //!Get total power flux
  void get_total_power_flux(std::vector<Point> h_point,
			    std::vector<RealGradient>& total_power_flux,
			    bool check_boundary);

 private:
 

    enum dd_var_kpart
    {
      CONDE,
      CONDH
    };


    struct model_options
   {
  
     bool particle_thermal_conductivity;

   };
    
   //! Current element 
   const Elem* _elem; 

  


   //!For particle solution
   std::set< ID >  dd_ID_kpart;

   std::vector<ID> ID_kpart;

   model_options model_opt;

  

   //!Iterator for thermal conductivity model
   typedef std::map<ID, ThermalConductivityInterface*>::iterator outer_conductivity_iterator;

   //!Iterator for heat source model
   typedef std::map<ID, HeatSourceInterface*>::iterator outer_source_iterator;

   //!Heat sources iterator within a specific heat source model
   typedef std::map<ID,double>::iterator inner_source_iterator;

   //!Thermal conductivity iterator within a specific thermal conductivity model
   typedef std::map<ID,Tensor2Sym>::iterator inner_conductivity_iterator;

   //!Heat flux source within a specific heat source model
   typedef std::map<ID,RealGradient>::iterator inner_flux_source_iterator;
 
   //!Update particle thermal conductivity
   void  update_particle_thermal_conductivity();

   //!Update lattice thermal conductivity
   void  update_lattice_thermal_conductivity(void); 
  
   //! Temperature 
   double _temperature;

   //! Lattice thermal conductivity
   Tensor2Sym _lattice_thermal_conductivity; 

   //! Electron thermal conductivity
   Tensor2Sym _electrons_thermal_conductivity; 
   
   //! Hole thermal conductivity
   Tensor2Sym _holes_thermal_conductivity; 
 
   //! Lattice thermal conductivity model
   LatticeThermalConductivity* kappa;

   //! Particle thermal conductivity model
   ParticleThermalConductivity* kappa_carrier;

   //! Pointer to a heat source model
   HeatSourceInterface* _heat_source_interface;
 
   //!Thermal conductivity model map
   std::map<ID, ThermalConductivityInterface*> _thermal_conductivity_models;

   //!Heat Source model map
   std::map<ID, HeatSourceInterface*> _heat_source_models;
 

   //! Add a heat source model 
   void add_heat_source_model(const std::string& model_name,
       			const ModelOptions& options = ModelOptions());

   //! Add a thermal conductivity model
    void add_thermal_conductivity_model(const std::string& model_name, 
  			       const ModelOptions& options);

   //! Clear all heat source models
   void clear_heat_sources(void);

    //! Clear thermal conductivity models
   void clear_thermal_conductivity(void);

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
  thermal_conductivity = _lattice_thermal_conductivity;
    // +  _electrons_thermal_conductivity
    // +  _holes_thermal_conductivity;
  
    
}



inline
void 
HeatModel::set_temperature(double temperature)
{
  _temperature = temperature;
}


inline
void 
HeatModel::set_element(const Elem* elem)
{

  _elem = elem;

}


inline
HeatSourceInterface*
HeatModel::get_heat_source_model(ID id)
{

HeatSourceInterface* heat_source_model = NULL;
outer_source_iterator it =  _heat_source_models.find(id);
if (it !=   _heat_source_models.end())
   heat_source_model = it->second;
  
   return heat_source_model;

   
 }



#endif

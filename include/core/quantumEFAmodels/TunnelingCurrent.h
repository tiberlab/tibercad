#ifndef _TUNNELINGCURRENT_H_ 
#define _TUNNELINGCURRENT_H_

#include "KspaceIntegration.h"



//! This class Calculates Tunneling current density, performing integration of the transmission  in k-space
class TunnelingCurrent: public KspaceIntegration
{

 public:


  //! Consructor
  TunnelingCurrent();

  //! Destructor
  ~TunnelingCurrent();


  //!creates a new object 
  static  TunnelingCurrent* create();


  struct options
  {
    double Efermi_left;
    double Efermi_right;
    double voltage_start;
    double voltage_stop; 
    unsigned int voltage_steps;

    unsigned int init_nodes_for_energy_int;
    double energy_int_tolerance;
    bool energy_int_refinement;
    bool energy_int_uniform_refinement;
    double  energy_int_refine_fraction;
    double  energy_int_zero_limit;
  };

 protected:

   //!calculates transmission (integrated on energy) at each  k_point and for a range of applied voltages 
   //virtual void calculate_at_each_k_point();

 
   void build_V_grid();

   //! Applied voltage mesh
   Mesh* Vmesh;


   virtual void do_plot (void);

   virtual void build_elemental_results(const std::set<std::string>& variables,
				   std::vector<double>& results, std::vector<std::string>& legend);


   virtual void do_solve(void);

   virtual void do_init(void);

   virtual void parse_options(void);


   virtual void calculate_density(void);

   //! we do not kelly here
   virtual void estimate_error_for_refinement(ErrorVector& error);

 private:

   //!used for all k-points
   const Elem* applied_voltage_elem;


   void k_space_output(void);

   //!calculates everything that is necessary at a single k-point
   double calculte_at_k_point(const Point& k);


   std::map<const Elem*, double> kspace_integral;


   std::map<const Elem*, double> volume;


   unsigned int how_many_elements_to_do();


   void calculate_volumes(void);


   //!performes integration over energy. May refine mesh, if required
   double integrate_over_energy(double k[3], double electric_potential);


   options opt;

   //!performes integration over energy on a fixed mesh
   double integrate_over_fix_energy(const Mesh* energy_mesh, double k[3], double electric_potential );


   //!fermi distribution factor
   /*!
     \param fermi_energy Fermi energy [eV]
     \param Energy [eV]
   */
   double thermal_probability(double fermi_energy, double Energy);



   std::map<const Elem*, double> energy_integral;


   void estimate_error_for_energy_refinement(ErrorVector& error);

   Mesh* Emesh;


};

//---------------------------------------------------------

inline TunnelingCurrent*  TunnelingCurrent::create()
{
  return (new TunnelingCurrent );
}



inline double TunnelingCurrent::thermal_probability(double Fermi_energy, double Energy)
{
  

  double T_EV = SimulationOptions::temperature * Constants::k_Boltzmann;
  double exp_arg =  (Energy - Fermi_energy)/T_EV;
  
  double el_fermi;

  if (exp_arg > 100 ) 
    el_fermi = 0.0;
  else
    el_fermi = 1.0/(  1 +  std::exp(exp_arg)  );


  return(el_fermi);
}



#endif

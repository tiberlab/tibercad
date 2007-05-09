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
   

 protected:

   //!calculates transmission (integrated on energy) at each  k_point and for a range of applied voltages 
   virtual void calculate_at_each_k_point();

 
   void build_V_grid();

   //! Applied voltage mesh
   Mesh* Vmesh;


   virtual void do_plot (void);

   virtual void build_elemental_results(const std::set<std::string>& variables,
				   std::vector<double>& results, std::vector<std::string>& legend);


   std::map<const Elem*, double> transmission_map;


};

//---------------------------------------------------------

inline TunnelingCurrent*  TunnelingCurrent::create()
{
  return (new TunnelingCurrent );
}


#endif

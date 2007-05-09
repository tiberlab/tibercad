#ifndef _TUNNELINGCURRENT_H_ 
#define _TUNNELINGCURRENT_H_

//! A class that calculates Tunneling current density, performing integration of the transmission 
// in k-space

// Basic include files needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "gmv_io.h"
#include "linear_implicit_system.h"
#include "nonlinear_implicit_system.h"
#include "equation_systems.h"

#include "getpot.h"
// For mesh refinement


#include "mesh_refinement.h"
#include "error_vector.h"
#include "kelly_error_estimator.h"

// Define the Finite Element object.
#include "fe.h"
#include "elem.h"
// Define Gauss quadrature rules.
#include "quadrature_gauss.h" 

// Define useful datatypes for finite element
// matrix and vector components.
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "dof_map.h"

#include "fe_interface.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <complex>


#include <fstream>
#include <iomanip>



#include <complex>
#include <vector>

#include <petsc_matrix.h>
#include "EFAbulkHamiltonian.h"
#include <algorithm>
#include <set>
#include <tecplot_io.h>
#include "mesh_data.h"
#include "DriftDiffusion.h"
#include "EnvelopFunctionApprox.h"
#include "SimulationInterface.h"
#include "tensor.h"

//! This class Calculates Tunneling current density, performing integration of the transmission  in k-space
class TunnelingCurrent  : public QuantumDensity
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

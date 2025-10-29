#ifndef _KMCINTERFACE_H_
#define _KMCINTERFACE_H_

#include "SimulationInterface.h"
#include <tr1/random>
#include "mesh.h"
#include "elem.h"
#include "TiberLinearSystem.h"
#include "Boundary.h"
#include "PotentialInterface.h"

#include "meshfree_interpolation.h"

//using namespace std;

class TBDLLOCAL KmcInterface : public SimulationInterface
{
 public:

  enum Solution
  {
    elDensity,
    hlDensity,
    xDensity,
    eCurrentDensity,
    hCurrentDensity,
    Generation,
    Recombination,
    Potential
  };
   
  virtual ~KmcInterface() {};
   
 protected:
 
  KmcInterface(const ModelOptions& opt) : SimulationInterface(opt){};
 
  virtual void do_init();
 
  //virtual void do_solve(){};
 
  //virtual void parse_options(){};
 
  virtual bool get_boundary_potentials(Boundary *bd, const Point& pt, double& phi, double& efermi, double& hfermi);
  
  SimulationInterface::SolutionProvider _pot_sol;
  SimulationInterface::SolutionProvider _mue_sol;
  SimulationInterface::SolutionProvider _muh_sol;

  //! this is used for storage and interpolation of data between grid and mesh 
  //libMesh::MeshfreeInterpolation* _griddata;

  std::vector<double> _eldensity;
  std::vector<double> _hldensity;
  std::vector<double> _xdensity;
  std::vector<double> _potential;

  //! Monte Carlo variables
  unsigned int _Nstep;

  unsigned int _IOstep;

};




#endif

#ifndef _STRAINSIMULATION_H_
#define _STRAINSIMULATION_H_
 
#include "SimulationInterface.h"
#include "tensor.h"

class  Device;

class StrainSimulation:  public SimulationInterface
{
 public:

  enum Variables
  {
    EPS_XX = 0,
    EPS_YY = 1,
    EPS_ZZ = 2,
    EPS_XY = 3,
    EPS_XZ = 4,
    EPS_YZ = 5,
    P_X = 6,
    P_Y = 7,
    P_Z = 8
  };



  /*! \copydoc SimulationInterface::create_physical_model() */
  virtual PhysicalModel*   create_physical_model(const ModelOptions& options,
			   const Material* mat) const  throw (ModelErrorException);



  virtual ID convert_variable_name_to_id(const std::string& variable_name) const;


  virtual void get_solution_secure(const Elem* elem,
				   const std::vector<Point>& p, const std::set<ID>& ids,
				   std::vector<std::map<ID, double> >& values);
  


 protected:
  
  StrainSimulation(const ModelOptions& options);

  
  //!map betwen the element and the result strain in crystal system
  /*!
    the map is created at the end of the method solve()
    it contains only the active elements for this strain simulation
   */
  std::map<const Elem*, Tensor2Sym> result_strain;
  
  static Device*   _device;

  virtual void do_init(void);

  //virtual void parse_options(void) ;
 
  //virtual void do_solve(void);

 private:


};

 
inline
StrainSimulation::StrainSimulation(const ModelOptions& options)
 : SimulationInterface(options)
 {
 }


#endif

#ifndef _CRACKSTRAIN_H_
#define _CRAINSTRAIN_H_
#include "StrainSimulation.h"
#include "Device.h"

class CrackStrain : public StrainSimulation
{
 public:

 
  // virtual void get_solution_secure(const Elem* elem,
  //				  const std::vector<Point>& p, const std::set<ID>& ids,
  //				   std::vector<std::map<ID, double> >& values) {} ;


  static CrackStrain* create(const ModelOptions& options);

 protected:

  CrackStrain(const ModelOptions& options);

  virtual void parse_options(void);
  
  virtual void do_init(void) ;
 
  virtual void do_solve(void);


 private:


  

  double _Ki;

  double _x0;

  double _y0;

  double _sigma_ys;

  //!calculate stress analitically
  void calculate_stress(Tensor2Sym& stress, const double x, const double y) const;


  /*! 
    \copydoc SimulationInterface::build_elemental_results()
    The variables are: "strain", "polarization"
    This means strain tensor components:
    \f$ \varepsilon_{xx}, \varepsilon_{yy},\varepsilon_{zz}, \varepsilon_{xy}, \varepsilon_{xz}, \varepsilon_{yz}\f$,
    and polarization vector \f$ \bf P \f$ components. 
    \f$ x, y,  z\f$ refer to calculation coordinate system.
   */
  virtual void build_elemental_results(const std::set<std::string>& variables,
				       std::vector<double>& results, std::vector<std::string>& legend) ;


  //! Preapare all 6 components of the strain tensor for output
  void prepare_strain_data_for_output( std::vector<std::string>& eps_names, std::vector<double>& eps_data );


  //! Preapare all 3 components of the polarization vector for output
  void prepare_polarization_data_for_output( std::vector<std::string>& polariz_names, std::vector<double>& polariz_data );

  //! Prepare stress tensor for output
  void prepare_stress_data_for_output(std::vector<std::string>& stress_names, std::vector<double>& stress_data);


};




inline
CrackStrain::CrackStrain(const ModelOptions& options)
 : StrainSimulation(options)
{
}



inline
CrackStrain* CrackStrain::create(const ModelOptions& options)
{
  return new CrackStrain(options);
}

#endif

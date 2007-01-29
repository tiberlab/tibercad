#ifndef _DDSEMICONDUCTOR_H_
#define _DDSEMICONDUCTOR_H_

//!  A class to provide all neccessary parameters for drift-diffusion calculation.
/*!
     The class can calculate information about the band structure, such as
     band edge energy, effective mass for the density of states calculation and
     degeneracy 
*/

#include "tensor.h"
#include <vector>
#include "PhysicalModelInterface.h"

#include "KPbulkHamiltonian.h"

class DDsemiconductor : public PhysicalModelInterface
{
 public:

  //! data structure for band extremum
  struct band_extremum
   {
     short degeneracy;/*!< degeneracy (including spin).*/
     double energy;   /*!< edge energy [eV].*/
     double mass_DOS; /*!< mass for density of states \f$ m = \left( \mathop {\rm det} \left(\frac{1}{m} \right)_{ij} \right)^{-1} \f$ */
  };


  void set_strain(const Tensor2Sym& strain);


  //! energy cut-off value
  /*!
    Conduction band: 
    \f$ E_i - E_i^{min} < E_{cut-off}    \f$

    Valence band:
    \f$ E_i^{max} - E_i > E_{cut-off}    \f$,
    where \f$ E_i \f$ are the band edge energies.
   */
  double energy_cutoff;
  

  //!Constructor
  DDsemiconductor(void);

  //!Desctructor
  virtual ~DDsemiconductor(void);


  //! returns information about conduction bands
  const std::vector<band_extremum>& get_conduction_band_energy_mass(void) const;
  
  //! returns information about valence bands
  const std::vector<band_extremum>& get_valence_band_energy_mass(void) const;



  //! calculate information about conduction bands
  virtual void  calculate_conduction_band_extremum(void) = 0;

  //! calculate information about valence bands
  virtual void  calculate_valence_band_extremum(void) = 0;


  //! calculates dispersion along a line in k-space
  /*!
    \param  k_i - initial k-point
    \param  k_f - final k-point
    \param  number_of_points number k-points between the initial and the final point 
   */
  std::vector<std::vector<double> >  get_valence_kp_dispersion(Tensor1 k_i, Tensor1 k_f, unsigned int number_of_points);

  std::vector< std::vector<double> > calculate_vb_bulk_states(const std::vector<Tensor1>& k_vector) ; 


   //! absolute value of k-vector for DOS mass calculation [a.u.]
  double k_max;

 
  
  //! creates new object
  static DDsemiconductor* create(const std::string& name,  const ModelOptions& options);
 

 private:

  KPbulkHamiltonian*  bulk_ham;

  
   

 

 protected:


  //! Hartree energy in eV
  static const double Hartree;

  //!semiconductor
  Semiconductor* semiconductor;

  //!  if \f$ ||\varepsilon_{ij}|| > 10^{-5} \f$, then true 
  bool strained ;

  //! strain tensor in crystal system
  Tensor2Sym   strain;


  //! information about conduction bands
  std::vector<band_extremum>  conduction_band;

  //! information about valence bands
  std::vector<band_extremum>  valence_band;


  virtual PhysicalModelInterface* create_new(void) const = 0;

  virtual void copy_from (const PhysicalModelInterface *rhs) ;

  virtual void do_init (void);

  virtual void read_database(void);

  virtual void read_bowing_parameters(void) {};
 
  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);

};




inline DDsemiconductor* DDsemiconductor::create(const std::string& name,
    const ModelOptions& options)
{
  return dynamic_cast<DDsemiconductor*>(
      PhysicalModelInterface::create("DDsemicond_" + name, options));
}
 




#endif

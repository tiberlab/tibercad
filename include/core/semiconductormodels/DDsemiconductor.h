#ifndef _DDSEMICONDUCTOR_H_
#define _DDSEMICONDUCTOR_H_

//!  A class to provide all neccessary parameters for drift-diffusion calculation.
/*!
     The class can calculate information about the band structure, such as
     band edge energy, effective mass for the density of states calculation and
     degeneracy 
*/

#include "PhysicalModelInterface.h"
#include "libMeshDefs.h"
#include "KPbulkHamiltonian.h"

#include "tensor.h"
#include "dense_vector.h"

#include <vector>

class TBDLEXPORT DDsemiconductor : public PhysicalModelInterface
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

  //!sets temperature for semiconductor object
  void set_temperature(const double T);

  //! sets temperature for semiconductor object based on coordinate
  void set_temperature(const Elem* element, const Point& point);


  //! energy cut-off value
  /*!
    Conduction band: 
    \f$ E_i - E_i^{min} < E_{cut-off}    \f$

    Valence band:
    \f$ E_i^{max} - E_i > E_{cut-off}    \f$,
    where \f$ E_i \f$ are the band edge energies.
   */
  double energy_cutoff;
  

  //!Desctructor
  virtual ~DDsemiconductor(void);


  //! returns information about conduction bands
  const std::vector<band_extremum>& get_conduction_band_energy_mass(void) const;
  
  //! returns information about valence bands
  const std::vector<band_extremum>& get_valence_band_energy_mass(void) const;



  //! calculate information about conduction bands
  void  calculate_conduction_band_extremum(void);

  //! calculate information about valence bands
  void  calculate_valence_band_extremum(void);



  //! calculate information about conduction bands
  void  calculate_conduction_band_extremum(const Elem* element, const Point& point);

  //! calculate information about valence bands
  void  calculate_valence_band_extremum(const Elem* element, const Point& point);

  //! Calculate valence band states
  void calculate_vb_bulk_states(const Tensor1& k_vector, std::vector<double>& eigenvalues, std::vector<libMesh::DenseVector<libMesh::Complex>>& eigenvectors);


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
  static DDsemiconductor* create(const Material* mat,  const ModelOptions& options);


 
 

  private:

  //! a pointer to an object that can calculate the valence band properties 
  KPbulkHamiltonian*  bulk_ham;
  

  //! Compute the inverse mass tensor for a given set of degenerate states
  /*!
   * This currently works only for hamiltonians without terms linear in k,
   * and only in standard directions in case of warping.
   */
  void calculate_inverse_mass(
      const std::vector<std::vector<KPbulkHamiltonian::MatrixElement> >& ham,
      const std::vector<libMesh::DenseVector<libMesh::Complex>>::const_iterator first,
      const std::vector<libMesh::DenseVector<libMesh::Complex>>::const_iterator last,
      std::map<ID, libMesh::RealTensor>& imasses);

   

 

 protected:

  //!Constructor
  DDsemiconductor(const ModelOptions& options);


  //!pointer to an object that can calculate the conduction band properties 
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

  virtual void do_init (void);

  virtual void prepare_submodels(void);

  virtual void read_database(void);


  //! calculate information about conduction bands
  virtual void  do_calculate_conduction_band_extremum(void) = 0;

  //! calculate information about valence bands
  virtual void  do_calculate_valence_band_extremum(void);

};





inline 
void  DDsemiconductor::calculate_conduction_band_extremum(const Elem* element, const Point& point)
{


  semiconductor->set_temperature( element, point);

  calculate_conduction_band_extremum();

}

inline
void  DDsemiconductor::calculate_valence_band_extremum(const Elem* element, const Point& point)
{


  semiconductor->set_temperature( element, point);
 

  calculate_valence_band_extremum();

}


inline void DDsemiconductor::calculate_conduction_band_extremum()
{
  do_calculate_conduction_band_extremum();
}


inline void  DDsemiconductor::calculate_valence_band_extremum()
{
  do_calculate_valence_band_extremum();
}


inline void DDsemiconductor::set_temperature(const double T)
{
  semiconductor->set_temperature(T);
 
}

inline void DDsemiconductor::set_temperature(const Elem* element, const Point& point)
{
  semiconductor->set_temperature( element, point);
}


inline const std::vector<DDsemiconductor::band_extremum>& DDsemiconductor::get_conduction_band_energy_mass(void) const
{
   return(conduction_band);
} 

inline  const std::vector<DDsemiconductor::band_extremum>& DDsemiconductor::get_valence_band_energy_mass(void) const
{
  return(valence_band);
}


#endif

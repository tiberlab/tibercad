/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file DDsemiconductor.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */

#ifndef TC_DDSEMICONDUCTOR_H
#define TC_DDSEMICONDUCTOR_H


#include "tibercad/physics/PhysicalModel.h"
#include "tibercad/base/libMeshDefs.h"
#include "tibercad/math/Tensor2.h"
#include "tibercad/physics/semiconductormodels/KPbulkHamiltonian.h"

#include "libmesh/dense_vector.h"

#include <vector>

//!  A class to provide all neccessary parameters for drift-diffusion calculation.
/*!
     The class can calculate information about the band structure, such as
     band edge energy, effective mass for the density of states calculation and
     degeneracy 
*/
class DDsemiconductor : public PhysicalModel
{
  public:



  //! data structure for band extremum
  struct band_extremum
   {
     short degeneracy = 1;/*!< degeneracy (including spin).*/
     double energy = 0;   /*!< edge energy [eV].*/
     double mass_DOS = 1; /*!< mass for density of states \f$ m = \left( \mathop {\rm det} \left(\frac{1}{m} \right)_{ij} \right)^{-1} \f$ */
  };


  void set_strain(const Tensor2& strain);

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
  std::vector<std::vector<double> >  get_valence_kp_dispersion(const Tensor1& k_i, const Tensor1& k_f, unsigned int number_of_points);

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
  Tensor2   strain;


  //! information about conduction bands
  std::vector<band_extremum>  conduction_band;

  //! information about valence bands
  std::vector<band_extremum>  valence_band;


  virtual PhysicalModel* create_new(void) const = 0;

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

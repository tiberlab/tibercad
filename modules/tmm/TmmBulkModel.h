/*  
 * This file is part of the tiberCAD module tmm.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file TmmBulkModel.h
 * \brief tiberCAD tmm module header.
 *
 * \note This file is part of module tmm.
 */


#ifndef TC_TMMBULKMODEL_H
#define TC_TMMBULKMODEL_H

#include "tibercad/physics/PhysicalModel.h"
#include "TmmDipoleSource.h"


//! This is the base class for the TMM bulk physical model
class TC_DLEXPORT TmmBulkModel : public PhysicalModel
{

public:


  //! Destructor
  virtual ~TmmBulkModel(void) {};

  //! Creator function
  static TmmBulkModel* create(const Material* mat, const ModelOptions& options);

  void  calculate(const libMesh::Elem* elem, const libMesh::Point& point, double lambda);


  //! Get the relative permittivity at a given wavelength
  libMesh::Complex get_permittivity(double lambda) const;


  //! Get the complex refractive index
  libMesh::Complex get_refractive_index(double lambda) const;

  double get_coherent_index(void) const;


  double get_emission_power(void) const;






protected:

  //! Constructor
  TmmBulkModel(const ModelOptions& options);

  /*!
   * \brief Reads complex refractive index from database
   *
   * Optical refractive index is provided in an extra ASCII file, and in the
   * \c permittivity block of the database this file is indicated by
   * \c optical_data, which has to be a relative path (relative to the directory where
   * the main database file resides)
   */
  virtual void read_database(void) override;

  /*!
   * \brief Setup the optical constants
   *
   * The data file is read, and the internal data structure containing complex
   * refractive index (n, k) is set up.
   */
  virtual void do_init(void) override;

  virtual void prepare_submodels(void) override;


  virtual void do_print_info(void) override;

  /*!
   * \brief Interpolate (n, k) to a given wavelength
   *
   * If \c wavelength is outside the data range, the nearest value will be taken
   */
  std::pair<double, double> interpolate(double wavelength) const;

private:

  //! The constructor method
  static TiberModelObject* _create(const ModelOptions& options);

  //! The destructor method
  static void _destroy(TiberModelObject* p);

  //! The file where we read optical data from
  std::string _datafile;

  //! The wavelengths
  std::vector<double> _wavelengths;

  //! n
  std::vector<double> _n_data;

  //! k
  std::vector<double> _k_data;

  //! An addressing array for faster access during interpolation
  std::vector<int> _addressing;

  double _incoherent_index{ 0 };

  std::vector<TmmDipoleSource* > _DS;

  double _emission_power;



};




#endif // TC_TMMBULKMODEL_H

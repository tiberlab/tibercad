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
 * \file TightBinding.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */

#ifndef TC_TIGHTBINDING_H
#define TC_TIGHTBINDING_H

//-----------------------------------------------------------------------------------------

#include "tibercad/physics/schroedinger/EigenvalueProblem.h"
#include "tibercad/module/SimulationEnvironment.h"
#include "tibercad/io/Database.h"
#include "tibercad/atomistic/Specie.h"


//forward declaration
class Device;
class DftbpWrapper;
class UptWrapper;

//!Main class for Atomistic Tight Binding simulation at equilibrium
class TightBinding : public EigenvalueProblem{


public:

  enum Unused
  {
    UNKNOWN = 0,
    CHARGE
  };

  enum Shell
  {
    NONE = 0,
    S = 1,
    P = 2,
    D = 3
  };


  //! Constructor
  TightBinding(const ModelOptions& options);

  //! Destructor
  ~TightBinding();

  // Create TightBinding object
  //static TightBinding* create(const ModelOptions& options);


  virtual void get_solution_secure(const Elem* elem,
      const std::set<ID>& ids, std::vector<std::map<ID, double> >& values);

  virtual void
  get_solution_secure(const Elem* elem, const std::vector<Point>& p,
      const std::set<ID>& ids, std::vector<std::map<ID, double> >& values);
 

  

private:


 
protected:

  virtual void do_init (void);

  virtual void do_solve (void);

  virtual void parse_options(void);

  //! Get hubbard parameters for species composing atomistic structure
  /*!
   *  (it fills a class member map _u_hub)
   */
  virtual void obtain_hubbard_parameters(void);


  //! Map of map containing hubbard parameters for any specie and any shell
  /*!
   * Usage: _u_hub[<specie>][shell] = hubbard_index
   */
  std::map<Specie, std::map<Shell, double> > _u_hub;

  //! Build charge density on given point
  double build_rho(const Elem* elem, const Point& r);

  // Build states density 
//  virtual void 
  //build_statedens(std::vector<double>& values, const Point& r){};

  //! Charge variation (Mulliken Analisys) on each atom
  std::vector<double> _mulliken_netcharges;

  //! Build a vector of potential projection over atom orbitals
  void project_potential(const std::string providing_model, const std::string mode);

  //!Vector for atom-projected potential shifts
  std::vector<double> _pot_shift;

  //! minimum potential
  double _pot_min;

  //!Vector for atom-projected electron chemical potential
  std::vector<double> _el_chem_pot;

  //!Vector for atom-projected hole chemical potential
  std::vector<double> _hl_chem_pot;


};


//inline
//TightBinding* TightBinding::create(const ModelOptions& options)
//{
//  return new  TightBinding(options);
//}



#endif

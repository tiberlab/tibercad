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
 * \file WzSemiconductor.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _WZSEMICONDUCTOR_H_
#define _WZSEMICONDUCTOR_H_


#include "tibercad/physics/semiconductormodels/Semiconductor.h"
#include "tibercad/physics/semiconductormodels/KPbulkHamiltonian.h"

#include <vector>
//! A class for wurtzite crystal
class  WzSemiconductor : public Semiconductor

{ 

 public:

  //! all the necessary parameters for wurtzite drift-diffusion calculation
  struct WzDDparameters
  {
    //!valence band energy (averaged)  [eV]  
    double       Ev; 
    
    //!band gap  [eV]
    double       EgGamma;
    
    //!conduction band mass along [0001] z direction [free electron mass]
    double       m_c_zz;

    //!conduction band mass along [1-210] x direction [free electron mass]
    double       m_c_xx; 

    //! kp parameter  \f$ A_1 \f$
    double       A1;

    //! kp parameter  \f$ A_2 \f$
    double       A2;

    //! kp parameter  \f$  A_3 \f$
    double       A3;

    //! kp parameter  \f$  A_4 \f$
    double       A4;

    //! kp parameter  \f$  A_5 \f$
    double       A5;

    //! kp parameter  \f$  A_6 \f$
    double       A6;

    //! kp parameter  \f$  A_7 \f$
    double       A7;

    //!hydrostatic deformation \f$ a_x \f$ potential of conduction band  [eV]
    double       a_x;
    
    //!hydrostatic deformation \f$ a_z \f$ potential of conduction band [eV]
    double       a_z; 

    //! deformation \f$ D_1 \f$ potential of valence  band [eV]
    double       D1;

    //! deformation \f$ D_2 \f$ potential of valence  band [eV]
    double       D2;

    //! deformation \f$ D_3 \f$ potential of valence  band  [eV]
    double       D3;
 
    //! deformation \f$ D_4 \f$ potential of valence  band  [eV]
    double       D4;


    //! deformation \f$ D_5 \f$ potential of valence  band [eV]
    double       D5;

    //!deformation \f$ D_6 \f$ potential of valence  band [eV]
    double       D6; 

    //!spin-orbit spliting energy  [eV]
    double       delta_s;

    //!crystal field spliting energy [eV]
    double  delta_cr; 
   
    //! optical matrix element \f$ 2\frac{\langle X |{\bf P}|S \rangle ^2}{m_0} \f$  [eV]
    double Ep_1 ;

    //! optical matrix element \f$ 2\frac{\langle Z |{\bf P}|S \rangle ^2}{m_0} \f$  [eV]
    double Ep_2 ;


    //! Varshni parameter alpha for Gamma valley
    double varshni_alpha_G;
    //!Varshni parameter beta for Gamma valley
    double varshni_beta_G; 


  };

 
 
  //! Get a writeable reference to the physical parameters
  const WzDDparameters& get_parameters(void);

  //! Get a  reference to the initial physical parameters
  const WzDDparameters& get_initial_parameters(void) const;


 
  static WzSemiconductor* create(const ModelOptions& options);
 
  //! apply varshni formulas
  virtual void apply_temperature(void);

 protected:

  //Constructor
  WzSemiconductor(const ModelOptions& options);

  virtual PhysicalModel* create_new(void) const;

  virtual void do_init(void);

  virtual void read_database(void);

  virtual void read_database_alloy(void);


  virtual void do_calculate_kp_params (KPparams& par);

 private:
  //!parameters that TiberCAD should use (e.g. for the actual temperature)
  WzDDparameters par; 

  //!initial parameters from the database (e.g. zero temperature)
  WzDDparameters  par_initial;


  //!bowing parameters
  WzDDparameters bow; 

  //! Calculates k.p parameters in atomic units for 6 band valence band calculation (see below)
  /*!
    \f$
    L_1 = \frac{1}{2} (A_5 + A_4 + A_2  - 1) \\
    L_2 = \frac{1}{2} (A_1 - 1.0)\\
    M_1 = \frac{1}{2} (A_4 + A_2 - A_5 - 1)\\
    M_2 = \frac{1}{2} (A_1 + A_3 - 1)\\
    M_3 = \frac{1}{2} (A_2 - 1.0) \\
    N_1 = A_5 \\
    N_2 = \frac{A_6}{\sqrt{2}}\\
    \mbox{} \\
    N_1^{yx} = M_1; N_1^{xy} = N_1 - N_1^{yx} \\
    N_2^{yx} = M_2; N_2^{xy} = N_2 - N_2^{xy} \\
    \mbox{}\\
    l_1 = D_5 + D_4 + D_2 \\
    l_2 = D_1 \\
    m_1 = D_4 + D_2 - D_5 \\
    m_2 = D_1 + D_3 \\
    m_3 = D_2 \\
    n_1 = 2 D_5\\
    n_2 = \sqrt{2} D_6 \\
    \mbox{}\\
    d_1 = \Delta_{cr} \\
    d_2 = \Delta_{so} \\
    d_3 = \Delta_{so} \\  
    \mbox{}\\
    P_1 = \sqrt(Ep_1/2)\\
    P_2 = \sqrt(Ep_2/2)\\
    \f$
    
  */
  void calculate_6x6_kp_params(KPparams& par);

  void calculate_2x2_kp_params(KPparams& par);

  void calculate_8x8_kp_params(KPparams& par);

  void calculate_14x14_kp_params(KPparams& par);
 
};

inline PhysicalModel* WzSemiconductor::create_new( ) const
{
  return ( new WzSemiconductor(get_options()) );
}

inline WzSemiconductor* WzSemiconductor::create(const ModelOptions& options)
{
  return  new WzSemiconductor(options) ;
}



inline
const WzSemiconductor::WzDDparameters& WzSemiconductor::get_parameters(void)
{
  if (_consider_temperature) apply_temperature();
  return par;
}

inline  
const WzSemiconductor::WzDDparameters& WzSemiconductor::get_initial_parameters() const
{
  return par_initial;
}


#endif 

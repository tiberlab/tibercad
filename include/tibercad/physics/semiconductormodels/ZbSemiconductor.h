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
 * \file ZbSemiconductor.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */

#ifndef TC_ZBSEMICONDUCTOR_H
#define TC_ZBSEMICONDUCTOR_H


#include "tibercad/physics/semiconductormodels/Semiconductor.h"
#include <vector>
#include <complex>
#include "tibercad/physics/semiconductormodels/KPbulkHamiltonian.h" 

//! A class for a zinc-blend (or diamond) crystal
class ZbSemiconductor  : public Semiconductor

{
 public:
  //! data structure  for zinc-blende structure parameters
  struct ZbDDparameters
  {
   

    double Ev;//!< valence band averaged energy \f$  E_v^{\Gamma} \f$ [eV]
    double EgGamma; //!< band gap \f$ E_c^{\Gamma} - E_v^{\Gamma} \f$  [eV]    
    double EgL;     //!< band gap \f$ E_c^{L} - E_v^{\Gamma} \f$  [eV]    
    double EgX;     //!< band gap \f$ E_c^{X} - E_v^{\Gamma} \f$   [eV]  
    double gamma1;  //!< Luttinger \f$ \gamma_1 \f$
    double gamma2;  //!< Luttinger \f$ \gamma_2 \f$
    double gamma3;  //!< Luttinger \f$ \gamma_3 \f$
    double m_t_L;   //!< transversal mass in L point of conduction band [free electron mass]
    double m_l_L;   //!< logitudinal mass in L point of conduction band [free electron mass]
    double m_t_X;   //!< transversal mass in X point of conduction band [free electron mass]
    double m_l_X;   //!< logitudinal mass in X point of conduction band [free electron mass]
    double m_G;     //!< mass in \f$ \Gamma \f$ minimum  point of conduction band 
    double a_c;     //!< conduction band deformation potential [eV]
    double a_v;     //!< valence band deformation potential (hydrostatic) [eV]
    double b;       //!< valence band deformation potential b (uniaxial) [eV] 
    double d;       //!< valence band deformation potential d (uniaxial) [eV] 
    double def_vol_X; //!<   volume deformation potential for X point \f$ \Xi_d + \frac{1}{3}\Xi_u  \f$ [eV]
    double def_uniax_X; //!< uniaxial deformation potential for X point \f$ \Xi_u \f$ [eV]
    double def_vol_L;  //!<  volume deformation potential for L point \f$ \Xi_d + \frac{1}{3}\Xi_u \f$ [eV]
    double def_uniax_L; //!< uniaxial deformation potential for L point \f$ \Xi_u \f$ [eV]
    double delta; //!< spin-orbit \f$ \Delta \f$ [eV]
    double Ep; //!< optical matrix element \f$ 2\frac{\langle X |{\bf P}|S \rangle ^2}{m_0} \f$  [eV]

    double varshni_alpha_G; //!< Varshni parameter alpha for Gamma valley
    double varshni_alpha_X; //!< Varshni parameter alpha for X valley
    double varshni_alpha_L; //!< Varshni parameter alpha for L valley

    double varshni_beta_G; //!< Varshni parameter beta for Gamma valley
    double varshni_beta_X; //!< Varshni parameter beta for X valley
    double varshni_beta_L; //!< Varshni parameter beta for L valley

    double Eg1;           //!< band gap between Gamma_8c and Gamma_7v (from valence to second conduction)
    double delta_c;       //!< spin-orbit in second conduction band
    double delta_cf;      //!< spin-orbit (crystal field?) between conduction-valence bands
    double Ep1;           //!< coupling terms between second conduction and conduction
    double Ep2;           //!< coupling terms between second conduction and valence
   
    double m_c2; 
    

  };

  

  

  //! Get a  reference to the physical parameters
  const ZbDDparameters& get_parameters(void) ;


  //! Get a  reference to the initial physical parameters
  const ZbDDparameters& get_initial_parameters(void) const;

  
  //! apply varshni formulas
  virtual void apply_temperature(void) ;

  static ZbSemiconductor* create(const ModelOptions& options);

 protected:

  //!Constructor
  ZbSemiconductor(const ModelOptions& options);

  virtual PhysicalModel* create_new(void) const;

  virtual void do_init(void);

  virtual void read_database(void);

  virtual void read_database_alloy(void);
 
  virtual void do_calculate_kp_params (KPparams& par);
 

 private:

  //-------------------------------------------------------------------------------//
  //material data block:

  //! parameters that TiberCAD should use (e.g. for the actual temperature)
  ZbDDparameters  par;

  //!initial parameters from the database (e.g. zero temperature)
  ZbDDparameters  par_initial;

  ZbDDparameters  bow;

  //! Calculates k.p parameters in atomic units for 2 band calculation
  void calculate_2x2_kp_params(KPparams& par);


  //! Calculates k.p parameters in atomic units for 6 band valence band calculation (see below)
  /*! Valence band k.p parameters:

      \f$

          L = \frac{1}{2} (-\gamma_1 - 4 \gamma_2 - 1)  ; \\

          M = \frac{1}{2} ( 2\gamma_2 - \gamma_1  - 1 ) ;  \\

          N = -3\gamma_3; \\

          N_{yx} = M; \\

          N_{xy} = N -  N_{yx};
          
      \f$

      Valence band deformation potential:

      \f$
        l  =  a_v + 2b; \\
        m  =  a_v  - b; \\
        n  =  \sqrt{3} d. \\
      \f$

     
      Averaged valence band energy:

      \f$
      \bar{E}_v = E_v - \frac{\Delta}{3};
      \f$

  */
   void calculate_6x6_kp_params(KPparams& par);

  //! Calculates k.p parameters in atomic units for 8 band valence band calculation (see below)
  /*!
    \f$ 
   
     S = 1, \mbox{i.e.}  H_{cc} = \frac{k^2}{2m_0}; \\
     E_c = E_v + \frac{\Delta}{3} + Eg; \\
     E_p = 2 \frac{P^2}{m_0};\\
     E_p = S \left(\frac{m_0}{m_c} - 1 \right) E_g \frac{E_g + \Delta}{E_g + 2/3 \Delta}\\
     \f$  

    \f$ 
    L^{8 \times 8} =  L^{6 \times 6} + \frac{P^2}{E_g + {\Delta}/{3}} ;\\
    N^{8 \times 8} =  N^{6 \times 6} + \frac{P^2}{E_g + {\Delta}/{3}} ;\\
    M^{8 \times 8} =  M^{6 \times 6}. \\
    \f$

  */
   void calculate_8x8_kp_params(KPparams& par);

  //! Calculates k.p parameters in atomic units for 14 band valence band calculation (see below)
  /*!
    \f$

     S =  \\
     E_c = E_v + \frac{\Delta}{3} + Eg \\

     \f$

    \f$
    L^{14 \times 14} =  L^{8 \times 8};\\
    N^{14 \times 14} =  N^{8 \times 8} + \frac{P^2}{E_g1 + {\Delta}/{3} + 2 {\Delta_c}/{3}} ;\\
    M^{14 \times 14} =  M^{8 \times 8} + \frac{P^2}{E_g1 + {\Delta}/{3} + 2 {\Delta_c}/{3}}. \\
    \f$

  */
  void calculate_14x14_kp_params(KPparams& par);
  //--------------------------------------------------------------------------------//


};



inline PhysicalModel* ZbSemiconductor::create_new(void) const
{
  return ( new ZbSemiconductor(get_options()) );
}

inline ZbSemiconductor* ZbSemiconductor::create(const ModelOptions& options)
{
  return new ZbSemiconductor(options);
}



inline const ZbSemiconductor::ZbDDparameters& ZbSemiconductor::get_parameters() 
{
  
  if (_consider_temperature) apply_temperature();

  return(par);
}


inline const ZbSemiconductor::ZbDDparameters& ZbSemiconductor::get_initial_parameters(void) const
{
  return (par_initial);
}


#endif 

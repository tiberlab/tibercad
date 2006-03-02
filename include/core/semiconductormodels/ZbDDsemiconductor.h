#ifndef _ZBDDSEMICONDUCTOR_H_
#define _ZBDDSEMICONDUCTOR_H_


#include "DDsemiconductor.h"
#include<vector>

//! A class to provide all neccessary parameters for drift-diffusion calculation for a zinc-blend (or diamond) crystal
class ZbDDsemiconductor  : public DDsemiconductor
/*!
  The class can calculate information about the band structure, such as
  band edge energy, effective mass for the density of states calculation and
  degeneracy.
  Conduction band masses do not depend on strain.
*/
{
 public:
  //! data structure  for zinc-blende structure parameters
  struct ZbDDparameters
  {
    double       Ev;//!< valence band top energy \f$  E_v^{\Gamma} \f$ [eV]
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
    double def_vol_X; //!<   volume deformation potential for X point \f$ \Xi_d \f$ [eV]
    double def_uniax_X; //!< uniaxial deformation potential for X point \f$ \Xi_u \f$ [eV]
    double def_vol_L;  //!<  volume deformation potential for L point \f$ \Xi_d \f$ [eV]
    double def_uniax_L; //!< uniaxial deformation potential for L point \f$ \Xi_u \f$ [eV]
    double delta; //!< spin-orbit \f$ \Delta \f$ [eV]
  };
  //!Constructor  
  ZbDDsemiconductor(void);
  
  //!Constructor with parameters
  /*!
    \param params all the necessary parameters
   */
  ZbDDsemiconductor(const ZbDDparameters& params);


  //! sets 3 bandgaps
  /*!
      \param EgGamma \f$ E_g^{\Gamma} =  E_c^{\Gamma} - E_v^{\Gamma} \f$
      \param EgL \f$ E_g^{L} = E_c^{L} - E_v^{\Gamma} \f$
      \param EgX \f$ E_g^{X} = E_c^{X} - E_v^{\Gamma} \f$ 
  */
  void set_Eg(double EgGamma, double EgL, double EgX);

  
  //! sets conduction band mass (isotropic)
  void set_mass_Gamma(double m);

  //! sets logitudinal and transversal mass for L point
  /*!
    \param m_t transversal mass [free electron mass]
    \param m_l logitudinal mass [free electron mass]
   */
  void set_masses_L(double m_t, double m_l);

  //! sets logitudinal and transversal mass for x point
  /*!
    \param m_t transversal mass [free electron mass]
    \param m_l logitudinal mass [free electron mass]
   */
  void set_masses_X(double m_t, double m_l);

  //! sets parameters for the valence band \f${\bf k \cdot p}\f$ description
  /*!
     \param gamma1 Luttinger \f$ \gamma_1 \f$
     \param gamma2 Luttinger \f$ \gamma_2 \f$
     \param gamma3 Luttinger \f$ \gamma_3 \f$
     \param delta  Spin-orbit spliting energy  \f$ \Delta \f$ [eV]  
   */
  void set_6x6kp_params(double gamma1, double gamma2, double gamma3, double delta);
  
  //! sets deformation potential for \f$ \Gamma \f$ point
  /*!
    \param a_c conduction band deformation potential [eV]
    \param a_v valence band volumic deformation potential [eV]
    \param b valence band deformation potential b [eV]
    \param d valence band deformation potential d [eV]
  */
  void set_deformation_parameters(double a_c, double a_v, double b, double d);

  //! calculates information about conduction bands
  /*!
    \f$ E_c^{\Gamma} = E_{c0}^{\Gamma} + a_c  \mathop{\rm Tr} (\varepsilon_{ij}) ;\f$

    \f$ E_c^{X} = E_{c0}^{X} + \Xi_d^{X}  \mathop{\rm Tr} (\varepsilon_{ij}) +
     \Xi_d^{X} k_i k_j \varepsilon_{ij} ;\f$

    \f$ E_c^{L} = E_{c0}^{L} + \Xi_d^{L}  \mathop{\rm Tr} (\varepsilon_{ij}) +
     \Xi_d^{L} k_i k_j \varepsilon_{ij} ,\f$

    where k is a unit vector from the Brillouin zone center to a minimum point 
  */
  virtual void  calculate_conduction_band_extremum(void);
  
  //! calculates information about valence bands
  /*!
    Uses 6 band Luttinger kp theory 
  */
  virtual void  calculate_valence_band_energy_extremum(void);

 private:

  //-------------------------------------------------------------------------------//
  //material data block:

  //! band gap \f$ E_c^{\Gamma} - E_v^{\Gamma} \f$  [eV] 
  double EgGamma; 

  //! band gap \f$ E_c^{L} - E_v^{\Gamma} \f$  [eV] 
  double EgL;
  
  //! band gap \f$ E_c^{X} - E_v^{\Gamma} \f$  [eV] 
  double EgX;

  //!Luttinger \f$ \gamma_1 \f$
  double gamma1;   

  //!Luttinger \f$ \gamma_2 \f$
  double gamma2; 

  //!Luttinger \f$ \gamma_3 \f$     
  double gamma3;         

  //!transversal mass in L point of conduction band [free electron mass]
  double m_t_L;

  //!logitudinal mass in L point of conduction band [free electron mass]
  double m_l_L;

  //!  transversal mass in X point of conduction band [free electron mass] 
  double m_t_X;

  //! logitudinal mass in X point of conduction band [free electron mass]
  double m_l_X;   
 
 // mass in \f$ \Gamma \f$ minimum  point of conduction band 
  double m_G;    

  // conduction band deformation potential [eV]
  double a_c;    

  // valence band deformation potential (hydrostatic) [eV]
  double a_v;    

  // valence band deformation potential b (uniaxial) [eV] 
  double b;       

  // valence band deformation potential d (uniaxial) [eV] 
  double d;      

  //   volume deformation potential for X point \f$ \Xi_d \f$ [eV]
  double def_vol_X; 

  // uniaxial deformation potential for X point \f$ \Xi_u \f$ [eV]
  double def_uniax_X; 

  //  volume deformation potential for L point \f$ \Xi_d \f$ [eV]
  double def_vol_L; 

  // uniaxial deformation potential for L point \f$ \Xi_u \f$ [eV]
  double def_uniax_L;

  // spin-orbit \f$ \Delta \f$ [eV]
  double delta;
  // end of material data block
 
};


#endif 

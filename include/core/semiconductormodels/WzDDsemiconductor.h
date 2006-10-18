#ifndef _WZDDSEMICONDUCTOR_H_
#define _WZDDSEMICONDUCTOR_H_


#include "DDsemiconductor.h"
#include "KPbulkHamiltonian.h"

#include<vector>
//! A class to provide all neccessary parameters for drift-diffusion calculation for a wurtzite  crystal.
class  WzDDsemiconductor : public DDsemiconductor
/*!
  The class can calculate information about the band structure, such as
  band edge energy, effective mass for the density of states calculation and
  degeneracy.
  Conduction band masses do not depend on strain.
  Only \f$ \Gamma \f$ minimum of conduction band is considered. 
*/
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
   
    //! optical matrix element \f$ \frac{\langle X|{\bf P}|S \rangle ^2}{2m_0} \f$  [eV]
    double Ep_1 ;

    //! optical matrix element \f$ \frac{\langle |{\bf P}|S \rangle ^2}{2m_0} \f$  [eV]
    double Ep_2 ;

  };

  //Constructor
  WzDDsemiconductor(void);

  //Constrauctor. Sets all necessary parameters
  WzDDsemiconductor(const WzDDparameters& params );

  //!sets valence band top energy
  void set_Ev(const double Ev);


  //!set band gap in \f$ \Gamma \f$ point
  void set_Eg(const double EgGamma);

  //!set conduction band effective mass \f$ m_{zz} \f$ and \f$ m_{xx} \f$
  void set_mass_Gamma(const double m_c_zz, const double m_c_xx);

   //!set conduction band deformation potential \f$ a_x \f$ and \f$ a_z \f$
  void set_deform_pot_cond(const double a_x, const double a_z);

  //! calculates information about conduction bands
  /*!
    \f$ E_c^{\Gamma} = E_{c0}^{\Gamma} + a_{x}  (\varepsilon_{xx} + \varepsilon_{yy}) + a_z \varepsilon_{zz} \f$
  */
  virtual void  calculate_conduction_band_extremum(void);

  //!calculates information about valence bands
  /*!
    Uses 6 band Luttinger kp theory 
  */
  virtual void  calculate_valence_band_extremum(void);

  //! Calculates k.p parameters in atomic units for 6 band valence band calculation (not coded yet!)
   virtual KPbulkHamiltonian::KPparams calculate_6x6_kp_params (void );

  //! Calculates k.p parameters in atomic units for 8 band valence band calculation (not coded yet!)
   virtual KPbulkHamiltonian::KPparams calculate_8x8_kp_params (void );

  //! Get a writeable reference to the physical parameters
  WzDDparameters& get_parameters(void);

  //! Set the physical parameters
  void set_parameters(const WzDDparameters& parameters);

  void read_database(const Dummy&);

  void build_alloy(const std::string& component2,
			   const std::string& bowing_params, double content);
 
 private:

  WzDDparameters par; 

  //! Hartree energy in eV
  static const double Hartree;
  
  


};


inline
WzDDsemiconductor::WzDDparameters&
WzDDsemiconductor::get_parameters(void)
{
  return par;
}

inline
void
WzDDsemiconductor::set_parameters(const WzDDparameters& parameters)
{
  par = parameters;
}



#endif 

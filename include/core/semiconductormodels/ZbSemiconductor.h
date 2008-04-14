#ifndef _ZBSEMICONDUCTOR_H_
#define _ZBSEMICONDUCTOR_H_


#include "Semiconductor.h"
#include<vector>
#include <complex>
#include "KPbulkHamiltonian.h" 

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

    

  };

  


  //!Constructor  
  ZbSemiconductor(void);
  

  //! Get a  reference to the physical parameters
  const ZbDDparameters& get_parameters(void) ;


  //! Get a  reference to the initial physical parameters
  const ZbDDparameters& get_initial_parameters(void) const;

  

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
 
  virtual KPparams do_calculate_6x6_kp_params (void );

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
  virtual KPparams do_calculate_8x8_kp_params (void );

  //! apply varshni formulas
  virtual void apply_temperature(void) ;

  static ZbSemiconductor* create();

 private:

  //-------------------------------------------------------------------------------//
  //material data block:

  //! parameters that TiberCAD should use (e.g. for the actual temperature)
  ZbDDparameters  par;

  //!initial parameters from the database (e.g. zero temperature)
  ZbDDparameters  par_initial;

  ZbDDparameters  bow;

  //---------------------------------------------------------------------------------//
  // k.p Hamiltonian section

  //! Hartree energy in eV
  static const double Hartree;



  //--------------------------------------------------------------------------------//


 protected:

  virtual PhysicalModelInterface* create_new(void) const;

  virtual void copy_from (const PhysicalModelInterface *rhs);

  virtual void do_init(void);

  virtual void read_database(void);

  virtual void read_bowing_parameters(void);
 
  virtual void do_calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);
  

 
};



inline PhysicalModelInterface* ZbSemiconductor::create_new(void) const
{
  return ( new ZbSemiconductor() );
}

inline ZbSemiconductor* ZbSemiconductor::create()
{
  return new ZbSemiconductor();
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

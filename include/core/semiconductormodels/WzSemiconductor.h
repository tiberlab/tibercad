#ifndef _WZSEMICONDUCTOR_H_
#define _WZSEMICONDUCTOR_H_


#include "Semiconductor.h"
#include "KPbulkHamiltonian.h"

#include<vector>
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
   
    //! optical matrix element \f$ \frac{\langle X|{\bf P}|S \rangle ^2}{2m_0} \f$  [eV]
    double Ep_1 ;

    //! optical matrix element \f$ \frac{\langle |{\bf P}|S \rangle ^2}{2m_0} \f$  [eV]
    double Ep_2 ;

  };

  //Constructor
  WzSemiconductor(void);

 
  

  //! Calculates k.p parameters in atomic units for 6 band valence band calculation (see below)
  /*!
    \f$
    L_1 = \dfrac{1}{2} (A_5 + A_4 + A_2  - 1); \\
    L_2 = \dfrac{1}{2} (A_1 - 1.0);\\
    M_1 = \dfrac{1}{2} (A_4 + A_2 - A_5 - 1);\\
    M_2 = \dfrac{1}{2} (A_1 + A_3 - 1);\\
    M_3 = \dfrac{1}{2} (A_2 - 1.0); \\
    N_1 = A_5; \\
    N_2 = \dfrac{A_6}{\sqrt{2}}\\;
    \mbox{}
    N_1^{yx} = M1; N_1^{xy} = N_1 - N_1^{yx} \\
    N_2^{yx} = N_1^{yx};N_2^{xy} = result.N1_xy
    
    \f$
  */
  virtual KPparams calculate_6x6_kp_params (void );

  //! Calculates k.p parameters in atomic units for 8 band valence band calculation (see below)
  /*
    
   */
  virtual KPparams calculate_8x8_kp_params (void );

  //! Get a writeable reference to the physical parameters
  WzDDparameters& get_parameters(void);

 
  static WzSemiconductor* create(void); 
 
 private:
  //!parameters
  WzDDparameters par; 

  //!bowing parameters
  WzDDparameters bow; 

  //! Hartree energy in eV
  static const double Hartree;
  
  

 protected:

  virtual PhysicalModelInterface* create_new(void) const;

  virtual void copy_from (const PhysicalModelInterface *rhs);

  virtual void do_init(void);

  virtual void read_database(void);

  virtual void read_bowing_parameters(void);
 
  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa);


};

inline PhysicalModelInterface* WzSemiconductor::create_new( ) const
{
  return ( new WzSemiconductor() );
}

inline WzSemiconductor* WzSemiconductor::create() 
{
  return  new WzSemiconductor() ;
}



inline
WzSemiconductor::WzDDparameters&
WzSemiconductor::get_parameters(void)
{
  return par;
}




#endif 

// $Id$

#ifndef _OPTICS_H_
#define _OPTICS_H_

#include "SimulationInterface.h"


//! A base class of optics calculation
/*!
 * The task of this class is to calculate the spectrum at a certain given
 * k-point
 */
class Optics : public SimulationInterface
{

  public:

    //! The constructor
    Optics(const ModelOptions& options) : SimulationInterface(options) { };

    //! The destructor
    virtual ~Optics(void) { };


  protected:

    //! do_solve() needs to be implemented by derived classes
    virtual void do_solve(void) = 0;

    //!calculate spectrum 
    /*!
      \f$
      
      P(\hbar \omega) = \sum_{i,j} \frac{1}{2\pi^2}  \frac{\omega^2_{ij} e^2 }{m^2 c^3}  |{\bf M_{i,j} e}|^2 f_i(E_i)(1 - f_j(E_j)) 
      \frac{\Gamma/2} {(\hbar \omega_{ij} - \hbar \omega)^2 + (\Gamma/2)^2} d\Omega
      \f$
      
      \param Energy energy grid [eV]
      \param spectrum calculated spectrum (atomic units)
      \param Gamma broadering parameter [eV]
      \param polariz polarization vector of a linearly polarized light (must be a normalized one, \f$ |{\bf e}| = 1 \f$)
      
    */
    //virtual void calculate_spectrum(const Mesh& Energy, double Gamma,const Tensor1& polariz, 
    //                      std::map<const Elem*, double>& spectrum);

  private:




};


#endif // _OPTICS_H_

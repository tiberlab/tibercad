/*
 * VariableRangeHopping.h
 *
 *  Created on: 25 Jan 2022
 *      Author: miesu
 */

#ifndef _VARIABLERANGEHOPPING_H_
#define _VARIABLERANGEHOPPING_H_

#include "MobilityModelInterface.h"

 /*!
  * \brief Implementation of variable range hopping mobility
  *
  * Formula of the mobility
  *
  * $$
  * \mu=\mu_{0} \exp \left
  * \{-\left(\frac{k_{B} T_{h}}{k_{B} T+q E r}\right)^{\gamma}\right\}
  * $$
  *
  * Chen, T., van Gelder, J., van de Ven, B. et al. Classification
  * with a disordered dopant-atom network in silicon. Nature 577, 341–345 (2020).
  * https://doi.org/10.1038/s41586-019-1901-0
  */
class VariableRangeHopping : public MobilityModelInterface
{
  public:


    //! Destructor
    virtual ~VariableRangeHopping(void);



    //! Create a VariableRangeHopping object
    static VariableRangeHopping* create(const ModelOptions& options);



    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void);



    //! \copydoc MobilityModelInterface::get_derivative_grad_potential()
    /*!
     *   Derivative of the mobility
     *   $$
     *   \frac{d \mu}{d \nabla \phi}=-\frac{d \mu}{d \bar{E}}=-A
     *   \cdot B^{\gamma} \cdot \exp \left[-B^{\gamma}\right]
     *   \cdot \frac{\bar{E}}{|E|}
     *   $$
     *   \bigskip
     *   $$
     *   A=\frac{\mu_{0} \cdot q \cdot r \cdot \gamma}{k_{B}
     *   \cdot T_{h} + q \cdot |E| \cdot r}
     *   $$
     *   $$
     *   B=\frac{K_{B} T_{h}}{K_{B} T+q|E| \cdot r}
     *   $$
     */
    virtual void get_derivative_grad_potential(libMesh::RealGradient& dm);



  protected:



    //! Constructor
    VariableRangeHopping(const ModelOptions& options);



    //! \copydoc MobilityModelInterface::do_init()
    virtual void do_init(void);



  private:

    //! prefactor of mobility
    double _mu_0;

    //! average hopping distance r
    double _hopping_distance;

    //! characteristic temperature of hopping conduction Th
    double _temp_hopping_conduction;

    //! exponent
    double _gamma;
};






inline
VariableRangeHopping::VariableRangeHopping(const ModelOptions& options)
 : MobilityModelInterface(options),
   _mu_0(100.0),
  _gamma(1.0/3.0),
  _hopping_distance(1.0e-6),
  _temp_hopping_conduction(1200.0)
{
}






inline
VariableRangeHopping*
VariableRangeHopping::create (const ModelOptions& options)
{
  return new VariableRangeHopping(options);
}






inline
VariableRangeHopping::~VariableRangeHopping(void)
{
}



#endif /*_VARIABLERANGEHOPPING_H_ */

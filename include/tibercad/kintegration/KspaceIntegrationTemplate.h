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
 * \file KspaceIntegrationTemplate.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */

#ifndef __KINTEGRATIONTEMPLATE__
#define __KINTEGRATIONTEMPLATE__

#include "tibercad/kintegration/KspaceIntegration.h"
#include "tibercad/base/ModelOptions.h"


//! Template class for KspaceIntegration.
/*!
 * The template is used to store a general pointer to an object invoking kintegration
 * The class must implement a method \c calculate_for_k_point(..)
 * With the same interface as described here.
 */
template <class T>
class KspaceIntegrationTemplate : public KspaceIntegration
{
  public:

    //! callback type for the simple integration
    typedef void (T::*Callback)(const Point&, DofField&, double&);

    //! callback type for integration in ordered mode
    /*!
     * The second Point is the last called k-point which is a neighor of the new one
     */
    typedef void (T::*Callback2)(const Point&, const Point&, DofField&, double&);

    virtual ~KspaceIntegrationTemplate(void){};


  protected:

    //! calculates density for each k-point
    /*!
     *
     */
    virtual void calculate_for_k_point(const Point& k_point,
        const Point& refpoint,
        DofField& density,
        double& estimator);

  private:

    //! Let base class be friend
    friend class KspaceIntegration;

    KspaceIntegrationTemplate(T* hook, 
                              const ModelOptions& opts,
                              const libMesh::Parallel::Communicator& k_comm);

    //! The object which knows how to calculate the density
    T* _hook;

    //! The method of object to be called
    Callback _callback;

    //! Another callback
    Callback2 _callback2;

};

template <class T>
inline
KspaceIntegrationTemplate<T>::KspaceIntegrationTemplate(T* hook, 
                              const ModelOptions& opt,
                              const libMesh::Parallel::Communicator& k_comm)
  : KspaceIntegration(opt, k_comm),
    _callback(0),
    _callback2(0)
{
   _hook = hook;
}





template <class T>
inline
void KspaceIntegrationTemplate<T>::calculate_for_k_point(const Point& k_point, 
    const Point& refpoint,
    DofField& density,
    double& integrated_quantity)
{
  if (_callback2 != 0)
    (_hook->*_callback2)(k_point, refpoint, density, integrated_quantity);
  else if (_callback != 0)
    (_hook->*_callback)(k_point, density, integrated_quantity);
}




#endif

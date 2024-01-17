#ifndef __KINTEGRATIONTEMPLATE__
#define __KINTEGRATIONTEMPLATE__

#include "KspaceIntegration.h"
#include "ModelOptions.h"


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

#ifndef __KINTEGRATIONTEMPLATE__
#define __KINTEGRATIONTEMPLATE__

#include "KspaceIntegration.h"
#include "ModelOptions.h"

//! Template class for KspaceIntegration.
//! The template is used to store a general pointer to an object invoking kintegration
//! The class must implement a method "calculate_for_k_point(..)"  
//! With the same interface as described here.
  
template <class T>
class KspaceIntegrationTemplate : public KspaceIntegration 
{
  public:

    KspaceIntegrationTemplate(T* hook, const ModelOptions& opts);

    virtual ~KspaceIntegrationTemplate(void){};

    static KspaceIntegration* create(T* hook, const ModelOptions& opts);
		 

  protected:

   //!calculates density that is necessary for eack k-point and a number that will be used for refinement 
   virtual void calculate_for_k_point(const Point& k_point, 
                                      DofField& density, 
                                      double& estimator);

  private:

    T* _hook;

};

template <class T>
inline
KspaceIntegrationTemplate<T>::KspaceIntegrationTemplate(T* hook, const ModelOptions& opt)
	: KspaceIntegration(opt)
{
   _hook = hook;
}

//template <class T>
//inline
//KspaceIntegrationTemplate<T>::~KspaceIntegrationTemplate
//{
//}

template <class T>
inline
void KspaceIntegrationTemplate<T>::calculate_for_k_point(const Point& k_point, 
        	         			     DofField& density, 
	         	         		     double& integrated_quantity)
{
    _hook->calculate_for_k_point(k_point, density, integrated_quantity);
}


template <class T>
inline
KspaceIntegration* KspaceIntegrationTemplate<T>::create(T* hook, const ModelOptions& opts)
{
  return new KspaceIntegrationTemplate<T>(hook, opts);
}

#endif

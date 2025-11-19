// $Id$

#ifndef _INITIALIZERBASE_H_
#define _INITIALIZERBASE_H_

class TiberModelObject;

//! The base class for initializer functors
template <typename T>
class InitializerBase
{

  public:

    //! Destructor
    virtual ~InitializerBase(void) {};

    //! The operator to be overloaded
    virtual void operator()(T& val) = 0;


  protected:


};

#endif /* _INITIALIZERBASE_H_ */

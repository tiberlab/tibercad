// $Id$

#ifndef _INITIALIZERBASE_H_
#define _INITIALIZERBASE_H_

class TiberModelObject;

//! The base class for initializer functors
class InitializerBase
{

  public:

    //! Destructor
    virtual ~InitializerBase(void) {};

    //! The operator to be overloaded
    virtual void operator()(void) {};


  protected:

    typedef void (TiberModelObject::*CastFunc)(void);

    //! Constructor
    InitializerBase(TiberModelObject& obj, CastFunc func)
      : _obj(&obj), _func(func) { };

    //! The object containing the function
    TiberModelObject* _obj;

    //! The function itself
    CastFunc _func;

};

#endif /* _INITIALIZERBASE_H_ */

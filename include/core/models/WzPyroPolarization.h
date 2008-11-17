// $Id$

#ifndef _WZPYROPOLARIZATION_H_
#define _WZPYROPOLARIZATION_H_

#include "PyroPolarization.h"


//! The wurtzite implementation of pyropolarization
class WzPyroPolarization : public PyroPolarization
{

  public:

    //! The create method
    static WzPyroPolarization* create(void);


  protected:

    //! Constructor
    WzPyroPolarization(void);

    //! Destructor
    virtual ~WzPyroPolarization(void);
  

    //! Create new PyroPolarization
    virtual PhysicalModelInterface* create_new(void) const;

    //! Read from database
    virtual void read_database(void);

    //! Read from alloy parameters database
    virtual void read_database_alloy(void);

    //! Initialize
    virtual void do_init(void);
  
    //! Calculate VCA for alloy
    virtual void do_init_alloy(const PhysicalModelInterface *comp_A,
        const PhysicalModelInterface *comp_B, double xa);


  private:

    //! The pyropolarization along the c-axis
    double _Pz;

    //! The bowing parameter for Pz
    double _Pz_bow;

};


//
// inline methods
// 

inline
WzPyroPolarization::WzPyroPolarization(void)
  : _Pz(0.0),
    _Pz_bow(0.0)
{
}


inline
WzPyroPolarization::~WzPyroPolarization(void)
{
}



#endif // _WZPYROPOLARIZATION_H_

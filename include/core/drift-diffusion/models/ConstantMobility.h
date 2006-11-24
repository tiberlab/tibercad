// $Id$

#ifndef _CONSTANTMOBILITY_H_
#define _CONSTANTMOBILITY_H_

#include "MobilityModelInterface.h"


class ConstantMobility : public MobilityModelInterface
{

  public:

    //! Destructor
    virtual ~ConstantMobility(void);

    //! Create a ConstantMobility object
    static ConstantMobility* create(void);

    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void);

    //! \copydoc MobilityModelInterface::get_mobility_derivatives()
    virtual void get_mobility_derivatives(std::vector<double>& dm);


  protected:

    //! constructor
    ConstantMobility(void);

    //! \copydoc MobilityModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc MobilityModelInterface::create_new()
    virtual PhysicalModelInterface* create_new(void) const;

    //! \copydoc MobilityModelInterface::copy_from()
    virtual void copy_from(const PhysicalModelInterface* rhs);

    /*! \copydoc MobilityModelInterface::calculate_VCA() */
    virtual void calculate_VCA(const PhysicalModelInterface* comp_A,
        const PhysicalModelInterface* comp_B, double xa);

  private:

    //! The (constant) mobility
    double _mu0;

};

//
// inline methods
// 

inline
ConstantMobility::ConstantMobility(void)
  :  _mu0(1000)
{
}


inline
ConstantMobility*
ConstantMobility::create(void)
{
  return new ConstantMobility();
}


inline
PhysicalModelInterface*
ConstantMobility::create_new(void) const
{
  return new ConstantMobility();
}


inline
void
ConstantMobility::copy_from(const PhysicalModelInterface* rhs)
{
  MobilityModelInterface::copy_from(rhs);

  const ConstantMobility* mod = dynamic_cast<const ConstantMobility*>(rhs);
  _mu0 = mod->_mu0;
}


inline
ConstantMobility::~ConstantMobility(void)
{
}

#endif // _CONSTANTMOBILITY_H_

// $Id: MobilityModelInterface.h 4184 2015-12-07 12:28:44Z maufder $

#ifndef _MOBILITYMODELINTERFACE_H_
#define _MOBILITYMODELINTERFACE_H_



#include "DriftDiffusionModelInterface.h"

#include "vector_value.h"

#include <vector>

//! The base class for mobility models
class TBDLEXPORT MobilityModelInterface : public DriftDiffusionModelInterface
{

  public:

    //! Destructor
    virtual ~MobilityModelInterface(void);

    //! Get the mobility
    virtual double get_mobility(void) = 0;

    //! Get the derivative of the mobility w.r.t. the electric potential
    virtual double get_derivative_potential(void);

    //! Get the derivative of the mobility w.r.t. the gradient of the electric potential
    virtual void get_derivative_grad_potential(libMesh::RealGradient& dm);

    //! Get the derivatives with respect to the fermi-level gradient
    virtual void get_derivative_grad_fermi(libMesh::RealGradient& dm);

    //! Creates a new named mobility model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     *
     * \param name the model name
     * \param options the options as given in the input file
     * \return a pointer to the newly created object
     */
    static MobilityModelInterface* create(const std::string& name,
        const Material* mat, const ModelOptions& options = ModelOptions());


    //! Set the carrier type
    /*!
     * \param type the carrier type, can be \c e or \c h
     */
    void set_carrier_type(char type);

    void set_carrier(ID carrier);

    //! Get the carrier type
    char get_carrier_type(void) const;

    ID get_carrier(void) const;


  protected:

    //! Create a mobility model
    MobilityModelInterface(const ModelOptions& options);



  private:

    //! The type of carriers this model is for
    /*!
     * Can be \c e or \c h
     */
    char _carrier_type;

    ID _carrier;

};



//
// inline methods
//


inline
MobilityModelInterface::MobilityModelInterface(const ModelOptions& options)
  : DriftDiffusionModelInterface(options)
{
  std::string p = get_options().get_option("carrier", "electron");
  p = get_options().get_option("particle", p);
  if (p == "electron")
    _carrier_type = 'e';
  else
    _carrier_type = 'h';
}


inline
MobilityModelInterface::~MobilityModelInterface(void)
{
}



inline
void
MobilityModelInterface::set_carrier_type(char type)
{
  _carrier_type = type;
}


inline
char
MobilityModelInterface::get_carrier_type(void) const
{
  return _carrier_type;
}

inline
void
MobilityModelInterface::set_carrier(ID carrier)
{
  _carrier = carrier;
}

inline
ID
MobilityModelInterface::get_carrier(void) const
{
  return _carrier;
}

inline
double
MobilityModelInterface::get_derivative_potential(void)
{
  return 0.0;
}

inline
void
MobilityModelInterface::get_derivative_grad_potential(libMesh::RealGradient& dm)
{
  dm.zero();
}

inline
void
MobilityModelInterface::get_derivative_grad_fermi(libMesh::RealGradient& dm)
{
  dm.zero();
}





#endif // _MOBILITYMODELINTERFACE_H_

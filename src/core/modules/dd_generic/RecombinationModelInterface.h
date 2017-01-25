// $Id: RecombinationModelInterface.h 3435 2012-11-15 15:06:50Z maufder $

#ifndef _RECOMBINATIONMODELINTERFACE_H_
#define _RECOMBINATIONMODELINTERFACE_H_



#include "DriftDiffusionModelInterface.h"

#include <vector>

class Boundary;


//! Base class for recombination models
/*!
 * This is the base class for recombination models. A new recombination model
 * can be implemented by deriving from this class.
 *
 * To simplify handling of the standard case of electrons and holes, the API
 * guarantees the order (electron, hole) whenever two particles with these names are
 * provided.
 *
 * Recombination models should be implemented as
 * \[ R_{12} = g(n_1, n_2)*\left( 1 - e^{(E_{F,2}-E_{F,1})/k_BT} \right) \]
 * where 1,2 are the two carriers. The net recombination for carrier 1 is then
 * \$R1 = -sign(q_1)R_{12}\$ and for carrier 2 \$R1 = sign(q_2)R_{12}\$.
 * The sign factor is provided as a weight, so that stoichiometric ratios can be
 * implemented. This handles both electron-hole recombination and transfers between
 * bands of the same type. The direction of transfer is given by the sign of the
 * Fermi level difference.
 *
 */
class TBDLEXPORT RecombinationModelInterface : public DriftDiffusionModelInterface
{

  public:

    //! Destructor
    virtual ~RecombinationModelInterface(void);


    //! Get the recombination rate and derivatives
    void get_net_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials);

    //! Get the weights for the different carriers
    const std::vector<double>& get_weights(void) const;


    //! Creates a new named recombination model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     *
     * \param name the model name
     * \param options the options as given in the input file
     * \return a pointer to the newly created object
     */
    static RecombinationModelInterface* create(const std::string& name,
        const PhysicalObject* owner, const ModelOptions& options = ModelOptions());


    //! Get the associated tunneling contact pointer, or NULL
    const Boundary* get_tunneling_contact(void);

    //! Get the IDs of the recombining carrers;
    const std::vector<ID>& get_carrier_ids(void) const;

    //! Get the names of the recombining carriers
    const std::vector<std::string>& get_carrier_names(void) const;


  protected:

    //! \copydoc DriftDiffusionProperties::DriftDiffusionProperties()
    RecombinationModelInterface(const ModelOptions& options);

    //! Set the associated contact for tunneling
    /*!
     * This method should be called if the recombination model
     * is modelling carrier tunneling to or from a contact.
     */
    void set_tunneling_contact(const Boundary* bd);

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;

    //! Calculate the recombination rate and its derivatives
    virtual void calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials);

    //! Reorder carriers according to a module specific order
    void reorder_carriers(const std::vector<ID>& new_order);

    //! Set the weights, if necessary
    void set_weights(std::vector<double>& weights);


  private:

    //! The associated tunneling contact, or NULL
    const Boundary* _tunneling_boundary;

    //! The names of the carriers as found in the input
    std::vector<std::string> _carriers;

    //! The global IDs for the carriers
    std::vector<ID> _carrier_ids;

    //! The weights for the different carrier
    std::vector<double> _weights;

};


//
// inline methods
//

inline
RecombinationModelInterface::RecombinationModelInterface(const ModelOptions& options)
 : DriftDiffusionModelInterface(options),
   _tunneling_boundary(NULL)
{

}

inline
RecombinationModelInterface::~RecombinationModelInterface(void)
{
}

inline
const Boundary*
RecombinationModelInterface::get_tunneling_contact(void)
{
  return _tunneling_boundary;
}

inline
void
RecombinationModelInterface::set_tunneling_contact(const Boundary* bd)
{
  _tunneling_boundary = bd;
}



inline
RecombinationModelInterface*
RecombinationModelInterface::create(const std::string& name,
    const PhysicalObject* owner, const ModelOptions& options)
{

  return dynamic_cast<RecombinationModelInterface*>(
      PhysicalModelInterface::create("recombination_" + name, owner, options));

}

inline
const std::vector<std::string>&
RecombinationModelInterface::get_carrier_names(void) const
{
  return(_carriers);
}

inline
const std::vector<ID>&
RecombinationModelInterface::get_carrier_ids(void) const
{
  return(_carrier_ids);
}

inline
const std::vector<double>&
RecombinationModelInterface::get_weights(void) const
{
  return(_weights);
}

inline
void
RecombinationModelInterface::set_weights(std::vector<double>& weights)
{
  _weights = weights;
}

#endif // _RECOMBINATIONMODELINTERFACE_H_

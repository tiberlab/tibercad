// $Id: SchottkyTunneling.h 3394 2012-08-25 22:21:45Z maufder $

#ifndef _SCHOTTKYTUNNELING_H_
#define _SCHOTTKYTUNNELING_H_

#include "RecombinationModelInterface.h"
#include "tibercad/base/HashMap.h"
#include "tibercad/base/TypeDefs.h"

#include "point.h"

class SimulationInterface;
class Elem;

//! Implementation of local Schottky Tunneling model
/*!
 * This class implements Schottky tunneling casted into
 * a local recombination model
 */
class TBDLLOCAL SchottkyTunneling : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~SchottkyTunneling(void) {};

    //! Create a ConstantMobility object
    static SchottkyTunneling* create(const ModelOptions& options);

    
  protected:

    //! Constructor
    SchottkyTunneling(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;

    //! Calculate the recombination rate and its derivatives
    virtual void calculate_rate_and_derivatives(std::vector<double>& R,
        std::vector<std::vector<double>>& dPotentials) override;


  private:


    //! The maximum tunneling length in nm
    double _max_tunnel_length;


    //! The name of the associated contact
    std::string _contact_name;


    //! The contact voltage
    double _contact_voltage;


    //! The band
    char _band;


    //! The Schottky barrier
    double _barrier;


    //! The effective mass
    double _mass;


    //! A map containing all elements inside the tunnel length
    HashMap<const Elem*, Point>::Type _elem_map;



};



//
// inline methods
// 



inline
SchottkyTunneling*
SchottkyTunneling::create(const ModelOptions& options)
{
  return new SchottkyTunneling(options);
}






#endif // _SCHOTTKYTUNNELING_H__

// $Id$

#ifndef _SOLUTIONPROVIDER_H_
#define _SOLUTIONPROVIDER_H_

#include "TypeDefs.h"

#include <utility>

class SimulationInterface;

/*!
 * \brief A class providing access to solutions in other modules
 */
class SolutionProvider : public std::pair<SimulationInterface*, ID>
{
  public:

    SolutionProvider()
      : std::pair<SimulationInterface*, ID>(nullptr, INVALID_ID)
    {}


    SolutionProvider(SimulationInterface* si, ID id)
      : std::pair<SimulationInterface*, ID>(si, id)
    {}


    explicit SolutionProvider(const std::pair<SimulationInterface*, ID>& in)
      : std::pair<SimulationInterface*, ID>(in)
    {}


    //! Get the pointer to the simulation object
    SimulationInterface* simulation(void) const
    {
      return(this->first);
    }

    //! Get the ID
    ID id(void) const
    {
      return(this->second);
    }

    //! Check if it is valid
    bool is_valid(void) const
    {
      return((this->first != nullptr) && (this->second != INVALID_ID));
    }

};


#endif // _SOLUTIONPROVIDER_H_

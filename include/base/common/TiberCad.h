#ifndef _TIBERCAD_H_
#define _TIBERCAD_H_

namespace TiberCad
{

  //! The statistics to be used
  enum Statistics
  {
    BOLTZMANN,     /*!< Boltzmann statistics */
    FERMIDIRAC,    /*!< Fermi-Dirac statistics */
    BOSEEINSTEIN   /*!< Bose-Einstein statistics */
  };


  //! Types of symmetries
  enum Symmetry
  {
    NONE,         /*!< no special symmetry */
    CYLINDRICAL   /*!< cylinder symmetry */
  };

};


#endif // _TIBERCAD_H_

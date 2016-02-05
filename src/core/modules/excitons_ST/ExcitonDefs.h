#ifndef _EXCITONDEFS_H_
#define _EXCITONDEFS_H_

//! A namespace for excitons specific definitions
/*!
 */
namespace ExcitonsDefs
{

  //! Specifies the type of coupling used
  /*!
   * The values are hardcoded to hex values so they can be easily
   * checked using logical operators
   */
  enum Coupling
  {
    SINGLET = 0x01,
    TRIPLET = 0x02,
    FULL = 0x03
  };

}


#endif // _EXCITONDEFS_H_

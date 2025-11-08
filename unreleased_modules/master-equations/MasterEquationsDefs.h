#ifndef _MASTEREQUATIONSDEFS_H_
#define _MASTEREQUATIONSDEFS_H_

//! A namespace for excitons specific definitions
/*!
 */
namespace MasterEquationsDefs
{

  //! Specifies the type of coupling used
  /*!
   * The values are hardcoded to hex values so they can be easily
   * checked using logical operators
   */

  enum Coupling
  {
    ELECTRONS = 0x01,
    HOLES = 0x02,
    BOTH = 0x03
  };


}


#endif // _MASTEREQUATIONSDEFS_H_

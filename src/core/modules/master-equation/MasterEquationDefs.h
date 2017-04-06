#ifndef _MASTEREQUATIONDEFS_H_
#define _MASTEREQUATIONDEFS_H_

//! A namespace for excitons specific definitions
/*!
 */
namespace MasterEquationDefs
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


#endif // _CONTINUOUSKINETICDEFS_H_

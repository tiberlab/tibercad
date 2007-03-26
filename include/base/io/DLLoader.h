// $Id$

#ifndef _DLLOADER_H_
#define _DLLOADER_H_

#include <string>

//! A wrapper class which tries to load a certain library
/*!
 * The given library is looked for in well defined paths. If it is found,
 * it is opened and well defined symbols are looked for.
 */
class DLLoader
{

  public:

    //! A structure which defines the library interface
    struct LibraryInterface
    {
      //! The library handle
      void* handle;

      //! The creation method
      void* create_fnc;

      //! The destroy method
      void* destroy_fnc;
    };

    
    //! Opens a library and gets the library interface
    /*!
     * \param name the name of the library (\em not the filename)
     * \param iface the LibraryInterface where the interface information will
     * be stored
     * \return true if the operation was successful
     */
    static bool open_library(const std::string& name, LibraryInterface& iface);

    
    //! Closes a library for the given handle
    static void close_library(void* handle);


  private:

    //! For static use only
    DLLoader(void);

    //! The library search path
    static std::string _libpath;

};


#endif // _DLLOADER_H_

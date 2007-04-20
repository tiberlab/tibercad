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

      LibraryInterface(void);
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


    //! Set the library search path
    static void set_library_path(const std::string& path);


    //! Prepend to library search path
    static void prepend_to_library_path(const std::string& path);



  private:

    //! For static use only
    DLLoader(void);

    //! The library search path
    static std::string _libpath;

};



//
// inline members
// 

inline
DLLoader::LibraryInterface::LibraryInterface(void)
  : handle(NULL),
    create_fnc(NULL),
    destroy_fnc(NULL)
{
}


#endif // _DLLOADER_H_

// $Id$

#ifndef _MATERIALBOUNDARY_H_
#define _MATERIALBOUNDARY_H_


#include "PhysicalObject.h"

class Material;


//! Description of a material boundary.
class MaterialBoundary : public PhysicalObject
{

  public:

    //! Destructor
    /*!
     * Deletes all \c PhysicalProperties objects
     */
    ~MaterialBoundary(void) {};


    //! Create a boundary between two given materials
    /*!
     * \param mat_A the material on one side of the boundary
     * \param mat_B the material on the other side of the boundary
     * \param options options for this boundary
     *
     * \c mat_A has to be non-NULL in any case.
     * \c mat_B is allowed to be NULL, which corresponds to an external
     * boundary
     */
    static MaterialBoundary* create(Material* mat_A, Material* mat_B,
        const ModelOptions& options);


    //! Return the material A
    Material* get_material_A(void) const;

    //! Return the material B
    Material* get_material_B(void) const;


    //! Return the name of component material A
    //const std::string& get_name_A(void) const;

    //! Return the name of component material B
    //const std::string& get_name_B(void) const;


  protected:

    //! Construct an  alloy material
    MaterialBoundary(Material* mat_A, Material* mat_B,
        const ModelOptions& options);


    //! \copydoc PhysicalObject::do_init()
    void do_init(void);


  private:

    //! The component A
    Material* _mat_A;


    //! The component B
    Material* _mat_B;

};



//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------


/*
inline
const std::string&
MaterialBoundary::get_name_A(void) const
{
  return _mat_A->get_name();
}


inline
const std::string&
MaterialBoundary::get_name_B(void) const
{
  return _mat_B->get_name();
}
*/


inline
Material*
MaterialBoundary::get_material_A(void) const
{
  return _mat_A;
}


inline
Material*
MaterialBoundary::get_material_B(void) const
{
  return _mat_B;
}



#endif /* _MATERIALBOUNDARY_H_ */

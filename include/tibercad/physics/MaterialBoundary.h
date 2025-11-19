// $Id$

#ifndef _MATERIALBOUNDARY_H_
#define _MATERIALBOUNDARY_H_


#include "tibercad/physics/PhysicalObject.h"

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
    static MaterialBoundary* create(ID id_A, Material* mat_A,
        ID id_B, Material* mat_B, const ModelOptions& options);


    //! Return the material A
    Material* get_material_A(void) const;

    //! Return the material B
    Material* get_material_B(void) const;


    //! Return the ID of material A
    ID get_id_A(void) const;

    //! Return the ID of material B
    ID get_id_B(void) const;


  protected:

    //! Construct an  alloy material
    MaterialBoundary(const ModelOptions& options);


    //! \copydoc PhysicalObject::do_init()
    void do_init(void);


  private:

    //! The ID of region A
    ID _id_A;


    //! The ID of region B
    ID _id_B;


    //! The component A
    Material* _mat_A;


    //! The component B
    Material* _mat_B;

};



//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------



inline
ID
MaterialBoundary::get_id_A(void) const
{
  return _id_A;
}


inline
ID
MaterialBoundary::get_id_B(void) const
{
  return _id_B;
}



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

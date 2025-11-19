// $Id$

#ifndef _PVMODULEBOUNDARYMODEL_H_
#define _PVMODULEBOUNDARYMODEL_H_

#include "tibercad/physics/PhysicalModel.h"
#include "tibercad/physics/MaterialBoundary.h"



class Elem;


//! This is the basic physical modle of PVModule
class PVModuleBoundaryModel : public PhysicalModel
{

  public:

    //! The type of contact
    enum ContactType
    {
      GND,  /*!< ground, will be node 0 */
      SRC,  /*!< Voltage source */
    };

    //! The layer of the contact
    enum ContactLayer
    {
      BOTTOM,  /*!< bottom layer */
      TOP,     /*!< top layer */
      BOTH,    /*!< contact on both layers*/
    };

    //! Destructor
    ~PVModuleBoundaryModel(void) {};

    //! Creator function
    static PVModuleBoundaryModel* create(const MaterialBoundary* boundary,
                                         const ModelOptions& options);

    //! Get contact type
    ContactType get_contact_type(void) const;

    //! Get contact layer
    ContactLayer get_contact_layer(void) const;

 

  protected:

    //! Constructor
    PVModuleBoundaryModel(const ModelOptions& options);

    virtual void do_init(void) final;


  private:

    //! The boundary type
    ContactType _contact_type {GND};

    //! The boundary layer
    ContactLayer _contact_layer {BOTTOM};


    //! The constructor method
    static TiberModelObject* _create(const ModelOptions& options, const void*);

    //! The destructor method
    static void _destroy(TiberModelObject* p);

};



inline
PVModuleBoundaryModel::ContactType
PVModuleBoundaryModel::get_contact_type(void) const
{
  return(_contact_type);
}

inline
PVModuleBoundaryModel::ContactLayer
PVModuleBoundaryModel::get_contact_layer(void) const
{
  return(_contact_layer);
}


#endif // _PVMODULEBOUNDARYMODEL_H_

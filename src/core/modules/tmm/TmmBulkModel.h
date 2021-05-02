// $Id: TmmBulkModel.h 4391 2017-04-07 11:16:58Z pecchia $

#ifndef _TMMBULKMODEL_H_
#define _TMMBULKMODEL_H_

#include "PhysicalModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"



//! This is the base class for the TMM bulk physical model
class TmmBulkModel : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~TmmBulkModel(void) {};

    //! Creator function
    static TmmBulkModel* create(const Material* mat, const ModelOptions& options);


    //! Get the relative permittivity at a given wavelength
    libMesh::Complex get_permittivity(double lambda) const;


  protected:

    //! Constructor
    TmmBulkModel(const ModelOptions& options);

    virtual void do_init(void) override;

    virtual void prepare_submodels(void) override;


  private:

    //! The constructor method
    static TiberModelObject* _create(const ModelOptions& options, const void*);

    //! The destructor method
    static void _destroy(TiberModelObject* p);

};




#endif // _TMMBULKMODEL_H_

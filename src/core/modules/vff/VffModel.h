#ifndef _VFFMODEL_H_
#define _VFFMODEL_H_


#include "PhysicalModel.h"
#include "tiber_dll.h"
#include "Messages.h"
#include "Database.h"



//! This is the base class for the Poisson physical model
class TBDLLOCAL VffModel : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~VffModel(void);


    //! Creator function
    static VffModel* create(const Material* mat, const ModelOptions& options);

    const double& get_alpha(void) const;

    const double& get_beta(void) const;

    const double& get_teta(void) const;

    const double& get_d(void) const;


  protected:

    //! Constructor
    VffModel(const ModelOptions& options);

    //! Read database
    void read_database(void);

    //! Initialize
    virtual void do_init(void);

    //virtual void prepare_submodels(void);

  private:


   static TiberModelObject* _create(const ModelOptions& options, const void*);

   static void  _destroy( TiberModelObject* p);

   double _c11;
   double _c12;
   double _c44;
   double _alpha;
   double _beta;
   double _a;
   double _teta;
   double _d;

};


inline
VffModel::VffModel(const ModelOptions& options) :
  PhysicalModel(options)
{
}

inline
VffModel::~VffModel()
{
}

inline
TiberModelObject*  VffModel::_create(const ModelOptions& options, const void*)
{

  return new VffModel(options);

}

inline
void  VffModel::_destroy( TiberModelObject* p)
{

  delete p;

}

inline
const double& VffModel::get_alpha(void) const
{
  return _alpha;
}

inline
const double& VffModel::get_beta(void) const
{
  return _beta;
}

inline
const double& VffModel::get_teta(void) const
{
  return _teta;
}

inline
const double& VffModel::get_d(void) const
{
  return _d;
}

#endif // _VFFMODEL_H_

#ifndef _VFFMODEL_H_
#define _VFFMODEL_H_


#include "PhysicalModel.h"
#include "tiber_dll.h"
#include "Messages.h"
#include "Database.h"
#include "Keating.h"
#include "Atom.h"

class Keating;


//! This is the base class for the Poisson physical model
class TBDLLOCAL VffModel : public PhysicalModel
{

public:

  //! Destructor
  virtual ~VffModel(void);


  //! Creator function
  static VffModel* create(const Material* mat, const ModelOptions& options);

  const double get_alpha(void) const;

  const double get_beta(void) const;

  const double get_teta(void) const;

  const double get_d(void) const;

  const double get_alpha(const Atom& atm1, const Atom& atm2) const;

  const double get_beta(const Atom& atm1, const Atom& atm2, const Atom& atm3) const;

  const double get_teta(const Atom& atm1, const Atom& atm2, const Atom& atm3) const;

  const double get_d(const Atom& atm1, const Atom& atm2) const;


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

  void prepare_submodels(void);

  double _c11;
  double _c12;
  double _c44;
  double _alpha;
  double _beta;
  double _a;
  double _teta;
  double _d;

  Keating* _keating;

  bool along_c(const Atom& atm1, const Atom& atm2) const;

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
const double VffModel::get_alpha(void) const
{
  //return _alpha;
  return _keating->get_alpha_0();
}

inline
const double VffModel::get_beta(void) const
{
  //return _beta;
  return _keating->get_beta_0();
}

inline
const double VffModel::get_teta(void) const
{
  return _keating->get_teta_0();
}

inline
const double VffModel::get_d(void) const
{
  return _keating->get_d_0();
}

#endif // _VFFMODEL_H_

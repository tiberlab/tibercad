/*
 * Negfmodel.h
 *
 *  Created on: Jan 18, 2012
 *      Author: fpalomba
 */

#ifndef NEGFMODEL_H_
#define NEGFMODEL_H_

#include "PhysicalModel.h"
#include "models/HamiltonianModel.h"

class TBDLLOCAL NegfModel : public PhysicalModel
{
  public:

    //! Destructor
    virtual ~NegfModel(void) {};

    //! Creator function
    static NegfModel* create(const Material* mat, const ModelOptions& options);

    const TensorValue<double>& get_inv_mass(ID id) const;

    unsigned int get_n_bands(void) const;

  protected:

    virtual void prepare_submodels(void);



  private:

    static TiberModelObject*  _create(const ModelOptions& options, const void*);

    static void  _destroy( TiberModelObject* p);

    NegfModel(const ModelOptions& options);

    std::vector<HamiltonianModel*> _ham_models;
};

inline
TiberModelObject*  NegfModel::_create(const ModelOptions& options, const void*)
{
  return new NegfModel(options);
}

inline
void  NegfModel::_destroy( TiberModelObject* p)
{
  delete p;
}


inline
const TensorValue<double>& NegfModel::get_inv_mass(ID id) const
{
   return _ham_models[id]->get_inv_mass();
}

inline
unsigned int NegfModel::get_n_bands(void) const
{
   return _ham_models.size();
}


#endif /* NEGFMODEL_H_ */

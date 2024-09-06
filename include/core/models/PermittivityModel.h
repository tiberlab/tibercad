// $Id$

#ifndef _PERMITTIVITYMODEL_H_
#define _PERMITTIVITYMODEL_H_

#include "PhysicalModel.h"
#include "Material.h"
#include "libMeshDefs.h"

#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"

USELIBMESHTYPE(RealTensor);


// Base class for charge density models
class TBDLEXPORT PermittivityModel : public PhysicalModel
{

  public:

  virtual ~PermittivityModel(void) {};
  
  const RealTensor& get_permittivity(void) const;
  
  virtual void calculate(const libMesh::Elem* elem, const libMesh::Point& point) = 0;
  
protected:
  
    PermittivityModel(const ModelOptions& options);

    void set_permittivity(const libMesh::RealVectorValue& permittivity);

  // void set_permittivity(const std::vector<double>& permittivity);
   
    void rotate(void);


  private:

   RealTensor _permittivity;

};


inline
PermittivityModel::PermittivityModel(const ModelOptions& options) :
  PhysicalModel(options)
{
  _permittivity = 0;
  _permittivity(0,0) = 1.0;
  _permittivity(1,1) = 1.0;
  _permittivity(2,2) = 1.0;
}

inline
const RealTensor&
PermittivityModel::get_permittivity() const
{
 
  return _permittivity;
}

inline
void
PermittivityModel::set_permittivity(const libMesh::RealVectorValue& permittivity_diag)
{
  
  _permittivity(0,0) =  permittivity_diag(0);
  _permittivity(1,1) =  permittivity_diag(1);
  _permittivity(2,2) =  permittivity_diag(2);
  
}

// inline
// void
// PermittivityModel::set_permittivity(const std::vector<double>& permittivity)
// {
//   _permittivity = 0;
//   _permittivity(0,0) =  permittivity[0];
//   _permittivity(1,1) =  permittivity[1];
//   _permittivity(2,2) =  permittivity[2];
  
// }

inline
void 
PermittivityModel::rotate(void)
{

  const Material* mat = get_material();
  const libMesh::RealTensor& rotm = mat->get_rotation_matrix();

  _permittivity = rotm * (_permittivity * rotm.transpose());

}

#endif // _POLARIZATIONMODEL_H_

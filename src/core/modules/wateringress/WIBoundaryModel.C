// $Id$

#include "WIBoundaryModel.h"
#include "MaterialBoundary.h"

using namespace std;

WIBoundaryModel*
WIBoundaryModel::create(const MaterialBoundary* boundary, const ModelOptions& options)
{
  std::string type = options.get_option("type", "pressure");
  WIBoundaryModel* mod = 
      PhysicalModel::create<WIBoundaryModel>("contact_" + type, boundary, options);

  if (mod == NULL)
  {
    ostringstream os;
    os << "water ingress boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;
}



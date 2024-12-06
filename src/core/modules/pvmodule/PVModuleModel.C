// $Id$

#include "PVModuleModel.h"
#include "Messages.h"

using namespace std;


PVModuleModel::PVModuleModel(const ModelOptions& options) :
  PhysicalModel(options)
{
}


TiberModelObject*
PVModuleModel::_create(const ModelOptions& options, const void*)
{
  return new PVModuleModel(options);
}


void
PVModuleModel::_destroy(TiberModelObject* p)
{
  delete p;
}


PVModuleModel*
PVModuleModel::create(const Material* mat, const ModelOptions& options)
{
  string type("default");
  options.get_option("type", type);

  PVModuleModel* pm = NULL;

  if (type == "default")
    // we create the default model from explicit creation method
    pm = PhysicalModel::create<PVModuleModel>(_create, _destroy, mat, options);
  else
  {
    // there is no such model, at the moment
    type = "pvstack_" + type;
    pm = PhysicalModel::create<PVModuleModel>(type, mat, options);
  }

  return(pm);
}



void
PVModuleModel::do_init(void)
{
  string type = get_option("region_type", "active");
  if (!has_option("region_type"))
  {
    // we try to guess it from the region name
    type = get_owner()->get_options().get_name();
  }
  if (type == "active") _region_type = ACTIVE;
  if (type == "P1") _region_type = P1;
  if (type == "P2") _region_type = P2;
  if (type == "P3") _region_type = P3;

  _top_rsheet = get_option("top_sheet_resistivity", _top_rsheet);
  _bottom_rsheet = get_option("bottom_sheet_resistivity", _bottom_rsheet);

  if (_region_type == P2)
    _connection_res = get_option("P2_resistivity", _connection_res);
}


void
PVModuleModel::prepare_submodels(void)
{
  
  //create_submodel(_charge_density_model, "charge_density");
  //create_submodels(_pm, "polarization");

}

void
PVModuleModel::do_print_info(void)
{
  Messages m;
  ostringstream os;
  os << "Type of region : ";
  switch (_region_type)
  {
    case ACTIVE:
      os << "active";
      break;

    case P1:
      os << "P1";
      break;

    case P2:
      os << "P2";
      break;

    case P3:
      os << "P3";
      break;

    default:
      os << "unknown";
      break;
  }
  os << "\n";

  if ((_region_type != P2) || (_region_type != P3))
  {
    os << "Top    sheet resistivity : " << _top_rsheet << " Ohms/sq\n";
    os << "Bottom sheet resistivity : " << _bottom_rsheet << " Ohms/sq\n";
  }
  m.info(os.str());
}




std::pair<double, double>
PVModuleModel::get_sheet_resistances(const Elem* elem,
                                     const Point& point) const
{
  return make_pair(_top_rsheet, _bottom_rsheet);
}



double
PVModuleModel::get_connection_resistance(const Elem* elem,
                                         const Point& point) const
{
  return _connection_res;
}

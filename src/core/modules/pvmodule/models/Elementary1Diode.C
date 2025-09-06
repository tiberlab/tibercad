#include "Elementary1Diode.h"
#include "Photocurrent.h"
#include "DegradationModel.h"

#include "TiberModule.h"


Elementary1Diode::Elementary1Diode(const ModelOptions& options) :
  ElementaryCell(options)
{
}


Elementary1Diode*
Elementary1Diode::create(const ModelOptions& options)
{
  Elementary1Diode* cd = new Elementary1Diode(options);

  return cd;
}


void
Elementary1Diode::prepare_submodels(void)
{
  create_submodel(_photocurr_model, "photocurrent");
  create_submodel(_degradation_model, "degradation");
}

void
Elementary1Diode::do_init(void)
{
  _rseries = get_option("series_resistance", _rseries);
  _rshunt = get_option("shunt_resistance", _rshunt);
  _photocurr = get_option("photocurrent", _photocurr);
  _isat = get_option("saturation_current", _isat);
  _ideality = get_option("ideality_factor", _ideality);
}


void
Elementary1Diode::do_write_netlist(unsigned int top_node, unsigned int bottom_node,
                                unsigned int& next_free,
                                double area,
                                const libMesh::Elem* elem,
                                const libMesh::Point& p,
                                std::ostream& os) const 
{
  double photocurr = _photocurr;
  if (_photocurr_model != nullptr)
    photocurr = _photocurr_model->get_photocurrent(elem, p);

  if (_degradation_model != nullptr)
  {
    DegradationModel::Parameters p;
  }

  photocurr *= area;

  os << "Rs" << top_node << " " << top_node << " " << next_free << " " << _rseries/area << "\n";
  os << "I" << top_node << " " << next_free << " " << bottom_node << " DC " << -photocurr << "\n";
  os << "D" << top_node << " " << next_free << " " << bottom_node << " DModel" << top_node << "\n";
  os << "Rsh" << top_node << " " << next_free << " " << bottom_node << " " << _rshunt/area << "\n";
  os << ".model DModel" << top_node << " D(IS=" << _isat*area << " N=" << _ideality << ")\n";

  ++next_free;
}

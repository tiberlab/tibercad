#include "Keating.h"
#include "Database.h"
#include "Messages.h"

#include <sstream>

using namespace std;


void
Keating::do_init(void)
{
  const Database& db = get_database();
  db.set_section("lattice");
  _a = db.get("a", 0.0, true);
  if (get_material()->get_structure() == "zb")
    {
      double d = _a * (sqrt(3.0) / 4.0) * 10.0;
      _d_0 = d; _d_1 = d;
      double teta = -1.0 / 3.0;
      _costeta_0 = teta; _costeta_1 = teta;
    }
  if (get_material()->get_structure() == "wz")
    {
      _c = db.get("c", 0.0, true);
      _u = db.get("u", 0.0, true);

      _d_1 = _c * _u * 10.0;
      double v = 1.0 - 2.0 * _u;
      double sq_3 = sqrt(3.0);
      _d_0 = (sqrt(3.0 * _c * _c * v * v + 4.0 * _a * _a) / (2.0 * sq_3)) * 10.0;
      //cerr << "d0 = " << _d_0 << " d1 = " << _d_1 << endl;

      _costeta_1 = (-1.0 * sq_3 * _c * v) / sqrt(3.0 * _c * _c * v * v + 4.0 * _a * _a);
      _costeta_0 = (3.0 * _c * _c * v * v - 2.0 * _a * _a) / (3.0 * _c * _c * v * v + 4.0 * _a * _a);
      //cerr << "costeta_0 = " << _costeta_0 << " costeta_1 = " << _costeta_1 << endl;
    }

}




void
Keating::do_print_info(void)
{
  if (get_material()->get_structure() == "zb")
  {
    ostringstream os;
    os << "Keating parameters: alpha = " <<
        get_alpha_0() << ", beta = " << get_beta_0();
    Messages::info(os.str());
  }
  else if (get_material()->get_structure() == "wz")
  {
    ostringstream os;
    os << "Keating parameters: alpha_0 = " <<
        get_alpha_0() << ", alpha_1 = " <<
        get_alpha_1() << ", beta_0 = " <<
        get_beta_0() << ", beta_1 = " << get_beta_1();
    Messages::info(os.str());
  }
}

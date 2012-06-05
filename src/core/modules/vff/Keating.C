#include "Keating.h"
#include "Database.h"

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
      double teta = -0.3333;
      _teta_0 = teta; _teta_1 = teta;
    }
  if (get_material()->get_structure() == "wz")
    {
      _c = db.get("c", 0.0, true);
      _u = db.get("u", 0.0, true);

      _d_1 = _c * _u * 10.0;
      double v = 1.0 - 2.0 * _u;
      double sq_3 = sqrt(3.0);
      _d_0 = (sqrt(3.0 * _c * _c * v * v + 4.0 * _a * _a) / (2.0 * sq_3)) * 10.0;

      _teta_1 = asin((-1.0 * sq_3 * _c * v) / sqrt(3.0 * _c * _c * v * v + 4.0 * _a * _a));
      _teta_0 = (3.0 * _c * _c * _a * _a - 2.0 * _a * _a) / (3.0 * _c * _c * _a * _a + 4.0 * _a * _a);
    }

}

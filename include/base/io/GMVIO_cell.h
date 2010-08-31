#ifndef _GMVIO_CELL_H
#define _GMVIO_CELL_H

#include <vector>

#include "gmv_io.h"
#include "tiber_dll.h"

class MeshBase;

class TBDLLOCAL GMVIO_cell : public GMVIO
{
 public:
  GMVIO_cell ( const MeshBase& mesh) : GMVIO(mesh){};

  void write_ascii_cell_data (const std::string& fname,
				const std::vector<Number>& soln,
				const std::vector<std::string>& names);
 private:
};

#endif // _GMVIO_CELL_H


 #include <iomanip>
 #include <fstream>
 #include <cstring> // for strcpy, memcpy
 #include <cstdio>  // for sprintf
 #include <vector>


 #include "libmesh_config.h"
 #include "gmv_io.h"
 #include "mesh_base.h"
 #include "elem.h"
 #include "equation_systems.h"

class GMVIO_cell : public GMVIO
{
 public:
  GMVIO_cell ( const MeshBase& mesh) : GMVIO(mesh){};

  void write_ascii_cell_data (const std::string& fname,
				const std::vector<Number>& soln,
				const std::vector<std::string>& names);
 private:
};


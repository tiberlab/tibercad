// $Id$

#include "tibercad/physics/BoundaryProperties.h"
#include "tibercad/geom/Boundary.h"
#include "tibercad/io/Messages.h"

using namespace std;


BoundaryProperties::~BoundaryProperties(void)
{
  Messages::debug("Deleted BoundaryProperties for boundary " +
      get_boundary()->get_name());
}



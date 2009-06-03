// $Id$

#include "BoundaryProperties.h"
#include "Boundary.h"
#include "Messages.h"

using namespace std;


BoundaryProperties::~BoundaryProperties(void)
{
  Messages::debug("Deleted BoundaryProperties for boundary " +
      get_boundary()->get_name());
}



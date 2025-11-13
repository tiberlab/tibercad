// $Id$

#include "XMonitor.h"
#include "petsc/TiberPetscUtils.h"


XMonitor::XMonitor(const std::string& title)
  : _title(title)
{
  PetscErrorCode ierr =
    PetscDrawCreate(PETSC_COMM_SELF, 0, title.c_str(), 0, 0, 500, 500, &_draw);
  TiberPetscUtils::checkerr(ierr);

  PetscDrawSetFromOptions(_draw);

  ierr = PetscDrawViewPortsCreate(_draw, 1, &_ports);
  TiberPetscUtils::checkerr(ierr);
  PetscDrawViewPortsSet(_ports, 0);

  PetscDrawLGCreate(_draw, 1, &_lg);
  // does not exist anymore
  //PetscDrawLGIndicateDataPoints(_lg);

  PetscDrawLGGetAxis(_lg, &_axis);
  PetscDrawAxisSetColors(_axis, PETSC_DRAW_BLACK, PETSC_DRAW_RED, PETSC_DRAW_BLUE);
}



XMonitor::~XMonitor(void)
{
  if (_draw != NULL)
  {
    PetscDrawFlush(_draw);
    PetscDrawViewPortsDestroy(_ports);
    PetscDrawLGDestroy(&_lg);
    PetscDrawDestroy(&_draw);
    _draw = NULL;
  }
}


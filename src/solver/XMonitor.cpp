/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file XMonitor.C
 * \brief Internal tiberCAD code.
 *
 * \internal
 */


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


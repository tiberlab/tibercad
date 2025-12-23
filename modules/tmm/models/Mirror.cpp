/*  
 * This file is part of the tiberCAD module tmm.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file Mirror.C
 * \brief tiberCAD tmm module implementation.
 *
 * \note This file is part of module tmm.
 */

/*
 * Mirror.cpp
 *
 *  Created on: 30 Sep 2021
 *      Author: pamiri
 */

#include "Mirror.h"

#include "tibercad/module/TiberModule.h"

void Mirror::do_init(void)
{
  write_type("Mirror");
  get_parameter("m00", _member00);
  get_parameter("m01", _member01);
  get_parameter("m10", _member10);
  get_parameter("m11", _member11);
}

void Mirror::Calculate_M_Matrix(void)
{
  set_elements(_member00, _member01, _member10, _member11);
}



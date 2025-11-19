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



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
 * \file Tensor4DSym.C
 * \brief tiberCAD API implementation.
 */

#include "tibercad/math/Tensor4DSym.h"
#include "tibercad/math/Tensor2.h"

Tensor4DSym::Tensor4DSym(double x)
{
  *this = x;
}

Tensor4DSym&
Tensor4DSym::operator=(double x)
{
  // 0.5*( delta_ik delta_jl + delta_jk delta_il )
  _data[0] = _data[1] = _data[2] = x;
  _data[3] = _data[4] = _data[5] = 0.5*x;
  for (unsigned int j = 6; j < 21; j++)
    _data[j] = 0.0;

  return *this;
}


void
Tensor4DSym::push_forward(const Tensor2& F)
{
  // calculate B_ijkl = F_ip F_jq F_kr F_ls A_pqrs

  Tensor4DSym B;

  for (unsigned int i = 1; i < 4; i++)
  {
    for (unsigned int j = 1; j <= i; j++)
    {
      int I;
      switch (i - j)
      {
      case 0:
        I = i - 1;
        break;
      case 1:
        I = i + 1;
        break;
      case 2:
        I = 5;
        break;
      }

      for (unsigned int k = 1; k < 4; k++)
      {
        for (unsigned int l = 1; l <= k; l++)
        {
          int J;
          switch (k - l)
          {
          case 0:
            J = k - 1;
            break;
          case 1:
            J = k + 1;
            break;
          case 2:
            J = 5;
            break;
          }

          if (J > I)
            continue;

          double value = 0.0;
          for (unsigned int p = 1; p < 4; p++)
            for (unsigned int q = 1; q < 4; q++)
              for (unsigned int r = 1; r < 4; r++)
                for (unsigned int s = 1; s < 4; s++)
                  value += F(i, p) * F(j, q) * F(k, r) * F(l, s) * (*this)(p, q, r, s);

          B(i, j, k, l) = value;
        }
      }
    }
  }

  *this = B;
}

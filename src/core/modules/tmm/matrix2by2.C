#include "matrix2by2.h"
#include<iostream>
#include <complex>


using namespace std;
matrix2by2::matrix2by2()
:_m00(0.0),_m01(0.0) ,_m10(0.0) ,_m11(0.0)
{
}

matrix2by2::~matrix2by2()
{
}

matrix2by2::matrix2by2(double a00, double a01, double a10, double a11)
:_m00(a00),_m01(a01) ,_m10(a10) ,_m11(a11)
{
}

void matrix2by2::print()
{
  cout << "m00 = " << _m00 << "  " << "m01 = "<< _m01 << "\r\n";
  cout << "m10 = " << _m10 << "  " << "m11 = "<< _m11 << "\r\n";
  cout << "                       " << "\r\n";

}

void matrix2by2::inv()
{
  complex<double> div (1,0);
  complex<double> m00;
  complex<double> m01;
  complex<double> m10;
  complex<double> m11;

  m00 = _m00;
  m01 = _m01;
  m10 = _m10;
  m11 = _m11;
  div = (m00 *m11 - m01 *m10);
  if (m11 == div)
    _m00 = 1;
  else
  _m00 = m11/div;

  if (m01 == div)
    _m01 = 1;
  else
    _m01 = -m01/div;;

  if (m10 == div)
    _m10 = 1;
  else
  _m10 = -m10/div;

  if (m00 == div)
    _m11 = 1;
  else
  _m11 = m00/div;
}

void matrix2by2::unit_matrix()
{
  _m00 = 1;
  _m01 = 0;
  _m10 = 0;
  _m11 = 1;
}

matrix2by2 matrix2by2::operator*(const matrix2by2 & old_matrix2by2)
{
  matrix2by2 new_Matrix_2by2;
  new_Matrix_2by2._m00 = (_m00 * old_matrix2by2._m00) + (_m01 * old_matrix2by2._m10);
  new_Matrix_2by2._m01 = (_m00 * old_matrix2by2._m01) + (_m01 * old_matrix2by2._m11);
  new_Matrix_2by2._m10 = (_m10 * old_matrix2by2._m00) + (_m11 * old_matrix2by2._m10);
  new_Matrix_2by2._m11 = (_m10 * old_matrix2by2._m01) + (_m11 * old_matrix2by2._m11);

  return(new_Matrix_2by2);
}



void matrix2by2::set(int elm , std::complex<double> a)
{
  switch(elm)
  {
    case 0:
      _m00 = a;
      break;
    case 1:
      _m01 = a;
      break;
    case 2:
      _m10 = a;
      break;
    case 3:
      _m11 = a;
      break;
  }
}

std::complex<double> matrix2by2::get(int elm)
{
  switch(elm)
  {
    case 0:
      return(_m00);
      break;
    case 1:
      return(_m01);
      break;
    case 2:
      return(_m10);
      break;
    case 3:
      return(_m11);
      break;
    default:
      return(0);
  }
}


#ifndef MATRIX2BY2_H
#define MATRIX2BY2_H
#include <complex>


class matrix2by2
{
	public:
		matrix2by2();
		matrix2by2(double, double,double,double);
		~matrix2by2();
		
		void inv();
		void unit_matrix();
		void print();
		matrix2by2 operator*(matrix2by2 const & old_matrix2by2);

		void set(int elm, std::complex<double> a);
		std::complex<double> get(int elm);
	
	
	private:
		std::complex<double> _m00;
		std::complex<double> _m01;
		std::complex<double> _m10;
		std::complex<double> _m11;
};
#endif
#ifndef TUTILS_H
#define TUTILS_H

// Tensor utilities

int QR_factor(double *A, const int *N, double *C, double *D);
int R_solve( const double *A, const int *M, const double *D, double *b );
int QR_solve(const double *A, const int *M, const double *c,
             const double *d, double *b);
int LU_Factor( double *A, const int *M, int *indx );
int LU_Solve( const double *A, const int *M, const int *indx, double *b );

void solveQR( double *A, const int *N, double *C, double *D, double *b );
void solveLU( double *A, const int *N, int *p, double *b );
void invFortranMatrix( double *m, int *N, int *p, double *hv );

#endif // TUTILS_H

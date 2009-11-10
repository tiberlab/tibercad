// Copyright (c) 2001 A.H. van den Boogaard

#ifndef XTENSOR_H
#define XTENSOR_H

#include "tensor.h"

Tensor2Sym solve( const Tensor4DSym &A, const Tensor2Sym& b );
Tensor2Sym solve( const Tensor4LSym &A, const Tensor2Sym& b );
Tensor2Sym deviatoric_solve( const Tensor4DSym &A, const Tensor2Sym& b );
Tensor2Sym deviatoric_solve( const Tensor4LSym &A, const Tensor2Sym& b );
Tensor4DSym inv( const Tensor4DSym& A );
Tensor4LSym inv( const Tensor4LSym& A );
Tensor4DSym push_forward( const Tensor4DSym& A, const Tensor2Gen& F );
Tensor4LSym push_forward( const Tensor4LSym& A, const Tensor2Gen& F );
Tensor2Gen expW( const Tensor2Gen& W );

#endif // XTENSOR_H

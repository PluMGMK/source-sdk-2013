// © 2015 Michael Keyes
// Function definitions for generalised Laguerre polynomials

#ifndef VIGOV_LAUGERRE_H
#define VIGOV_LAGUERRE_H

#include "tier0/basetypes.h" // For vec_t

vec_t NchooseR(int n, int r); // Binomial coefficient
vec_t GeneralisedLaguerreInt(int n, int alpha, vec_t x);
vec_t GeneralisedLaguerre(int n, vec_t alpha, vec_t x);

#endif //ndef VIGOV_LAGUERRE_H
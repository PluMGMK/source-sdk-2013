// © 2015 Michael Keyes
// Functions for generalised Laguerre polynomials

#include <math.h>
#include <float.h>	// Needed for FLT_EPSILON
#include "basetypes.h"
#include <memory.h>
#include "tier0/dbg.h"
#include "mathlib/mathlib.h"

// memdbgon must be the last include file in a .cpp file!!!
#include "tier0/memdbgon.h"

extern float s_flFactorials[]; // in spherical.cpp

vec_t NchooseR(int n, int r) {
	if( r > n )
		return -1;
	
	return s_flFactorials[n] / (s_flFactorials[n-r] * s_flFactorials[r]);
}

vec_t GeneralisedLaguerreInt(int n, int alpha, vec_t x) {
	vec_t result = 0.;
	
	for(int i=0; i<=n; i++) {
		vec_t increment = NchooseR(n+alpha, n-i) * powf(x, i) / s_flFactorials[i];
		if(i%2) //Odd
			result -= increment;
		else //Even
			result += increment;
	}
	
	return result;
}

vec_t GeneralisedLaguerre(int n, vec_t alpha, vec_t x) {
	if( floorf(alpha) == alpha ) // Alpha is an integer
		return GeneralisedLaguerreInt(n, (int)alpha, x);
	
	// Otherwise we must use a recurrence relation
	if( n == 0 )
		return 1;
	if( n == 1 )
		return 1 + alpha - x;
	
	int k = n - 1;
	return ( (2*k + 1 + alpha - x)*GeneralisedLaguerre(k, alpha, x) - (k + alpha)*GeneralisedLaguerre(k-1, alpha, x) ) / (vec_t)n;
}
// © 2014 Michael Keyes
// An extension of Vector2D to allow it to be used for complex number maths.
// Hopefully this can let me replace the GSL and play nicely with the rest of Source…

#ifndef VIGOV_COMPLEX_H
#define VIGOV_COMPLEX_H

#include "mathlib/vector2d.h"

//This can surely work as a macro:
#define COMPLEX_SQRT(n) ((n<0.0f)?ComplexNumber(0,sqrt(-n)):ComplexNumber(sqrt(n)))

class ComplexNumber : public Vector2D {
public:
	// Construction/destruction
	ComplexNumber(void);
	ComplexNumber(vec_t X); //Construct from real
	ComplexNumber(vec_t X, vec_t Y);
	ComplexNumber(vec_t r, vec_t theta, bool polar);
	ComplexNumber(const float *pFloat);
	// copy constructor
	ComplexNumber(const ComplexNumber &vOther);

	//Initialisation (polar form)
	void InitPolar( vec_t r, vec_t theta );

	// negate the components
	void	Negate(); 

	// equality
	bool operator==(const ComplexNumber& z) const;
	bool operator!=(const ComplexNumber& z) const;	
	// a complex can be equal to a real
	bool operator==(const vec_t X) const;
	bool operator!=(const vec_t X) const;	
	
	//arithmetic operations
	ComplexNumber&	operator+=(const ComplexNumber &z);			
	ComplexNumber&	operator-=(const ComplexNumber &z);		
	ComplexNumber&	operator*=(const ComplexNumber &z);
	ComplexNumber&	operator*=(float s);
	ComplexNumber&	operator/=(const ComplexNumber &z);
	ComplexNumber&	operator/=(float s);		

	// assignment
	ComplexNumber& operator=(const ComplexNumber &vOther);
	ComplexNumber& operator=(const vec_t &X);
	
#ifndef VECTOR_NO_SLOW_OPERATIONS
	// arithmetic 
	ComplexNumber	operator-(void) const;
				
	ComplexNumber	operator+(const ComplexNumber& z) const;	
	ComplexNumber	operator-(const ComplexNumber& z) const;	
	ComplexNumber	operator*(const ComplexNumber& z) const;
	ComplexNumber	operator/(const ComplexNumber& z) const;
	ComplexNumber	operator*(float fl) const;
	ComplexNumber	operator/(float fl) const;
#endif
};

// ComplexNumber operations
void ComplexAdd( const ComplexNumber& a, const ComplexNumber& b, ComplexNumber& result );
void ComplexSubtract( const ComplexNumber& a, const ComplexNumber& b, ComplexNumber& result );
void ComplexMultiply( const ComplexNumber& a, vec_t b, ComplexNumber& result );
void ComplexMultiply( const ComplexNumber& a, const ComplexNumber& b, ComplexNumber& result );
void ComplexDivide( const ComplexNumber& a, vec_t b, ComplexNumber& result );
void ComplexDivide( const ComplexNumber& a, const ComplexNumber& b, ComplexNumber& result );
void PolarForm( ComplexNumber& z, vec_t& r, vec_t& theta );

//-----------------------------------------------------------------------------
// constructors
//-----------------------------------------------------------------------------

inline ComplexNumber::ComplexNumber(void)									
{ 
#ifdef _DEBUG
	// Initialize to NAN to catch errors
	x = y = VEC_T_NAN;
#endif
}

inline ComplexNumber::ComplexNumber(vec_t X)						
{ 
	x = X; y = 0;
	Assert( IsValid() );
}

inline ComplexNumber::ComplexNumber(vec_t X, vec_t Y)						
{ 
	x = X; y = Y;
	Assert( IsValid() );
}

inline ComplexNumber::ComplexNumber(vec_t r, vec_t theta, bool polar)						
{ 
	if(polar) {
		InitPolar(r, theta);
	} else {
		x = r; y = theta; //Just treat it as normal.
		Assert( IsValid() );
	}
}

inline ComplexNumber::ComplexNumber(const float *pFloat)					
{
	Assert( pFloat );
	x = pFloat[0]; y = pFloat[1];	
	Assert( IsValid() );
}


//-----------------------------------------------------------------------------
// copy constructor
//-----------------------------------------------------------------------------

inline ComplexNumber::ComplexNumber(const ComplexNumber &vOther)					
{ 
	Assert( vOther.IsValid() );
	x = vOther.x; y = vOther.y;
}

//-----------------------------------------------------------------------------
// initialisation
//-----------------------------------------------------------------------------

inline void ComplexNumber::InitPolar( vec_t r, vec_t theta )
{ 
	FastSinCos(theta, &y, &x);
	x *= r;
	y *= r;
	Assert( IsValid() );
}

//-----------------------------------------------------------------------------
// assignment
//-----------------------------------------------------------------------------

inline ComplexNumber& ComplexNumber::operator=(const ComplexNumber &vOther)	
{
	Assert( vOther.IsValid() );
	x=vOther.x; y=vOther.y;
	return *this; 
}

inline ComplexNumber& ComplexNumber::operator=(const vec_t &X)	
{
	x=X; y=0;
	return *this; 
}

//-----------------------------------------------------------------------------
// comparison
//-----------------------------------------------------------------------------

inline bool ComplexNumber::operator==( const ComplexNumber& z ) const
{
	Assert( z.IsValid() && IsValid() );
	return (z.x == x) && (z.y == y);
}

inline bool ComplexNumber::operator!=( const ComplexNumber& z ) const
{
	Assert( z.IsValid() && IsValid() );
	return (z.x != x) || (z.y != y);
}

inline bool ComplexNumber::operator==( const vec_t X ) const
{
	Assert( IsValid() );
	return (y != 0.0f) && (X == x);
}

inline bool ComplexNumber::operator!=( const vec_t X ) const
{
	Assert( src.IsValid() && IsValid() );
	return (y == 0.0f) || (X != x);
}

//-----------------------------------------------------------------------------
// standard math operations
//-----------------------------------------------------------------------------

inline void ComplexNumber::Negate()
{ 
	Assert( IsValid() );
	x = -x; y = -y;
} 

inline ComplexNumber& ComplexNumber::operator+=(const ComplexNumber& v)	
{ 
	Assert( IsValid() && v.IsValid() );
	x+=v.x; y+=v.y;	
	return *this;
}

inline ComplexNumber& ComplexNumber::operator-=(const ComplexNumber& v)	
{ 
	Assert( IsValid() && v.IsValid() );
	x-=v.x; y-=v.y;	
	return *this;
}

inline ComplexNumber& ComplexNumber::operator*=(float fl)	
{
	x *= fl;
	y *= fl;
	Assert( IsValid() );
	return *this;
}

inline ComplexNumber& ComplexNumber::operator*=(const ComplexNumber& z)	
{ 
	//Standard binomial multiplication
	x = x * z.x - y * z.y; //Minus because i² is of course -1.
	y = x * z.y + y * z.x;
	Assert( IsValid() );
	return *this;
}

inline ComplexNumber& ComplexNumber::operator/=(float fl)	
{
	Assert( fl != 0.0f );
	float oofl = 1.0f / fl;
	x *= oofl;
	y *= oofl;
	Assert( IsValid() );
	return *this;
}

inline ComplexNumber& ComplexNumber::operator/=(const ComplexNumber& z)	
{ 
	Assert( v.x != 0.0f && v.y != 0.0f );
	//Multiply above and below by conjugate.
	//Above:
	x = x * z.x + y * z.y; //Plus because it's the conjugate (change y sign).
	y = x * (-z.y) + y * z.x;
	//Below:
	vec_t denominator = z.LengthSqr(); //Difference of two squares - veery clever… :)
	x /= denominator;
	y /= denominator;
	Assert( IsValid() );
	return *this;
}

inline void ComplexAdd( const ComplexNumber& a, const ComplexNumber& b, ComplexNumber& c )
{
	Assert( a.IsValid() && b.IsValid() );
	c.x = a.x + b.x;
	c.y = a.y + b.y;
}

inline void ComplexSubtract( const ComplexNumber& a, const ComplexNumber& b, ComplexNumber& c )
{
	Assert( a.IsValid() && b.IsValid() );
	c.x = a.x - b.x;
	c.y = a.y - b.y;
}

inline void ComplexMultiply( const ComplexNumber& a, vec_t b, ComplexNumber& c )
{
	Assert( a.IsValid() && IsFinite(b) );
	c.x = a.x * b;
	c.y = a.y * b;
}

inline void ComplexMultiply( const ComplexNumber& a, const ComplexNumber& b, ComplexNumber& c )
{
	Assert( a.IsValid() && IsFinite(b) );
	//Standard binomial multiplication
	c.x = a.x * b.x - a.y * b.y; //Minus because i² is of course -1.
	c.y = a.x * b.y + a.y * b.x;
}


inline void ComplexDivide( const ComplexNumber& a, vec_t b, ComplexNumber& c )
{
	Assert( a.IsValid() );
	Assert( b != 0.0f );
	vec_t oob = 1.0f / b;
	c.x = a.x * oob;
	c.y = a.y * oob;
}

inline void ComplexDivide( const ComplexNumber& a, const ComplexNumber& b, ComplexNumber& c )
{
	Assert( a.IsValid() );
	Assert( (b.x != 0.0f) && (b.y != 0.0f) );
	//Multiply above and below by conjugate.
	ComplexNumber conjugate(b.x, -b.y);
	//Above:
	ComplexNumber numerator;
	ComplexMultiply(a, conjugate, numerator);
	//Below:
	vec_t denominator = b.LengthSqr(); //Difference of two squares - veery clever… :)
	
	c = numerator / denominator;
}

#ifndef VECTOR_NO_SLOW_OPERATIONS

//-----------------------------------------------------------------------------
// arithmetic operations
//-----------------------------------------------------------------------------

inline ComplexNumber ComplexNumber::operator-(void) const
{ 
	return ComplexNumber(-x,-y);				
}

inline ComplexNumber ComplexNumber::operator+(const ComplexNumber& z) const	
{ 
	ComplexNumber res;
	ComplexAdd( *this, z, res );
	return res;	
}

inline ComplexNumber ComplexNumber::operator-(const ComplexNumber& z) const	
{ 
	ComplexNumber res;
	ComplexSubtract( *this, z, res );
	return res;	
}

inline ComplexNumber ComplexNumber::operator*(float fl) const	
{ 
	ComplexNumber res;
	ComplexMultiply( *this, fl, res );
	return res;	
}

inline ComplexNumber ComplexNumber::operator*(const ComplexNumber& z) const	
{ 
	ComplexNumber res;
	ComplexMultiply( *this, z, res );
	return res;	
}

inline ComplexNumber ComplexNumber::operator/(float fl) const	
{ 
	ComplexNumber res;
	ComplexDivide( *this, fl, res );
	return res;	
}

inline ComplexNumber ComplexNumber::operator/(const ComplexNumber& z) const	
{ 
	ComplexNumber res;
	ComplexDivide( *this, z, res );
	return res;	
}

inline ComplexNumber operator*(float fl, const ComplexNumber& z)	
{ 
	return z * fl; 
}

#endif //slow

//-----------------------------------------------------------------------------
// get polar form
//-----------------------------------------------------------------------------
inline void PolarForm( ComplexNumber& z, vec_t& r, vec_t& theta )
{
	Assert( z.IsValid() );
	vec_t l = z.Length();
	r = l;
	if (l != 0.0f)
	{
		vec_t arg = atan(z.y/z.x); //This gives us a quantity <pi/2, corresponding to 1st/4th quadrant.
		if(z.x<0.0f && z.y>0.0f) {
			//2nd quadrant
			arg += M_PI; //We have a minus value, so it can just be added to pi to rotate it back into the right quadrant.
		} else if(z.x<0.0f && z.y<0.0f) {
			//3rd quadrant
			arg -= M_PI; //We have a plus value, so add it to -pi to rotate it into the right quadrant.
		}
		theta = arg;
	}
	else
	{
		theta = 0.0f;
	}
}

#endif //ndef VIGOV_COMPLEX_H
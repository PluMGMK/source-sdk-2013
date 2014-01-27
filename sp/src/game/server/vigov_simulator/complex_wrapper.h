//////////////////////////////////////////
// Class to simplify GSL representation //
// of complex numbers                   //
//                             J Walshe //
//                             Mar 2013 //
//////////////////////////////////////////

#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>

//#include<iostream>
//#include<fstream>
//#include<string>
#include<cmath>

using namespace std;

class complex {
    private:
        gsl_complex z;
    public:
        // Constructors
        complex(double realin=0, double imin=0);
        complex(gsl_complex zin);
        
        void polar(double magin, double argin); // Set value using polar coordinates
        gsl_complex gsl() const; // Return the gsl_complex variable - changed to const to avoid errors with const operators in particle simulator
        
        // Return components
        double real();
        double im();
        double abs();
        double arg();
        
        // Operators - changed to const for particle simulator
        complex operator+(const complex& zin);
        complex operator-(const complex& zin);
        complex operator*(const complex& zin);
        complex operator/(const complex& zin);
		//Particle simulator - Real operator prototypes
        complex operator+(double in);
        complex operator-(double in);
        complex operator*(double in);
        complex operator/(double in);

};

// Member definitions ================================================================
complex::complex(double realin, double imin){ z = gsl_complex_rect(realin,imin); }
complex::complex(gsl_complex zin){z=zin;}

void complex::polar(double magin, double argin){
    z = gsl_complex_polar(magin,argin);
}
gsl_complex complex::gsl() const{return z;}

double complex::real(){return GSL_REAL(z);}
double complex::im(){return GSL_IMAG(z);}
double complex::abs(){return gsl_complex_abs(z);}
double complex::arg(){return gsl_complex_arg(z);}

complex complex::operator+ (const complex& zin){ return complex( gsl_complex_add(z,zin.gsl()) ); }
complex complex::operator- (const complex& zin){ return complex( gsl_complex_sub(z,zin.gsl()) ); }
complex complex::operator* (const complex& zin){ return complex( gsl_complex_mul(z,zin.gsl()) ); }
complex complex::operator/ (const complex& zin){ return complex( gsl_complex_div(z,zin.gsl()) ); }
//Particle simulator - Real operators
complex complex::operator+ (double in){ return complex( gsl_complex_add_real(z,in) ); }
complex complex::operator- (double in){ return complex( gsl_complex_sub_real(z,in) ); }
complex complex::operator* (double in){ return complex( gsl_complex_mul_real(z,in) ); }
complex complex::operator/ (double in){ return complex( gsl_complex_div_real(z,in) ); }
//Particle simulator - Real operators (other direction)
complex operator+ (double in, const complex& zin){ return complex( gsl_complex_add_real(zin.gsl(),in) ); }
complex operator- (double in, const complex& zin){ return complex( gsl_complex_sub_real(zin.gsl(),in) ); }
complex operator* (double in, const complex& zin){ return complex( gsl_complex_mul_real(zin.gsl(),in) ); }
complex operator/ (double in, const complex& zin){ return complex( gsl_complex_div_real(zin.gsl(),in) ); }

// Non-member functions =======================================================
complex conjugate(complex zin){ return complex(gsl_complex_conjugate(zin.gsl())); }
complex inverse(complex zin){ return complex(gsl_complex_inverse(zin.gsl())); }
complex sqrt(complex zin){ return complex(gsl_complex_sqrt(zin.gsl())); }
complex complex_sqrt(double xin){ return complex(gsl_complex_sqrt_real(xin)); } // Can't redefine sqrt(double)
complex pow(complex zin,complex power){ return complex(gsl_complex_pow(zin.gsl(),power.gsl())); }
complex pow(complex zin,double power){ return complex(gsl_complex_pow_real(zin.gsl(),power)); }
complex exp(complex zin){ return complex(gsl_complex_exp(zin.gsl())); }
complex ln(complex zin){ return complex(gsl_complex_log(zin.gsl())); }
complex log(complex zin, complex base){ return complex(gsl_complex_log_b(zin.gsl(),base.gsl())); }

// Comparison operators =======================================================
inline bool operator==(complex z1, complex z2){ return z1.real()==z2.real()&&z1.im()==z2.im(); }
inline bool operator!=(complex z1, complex z2){ return !operator==(z1,z2); }
inline bool operator<(complex z1, complex z2){ return z1.abs()<z2.abs(); }
inline bool operator>(complex z1, complex z2){ return operator<(z2,z1); }
inline bool operator<=(complex z1, complex z2){ return !operator>(z1,z2); }
inline bool operator>=(complex z1, complex z2){ return !operator<(z1,z2); }

//Particle simulator - commented out - not needed
/*ostream& operator<<(ostream& os, complex zin){
    string sign = zin.im()>=0 ? "+" : "-";
    os << zin.real() << sign << fabs(zin.im()) << "i";
    return os;
}*/

//=============================================================================









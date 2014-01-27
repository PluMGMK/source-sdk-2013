// This is part of the Subatomic Particle Simulator
// - Created by VIGoV Interactive
// - Copyright © 2013 Michael Keyes

#define CULL_SPRITES 1

#define MASS_SCALE 1e30
#define DISTANCE_SCALE 1e13
#define TIME_SCALE 1e16 //We're dealing with huge velocities here, which we want to scale down so something humanly interpretable is produced

#define WEAK_BOSON_LIFETIME 2.0f //VERY unrealistic, but the concept needs to be clearly demonstrated

//Important physics constants
#define PLANCK 6.62606957e-34 //h
#define DIRAC 1.054571726e-34 //h-bar - equal to h/(2*pi)
#define COULOMB 8.98755178e9 //ke
//#define ALPHA 7.29735257e-3 //Fine structure constant
#define WEAK_COUPLING_CONSTANT 5e-7 //g - its value is said to be between 10^-7 and 10^-6, so we'll go half-way!

//masses are scaled by an order of magnitude of 30, because VPhysics won't take anything under 0.1 or above 50,000
#define ELECTRON_MASS 9.10938291e-1 //e-31
#define AMU 1.660538921e+3 //e-27
#define PROTON_MASS 1.67262178e+3 //e-27
#define NEUTRON_MASS 1.674927351e+3 //e-27
#define WEAK_BOSON_MASS 1.435e+5 //e-25 - this may pose a problem... Weak bosons don't stick around very long though, so we'll probably be fine
#define NEUTRAL_PION_MASS 2.406176348e+2 //e-28

#define C_SQ 8.98755179e16
#define SPD_LGT 299792458.00f

//eV constants
#define ELEMENTARY_CHARGE 1.6e-19
#define EV_MASS 1.783e-36
#define EV_MOMENTUM 5.344286e-28

//Planck Units
#define PLANCK_MASS 2.17651e-8
#define PLANCK_ENERGY 1.9561e9
#define PLANCK_LENGTH 1.616199e-35
#define PLANCK_CHARGE 1.875545956e-18
#define PLANCK_MOMENTUM 6.52485f

#define PARTICLE_GLOW_SPRITE	"sprites/glow1.vmt"

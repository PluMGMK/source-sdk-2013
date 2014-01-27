// This is part of the Subatomic Particle Simulator
// - Created by VIGoV Interactive
// - Copyright © 2013 Michael Keyes

#include "cbase.h"
#include "vigov_sim_defs.h"

#if VIGOV_SIMULATOR
#define _USE_MATH_DEFINES //Required on Windows to get M_E constant

//Valve headers:
#include "Sprite.h"
#include "mathlib/mathlib.h"
#include "mathlib/spherical_geometry.h"
#include "utlvector.h"
#include "saverestore_utlvector.h"
#include "beam_shared.h"
#include "fmtstr.h"

//GNU Scientific Library:
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_complex_math.h>
//#include <gsl/gsl_sf_legendre.h>

//Obtained from http://joewalshe.net/2013/03/a-wrapper-for-complex-numbers-with-gsl/
#include "complex_wrapper.h"

// memdbgon must be the last include file in a .cpp file!!!
#include "tier0/memdbgon.h"

static double ALPHA = COULOMB * Square<double>(ELEMENTARY_CHARGE) / (DIRAC * SPD_LGT); //Fine structure constant
static double A0 = DISTANCE_SCALE * DIRAC / ((ELECTRON_MASS/MASS_SCALE) * SPD_LGT * ALPHA); //The Bohr radius in game units as defined in vigov_sim_defs.h
ConVar draw_bohr_radius( "draw_bohr_radius", "0" );
ConVar trace_orbitals( "trace_orbitals", "0" );
ConVar test_p_orbitals( "test_p_orbitals", "0", FCVAR_CHEAT );

ConVar debug_general( "debug_general", "0", FCVAR_CHEAT );
ConVar debug_wavefunctions( "debug_wavefunctions", "0", FCVAR_CHEAT );
ConVar debug_electrons( "debug_electrons", "0", FCVAR_CHEAT );
ConVar debug_nucleus( "debug_nucleus", "0", FCVAR_CHEAT );

//List of elements
struct Element {
	int m_iAtomicNo;
	float m_flRelMass;
	char *m_szSymbol;
	char *m_szName;
	
	bool NameIs(const char *szName) {
		return !Q_stricmp(szName, m_szName);
	}
	bool SymbolIs(const char *szSymbol) {
		return !Q_stricmp(szSymbol, m_szSymbol);
	}
	
	Element(int iAtomicNo, float flRelMass, char *szSymbol, char *szName) :
		m_iAtomicNo(iAtomicNo),
		m_flRelMass(flRelMass),
		m_szSymbol(szSymbol),
		m_szName(szName) {}
		
	Element(const Element &other) :
		m_iAtomicNo(other.m_iAtomicNo),
		m_flRelMass(other.m_flRelMass),
		m_szSymbol(other.m_szSymbol),
		m_szName(other.m_szName) {}
};

//Info from http://www.webelements.com/
static Element g_pElements [118] = {Element(1, 1.008f, "H", "Hydrogen"),
										Element(2, 4.0026f, "He", "Helium"),
										Element(3, 6.94f, "Li", "Lithium"),
										Element(4, 9.0122f, "Be", "Beryllium"),
										Element(5, 10.81f, "B", "Boron"),
										Element(6, 12.011f, "C", "Carbon"),
										Element(7, 14.007f, "N", "Nitrogen"),
										Element(8, 15.999f, "O", "Oxygen"),
										Element(9, 18.998f, "F", "Fluorine"),
										Element(10, 20.18f, "Ne", "Neon"),
										Element(11, 22.99f, "Na", "Sodium"),
										Element(12, 24.305f, "Mg", "Magnesium"),
										Element(13, 26.982f, "Al", "Aluminium"),
										Element(14, 28.085f, "Si", "Silicon"),
										Element(15, 30.974f, "P", "Phosphorus"),
										Element(16, 32.06f, "S", "Sulphur"),
										Element(17, 35.45f, "Cl", "Chlorine"),
										Element(18, 39.948f, "Ar", "Argon"),
										Element(19, 39.098f, "K", "Potassium"),
										Element(20, 40.078f, "Ca", "Calcium"),
										Element(21, 44.956f, "Sc", "Scandium"),
										Element(22, 47.867f, "Ti", "Titanium"),
										Element(23, 50.942f, "V", "Vanadium"),
										Element(24, 51.996f, "Cr", "Chromium"),
										Element(25, 54.938f, "Mn", "Manganese"),
										Element(26, 55.845f, "Fe", "Iron"),
										Element(27, 58.933f, "Co", "Cobalt"),
										Element(28, 58.693f, "Ni", "Nickel"),
										Element(29, 63.546f, "Cu", "Copper"),
										Element(30, 65.38f, "Zn", "Zinc"),
										Element(31, 69.723f, "Ga", "Gallium"),
										Element(32, 72.63f, "Ge", "Germanium"),
										Element(33, 74.922f, "As", "Arsenic"),
										Element(34, 78.96f, "Se", "Selenium"),
										Element(35, 79.904f, "Br", "Bromine"),
										Element(36, 83.798f, "Kr", "Krypton"),
										Element(37, 85.438f, "Rb", "Rubidium"),
										Element(38, 87.62f, "Sr", "Strontium"),
										Element(39, 88.906f, "Y", "Yttrium"),
										Element(40, 91.224f, "Zr", "Zirconium"),
										Element(41, 92.906f, "Nb", "Niobium"),
										Element(42, 95.96f, "Mo", "Molybdenum"),
										Element(43, 97.91f, "Tc", "Technetium"),
										Element(44, 101.07f, "Ru", "Ruthenium"),
										Element(45, 102.91f, "Rh", "Rhodium"),
										Element(46, 106.42f, "Pd", "Palladium"),
										Element(47, 107.87f, "Ag", "Silver"),
										Element(48, 112.41f, "Cd", "Cadmium"),
										Element(49, 114.82f, "In", "Indium"),
										Element(50, 118.71f, "Sn", "Tin"),
										Element(51, 121.76f, "Sb", "Antimony"),
										Element(52, 127.6f, "Te", "Tellurium"),
										Element(53, 126.9f, "I", "Iodine"),
										Element(54, 131.29f, "Xe", "Xenon"),
										Element(55, 132.91f, "Cs", "Caesium"),
										Element(56, 137.33f, "Ba", "Barium"),
										Element(57, 138.91f, "La", "Lanthanum"),
										Element(58, 140.12f, "Ce", "Cerium"),
										Element(59, 140.91f, "Pr", "Praseodymium"),
										Element(60, 144.24f, "Nd", "Neodymium"),
										Element(61, 144.91f, "Pm", "Promethium"),
										Element(62, 150.36f, "Sm", "Samarium"),
										Element(63, 151.96f, "Eu", "Europium"),
										Element(64, 157.25f, "Gd", "Gadolinium"),
										Element(65, 158.93f, "Tb", "Terbium"),
										Element(66, 162.5f, "Dy", "Dysprosium"),
										Element(67, 164.93f, "Ho", "Holmium"),
										Element(68, 167.26f, "Er", "Erbium"),
										Element(69, 168.93f, "Tm", "Thulium"),
										Element(70, 173.05f, "Yb", "Ytterbium"),
										Element(71, 174.97f, "Lu", "Lutetium"),
										Element(72, 178.49f, "Hf", "Hafnium"),
										Element(73, 180.95f, "Ta", "Tantalum"),
										Element(74, 183.84f, "W", "Tungsten"),
										Element(75, 186.21f, "Re", "Rhenium"),
										Element(76, 190.23f, "Os", "Osmium"),
										Element(77, 192.22f, "Ir", "Iridium"),
										Element(78, 195.08f, "Pt", "Platinum"),
										Element(79, 196.97f, "Au", "Gold"),
										Element(80, 200.59f, "Hg", "Mercury"),
										Element(81, 204.38f, "Tl", "Thallium"),
										Element(82, 207.2f, "Pb", "Lead"),
										Element(83, 208.98f, "Bi", "Bismuth"),
										Element(84, 208.98f, "Po", "Polonium"),
										Element(85, 209.99f, "At", "Astatine"),
										Element(86, 222.02f, "Rn", "Radon"),
										Element(87, 223.02f, "Fr", "Francium"),
										Element(88, 226.03f, "Ra", "Radium"),
										Element(89, 227.03f, "Ac", "Actinium"),
										Element(90, 232.04f, "Th", "Thorium"),
										Element(91, 231.04f, "Pa", "Protactinium"),
										Element(92, 238.03f, "U", "Uranium"),
										Element(93, 237.05f, "Np", "Neptunium"),
										Element(94, 244.06f, "Pu", "Plutonium"),
										Element(95, 243.06f, "Am", "Americium"),
										Element(96, 247.07f, "Cm", "Curium"),
										Element(97, 247.07f, "Bk", "Berkelium"),
										Element(98, 251.08f, "Cf", "Californium"),
										Element(99, 252.08f, "Es", "Einsteinium"),
										Element(100, 257.1f, "Fm", "Fermium"),
										Element(101, 258.1f, "Md", "Mendelevium"),
										Element(102, 259.1f, "No", "Nobelium"),
										Element(103, 262.11f, "Lr", "Lawrencium"),
										Element(104, 265.12f, "Rf", "Rutherfordium"),
										Element(105, 268.13f, "Db", "Dubnium"),
										Element(106, 271.13f, "Sg", "Seaborgium"),
										Element(107, 270, "Bh", "Bohrium"),
										Element(108, 277.15f, "Hs", "Hassium"),
										Element(109, 276.15f, "Mt", "Meitnerium"),
										Element(110, 281.16f, "Ds", "Darmstadium"),
										Element(111, 280.16f, "Rg", "Roentgenium"),
										Element(112, 285.17f, "Cn", "Copernicum"),
										Element(113, 284.18f, "Uut", "Ununtrium"),
										Element(114, 289.19f, "Fl", "Flerovium"),
										Element(115, 288.19f, "Uup", "Ununpentium"),
										Element(116, 293, "Lv", "Livermorium"),
										Element(117, 294, "Uus", "Ununseptium"),
										Element(118, 294, "Uuo", "Ununoctium")
};

enum ColourCharges {
	CHARGE_ANTIBLUE = -3,
	CHARGE_ANTIGREEN = -2,
	CHARGE_ANTIRED = -1,
	CHARGE_NEUTRAL,
	CHARGE_RED,
	CHARGE_GREEN,
	CHARGE_BLUE
};
static Vector g_vecNuclearOrigin;

class CBaseParticle : public CBaseAnimating {
	DECLARE_CLASS(CBaseParticle, CBaseAnimating);
	
public:
	
	void ParticleThink();

	int m_iPrincipal; 				//quantum number n - quantised energy
	int m_iAzimuthal; 				//quantum number l - allowed values for angular momentum
	int m_iMagnetic;					//quantum number ml - quantised angular momentum
	float m_flSpinProjection;		//quantum number ms - direction of angular momentum
	
	bool m_bFree; //Is this a free particle? Used to make beta particles ignore quantised energy from a time
	float m_flFreeUntil; //Time to become quantised again - not really physically correct, but we'll use it for the moment and see how it goes
	void SetFreeFor(float flTime) {m_flFreeUntil = gpGlobals->curtime + flTime; m_bFree = true;}
	
	//Energies:
	float m_flPotential;
	float m_flKinetic;
	float m_flTotalEnergy;
	float m_flRestEnergy;
	
	//Momenta:
	Vector m_vecMomentum;
	Vector m_vecAngularMomentum;
	
	float m_flAppropriateMass; //Mass of the particle in question (kg) - scaled
	float m_flCharge; //Electric charge (C)
	float m_flWavelength; //Wavelength of the wave-function
	
	bool m_bAntiParticle; //Very important...
	bool m_bUsesVPhysics;
	bool m_bAnnihilating;
	
	//Can't be wasting time every think for these
	bool m_bScaleUpToDate;
	bool m_bColourUpToDate;
	
	bool m_bStrongInt; //Is this a strong-interacting particle?
	int m_iColour;
	
	CHandle<CSprite> m_hSprite; //Handle to our sprite
	
	DECLARE_DATADESC();
	
	CBaseParticle();
	~CBaseParticle();
	
	void Precache( void );
	void Spawn( void );
	void UpdateOnRemove();
	virtual bool CreateVPhysics();
	virtual void ParticleUpdateSprite();
	virtual void ApplyParticlePhysics();
	virtual bool AcceptEnergy(float flEnergy, float flTolerancePortion);
	
	virtual float GetCharge(float *pflChargeInUnits = NULL);
	float GetParticleMass(double *pdKG = NULL, double *pdEV = NULL, double *pdPlanck = NULL);
	double ReducedMass();
	
	void Annihilate(){m_bAnnihilating = true;}
	
	//Wavefunction components
	complex g(int n, int k, double r, float flPotential);
	complex f(int n, int k, double r, float flPotential);
	
	complex Y(double a, double b, Vector &vecDirection); //Spherical harmonic function
	gsl_complex Psi(Vector vecLocation); //Particle wavefunction
	double Psi2(Vector vecLocation); //|Psi|^2
	double DiagnosticIntegral(int iStart, int iEnd, int iCoarseness);
	
	virtual bool IsElectron(){return false;}
	virtual bool IsProton(){return false;}
	virtual bool IsQuark(){return false;}
};

static CUtlVector<CBaseParticle *> g_vpParticles; //Just use normal pointers because this isn't saved/loaded anyway.
static float g_flNextAufbauTime = 0; //Next time to go through the list of particles and figure out their ground-states

static void Simulate ( const CCommand &args )
{
	if(args.ArgC() < 1) {
		Msg("Usage: simulate_element <name or symbol> <mass no. (optional)> <charge/oxidation number (optional)>\n");
	} else {
		int Z = 0;
		int A = 0;
		for(int i=0; i<( sizeof(g_pElements)/sizeof(g_pElements[0]) ); i++) {
			if(g_pElements[i].NameIs(args.Arg(1)) || g_pElements[i].SymbolIs(args.Arg(1))) {
				Z = g_pElements[i].m_iAtomicNo;
				A = (args.ArgC() > 2) ? atoi(args.Arg(2)) : (int)(g_pElements[i].m_flRelMass + 0.5f); //This rounds up above .5, and down below .5
			}
		}
		if(Z == 0) {
			Msg("Invalid element specified\n");
		} else if (A < Z) {
			Msg("Invalid mass number specified: %i\n", A);
		} else {
			int iElectrons = Z;
			if(args.ArgC() > 3) {
				iElectrons -= atoi(args.Arg(2)); //Take away the oxidation number if specified
			}
			
			//Get rid of the old particles...
			int iParticleCount = g_vpParticles.Count();
			for(int j=0; j<iParticleCount; j++) {
				if(!g_vpParticles[j]->IsQuark()){
					Msg("Deleting %i (particle %i of %i)\n", g_vpParticles[j]->entindex(), j+1, iParticleCount);
					UTIL_Remove(g_vpParticles[j]);
					if(g_vpParticles[j]->m_hSprite) {
						UTIL_Remove(g_vpParticles[j]->m_hSprite);
					}
				} else {
					Msg("Not deleting %i (particle %i of %i)\n", g_vpParticles[j]->entindex(), j+1, iParticleCount);
				}
			}
			g_vecNuclearOrigin = vec3_origin;
			
			//Spawn the new particles!
			int k=0;
			int iAngleIncrement = 360;
			int iRadius = 0;
			while(k<A) {
				//Arrange them in a sphere - adapted from code proposed at http://uclue.com/?xq=1288
				for(float angXY=-180; angXY<180; angXY+=iAngleIncrement) {
					for(float angYZ=-180; angYZ<180; angYZ+=iAngleIncrement) {
						QAngle angEuler = QAngle(angXY, angYZ, 0);
						Vector vecPosition = vec3_origin;
						AngleVectors(angEuler, &vecPosition);
						vecPosition *= iRadius;
						CBaseParticle *pParticle = (CBaseParticle *)CreateEntityByName(k>(Z-1)?"particle_neutron":"particle_proton"); //Protons first, then neutrons
						pParticle->SetAbsOrigin(vecPosition);
						pParticle->Spawn();
						
						k++;
						if(k>=A){break;}
					}
					if(k>=A){break;}
				}
				iAngleIncrement /= 3;
				iRadius += 128;
			}
			for(int l=0; l<iElectrons; l++) {
				Vector vecPosition = Vector(5,0,-50*iElectrons+1+50*l); //This arranges them in a vertical line centred on the zero-point, but offset so they're not inside nucleons
				CBaseParticle *pParticle = (CBaseParticle *)CreateEntityByName("particle_electron");
				pParticle->SetAbsOrigin(vecPosition);
				pParticle->Spawn();
			}
		}
	}
}

static int ElementAutoComp ( char const *partial, 
char commands[ COMMAND_COMPLETION_MAXITEMS ][ COMMAND_COMPLETION_ITEM_LENGTH ] )
{
	//This code is adapted from the ent_fire completion code, found in baseentity.cpp starting from line 5385.s
	const char *cmdname = "simulate_element";

	char *substring = (char *)partial;
	if ( Q_strstr( partial, cmdname ) )
	{
		substring = (char *)partial + strlen( cmdname ) + 1;
	}

	int checklen = 0;
	char *space = Q_strstr( substring, " " );
	/*if ( space )
	{
		//A name or symbol has already been specified
		return 0;
	}
	else*/
	{
		checklen = Q_strlen( substring );
	}
	
	CUtlVector< char * > symbols;

	int iElementIndex = -1;
	int iElementCount = sizeof(g_pElements)/sizeof(g_pElements[0]);
	while ( iElementIndex<iElementCount )
	{
		iElementIndex++;

		// Check target name against partial string
		if ( Q_strnicmp( g_pElements[iElementIndex].m_szName, substring, checklen ) )
			continue;

		symbols.AddToTail( g_pElements[iElementIndex].m_szName );
		
		// Too many
		if ( symbols.Count() >= COMMAND_COMPLETION_MAXITEMS )
			break;
	}

	// Now fill in the results
	for ( int i = 0; i < symbols.Count(); i++ )
	{
		const char *name = symbols[ i ];

		char buf[ 512 ];
		Q_strncpy( buf, name, sizeof( buf ) );

		CUtlString command;
		command = CFmtStr( "%s %s", cmdname, buf );
		Q_strncpy( commands[i], command.Get(), min(strlen( command.Get() ) + 1, COMMAND_COMPLETION_ITEM_LENGTH) );
	}

	return symbols.Count();
}

static ConCommand simulate_element("simulate_element", Simulate, "Simulate a named element (or give its symbol) with optional control of mass no. and charge", 0, ElementAutoComp);

static void Aufbau(void) {
	g_flNextAufbauTime = gpGlobals->curtime + 5.0f; //Run this every five seconds to take into account new particles added to the system
	
	int n, l, ml, i;
	n=1; l=ml=i=0; //Initialisation
	float ms = 0.5f; //Spin is a half-integer, so this is awkward!
	
	//Diagnostic: This starts with p orbitals instead:
	if(test_p_orbitals.GetBool()) {
		n = 2;
		l = 1;
	}
	
	CBaseParticle *pParticle;
	while(i < g_vpParticles.Count()) {
		pParticle=g_vpParticles[i];
		if(pParticle->IsElectron()) {
			//Assign the next available set of quantum numbers
			pParticle->m_iPrincipal = n;
			pParticle->m_iAzimuthal = l;
			pParticle->m_iMagnetic = ml;
			pParticle->m_flSpinProjection = ms;
			if(debug_electrons.GetBool()) {
				Msg("Quantum numbers on electron %i are (%i, %i, %i, %.1f)\n", pParticle->entindex(), n, l, ml, ms);
			}
			
			//"Roll over" quantum numbers as necessary
			ml++;
			if(ml>l) { //-l<=ml<=l
				if(ms > -0.5f) {
					ms -= 1.0f; //Opposite spin (obeys Hund's Rule of Maximum Multiplicity)
				} else {
					l++; //Next Azimuthal (sublevel)
				}
				ml=-l;
				if(l==n) { //0<=l<n
					n++;
					l=0;
				}
			}
		}
		i++; //Next particle
	}
	
	n=1; l=ml=i=0; ms=0.5f; //Reset
	while(i < g_vpParticles.Count()) {
		pParticle=g_vpParticles[i];
		if(pParticle->IsProton()) {
			//Assign the next available set of quantum numbers
			pParticle->m_iPrincipal = n;
			pParticle->m_iAzimuthal = l;
			pParticle->m_iMagnetic = ml;
			pParticle->m_flSpinProjection = ms;
			if(debug_nucleus.GetBool()) {
				Msg("Quantum numbers on proton %i are (%i, %i, %i, %.1f)\n", pParticle->entindex(), n, l, ml, ms);
			}
			
			//"Roll over" quantum numbers as necessary
			ml++;
			if(ml>l) { //-l<=ml<=l
				if(ms > 0.5f) {
					ms -= 1.0f; //Opposite spin (obeys Hund's Rule of Maximum Multiplicity)
				} else {
					l++; //Next Azimuthal (sublevel)
				}
				ml=-l;
				if(l==n) { //0<=l<n
					n++;
					l=0;
				}
			}
		}
		i++; //Next particle
	}
}

BEGIN_DATADESC(CBaseParticle)
	DEFINE_THINKFUNC(ParticleThink),
	
	DEFINE_FIELD(m_hSprite, FIELD_EHANDLE),
	
	DEFINE_FIELD(m_vecMomentum, FIELD_VECTOR),
	DEFINE_FIELD(m_vecAngularMomentum, FIELD_VECTOR),
	
	DEFINE_FIELD(m_iPrincipal, FIELD_INTEGER),
	DEFINE_FIELD(m_iAzimuthal, FIELD_INTEGER),
	DEFINE_FIELD(m_iMagnetic, FIELD_INTEGER),
	DEFINE_FIELD(m_flSpinProjection, FIELD_FLOAT),
	
	DEFINE_FIELD(m_flPotential, FIELD_FLOAT),
	DEFINE_FIELD(m_flKinetic, FIELD_FLOAT),
	DEFINE_FIELD(m_flRestEnergy, FIELD_FLOAT),
	DEFINE_FIELD(m_flAppropriateMass, FIELD_FLOAT),
	DEFINE_FIELD(m_flCharge, FIELD_FLOAT),	
	DEFINE_FIELD(m_flWavelength, FIELD_FLOAT),
	
	DEFINE_FIELD(m_bAntiParticle, FIELD_BOOLEAN),
	DEFINE_FIELD(m_bAnnihilating, FIELD_BOOLEAN),
	DEFINE_FIELD(m_bUsesVPhysics, FIELD_BOOLEAN),
	
	DEFINE_FIELD(m_bScaleUpToDate, FIELD_BOOLEAN),
	DEFINE_FIELD(m_bColourUpToDate, FIELD_BOOLEAN),
	
	DEFINE_FIELD(m_bStrongInt, FIELD_BOOLEAN),
	DEFINE_FIELD(m_iColour, FIELD_INTEGER),
END_DATADESC()

CBaseParticle::CBaseParticle() {
	m_bFree = false;
	m_flFreeUntil = 0.0f;
	m_hSprite = NULL;
	m_flPotential = 0;
	m_flKinetic = 0;
	m_flRestEnergy = 0;
	m_flTotalEnergy = 0;
	m_vecMomentum = vec3_origin;
	m_vecAngularMomentum = vec3_origin;
	m_iPrincipal = 1;
	m_iAzimuthal = 0;
	m_iMagnetic = 0;
	m_flSpinProjection = -0.5;
	m_flAppropriateMass = ELECTRON_MASS;
	m_flCharge = 0;
	m_flWavelength = 0;
	m_bAntiParticle = false;
	m_bAnnihilating = false;
	m_bUsesVPhysics = false;
	m_bScaleUpToDate = false;
	m_bColourUpToDate = false;
	m_bStrongInt = false;
	m_iColour = CHARGE_NEUTRAL;
	
	g_vpParticles.AddToHead(this);
}

CBaseParticle::~CBaseParticle() {
	g_vpParticles.FindAndRemove(this);
}

static double TotalNuclearEnergy(int iProtons = -1); //Prototype

double CBaseParticle::ReducedMass() {
	double dIndividualMass = 0;
	GetParticleMass(&dIndividualMass);
	double dNuclearMass = TotalNuclearEnergy() / C_SQ; //Treat the nucleus as a single body, with E=mc^2
	return (dIndividualMass*dNuclearMass) / (dIndividualMass + dNuclearMass);
}

void CBaseParticle::UpdateOnRemove() {
	g_vpParticles.FindAndRemove(this);
	if(m_bUsesVPhysics) {
		VPhysicsDestroyObject();
	}
	BaseClass::UpdateOnRemove();
}

void CBaseParticle::Precache(void) {
	PrecacheModel( PARTICLE_GLOW_SPRITE );
	PrecacheModel( "models/roller.mdl" );
	
	BaseClass::Precache();
}

void CBaseParticle::Spawn() {
	Precache();
	
	BaseClass::Spawn();
	
	m_bUsesVPhysics = CreateVPhysics();

	RegisterThinkContext( "ParticleContext" ); //This allows us to run a think function for the particle
	SetContextThink( &CBaseParticle::ParticleThink, gpGlobals->curtime, "ParticleContext" );
}

static float g_flNextNuclearOriginTime = 0; //Next time to go through the list of nucleons
static void RecalculateNuclearOrigin(); //Prototype

void CBaseParticle::ParticleThink() {
	if( g_flNextAufbauTime < gpGlobals->curtime ) {
		Aufbau();
	}
	if(g_flNextNuclearOriginTime < gpGlobals->curtime) {
		RecalculateNuclearOrigin();
	}
	
	//These are virtual, so derived classes don't need to define their own think function.
	ParticleUpdateSprite();
	ApplyParticlePhysics();
	
	SetNextThink(gpGlobals->curtime + 0.04f, "ParticleContext"); //25 times a second!
}

//Basically the same as the prop_sphere, only smaller
bool CBaseParticle::CreateVPhysics()
{
	SetSolid( SOLID_BBOX );
	SetCollisionBounds( -Vector(1,1,1), Vector(1,1,1) );
	objectparams_t params = g_PhysDefaultObjectParams;
	params.mass = m_flAppropriateMass;
	params.pGameData = static_cast<void *>(this);
	
	IPhysicsObject *pPhysicsObject = physenv->CreateSphereObject( 1, 0, GetAbsOrigin(), GetAbsAngles(), &params, false );
	if ( pPhysicsObject )
	{
		VPhysicsSetObject( pPhysicsObject );
		SetMoveType( MOVETYPE_VPHYSICS );
		pPhysicsObject->Wake();
		
		pPhysicsObject->EnableGravity( false );
		pPhysicsObject->EnableDrag( false );
	} else {
		return false;
	}
	
	pPhysicsObject->ApplyForceCenter(Vector(10,10,0)); //Start us off with some momentum, because we need it to solve the Schrödinger equation

	return true;
}

void CBaseParticle::ParticleUpdateSprite() {
	// Must have a valid sprite to affect
	if ( !m_hSprite )
	{
		// Create our sprite
		m_hSprite = CSprite::SpriteCreate( PARTICLE_GLOW_SPRITE, GetAbsOrigin(), true );
		if ( !m_hSprite )
			return;

		m_hSprite->SetTransparency( kRenderWorldGlow, 255, 0, 0, 255, kRenderFxNoDissipation );
		m_hSprite->SetParent(this);
		
		m_bScaleUpToDate = false;
		m_bColourUpToDate = false;
	}
	
	if( m_hSprite->GetLocalOrigin().LengthSqr() > 4 ) {
		//More than 2 units out of position
		m_hSprite->SetLocalOrigin(vec3_origin); //Put it back!
	}
	
	if( !m_bScaleUpToDate )
	{
		//Scaling in relation to the mass of the electron seems like a good reference - using ln because mass sizes are so different
		float flScale = Clamp<float>(log(GetParticleMass() / ELECTRON_MASS), 0.0f, 25.0f);
		if(debug_general.GetBool()) {Msg("Scaling particle %i to %.2f\n", entindex(), flScale);}
		m_hSprite->SetScale(flScale);
		m_bScaleUpToDate = true;
	}
	
	if( !m_bColourUpToDate && m_bStrongInt )
	{
		switch (m_iColour)
		{
			case CHARGE_ANTIBLUE:
				m_hSprite->SetColor(255,255,0);
				break;
			case CHARGE_ANTIGREEN:
				m_hSprite->SetColor(255,0,255);
				break;
			case CHARGE_ANTIRED:
				m_hSprite->SetColor(0,255,255);
				break;
			case CHARGE_NEUTRAL:
				break; //Let the derived class take care of it
			case CHARGE_RED:
				m_hSprite->SetColor(255,0,0);
				break;
			case CHARGE_GREEN:
				m_hSprite->SetColor(0,255,0);
				break;
			case CHARGE_BLUE:
				m_hSprite->SetColor(0,0,255);
				break;
		}
		m_bColourUpToDate = true;
	}
}

void CBaseParticle::ApplyParticlePhysics() {
	//Coulombic attraction/repulsion:
	Vector vecForce;
	int iParticleCount = g_vpParticles.Count();
	for(int i=0; i<iParticleCount; i++) {
		CBaseParticle *pParticle = g_vpParticles[i];
		if(pParticle && pParticle != this && pParticle->GetLocalOrigin() != GetLocalOrigin()) {
			Vector vecAddition = GetLocalOrigin() - pParticle->GetLocalOrigin();
			float flDistance = VectorNormalize(vecAddition);
			vecAddition *= pParticle->GetCharge();
			vecAddition /= Square<float>(flDistance / DISTANCE_SCALE);
			vecForce += vecAddition;
		}
	}
	vecForce *= COULOMB * GetCharge();
	VPhysicsGetObject()->ApplyForceCenter(vecForce);
}

bool CBaseParticle::AcceptEnergy(float flEnergy, float flTolerancePortion) {
	float flNewEnergy = m_flTotalEnergy + flEnergy;
	if(flNewEnergy > m_flRestEnergy) {
		m_flTotalEnergy = flNewEnergy;
		m_bFree = true;
		return true;
	}
	return false;
}

float CBaseParticle::GetParticleMass(double *pdKG, double *pdEV, double *pdPlanck) {
	float flMass = m_bUsesVPhysics ? VPhysicsGetObject()->GetMass() : m_flAppropriateMass; //Quarks and photons are not simulated using VPhysics
	
	if(pdKG) { //if we get a pointer to a variable for the real mass in kg, return that.
		(*pdKG) = flMass / MASS_SCALE;
	}
	if(pdEV) { //if we get a pointer to a variable for the real mass in eV/c^2, return that.
		(*pdEV) = (flMass / MASS_SCALE) / EV_MASS;
	}
	if(pdPlanck) { //if we get a pointer to a variable for the real mass in Planck units, return that.
		(*pdPlanck) = (flMass / MASS_SCALE) / PLANCK_MASS;
	}
	
	return flMass; //return the VPhys mass (scaled by 1e30) in kg
}

float CBaseParticle::GetCharge(float *pflChargeInUnits) {
	float flCharge = m_bAntiParticle ? -m_flCharge : m_flCharge;
	
	if(pflChargeInUnits) { //if we get a pointer to a variable for the charge in atomic charge units, return that.
		(*pflChargeInUnits) = flCharge / ELEMENTARY_CHARGE;
	}
	
	return flCharge; //return the electric charge in coulombs
}

/*complex CBaseParticle::g(int n, int k, double r, float flPotential) {
	//C has dimensions of metre, so probably should be DISTANCE_SCALEd too.
	double C = sqrt(Square<double>(m_flRestEnergy) - Square<double>(m_flTotalEnergy)) / (DIRAC * SPD_LGT);
	double ZALPHA = -flPotential * (r / (DIRAC*SPD_LGT));
	if( k==-n ) {
		C = ZALPHA / n * m_flRestEnergy / (DIRAC * SPD_LGT * DISTANCE_SCALE);
	}
	if(debug_wavefunctions.GetBool()){Msg("g on %i: C is %.20f == sqrt(%.40f - %.40f) / %.40f\n", entindex(), C, Square<double>(m_flRestEnergy), Square<float>(m_flTotalEnergy), DIRAC * SPD_LGT);}
	double rho = 2*C*r; //"Scaled radius"
	if(debug_wavefunctions.GetBool()){Msg("g on %i: rho is %.20f\n", entindex(), rho);}
	double gamma = sqrt(Square<int>(k) - Square<double>(ZALPHA));
	if(debug_wavefunctions.GetBool()){Msg("g on %i: gamma is %.20f\n", entindex(), gamma);}
	if( k==-n ) {
		double A = (1 / sqrt(2*n*(n+gamma))) * sqrt(C / (gamma * gsl_sf_gamma(2*gamma)));
		//double A = 1 / (sqrt(2*n*(n+gamma)) * sqrt(C / (gamma * gsl_sf_gamma(2*gamma))));
		if(debug_wavefunctions.GetBool()){Msg("g on %i: A is %.20f\n", entindex(), A);}
		if(debug_wavefunctions.GetBool()){Msg("g on %i: n+gamma is %.20f, rho^gamma is %.20f, e^(-rho/2) is %.20f\n", entindex(), n+gamma, pow(rho, gamma), pow(M_E, -rho/2));}
		return A * (n+gamma) * pow(rho, gamma) * pow(M_E, -rho/2);
	}
	int absk = k>=0 ? k : -k;
	if(debug_wavefunctions.GetBool()){Msg("g on %i: |k| is %i\n", entindex(), absk);}
	double A = (1 / sqrt(2*k*(k-gamma))) * sqrt( (C / (n-absk+gamma)) * (gsl_sf_fact(n-absk-1) / gsl_sf_gamma(n - absk + 2*gamma + 1)) * 0.5f*(Square<double>(m_flTotalEnergy*k / (gamma * m_flRestEnergy)) + (m_flTotalEnergy*k / (gamma * m_flRestEnergy))) );
	//double A = 1 / (sqrt(2*k*(k-gamma)) * sqrt( (C / (n-absk+gamma)) * (gsl_sf_fact(n-absk-1) / gsl_sf_gamma(n - absk + 2*gamma + 1)) * 0.5f*(Square<double>(m_flTotalEnergy*k / (gamma * m_flRestEnergy)) + (m_flTotalEnergy*k / (gamma * m_flRestEnergy))) ));
	if(debug_wavefunctions.GetBool()){Msg("g on %i: A is %.20f\n", entindex(), A);}
	return A * pow(r, gamma) * pow(M_E, -rho/2) * ((ZALPHA * rho * gsl_sf_laguerre_n(n-absk-1, 2*gamma + 1, rho)) + ((gamma-k) * ((gamma*m_flRestEnergy - k*m_flTotalEnergy)/(DIRAC*SPD_LGT * DISTANCE_SCALE*C)) * gsl_sf_laguerre_n(n-absk, 2*gamma + 1, rho)));
}*/

//In Planck units
complex CBaseParticle::g(int n, int k, double r, float flPotential) {
	double C = sqrt(Square<double>(m_flRestEnergy / PLANCK_ENERGY) - Square<double>(m_flTotalEnergy / PLANCK_ENERGY));
	double ZALPHA = -(flPotential / PLANCK_ENERGY) * (r / PLANCK_LENGTH);
	if(debug_wavefunctions.GetBool()){Msg("g on %i: Z*alpha is %.20f\n", entindex(), ZALPHA);}
	if( k==-n ) {
		C = (ZALPHA / n) * (m_flRestEnergy / PLANCK_ENERGY);
	}
	if(debug_wavefunctions.GetBool()){Msg("g on %i: C is %.20f\n", entindex(), C);}
	double rho = 2*C*(r/PLANCK_LENGTH); //"Scaled radius"
	if(debug_wavefunctions.GetBool()){Msg("g on %i: r is %.20f\n", entindex(), r);}
	if(debug_wavefunctions.GetBool()){Msg("g on %i: rho is %.20f\n", entindex(), rho);}
	double gamma = sqrt(Square<int>(k) - Square<double>(ZALPHA));
	if(debug_wavefunctions.GetBool()){Msg("g on %i: gamma is %.20f\n", entindex(), gamma);}
	if( k==-n ) {
		double A = (1 / sqrt(2*n*(n+gamma))) * sqrt(C / (gamma * gsl_sf_gamma(2*gamma)));
		if(debug_wavefunctions.GetBool()){Msg("g on %i: A is %.20f\n", entindex(), A);}
		if(debug_wavefunctions.GetBool()){Msg("g on %i: n+gamma is %.20f, rho^gamma is %.20f, e^(-rho/2) is %.20f\n", entindex(), n+gamma, pow(rho, gamma), pow(M_E, -rho/2));}
		return A * (n+gamma) * pow(rho, gamma) * pow(M_E, -rho/2);
	}
	int absk = k>=0 ? k : -k;
	if(debug_wavefunctions.GetBool()){Msg("g on %i: |k| is %i\n", entindex(), absk);}
	double A = (1 / sqrt(2*k*(k-gamma))) * sqrt( (C / (n-absk+gamma)) * (gsl_sf_fact(n-absk-1) / gsl_sf_gamma(n - absk + 2*gamma + 1)) * 0.5f*(Square<double>(m_flTotalEnergy/PLANCK_ENERGY*k / (gamma * m_flRestEnergy/PLANCK_ENERGY)) + (m_flTotalEnergy/PLANCK_ENERGY*k / (gamma * m_flRestEnergy/PLANCK_ENERGY))) );
	//double A = 1 / (sqrt(2*k*(k-gamma)) * sqrt( (C / (n-absk+gamma)) * (gsl_sf_fact(n-absk-1) / gsl_sf_gamma(n - absk + 2*gamma + 1)) * 0.5f*(Square<double>(m_flTotalEnergy*k / (gamma * m_flRestEnergy)) + (m_flTotalEnergy*k / (gamma * m_flRestEnergy))) ));
	if(debug_wavefunctions.GetBool()){Msg("g on %i: A is %.20f\n", entindex(), A);}
	return A * pow(r, gamma) * pow(M_E, -rho/2) * ((ZALPHA * rho * gsl_sf_laguerre_n(n-absk-1, 2*gamma + 1, rho)) + ((gamma-k) * ((gamma*m_flRestEnergy - k*m_flTotalEnergy)/(PLANCK_ENERGY*C)) * gsl_sf_laguerre_n(n-absk, 2*gamma + 1, rho)));
}

/*complex CBaseParticle::f(int n, int k, double r, float flPotential) {
	double C = sqrt(Square<double>(m_flRestEnergy) - Square<float>(m_flTotalEnergy)) / (DIRAC * SPD_LGT * DISTANCE_SCALE);
	double ZALPHA = -flPotential * (r / (DIRAC*SPD_LGT));
	if( k==-n ) {
		C = ZALPHA / n * m_flRestEnergy / (DIRAC * SPD_LGT * DISTANCE_SCALE);
	}
	double rho = 2*C*r; //"Scaled radius"
	double gamma = sqrt(Square<int>(k) - Square<double>(ZALPHA));
	if( k==-n ) {
		double A = (1 / sqrt(2*n*(n+gamma))) * sqrt(C / (gamma * gsl_sf_gamma(2*gamma)));
		//double A = 1 / (sqrt(2*n*(n+gamma)) * sqrt(C / (gamma * gsl_sf_gamma(2*gamma))));
		if(debug_wavefunctions.GetBool()){Msg("f on %i: A is %.20f\n", entindex(), A);}
		return A * ZALPHA * pow(rho, gamma) * pow(M_E, -rho/2);
	}
	int absk = k>=0 ? k : -k;
	double A = (1 / sqrt(2*k*(k-gamma))) * sqrt( (C / (n-absk+gamma)) * (gsl_sf_fact(n-absk-1) / gsl_sf_gamma(n - absk + 2*gamma + 1)) * 0.5f*(Square<double>(m_flTotalEnergy*k / (gamma * m_flRestEnergy)) + (m_flTotalEnergy*k / (gamma * m_flRestEnergy))) );
	//double A = 1 / (sqrt(2*k*(k-gamma)) * sqrt( (C / (n-absk+gamma)) * (gsl_sf_fact(n-absk-1) / gsl_sf_gamma(n - absk + 2*gamma + 1)) * 0.5f*(Square<double>(m_flTotalEnergy*k / (gamma * m_flRestEnergy)) + (m_flTotalEnergy*k / (gamma * m_flRestEnergy))) ));
	return A * pow(rho, gamma) * pow(M_E, -rho/2) * ((gamma-k) * rho * gsl_sf_laguerre_n(n-absk-1, 2*gamma + 1, rho) + ZALPHA * ((gamma*ReducedMass()*C_SQ - k*m_flTotalEnergy)/(DIRAC*SPD_LGT * DISTANCE_SCALE*C)) * gsl_sf_laguerre_n(n-absk, 2*gamma + 1, rho));
}*/

//In Planck units
complex CBaseParticle::f(int n, int k, double r, float flPotential) {
	double C = sqrt(Square<double>(m_flRestEnergy / PLANCK_ENERGY) - Square<double>(m_flTotalEnergy / PLANCK_ENERGY));
	double ZALPHA = -(flPotential / PLANCK_ENERGY) * (r / PLANCK_LENGTH);
	if(debug_wavefunctions.GetBool()){Msg("f on %i: Z*alpha is %.20f\n", entindex(), ZALPHA);}
	if( k==-n ) {
		C = (ZALPHA / n) * (m_flRestEnergy / PLANCK_ENERGY);
	}
	if(debug_wavefunctions.GetBool()){Msg("f on %i: C is %.20f\n", entindex(), C);}
	double rho = 2*C*(r/PLANCK_LENGTH); //"Scaled radius"
	if(debug_wavefunctions.GetBool()){Msg("f on %i: r is %.20f\n", entindex(), r);}
	if(debug_wavefunctions.GetBool()){Msg("f on %i: rho is %.20f\n", entindex(), rho);}
	double gamma = sqrt(Square<int>(k) - Square<double>(ZALPHA));
	if(debug_wavefunctions.GetBool()){Msg("f on %i: gamma is %.20f\n", entindex(), gamma);}
	if( k==-n ) {
		double A = (1 / sqrt(2*n*(n+gamma))) * sqrt(C / (gamma * gsl_sf_gamma(2*gamma)));
		if(debug_wavefunctions.GetBool()){Msg("f on %i: A is %.20f\n", entindex(), A);}
		if(debug_wavefunctions.GetBool()){Msg("f on %i: n+gamma is %.20f, rho^gamma is %.20f, e^(-rho/2) is %.20f\n", entindex(), n+gamma, pow(rho, gamma), pow(M_E, -rho/2));}
		return A * ZALPHA * pow(rho, gamma) * pow(M_E, -rho/2);
	}
	int absk = k>=0 ? k : -k;
	double A = (1 / sqrt(2*k*(k-gamma))) * sqrt( (C / (n-absk+gamma)) * (gsl_sf_fact(n-absk-1) / gsl_sf_gamma(n - absk + 2*gamma + 1)) * 0.5f*(Square<double>(m_flTotalEnergy*k / (gamma * m_flRestEnergy)) + (m_flTotalEnergy*k / (gamma * m_flRestEnergy))) );
	//double A = 1 / (sqrt(2*k*(k-gamma)) * sqrt( (C / (n-absk+gamma)) * (gsl_sf_fact(n-absk-1) / gsl_sf_gamma(n - absk + 2*gamma + 1)) * 0.5f*(Square<double>(m_flTotalEnergy*k / (gamma * m_flRestEnergy)) + (m_flTotalEnergy*k / (gamma * m_flRestEnergy))) ));
	return A * pow(rho, gamma) * pow(M_E, -rho/2) * ((gamma-k) * rho * gsl_sf_laguerre_n(n-absk-1, 2*gamma + 1, rho) + ZALPHA * ((gamma*ReducedMass()*C_SQ - k*m_flTotalEnergy)/(DIRAC*SPD_LGT * DISTANCE_SCALE*C)) * gsl_sf_laguerre_n(n-absk, 2*gamma + 1, rho));
}

complex CBaseParticle::Y(double a, double b, Vector &vecDirection) {
	if (a<0){
		a = -a-1;
	}
	if(debug_wavefunctions.GetBool()){Msg("Y on %i: a is %.2f, b is %.2f\n", entindex(), a, b);}
	
	//This is rather complicated. This function is supposed to be complex, yet Valve's function is real.
	//On investigation, it turns out that Valve's function will return the real part when given a positive b-value, and the imaginary part when given a negative b-value.
	//The actual function involves an e^i term, but Valve's has cos for positive b and sin for negative b.
	//It also uses b in the cos but -b in the sin, confirming that this behaviour is designed to give real and imaginary parts.
	//Unfortunately, the actual b-values we use can be either positive or negative.
	//Therefore, we can first check if b is negative. Then, we can minus it as appropriate.
	//There is also a (-1)^b term in the quantum mechanics function, which is not present in the Valve function
	if(b < 0) {
		//Simply pass -b to the real part, since cos(-x) = cos(x)
		//Pass b to the imaginary part, but minus the answer, since sin(-x) = -sin(x)
		return pow(-1, b) * complex(SphericalHarmonic(a, -b, vecDirection), -SphericalHarmonic(a, b, vecDirection));
	} else if(b == 0) {
		//The e^i term equals 1 in this case, so there is only a real component.
		return pow(-1, b) * complex(SphericalHarmonic(a, b, vecDirection), 0);
	} else {
		//Pass b to the real part and -b to the imaginary part.
		return pow(-1, b) * complex(SphericalHarmonic(a, b, vecDirection), SphericalHarmonic(a, -b, vecDirection));
	}
}

gsl_complex CBaseParticle::Psi(Vector vecLocation) {
	if(m_bFree) {
		//Solution to Dirac equation for free particle (e^(ip.x/h-bar)) - http://www.nyu.edu/classes/tuckerman/quant.mech/lectures/lecture_7/node1.html
		return (pow(M_E, complex(0, DotProduct(m_vecMomentum, vecLocation) / DIRAC) )).gsl();
	}
	
	//Calculate Potential Energy at this position:
	int iParticleCount = g_vpParticles.Count();
	float flPotential = 0;
	for(int i=0; i<iParticleCount; i++) {
		CBaseParticle *pParticle = g_vpParticles[i];
		if(pParticle && pParticle != this) {
			if(pParticle->GetLocalOrigin() != vecLocation) {
				//This is the Q/r part of the formula, i.e. charge over distance
				double dAddition = pParticle->GetCharge();
				dAddition /= ((pParticle->GetLocalOrigin() - vecLocation).Length()/DISTANCE_SCALE);
				flPotential += dAddition;
			}
		}
	}
	flPotential *= m_flCharge * COULOMB;
	float flKinetic = m_flTotalEnergy - m_flRestEnergy - flPotential;
	
	if(debug_wavefunctions.GetBool()){Msg("Psi on %i: ", entindex());}
	
	Vector vecDirection = vecLocation - g_vecNuclearOrigin;
	double r = VectorNormalize(vecDirection) / DISTANCE_SCALE; //Small radius - g and f will scale it back up
	
	if(m_flSpinProjection == 0) {
		if(debug_wavefunctions.GetBool()){Msg("Psi on %i: 0 Spin projection???\n", entindex());}
		Aufbau();
	}
	if(debug_wavefunctions.GetBool()){Msg("Psi on %i: Spin projection is %.1f\n", entindex(), m_flSpinProjection);}
	
	float j = m_flSpinProjection + m_iAzimuthal;
	if(j<0){j=-j;}
	float m = (j-m_iAzimuthal) + m_iMagnetic;
	float k = (j==0.5f+m_iAzimuthal)?(-j - 0.5f):(j + 0.5f);
	float absk = k>=0 ? k : -k;
	int sgnk = k==0 ? 0 : (k>0 ? 1 : -1);
	if(debug_wavefunctions.GetBool()){Msg("Psi on %i: j is %.1f, mj is %.1f\n", entindex(), j, m);}
	if(debug_wavefunctions.GetBool()){Msg("Psi on %i: k is %.1f, |k| is %.1f\n", entindex(), k, absk);}
	
	complex Answer;
	//Determine which one to use
	if(!m_bAntiParticle) { //Positive energy is matter
		if (m_flSpinProjection > 0){
			if(debug_wavefunctions.GetBool()){Msg("Psi on %i: Spin-up, positive energy\n", entindex());}
			//Spin-up, positive energy
			Answer = g(m_iPrincipal,k,r,flPotential)/r * complex_sqrt((k+0.5f-m)/(2*k+1)) * Y(k, m-0.5f, vecDirection);
		} else {
			if(debug_wavefunctions.GetBool()){Msg("Psi on %i: Spin-down, positive energy\n", entindex());}
			//Spin-down, positive energy
			Answer = -1.0f * g(m_iPrincipal,k,r,flPotential)/r * sgnk * complex_sqrt((k+0.5f+m)/(2*k+1)) * Y(k, m+0.5f, vecDirection);
		}
	} else { //Negative energy is antimatter
		if (m_flSpinProjection > 0){
			if(debug_wavefunctions.GetBool()){Msg("Psi on %i: Spin-up, negative energy\n", entindex());}
			//Spin-up, negative energy
			Answer = complex(0,1) * f(m_iPrincipal,k,r,flPotential)/r * complex_sqrt((-k+0.5f-m)/(-2*k+1)) * Y(-k, m-0.5f, vecDirection);
		} else {
			if(debug_wavefunctions.GetBool()){Msg("Psi on %i: Spin-down, negative energy\n", entindex());}
			//Spin-down, negative energy
			Answer = complex(0,-1) * f(m_iPrincipal,k,r,flPotential)/r * sgnk * complex_sqrt((-k+0.5f-m)/(-2*k+1)) * Y(-k, m+0.5f, vecDirection);
		}
	}
	/*complex supe = g(m_iPrincipal,k,r,flPotential)/r * complex_sqrt((k+0.5f-m)/(2*k+1)) * Y(k, m-0.5f, vecDirection);
	complex sdpe = -1.0f * g(m_iPrincipal,k,r,flPotential)/r * sgnk * complex_sqrt((k+0.5f+m)/(2*k+1)) * Y(k, m+0.5f, vecDirection);
	complex sune = complex(0,1) * f(m_iPrincipal,k,r,flPotential)/r * complex_sqrt((-k+0.5f-m)/(-2*k+1)) * Y(-k, m-0.5f, vecDirection);
	complex sdne = complex(0,-1) * f(m_iPrincipal,k,r,flPotential)/r * sgnk * complex_sqrt((-k+0.5f-m)/(-2*k+1)) * Y(-k, m+0.5f, vecDirection);
	if(!m_bAntiParticle) { //Positive energy is matter
		if (m_flSpinProjection > 0){
			if(debug_wavefunctions.GetBool()){Msg("Psi on %i: Spin-up, positive energy\n", entindex());}
			//Spin-up, positive energy
			Answer = supe - sune;//(sdpe + sune + sdne);
		} else {
			if(debug_wavefunctions.GetBool()){Msg("Psi on %i: Spin-down, positive energy\n", entindex());}
			//Spin-down, positive energy
			Answer = sdpe - sdne;//(supe + sune + sdne);
		}
	} else { //Negative energy is antimatter
		if (m_flSpinProjection > 0){
			if(debug_wavefunctions.GetBool()){Msg("Psi on %i: Spin-up, negative energy\n", entindex());}
			//Spin-up, negative energy
			Answer = sune - supe;//(supe + sdpe + sdne);
		} else {
			if(debug_wavefunctions.GetBool()){Msg("Psi on %i: Spin-down, negative energy\n", entindex());}
			//Spin-down, negative energy
			Answer = sdne - sdpe;//(supe + sdpe + sune);
		}
	}*/
	if(debug_wavefunctions.GetBool()){Msg("Psi on %i: About to return %.20f + (%.20f i)\n", entindex(), Answer.real(), Answer.im());}
	
	return Answer.gsl();
}

double CBaseParticle::Psi2(Vector vecLocation) {
	return gsl_complex_abs2(Psi(vecLocation));// / pow(DISTANCE_SCALE/2, 3);
}

//Integrate |Psi|^2 over a region using the trapezoidal rule
double CBaseParticle::DiagnosticIntegral(int iStart, int iEnd, int iCoarseness) {
	double dResult = 0.0f;
	double dPrevX = 0.0f;
	double dCurX = 0.0f;
	for(int x=iStart; x<=iEnd; x+=iCoarseness) {
		double dPrevY = 0.0f;
		double dCurY = 0.0f;
		for(int y=iStart; y<=iEnd; y+=iCoarseness) {
			double dPrevZ = 0.0f;
			double dCurZ = 0.0f;
			for(int z=iStart; z<=iEnd; z+=iCoarseness) {
				dCurZ = Psi2(Vector(x,y,z));
				dCurY += iCoarseness * (dCurZ + dPrevZ) / 2;
				dPrevZ = dCurZ;
			}
			//Msg("Finished a line\n");
			dCurX += iCoarseness * (dCurY + dPrevY) / 2;
			dPrevY = dCurY;
			dCurY = 0.0f;
		}
		//Msg("Finished a plane\n");
		dResult += iCoarseness * (dCurX + dPrevX) / 2;
		dPrevX = dCurX;
		dCurX = 0.0f;
	}
	return dResult;
}

static void GlobalDiagnosticIntegral(const CCommand &args) {
	int iStart = -4096;
	int iEnd = 4096;
	if(args.ArgC() > 1 && atoi(args.Arg(1))) {
		iStart = -(iEnd = atoi(args.Arg(1)));
		if(args.ArgC() > 2 && atoi(args.Arg(2))) {
			iEnd = atoi(args.Arg(2));
		}
	}
	int iSize = iEnd - iStart;
	int iCoarseness;
	if(iSize < 0) {
		return;
	} else if (iSize < 150) {
		iCoarseness = 1;
	} else if (iSize < 300) {
		iCoarseness = 2;
	} else if (iSize < 500) {
		iCoarseness = 4;
	} else if (iSize < 1000) {
		iCoarseness = 8;
	} else if (iSize < 2000) {
		iCoarseness = 16;
	} else if (iSize < 4000) {
		iCoarseness = 32;
	} else {
		iCoarseness = 64;
	}
	CBaseParticle *pElectron = NULL;
	int iParticleIndex = 0;
	while(!pElectron && iParticleIndex < g_vpParticles.Count()) {
		if(g_vpParticles[iParticleIndex] && g_vpParticles[iParticleIndex]->IsElectron()) {
			pElectron = g_vpParticles[iParticleIndex];
		} else {
			iParticleIndex++;
		}
	}
	Msg("Diagnostic Integral of electron %i returned %.20f\n", pElectron->entindex(), pElectron->DiagnosticIntegral(iStart, iEnd, iCoarseness));
}
static ConCommand diagnostic_integral("diagnostic_integral", GlobalDiagnosticIntegral, "Integrates |Psi|^2 of the first electron in the global vector over space specified");

class CPhoton : public CBaseParticle { //Massless particle
	DECLARE_CLASS(CPhoton, CBaseParticle);
	
public:
	CPhoton();
	bool m_bNeutrino; //The massless particle code for a photon works just as well for a neutrino
	DECLARE_DATADESC();
	
	virtual void ParticleUpdateSprite();
	virtual void ApplyParticlePhysics() {} //No special behaviour - just keep moving 'till we hit something!
	
	virtual void Touch(CBaseEntity *pOther);
	virtual bool CreateVPhysics();
	
	static CPhoton *CreatePhoton(float flEnergy, Vector vecPosition, QAngle angDirection, CBaseParticle *pOwner);
	static CPhoton *CreateNeutrino(float flEnergy, Vector vecPosition, QAngle angDirection, CBaseParticle *pOwner);
};
LINK_ENTITY_TO_CLASS(particle_photon, CPhoton);

BEGIN_DATADESC(CPhoton)
	DEFINE_FIELD(m_bNeutrino, FIELD_BOOLEAN),
END_DATADESC()

CPhoton *CPhoton::CreatePhoton(float flEnergy, Vector vecPosition, QAngle angDirection, CBaseParticle *pOwner) {
	if(flEnergy < 0.0f) {
		Msg("Attempted to create photon with negative energy!\n");
		return NULL;
	}
	
	CPhoton *pPhoton = (CPhoton *)CreateEntityByName("particle_photon");
	if(!pPhoton) {
		return NULL;
	}
	
	pPhoton->m_flTotalEnergy = flEnergy;
	pPhoton->SetOwnerEntity(pOwner);
	pPhoton->SetAbsOrigin(vecPosition);
	pPhoton->SetAbsAngles(angDirection);
	pPhoton->Spawn();
	
	return pPhoton;
}

CPhoton *CPhoton::CreateNeutrino(float flEnergy, Vector vecPosition, QAngle angDirection, CBaseParticle *pOwner) {
	CPhoton *pNeutrino = CreatePhoton(flEnergy, vecPosition, angDirection, pOwner);
	if(pNeutrino) {
		pNeutrino->m_bNeutrino = true;
	}
	return pNeutrino;
}

CPhoton::CPhoton() {
	m_bNeutrino = false;
	m_flAppropriateMass = 0.0f;
	m_flCharge = 0.0f;
}

void CPhoton::ParticleUpdateSprite() {
	// Must have a valid sprite to affect
	if ( !m_hSprite )
	{
		// Create our sprite
		m_hSprite = CSprite::SpriteCreate( PARTICLE_GLOW_SPRITE, GetLocalOrigin(), false );
		if ( !m_hSprite )
			return;

		m_hSprite->SetTransparency( kRenderWorldGlow, 255, 0, 0, 255, kRenderFxNoDissipation );
		m_hSprite->SetParent(this);
		
		m_bScaleUpToDate = false;
		m_bColourUpToDate = false;
	}
	
	if( !m_bScaleUpToDate )
	{
		m_hSprite->SetScale(m_bNeutrino ? 0.5f : 3.0f);
		m_bScaleUpToDate = true;
	}
	
	if( !m_bColourUpToDate )
	{
		if(m_bNeutrino) {
			m_hSprite->SetColor(255,255,255); //White
		} else {
			float flFrequency = m_flTotalEnergy / PLANCK; //E=hf
			m_flWavelength = SPD_LGT / flFrequency; //c=f(lambda)
			Msg("Wavelength of photon %i is %.20f\n", entindex(), m_flWavelength);
			if( m_flWavelength < 380e-9 ) {
				//Ultraviolet
				m_hSprite->SetColor(128,64,255); //Pinky?
			} else if( m_flWavelength > 780e-9 ) {
				//Infrared
				m_hSprite->SetColor(128,0,0); //Brown
				//The following approximate conversions are taken from http://stackoverflow.com/questions/3407942/rgb-values-of-visible-spectrum
				//Since intensity isn't an issue they should be okay
			} else if( m_flWavelength >= 645e-9 ) {
				//Red
				m_hSprite->SetColor(255,0,0);
			} else if( m_flWavelength >= 580e-9 ) {
				//Orange
				float flGreen = -(m_flWavelength-645e-9) / (645e-9-580e-9);
				flGreen *= 256;
				int iGreen = (int)clamp<float>(flGreen, 0, 255);
				m_hSprite->SetColor(255,iGreen,0);
			} else if( m_flWavelength >= 510e-9 ) {
				//Yellow
				float flRed = -(m_flWavelength-580e-9) / (580e-9-510e-9);
				flRed *= 256;
				int iRed = (int)clamp<float>(flRed, 0, 255);
				m_hSprite->SetColor(iRed,255,0);
			} else if( m_flWavelength >= 490e-9 ) {
				//Green
				float flBlue = -(m_flWavelength-510e-9) / (510e-9-490e-9);
				flBlue *= 256;
				int iBlue = (int)clamp<float>(flBlue, 0, 255);
				m_hSprite->SetColor(0,255,iBlue);
			} else if( m_flWavelength >= 440e-9 ) {
				//Blue
				float flGreen = -(m_flWavelength-490e-9) / (490e-9-440e-9);
				flGreen *= 256;
				int iGreen = (int)clamp<float>(flGreen, 0, 255);
				m_hSprite->SetColor(0,iGreen,255);
			} else if( m_flWavelength >= 380e-9 ) {
				//Purple
				float flRed = -(m_flWavelength-440e-9) / (440e-9-380e-9);
				flRed *= 256;
				int iRed = (int)clamp<float>(flRed, 0, 255);
				m_hSprite->SetColor(iRed,0,255);
			}
		}
		
		m_bColourUpToDate = true;
	}
}

void CPhoton::Touch(CBaseEntity *pOther) {
	if(!pOther || pOther == GetOwnerEntity()) {
		return;
	}
	Msg("Photon %i Touch()ed %i!\n", entindex(), pOther->entindex());
	if(pOther->IsWorld()) {
		//Touched the walls...
		UTIL_Remove(this);
	} else {
		CBaseParticle *pParticle = dynamic_cast<CBaseParticle *>(pOther);
		if(pParticle && pParticle->AcceptEnergy(m_flTotalEnergy, 1e-15)) {
			UTIL_Remove(this);
		}
	}
	BaseClass::Touch(pOther);
}

bool CPhoton::CreateVPhysics() {
	//The photon is massless, so there's no real point in using VPhysics
	SetSolid( SOLID_BBOX );
	SetCollisionBounds( -Vector(3,3,3), Vector(3,3,3) );
	AddSolidFlags(FSOLID_NOT_SOLID | FSOLID_TRIGGER); //We don't actually want to collide with particles, just know when they enter our bounding box so we can transfer energy
	
	SetMoveType( MOVETYPE_FLY );
	Vector vecVelocity = vec3_origin;
	AngleVectors(GetAbsAngles(), &vecVelocity);
	vecVelocity *= SPD_LGT * (DISTANCE_SCALE/TIME_SCALE);
	SetAbsVelocity(vecVelocity);

	return false;
}

class CElectron : public CBaseParticle {
	DECLARE_CLASS(CElectron, CBaseParticle);
	
public:
	//DECLARE_DATADESC();
	
	CElectron();
	
	virtual void ParticleUpdateSprite();
	virtual void ApplyParticlePhysics();
	virtual bool AcceptEnergy(float flEnergy, float flTolerancePortion);
	virtual bool IsElectron(){return true;}
};

static CUtlVector<CElectron *> g_vpElectrons;

LINK_ENTITY_TO_CLASS(particle_electron, CElectron);

CElectron::CElectron() {
	m_flCharge = -ELEMENTARY_CHARGE;
	
	g_vpElectrons.AddToHead(this);
}

void CElectron::ParticleUpdateSprite() {
	BaseClass::ParticleUpdateSprite();
	
	if( !m_bColourUpToDate )
	{
		if(m_bAntiParticle) {
			m_hSprite->SetColor(0,0,128); //Opposite of pale yellow below
		} else {
			m_hSprite->SetColor(255,255,128); //Pale yellow, so it's not mixed up with antiblue
		}
		m_bColourUpToDate = true;
	}
}

void CElectron::ApplyParticlePhysics() {
	//Calculate our potential energy
	m_flPotential = 0; //Starting from a base of 0
	//Coulombic:
	int iParticleCount = g_vpParticles.Count();
	float Z = 0; //Effective nuclear charge
	if(!m_bAnnihilating) {
		for(int i=0; i<iParticleCount; i++) {
			CBaseParticle *pParticle = g_vpParticles[i];
			if(pParticle && pParticle != this) {
				if(pParticle->IsElectron() && pParticle->m_bAntiParticle != m_bAntiParticle && (pParticle->GetLocalOrigin() - GetLocalOrigin()).LengthSqr() <= 400) {
					//Positron less than 20 units away
					m_bAnnihilating = true;
					pParticle->Annihilate(); //They don't seem to do this independently, so ensure mutual agreement
					break;
				} else if(pParticle->GetLocalOrigin() != GetLocalOrigin()) {
					//This is the Q/r part of the formula, i.e. charge over distance
					double dAddition = pParticle->GetCharge();
					if(pParticle->IsProton()) {
						Z += 1;
					} else if(pParticle->IsElectron() && pParticle->m_iPrincipal <= m_iPrincipal) {
						//This involves using Slater's rules for the screening effect to determine an effective nuclear charge
						if(m_iPrincipal==1) {
							//1s-1s shielding
							dAddition *= 0.3f; //Effectively, the electron has only 0.3 times the elementary charge in its dealings with us...
							Z -= 0.3f;
						} else if(m_iAzimuthal <= 1) {
							//s or p
							if(pParticle->m_iPrincipal == m_iPrincipal) {
								//Same group
								dAddition *= 0.35f;
								Z -= 0.35f;
							} else if(pParticle->m_iPrincipal == m_iPrincipal-1) {
								//Group just below
								dAddition *= 0.85f;
								Z -= 0.85f;
							} else if(pParticle->m_iPrincipal < m_iPrincipal-1) {
								//Groups further down
								Z -= 1;
							}
						} else if(m_iAzimuthal > 1) {
							//d or f (or higher)
							if(pParticle->m_iPrincipal == m_iPrincipal) {
								//Same subgroup
								dAddition *= 0.35f;
								Z -= 0.35f;
							} else if(pParticle->m_iAzimuthal < m_iAzimuthal || pParticle->m_iPrincipal < m_iPrincipal) {
								//(sub)groups below
								Z -= 1;
							}
						}
					}
					dAddition /= ((pParticle->GetLocalOrigin() - GetLocalOrigin()).Length()/DISTANCE_SCALE);
					m_flPotential += dAddition;
				}
			}
		}
	}
	m_flPotential *= GetCharge() * COULOMB;
	
	if(debug_electrons.GetBool()){Msg("Coulombic potential energy of %i is %.30f\n", entindex(), m_flPotential);}
	
	double dMass = 0;
	GetParticleMass(&dMass); //Real kg
		
	VPhysicsGetObject()->GetVelocity(&m_vecMomentum, NULL);
	m_vecMomentum *= VPhysicsGetObject()->GetMass();
	if(debug_electrons.GetBool()){Msg("Momentum of %i is (%.20f, %.20f, %.20f)\n", entindex(), m_vecMomentum.x, m_vecMomentum.y, m_vecMomentum.z);}
	
	m_flRestEnergy = ReducedMass() * C_SQ;
	float flTotalEnergy = 0;
	if(!m_bAnnihilating) {
		//Extra quantum numbers required for kinetic energy
		float j = m_iAzimuthal + m_flSpinProjection;
		if(j<0){j=-j;}
		float m = m_iMagnetic + (j-m_iAzimuthal);
		float k = (j==m_iAzimuthal+0.5f)?(-j - 0.5f):(j + 0.5f);
		float absk = k>=0 ? k : -k;
		
		//This method of calculating kinetic energy is outlined on paper
		//double ZALPHA = -m_flPotential * ((rho/DISTANCE_SCALE) / (DIRAC*SPD_LGT));
		double ZALPHA = Z * ALPHA;
		if(debug_electrons.GetBool()) {Msg("Z*alpha on %i is %.20f\n", entindex(), ZALPHA);}
		double gamma = sqrt(Square<float>(absk) - Square<double>(ZALPHA));
		if(debug_electrons.GetBool()) {Msg("gamma on %i is %.20f\n", entindex(), gamma);}
		double dSquareInSqrt = Square<double>(ZALPHA / (m_iPrincipal - absk + gamma));
		if(debug_electrons.GetBool()) {Msg("Square in square root in denominator of rest energy coefficient on %i is %.20f\n", entindex(), dSquareInSqrt);}
		double dRestCoeff = (1 / sqrt(1 + dSquareInSqrt));// - 1;
		if(debug_electrons.GetBool()) {Msg("Rest energy coefficient on %i is %.20f\n", entindex(), dRestCoeff);}
		
		flTotalEnergy = m_flRestEnergy*dRestCoeff;
		m_flKinetic = flTotalEnergy - (m_flRestEnergy + m_flPotential);
		//flTotalEnergy = m_flRestEnergy + m_flPotential + m_flKinetic;
	}
	
	if(debug_electrons.GetBool()){Msg("Total energy of %i is %.30f\n", entindex(), flTotalEnergy);}
	if(m_flTotalEnergy != flTotalEnergy) {
		float flDiff = m_flTotalEnergy - flTotalEnergy;
		if(flDiff > 0) {
			QAngle angDirection = vec3_angle;
			VectorAngles(m_vecMomentum, angDirection);
			int iDimension = random->RandomInt(1,3);
			switch(iDimension) {
				case 1:
					angDirection += QAngle(90,0,0);
					break;
				case 2:
					angDirection += QAngle(0,90,0);
					break;
				case 3:
					angDirection += QAngle(0,0,90);
			}
			CPhoton::CreatePhoton(flDiff, GetLocalOrigin(), angDirection, this);
		}
	}
	if(m_bAnnihilating) {
		UTIL_Remove(this);
		return;
	}
	m_flTotalEnergy = flTotalEnergy;
	if(debug_electrons.GetBool()){Msg("Potential energy of %i is %.30f\n", entindex(), m_flPotential);}
	if(debug_electrons.GetBool()){Msg("Rest energy of %i is %.30f\n", entindex(), m_flRestEnergy);}
	if(debug_electrons.GetBool()){Msg("Kinetic energy of %i is %.30f\n", entindex(), m_flKinetic);}
	if(m_flKinetic >= 1e-20) {
		//The MASS_SCALEs balance out here and are used to keep the number big, so it is not neglected by the datatype
		double dPermittedMomentum = (sqrtf(m_flKinetic * 2 * dMass * MASS_SCALE) / sqrtf(MASS_SCALE)); //Ek=p^2/2m
		dPermittedMomentum *= (DISTANCE_SCALE/TIME_SCALE);
		double dPermittedVelocity = dPermittedMomentum / dMass; //p=mv
		if(debug_electrons.GetBool()){Msg("Permitted momentum and velocity of %i are %.90f and %.90f\n", entindex(), dPermittedMomentum, dPermittedVelocity);}
		
		Vector vecForce = -m_vecMomentum; //Required impulse to stop us in our tracks...
		double dFinalProbability; //For brightness change calculation below
		if(m_bFree) {
			//Just keep propagating forwards, changing brightness to represent phases of wavefunction
			VectorNormalize(m_vecMomentum);
			Vector vNextPosition = m_vecMomentum * dPermittedVelocity / 25;
			
			m_vecMomentum *= dPermittedMomentum; //This is the value we'll give to the wave function - same direction, but updated magnitude
			vecForce += m_vecMomentum * MASS_SCALE;
			
			dFinalProbability = Psi2(vNextPosition + GetLocalOrigin());
		} else {
			//What we have to do here is:
			//-pick a few random points to go
			//-square the wave-function for each of these points (Psi^2)
			//-pick a random number between 0 and 1
			//-go to the point for which Psi^2 most closely matches the random number
			//Also, we might scale the sprite's brightness, based on how likely the particle was to go here
			//(In other words, how close the random number is to 1)
			
			QAngle aRandom = QAngle(random->RandomFloat(-180.0f, 180.0f), random->RandomFloat(-180.0f, 180.0f), random->RandomFloat(-180.0f, 180.0f));
			Vector vRandForward, vRandRight, vRandUp;
			AngleVectors(aRandom, &vRandForward, &vRandRight, &vRandUp);
			Vector vRandBack = -vRandForward;
			Vector vRandLeft = -vRandRight;
			Vector vRandDown = -vRandUp;
			//Now we have six directions to go, based on a random angle
			//Whittle it down to three, which don't make it look like it's reversing
			Vector vCurDirection = m_vecMomentum;
			if( vCurDirection == vec3_origin) {
				//This probably means that the player has grabbed us with the physcannon.
				//Whatever it means, it will cause a crash. We could just stop the code from executing, but then we run the risk of stopping forever.
				//Instead, put a fake direction vector here to keep it going.
				//It won't intefere with physcannon behaviour anyway, and will ensure continued movement as soon as the player lets go.
				vCurDirection = Vector(1,0,0);
			} else {
				VectorNormalize(vCurDirection);
			}
			Vector *ppDirections[6] = {&vRandForward, &vRandRight, &vRandUp, &vRandBack, &vRandLeft, &vRandDown};
			Vector *ppAcceptedDirections[3] = {NULL, NULL, NULL};
			for(int i=0;i<6;i++) {
				if(DotProduct(*ppDirections[i],vCurDirection) > 0) { //No acute angles!
					for(int j=0;j<3;j++) {
						if(ppAcceptedDirections[j] == NULL) {
							ppAcceptedDirections[j] = ppDirections[i];
							if(debug_electrons.GetBool()){Msg("Assigned ppDirections[%i] to ppAcceptedDirections[%i]\n", i, j);}
							(*ppAcceptedDirections[j]) *= dPermittedVelocity;
							break;
						}
					}
				}
			}
			
			VectorNormalize(m_vecMomentum);
			m_vecMomentum *= dPermittedMomentum; //This is the value we'll give to the wave function - same direction, but updated magnitude according to quantised energy
			
			double pdProbabilities[3] = {Psi2(*ppAcceptedDirections[0]/25 + GetLocalOrigin()), //The division by 25 gives us the location it'd be at that velocity by our next think
											Psi2(*ppAcceptedDirections[1]/25 + GetLocalOrigin()),
											Psi2(*ppAcceptedDirections[2]/25 + GetLocalOrigin())};
			
			double dRandom = random->RandomFloat(0.0f, 1e-4); //Random number between 0 and 0.0001
			int iPossibility = -1;
			//Determine which direction to take (which has probability closest to random number)
			double pdProbDifferences[3] =	{pdProbabilities[0] - dRandom,
												pdProbabilities[1] - dRandom,
												pdProbabilities[2] - dRandom};
			for(int i=0; i<3; i++) {
				if(pdProbDifferences[i] < 0.0f) {
					//This way we get the absolute difference between each probability and the random number
					pdProbDifferences[i] *= -1.0f;
				}
			}
			//Select the least of the 3 differences
			if(pdProbDifferences[0] < pdProbDifferences[1]) {
				if(pdProbDifferences[0] < pdProbDifferences[2]) {
					//0 is the least
					iPossibility = 0;
				} else {
					//2 is the least
					iPossibility = 2;
				}
			} else {
				//1 is the least
				iPossibility = 1;
			}
			
			if(iPossibility >= 0) {
				vecForce += (*ppAcceptedDirections[iPossibility]) * GetParticleMass();
				if(trace_orbitals.GetBool()) {
					NDebugOverlay::Line(GetLocalOrigin(), GetLocalOrigin() + (*ppAcceptedDirections[iPossibility])/25.0f, 255, 255, 128, 0, 180.0f); //Keep the line there for 3 minutes?
				}
			} else {
				iPossibility = 0; //Avoid null-pointer on brightness setting below
				vecForce = vec3_origin; //We'll just go nowhere if we let it apply a non-zero force, so just keep getting readings until we go somewhere...
			}
			
			dFinalProbability = pdProbabilities[iPossibility];
		}
		VPhysicsGetObject()->ApplyForceCenter(vecForce);
		if(m_hSprite) {
			float flBrightnessAddition = 64 * clamp<float>(dFinalProbability*1e4,0,1);
			if(debug_electrons.GetBool()){Msg("%i: flBrightnessAddition is %.20f or %i, coming from probability of %.20f\n", entindex(), flBrightnessAddition, (int)flBrightnessAddition, dFinalProbability);}
			m_hSprite->SetBrightness(clamp<int>(192 + (int)flBrightnessAddition,0,255), 0.04f);
		}
	} else {
		//Negative or negligible kinetic energy - revert to classical Coulomb force
		Vector vecForce = vec3_origin;
		for(int i=0; i<iParticleCount; i++) {
			CBaseParticle *pParticle = g_vpParticles[i];
			if(pParticle && pParticle != this) {
				if(pParticle->GetLocalOrigin() != GetLocalOrigin()) {
					//This is the Q/r part of the formula, i.e. charge over distance
					Vector vecAddition = pParticle->GetLocalOrigin() - GetLocalOrigin();
					float flDistance = VectorNormalize(vecAddition) / DISTANCE_SCALE;
					vecAddition *= pParticle->GetCharge() / Square<float>(flDistance);
					vecForce += vecAddition;
				}
			}
		}
		//Force is the change in momentum in one second, but we deal in 1/25 of a second, so adjust the impulse accordingly
		vecForce /= 25;
		//Scale it
		vecForce *= MASS_SCALE * DISTANCE_SCALE / Square<float>(TIME_SCALE);
		VPhysicsGetObject()->ApplyForceCenter(vecForce);
		if(m_hSprite){m_hSprite->SetBrightness(192);}
	}
}

bool CElectron::AcceptEnergy(float flEnergy, float flTolerancePortion) {
	int iParticleCount = g_vpParticles.Count();
	float Z = 0; //Effective nuclear charge
	for(int i=0; i<iParticleCount; i++) {
		CBaseParticle *pParticle = g_vpParticles[i];
		if(pParticle && pParticle != this) {
			if(pParticle->GetLocalOrigin() != GetLocalOrigin()) {
				if(pParticle->IsProton()) {
					Z += 1;
				} else if(pParticle->IsElectron() && pParticle->m_iPrincipal <= m_iPrincipal) {
					//This involves using Slater's rules for the screening effect to determine an effective nuclear charge
					if(m_iPrincipal==1) {
						//1s-1s shielding
						Z -= 0.3f;
					} else if(m_iAzimuthal <= 1) {
						//s or p
						if(pParticle->m_iPrincipal == m_iPrincipal) {
							//Same group
							Z -= 0.35f;
						} else if(pParticle->m_iPrincipal == m_iPrincipal-1) {
							//Group just below
							Z -= 0.85f;
						} else if(pParticle->m_iPrincipal < m_iPrincipal-1) {
							//Groups further down
							Z -= 1;
						}
					} else if(m_iAzimuthal > 1) {
						//d or f (or higher)
						if(pParticle->m_iPrincipal == m_iPrincipal) {
							//Same subgroup
							Z -= 0.35f;
						} else if(pParticle->m_iAzimuthal < m_iAzimuthal || pParticle->m_iPrincipal < m_iPrincipal) {
							//(sub)groups below
							Z -= 1;
						}
					}
				}
			}
		}
	}
	float flNewEnergy = m_flTotalEnergy + flEnergy;
	int n, l, ml;
	n=1; l=ml=0; //Initialisation
	float ms = 0.5f;
	while (n < 10) { //Unrealistically high...
		//"Roll over" quantum numbers as necessary
		ml++;
		if(ml>l) { //-l<=ml<=l
			if(ms > -0.5f) {
				ms -= 1.0f; //Opposite spin (obeys Hund's Rule of Maximum Multiplicity)
			} else {
				l++; //Next Azimuthal (sublevel)
			}
			ml=-l;
			if(l==n) { //0<=l<n
				n++;
				l=0;
			}
		}
		
		float j = m_iAzimuthal + m_flSpinProjection;
		if(j<0){j=-j;}
		float m = m_iMagnetic + (j-m_iAzimuthal);
		float k = (j==m_iAzimuthal+0.5f)?(-j - 0.5f):(j + 0.5f);
		float absk = k>=0 ? k : -k;
		
		//This method of calculating kinetic energy is outlined on paper
		//double ZALPHA = -m_flPotential * ((rho/DISTANCE_SCALE) / (DIRAC*SPD_LGT));
		double ZALPHA = Z * ALPHA;
		double gamma = sqrt(Square<float>(absk) - Square<double>(ZALPHA));
		double dSquareInSqrt = Square<double>(ZALPHA / (m_iPrincipal - absk + gamma));
		double dRestCoeff = (1 / sqrt(1 + dSquareInSqrt));
		
		float flTestEnergy = m_flRestEnergy*dRestCoeff;
		if(flNewEnergy >= flTestEnergy*(1-flTolerancePortion) && flNewEnergy <= flTestEnergy*(1+flTolerancePortion)) {
			//Within tolerance threshold
			m_iPrincipal = n;
			m_iAzimuthal = l;
			m_iMagnetic = ml;
			m_flSpinProjection = ms;
			return true;
		}
	}
	return BaseClass::AcceptEnergy(flEnergy, flTolerancePortion);
}

class CWeakBoson : public CBaseParticle {
	DECLARE_CLASS(CWeakBoson, CBaseParticle);
	
public:
	float m_flDecayTime; //When to decay into a lepton and neutrino

	virtual void ApplyParticlePhysics();
	virtual void ParticleUpdateSprite();
	CWeakBoson();
	
	static CWeakBoson *CreateWeak(float flEnergy, int iCharge, Vector vecPosition, QAngle angDirection, CBaseParticle *pOwner);
};
LINK_ENTITY_TO_CLASS(weak_boson, CWeakBoson);

CWeakBoson::CWeakBoson() {
	m_flAppropriateMass = WEAK_BOSON_MASS;
	m_flCharge = ELEMENTARY_CHARGE;
	m_flDecayTime = 0.0f;
}

CWeakBoson *CWeakBoson::CreateWeak(float flEnergy, int iCharge, Vector vecPosition, QAngle angDirection, CBaseParticle *pOwner) {
	if(flEnergy < 0.0f) {
		Msg("Attempted to create weak boson with negative energy!\n");
		return NULL;
	}
	
	CWeakBoson *pWeak = (CWeakBoson *)CreateEntityByName("weak_boson");
	if(!pWeak) {
		return NULL;
	}
	
	pWeak->m_flCharge *= iCharge;
	pWeak->m_flKinetic = flEnergy; //The rest energy of a W boson is independent of this - it comes out of nowhere, or more specifically, the Higgs field!!!
	pWeak->SetOwnerEntity(pOwner);
	pWeak->SetAbsOrigin(vecPosition);
	pWeak->SetAbsAngles(angDirection);
	pWeak->m_flDecayTime = gpGlobals->curtime + WEAK_BOSON_LIFETIME;
	pWeak->Spawn();
	
	return pWeak;
}

void CWeakBoson::ParticleUpdateSprite() {
	BaseClass::ParticleUpdateSprite();

	if( !m_bColourUpToDate )
	{
		m_hSprite->SetColor(64,64,64); //Dark grey
		m_bColourUpToDate = true;
	}
}

void CWeakBoson::ApplyParticlePhysics() {
	//Calculate our potential energy
	//Coulombic:
	/*int iParticleCount = g_vpParticles.Count();
	int iProtons = 0;
	for(int i=0; i<iParticleCount; i++) {
		CBaseParticle *pParticle = g_vpParticles[i];
		if(pParticle && pParticle != this && pParticle->GetLocalOrigin() != GetLocalOrigin()) {
			//This is the Q/r part of the formula, i.e. charge over distance
			double dAddition = pParticle->GetCharge();
			if(dAddition>0) {iProtons++;} //Condensed onto one line for tidiness
			dAddition /= ((pParticle->GetLocalOrigin() - GetLocalOrigin()).Length()/DISTANCE_SCALE);
			m_flPotential += dAddition;
		}
	}
	m_flPotential *= GetCharge() * COULOMB;
	
	double dMass = m_flAppropriateMass / MASS_SCALE;
	float flRestEnergy = dMass * C_SQ;
	m_flPotential += flRestEnergy;
	
	m_flKinetic = m_flTotalEnergy - m_flPotential;
	if( m_flKinetic < 0.0f ) {
		Msg("Weak boson %i has negative kinetic energy! What's the deal with that?", entindex());
	} else {*/
	double dMass = m_flAppropriateMass / MASS_SCALE;
		Vector vecVelocity;
		VPhysicsGetObject()->GetVelocity(&vecVelocity, NULL);
		if(vecVelocity == vec3_origin) {
			AngleVectors(GetAbsAngles(), &vecVelocity);
			vecVelocity *= sqrt(m_flKinetic / (dMass/2)); //Ek = (m/2)(v^2)
			vecVelocity *= (DISTANCE_SCALE/TIME_SCALE);
			VPhysicsGetObject()->SetVelocity(&vecVelocity, &vec3_origin);
		}
	//}
	
	if(m_flDecayTime < gpGlobals->curtime) {
		float flBetaEnergy = random->RandomFloat(0,1) * m_flKinetic; //Random portion of decay energy to be imparted to electron/positron (rather than neutrino)
		float flNeutrinoEnergy = m_flKinetic - flBetaEnergy;
		
		CElectron *pBeta = (CElectron *)CreateEntityByName("particle_electron");
		pBeta->m_bAntiParticle = (GetCharge() > 0);
		pBeta->SetAbsOrigin(GetAbsOrigin());
		pBeta->m_vecMomentum = vecVelocity; //Keep it going in the same direction
		pBeta->Spawn();
		pBeta->SetFreeFor(15.0f);
		pBeta->m_flTotalEnergy = flBetaEnergy;
		
		QAngle angNeutrinoDirection = vec3_angle;
		VectorAngles(vecVelocity, angNeutrinoDirection);
		angNeutrinoDirection.y += 90;
		CPhoton::CreateNeutrino(flNeutrinoEnergy, GetAbsOrigin(), angNeutrinoDirection, pBeta);
		
		UTIL_Remove(this);
	}
}

static ConVar hadron_mode("hadron_mode", "1", 0, "How to display hadrons: 0 = only outlines; 1 = outlines and quarks; 2 = only quarks");

class CQuark : public CBaseParticle {
	DECLARE_CLASS(CQuark, CBaseParticle);
	
public:
	DECLARE_DATADESC();
	
	//Flavour quantum numbers
	float m_flIsospin;
	float m_flBaryonNo;
	int m_iStrangeness;
	int m_iCharm;
	int m_iBeauty; //Bottom quark
	int m_iTruth; //Top quark

	CQuark();
	
	virtual bool CreateVPhysics();
	virtual void ApplyParticlePhysics();
	virtual void ParticleUpdateSprite();
	
	virtual bool IsQuark(){return true;}
};
LINK_ENTITY_TO_CLASS(particle_quark, CQuark);

BEGIN_DATADESC(CQuark)
	DEFINE_FIELD(m_flIsospin, FIELD_FLOAT),
	DEFINE_FIELD(m_flBaryonNo, FIELD_FLOAT),
	DEFINE_FIELD(m_iStrangeness, FIELD_INTEGER),
	DEFINE_FIELD(m_iCharm, FIELD_INTEGER),
	DEFINE_FIELD(m_iBeauty, FIELD_INTEGER),
	DEFINE_FIELD(m_iTruth, FIELD_INTEGER),
END_DATADESC()

CQuark::CQuark() {
	m_bStrongInt = true;
	
	m_flIsospin = 0.5f;
	m_flBaryonNo = 1.0f/3.0f;
	m_iStrangeness = 0;
	m_iCharm = 0;
	m_iBeauty = 0;
	m_iTruth = 0;
}

bool CQuark::CreateVPhysics() {
	SetSolid(SOLID_NONE);
	SetMoveType(MOVETYPE_NONE);
	return false;
}

void CQuark::ApplyParticlePhysics() {
	//Quarks in this code let their parent hadrons take care of most of their physics.
	//However, it is possible for quarks to change into antiquarks and vice-versa.
	//This updates the quantum numbers if that happens.
	bool bWasAnti = m_bAntiParticle;
	m_bAntiParticle = (m_iColour < 0);
	if(bWasAnti != m_bAntiParticle) {
		//Negate quantum numbers
		m_flIsospin *= -1.0f;
		m_flBaryonNo *= -1.0f;
		m_iStrangeness *= -1;
		m_iCharm *= -1;
		m_iBeauty *= -1;
		m_iTruth *= -1;
	}
	
	//Emergency situations can arise in which a quark floats off into space...
	if(GetLocalOrigin().LengthSqr() > 45000.0f) {
		SetAbsVelocity(vec3_origin);
		SetLocalVelocity(vec3_origin);
		Vector vecOrigin = GetLocalOrigin();
		VectorNormalize(vecOrigin);
		vecOrigin *= 40.0f;
		SetLocalOrigin(vecOrigin);
	}
}

void CQuark::ParticleUpdateSprite() {
#if CULL_SPRITES
	//Instead of turning the sprite off when the user doesn't want it, we delete it - this saves on edicts, important for bigger atoms.
	if(m_hSprite && hadron_mode.GetInt() < 1) {
		UTIL_Remove(m_hSprite);
	} else if (hadron_mode.GetInt() > 0) {
		BaseClass::ParticleUpdateSprite();
	}
#else
	if(m_hSprite) {
		m_hSprite->SetBrightness( (hadron_mode.GetInt() < 1) ? 0 : 255 );
	}
	BaseClass::ParticleUpdateSprite();
#endif
}

class CBaseHadron : public CBaseParticle { //Particles composed of quarks
	DECLARE_CLASS(CBaseHadron, CBaseParticle);
	
public:
	CBaseHadron();
	
	virtual float GetCharge(float *pflChargeInUnits = NULL);
	void UpdateOnRemove();
	virtual void ApplyParticlePhysics();
	virtual void ParticleUpdateSprite();
	virtual void QuarkInteraction();
	CUtlVector< CHandle<CQuark> > m_vhQuarks;
	bool m_bQuarkIntDisabled;
	bool m_bTargettedByMeson;
	float m_flNextQuarkIntTime;
	int m_iQuarkToStop; //Allow us to move a quark a certain distance in a set amount of time
	
	DECLARE_DATADESC();
};

BEGIN_DATADESC(CBaseHadron)
	DEFINE_UTLVECTOR(m_vhQuarks, FIELD_EHANDLE),
	DEFINE_FIELD(m_bQuarkIntDisabled, FIELD_BOOLEAN),
	DEFINE_FIELD(m_bTargettedByMeson, FIELD_BOOLEAN),
	DEFINE_FIELD(m_flNextQuarkIntTime, FIELD_TIME),
	DEFINE_FIELD(m_iQuarkToStop, FIELD_TIME),
END_DATADESC()

CBaseHadron::CBaseHadron() {
	m_flNextQuarkIntTime = 0.0f;
	m_bQuarkIntDisabled = false;
	m_bTargettedByMeson = false;
	m_iQuarkToStop = -1;
}

float CBaseHadron::GetCharge(float *pflChargeInUnits){
	//Gell-Mann–Nishijima formula
	//Msg("Hadron %i: Starting Gell-Mann–Nishijima formula\n", entindex());
	float I=0, B=0;
	int S=0, C=0, B1=0, T=0;
	for(int i=0;i<m_vhQuarks.Count();i++) {
		CQuark *pQuark = m_vhQuarks[i].Get();
		I += pQuark->m_flIsospin;
		B += pQuark->m_flBaryonNo;
		S += pQuark->m_iStrangeness;
		C += pQuark->m_iCharm;
		B1+= pQuark->m_iBeauty;
		T += pQuark->m_iTruth;
	}
	float Q = I + (B+S+C+B1+T) / 2;
	if( pflChargeInUnits ) {
		*pflChargeInUnits = Q;
	}
	//if(debug_nucleus.GetBool()){Msg("Hadron %i: Charge number is %.1f == %.1f + (%.1f + %i + %i + %i + %i)/2\n", entindex(), Q, I, B, S, C, B1, T);}
	return Q * ELEMENTARY_CHARGE;
}

void CBaseHadron::UpdateOnRemove() {
	for(int i=0;i<m_vhQuarks.Count();i++) {
		UTIL_Remove(m_vhQuarks[i]);
	}
	m_vhQuarks.RemoveAll();
	m_vhQuarks.Purge();
}

void CBaseHadron::ApplyParticlePhysics() {
	if(hadron_mode.GetInt() >= 1) {
		//No point in simulating quark interactions if the quarks themselves aren't being rendered!
		if(gpGlobals->curtime > m_flNextQuarkIntTime && !m_bQuarkIntDisabled) {
			QuarkInteraction();
		} else if (gpGlobals->curtime > m_flNextQuarkIntTime - 1.5f) {
			//This happens 0.5 seconds after the last interaction, i.e. half-way through
			for(int i=0;i<m_vhQuarks.Count();i++) {
				//Tell each quark to update its sprite's colour
				m_vhQuarks[i]->m_bColourUpToDate = false;
			}
		}
	}
}

void CBaseHadron::ParticleUpdateSprite() {
#if CULL_SPRITES
	//Instead of turning the sprite off when the user doesn't want it, we delete it - this saves on edicts, important for bigger atoms.
	if(m_hSprite && hadron_mode.GetInt() > 1) {
		UTIL_Remove(m_hSprite);
	} else if (hadron_mode.GetInt() < 2) {
		BaseClass::ParticleUpdateSprite();
	}
#else
	if(m_hSprite) {
		m_hSprite->SetBrightness( (hadron_mode.GetInt() > 1) ? 0 : 255 );
	}
	BaseClass::ParticleUpdateSprite();
#endif
}

void CBaseHadron::QuarkInteraction() {
	if(m_iQuarkToStop >= 0) {
		m_vhQuarks[m_iQuarkToStop]->SetLocalVelocity(vec3_origin);
		m_vhQuarks[m_iQuarkToStop]->SetMoveType(MOVETYPE_NONE);
		m_iQuarkToStop = -1;
	}
	
	CQuark *pFirstQuark;
	CQuark *pSecondQuark;
	if(m_vhQuarks.Count() <= 2 ){
		//Meson
		pFirstQuark = m_vhQuarks[0];
		pSecondQuark = m_vhQuarks[1];
	} else {
		//Baryon
		//Randomly select which quark will *not* be taking part in the interaction
		switch(random->RandomInt(0, 2)) {
			case 0:
				pFirstQuark = m_vhQuarks[1];
				pSecondQuark = m_vhQuarks[2];
				break;
			case 1:
				pFirstQuark = m_vhQuarks[0];
				pSecondQuark = m_vhQuarks[2];
				break;
			case 2:
				pFirstQuark = m_vhQuarks[0];
				pSecondQuark = m_vhQuarks[1];
		}
	}
	
	int iGluonColour1 = pFirstQuark->m_iColour;
	int iGluonColour2 = pSecondQuark->m_iColour;
	if(iGluonColour1*iGluonColour2 > 0) {
		//If the product is positive, it means they have the same sign...
		//Randomly select which will be the gluon's anticolour
		switch(random->RandomInt(1, 2)) {
			case 1:
				iGluonColour1 *= -1;
				break;
			case 2:
				iGluonColour2 *= -1;
		}
	}
	
	Vector vecBeamOffset = Vector(0, 10, 0);
	CBeam *pGluon1 = CBeam::BeamCreate( "sprites/glow_test02.vmt", 3 );
	CBeam *pGluon2 = CBeam::BeamCreate( "sprites/glow_test02.vmt", 3 );
	pGluon1->SetEndWidth( pGluon1->GetWidth() );
	pGluon2->SetEndWidth( pGluon2->GetWidth() );
	switch (iGluonColour1)
	{
		case CHARGE_ANTIBLUE:
			pGluon1->SetColor(255,255,0);
			break;
		case CHARGE_ANTIGREEN:
			pGluon1->SetColor(255,0,255);
			break;
		case CHARGE_ANTIRED:
			pGluon1->SetColor(0,255,255);
			break;
		case CHARGE_RED:
			pGluon1->SetColor(255,0,0);
			break;
		case CHARGE_GREEN:
			pGluon1->SetColor(0,255,0);
			break;
		case CHARGE_BLUE:
			pGluon1->SetColor(0,0,255);
			break;
	}
	switch (iGluonColour2)
	{
		case CHARGE_ANTIBLUE:
			pGluon2->SetColor(255,255,0);
			break;
		case CHARGE_ANTIGREEN:
			pGluon2->SetColor(255,0,255);
			break;
		case CHARGE_ANTIRED:
			pGluon2->SetColor(0,255,255);
			break;
		case CHARGE_RED:
			pGluon2->SetColor(255,0,0);
			break;
		case CHARGE_GREEN:
			pGluon2->SetColor(0,255,0);
			break;
		case CHARGE_BLUE:
			pGluon2->SetColor(0,0,255);
			break;
	}
	pGluon1->PointsInit(pFirstQuark->GetAbsOrigin() + vecBeamOffset, pSecondQuark->GetAbsOrigin() + vecBeamOffset);
	pGluon2->PointsInit(pFirstQuark->GetAbsOrigin() - vecBeamOffset, pSecondQuark->GetAbsOrigin() - vecBeamOffset);
	/*pGluon1->EntsInit(pFirstQuark, pSecondQuark);
	pGluon2->EntsInit(pFirstQuark, pSecondQuark);
	//Offset them a bit so they are distinct
	pGluon1->SetStartPos(vecBeamOffset);
	pGluon1->SetEndPos(vecBeamOffset);
	pGluon1->RelinkBeam();
	pGluon2->SetStartPos(-vecBeamOffset);
	pGluon2->SetEndPos(-vecBeamOffset);
	pGluon2->RelinkBeam();*/
	pGluon1->LiveForTime(1);
	pGluon2->LiveForTime(1);
	
	if(m_vhQuarks.Count() > 2 ) {
		//Baryon - we don't swap the colours in a meson because that would swap the quark and antiquark, introducing excessive confusion.
		int iQuark1Colour = pFirstQuark->m_iColour;
		pFirstQuark->m_iColour = pSecondQuark->m_iColour;
		pSecondQuark->m_iColour = iQuark1Colour;
	}
	
	m_flNextQuarkIntTime = gpGlobals->curtime + 2.0f;
}

class CMeson : public CBaseHadron {
	DECLARE_CLASS(CMeson, CBaseHadron);
	
public:
	DECLARE_DATADESC();
	
	CHandle<CBaseHadron> m_hTarget;
	CHandle<CQuark> m_hQuarkTarget;
	CHandle<CQuark> m_hMyQuark;
	CHandle<CQuark> m_hAntiquark;
	bool m_bHitTarget;
	
	float m_flTimeToMove;
	
	Vector m_vecQuarkTargetOrigin;
	
	CMeson();
	
	void Spawn();
	void UpdateOnRemove();
	static CMeson *CreateMeson(CQuark *pFirstQuark, CQuark *pSecondQuark, CBaseHadron *pTarget, Vector vecLocation, float flMovementDelay = 1.0f);
	
	virtual void ApplyParticlePhysics();
	virtual void ParticleUpdateSprite();
};

LINK_ENTITY_TO_CLASS(particle_pion, CMeson);

BEGIN_DATADESC(CMeson)
	DEFINE_FIELD(m_hTarget, FIELD_EHANDLE),
	DEFINE_FIELD(m_hQuarkTarget, FIELD_EHANDLE),
	DEFINE_FIELD(m_hAntiquark, FIELD_EHANDLE),
	DEFINE_FIELD(m_bHitTarget, FIELD_BOOLEAN),
	DEFINE_FIELD(m_flTimeToMove, FIELD_TIME),
	DEFINE_FIELD(m_vecQuarkTargetOrigin, FIELD_VECTOR),
END_DATADESC()

CMeson::CMeson() {
	m_hTarget.Set(NULL);
	m_hQuarkTarget.Set(NULL);
	m_hMyQuark.Set(NULL);
	m_hAntiquark.Set(NULL);
	m_bHitTarget = false;
	m_vecQuarkTargetOrigin = vec3_origin;
	m_flAppropriateMass = NEUTRAL_PION_MASS;
}

CMeson *CMeson::CreateMeson(CQuark *pFirstQuark, CQuark *pSecondQuark, CBaseHadron *pTarget, Vector vecLocation, float flMovementDelay) {
	CMeson *pMeson = (CMeson *)CreateEntityByName("particle_pion");
	if(!pMeson) {
		return NULL;
	}
	pTarget->m_bTargettedByMeson = true;
	
	pMeson->m_vhQuarks.AddToHead(pFirstQuark);
	pMeson->m_vhQuarks.AddToHead(pSecondQuark);
	pMeson->m_hTarget = pTarget;
	pMeson->m_flNextQuarkIntTime = pMeson->m_flTimeToMove = gpGlobals->curtime + flMovementDelay; //Block quark interactions too - otherwise they will immediately change places and mess up the visualisation
	pMeson->SetAbsOrigin(vecLocation);
	if(debug_nucleus.GetBool()){Msg("Created meson %i at (%.2f, %.2f, %.2f)\n", pMeson->entindex(), vecLocation.x, vecLocation.y, vecLocation.z);}
	pMeson->Spawn();
	
	return pMeson;
}

void CMeson::Spawn() {
	if(!m_vhQuarks.Count()) {
		Msg("WARNING: Pion created with no quarks!\n");
		UTIL_Remove(this);
		return;
	}
	for(int i=0;i<m_vhQuarks.Count();i++) {
		CQuark *pQuark = m_vhQuarks[i].Get();
		pQuark->SetParent(this);
		pQuark->SetOwnerEntity(this);
		pQuark->SetMoveType(MOVETYPE_NOCLIP); //We will be moving them around.
		if(pQuark->m_hSprite) {
			pQuark->m_hSprite->SetColor(255,255,255);
		}
		pQuark->m_bColourUpToDate = false; //This seems to be the only way to make the sprite move properly.
		if(debug_nucleus.GetBool()){Msg("Quark %i attached to meson %i - it is now located at (%.2f, %.2f, %.2f)\n", pQuark->entindex(), entindex(), pQuark->GetAbsOrigin().x, pQuark->GetAbsOrigin().y, pQuark->GetAbsOrigin().z);}
	}
	
	BaseClass::Spawn();
	
	if(VPhysicsGetObject()) {
		VPhysicsGetObject()->EnableCollisions(false);
	} else {
		if(debug_nucleus.GetBool()) {Msg("WARNING: Meson %i failed to create physics object!\n", entindex());}
	}
}

void CMeson::UpdateOnRemove() {
	if(m_hTarget) {
		m_hTarget->m_bTargettedByMeson = false;
	}
	if(m_hSprite) {
		UTIL_Remove(m_hSprite);
	}
	BaseClass::UpdateOnRemove();
}

static ConVar meson_speed("meson_speed", "250.0");

void CMeson::ApplyParticlePhysics() {
	if(!(m_vhQuarks.Count())) {
		if(!m_hMyQuark) {
			UTIL_Remove(this);
			return;
		}
		//This indicates we're already finished, but we still have to move that last quark into position.
		Vector vecQuarkDistance = m_vecQuarkTargetOrigin - m_hMyQuark->GetLocalOrigin();
		if(vecQuarkDistance.LengthSqr() > 1) {
			//More than 1 unit out of position
			if(m_hMyQuark->GetLocalVelocity() == vec3_origin || DotProduct(m_hMyQuark->GetLocalVelocity(), vecQuarkDistance) < 0.9f) {
				//It's going the wrong way - redirect it!
				VectorNormalize(vecQuarkDistance);
				m_hMyQuark->SetLocalVelocity(vecQuarkDistance * 15.0f);
			}
		} else {
			//STOP!
			m_hMyQuark->SetLocalVelocity(vec3_origin);
			m_hMyQuark->SetMoveType(MOVETYPE_NONE);
			UTIL_Remove(this); //We're done here!
		}
		return;
	}
	BaseClass::ApplyParticlePhysics(); //This handles quark interaction scheduling
	if(!m_bHitTarget) {
		if(!m_hTarget) {
			Msg("WARNING: Meson %i's target has disappeared", entindex());
			UTIL_Remove(this);
			return;
		}
		if(m_flTimeToMove < gpGlobals->curtime) {
			/*for(int i=0;i<m_vhQuarks.Count();i++) {
				//m_vhQuarks[i]->m_bScaleUpToDate = false; //This seems to be the only way to make the sprite move properly.
				if(m_vhQuarks[i]->m_hSprite)
					{m_vhQuarks[i]->m_hSprite->SetTransmitState( FL_EDICT_ALWAYS );} //This seems to be the only way to make the sprite move properly.
				if(debug_nucleus.GetBool())
					{Msg("Quark %i belongs to meson %i and is located at (%.2f, %.2f, %.2f)\n", m_vhQuarks[i]->entindex(), entindex(), m_vhQuarks[i]->GetAbsOrigin().x, m_vhQuarks[i]->GetAbsOrigin().y, m_vhQuarks[i]->GetAbsOrigin().z);}
			}*/
			//Move towards my target
			Vector vecDirection = (m_hTarget->GetAbsOrigin() - GetAbsOrigin());
			if(vecDirection.LengthSqr() <= 100) {
				//We are within 10 units of each other!
				m_bHitTarget = true;
			} else {
				VectorNormalize(vecDirection);
				vecDirection *= meson_speed.GetFloat();
				if(m_bUsesVPhysics) {
					VPhysicsGetObject()->SetVelocity(&vecDirection, NULL);
				} else {
					if(debug_nucleus.GetBool()) {Msg("Meson %i not using VPhysics for some reason\n", entindex());}
					SetAbsVelocity(vecDirection);
				}
			}
		}
	} else {
		if(m_bUsesVPhysics) {
			//Initial code when we hit our target - turn off our physics, but also do a couple of other things
			m_bUsesVPhysics = false;
			SetMoveType(MOVETYPE_NONE);
			SetSolid(SOLID_NONE);
			VPhysicsDestroyObject();
			
			if(hadron_mode.GetInt() < 1) {
				UTIL_Remove(this);
				return;
			}
			if(!m_hTarget) {
				Msg("WARNING: Meson %i's target has disappeared", entindex());
				UTIL_Remove(this);
				return;
			}
			
			//Turn off quark interactions in all particles concerned
			m_hTarget->m_bQuarkIntDisabled = m_bQuarkIntDisabled = true;
			
			//Select a quark to annihilate in the target hadron
			for(int i=0;i<m_vhQuarks.Count();i++) {
				m_vhQuarks[i]->SetParent(m_hTarget);
				m_vhQuarks[i]->SetOwnerEntity(m_hTarget);
				if(m_vhQuarks[i]->m_bAntiParticle != m_hTarget->m_bAntiParticle) {
					m_hAntiquark = m_vhQuarks[i];
				}
			}
			if(!m_hAntiquark) {
				if(debug_nucleus.GetBool()) {Msg("WARNING: Meson %i can't find its antiquark!!!\n", entindex());}
				return;
			}
			
			for(int j=0;j<m_hTarget->m_vhQuarks.Count();j++) {
				CQuark *pQuark = m_hTarget->m_vhQuarks[j];
				if(pQuark->m_iColour + m_hAntiquark->m_iColour != 0 && pQuark->m_flIsospin + m_hAntiquark->m_flIsospin == 0.0f) {
					//The colours don't sum to zero, meaning the quark's colour is different from the antiquark's anticolour.
					//The isospins *do* sum to zero, meaning the quark's flavour is *the same* as the antiquark's antiflavour.
					m_hQuarkTarget = pQuark;
				}
			}
			
			if(!m_hQuarkTarget) {
				//We've failed to select a quark to annihilate, so do an emergency swap of our quark and antiquark, then try again.
				QuarkInteraction();
				for(int i=0;i<m_vhQuarks.Count();i++) {
					m_vhQuarks[i]->m_bColourUpToDate = false; //Update colours straight away, since it's an emergency swap
					m_vhQuarks[i]->ApplyParticlePhysics(); //Update quark/antiquark status straight away!
					if(m_vhQuarks[i]->m_bAntiParticle != m_hTarget->m_bAntiParticle) {
						m_hAntiquark = m_vhQuarks[i];
					}
				}
				for(int j=0;j<m_hTarget->m_vhQuarks.Count();j++) {
					CQuark *pQuark = m_hTarget->m_vhQuarks[j];
					if(pQuark->m_iColour + m_hAntiquark->m_iColour != 0 && pQuark->m_flIsospin + m_hAntiquark->m_flIsospin == 0.0f) {
						//The colours don't sum to zero, meaning the quark's colour is different from the antiquark's anticolour.
						//The isospins *do* sum to zero, meaning the quark's flavour is *the same* as the antiquark's antiflavour.
						m_hQuarkTarget = pQuark;
					}
				}
				
				if(!m_hQuarkTarget) {
					//Bottle it!
					Msg("WARNING: Meson %i couldn't find a quark to annihilate in hadron %i!\n", entindex(), m_hTarget->entindex());
					/*for(int i=0;i<m_vhQuarks.Count();i++) {
						UTIL_Remove(m_vhQuarks[i]);
					}*/
					m_hTarget->m_bQuarkIntDisabled = false;
					UTIL_Remove(this);
					return;
				}
			}
		} else if (hadron_mode.GetInt() < 1) {
			UTIL_Remove(this);
			return;
		}
		
		//Move our antiquark towards the quark to be annihilated
		Vector vecQuarkDistance = m_hQuarkTarget->GetAbsOrigin() - m_hAntiquark->GetAbsOrigin();
		if(vecQuarkDistance.LengthSqr() > 25) {
			//More than 5 units apart
			if(m_hAntiquark->GetAbsVelocity() == vec3_origin || DotProduct(m_hAntiquark->GetAbsVelocity(), vecQuarkDistance) < 0.9f) {
				//It's going the wrong way - redirect it!
				VectorNormalize(vecQuarkDistance);
				m_hAntiquark->SetAbsVelocity(vecQuarkDistance * 25.0f); //Move at 25 units per second
			}
		} else {
			CQuark *pChangeColour;
			for(int i=0;i<m_hTarget->m_vhQuarks.Count();i++) {
				CQuark *pQuark = m_hTarget->m_vhQuarks[i];
				if(pQuark->m_iColour + m_hAntiquark->m_iColour == 0) {
					//The colours sum to zero, so we can send a gluon to this quark to change its colour
					pChangeColour = pQuark;
				}
			}
			if(hadron_mode.GetInt() >= 1) {
				CBeam *pGluon1 = CBeam::BeamCreate( "sprites/glow_test02.vmt", 3 );
				CBeam *pGluon2 = CBeam::BeamCreate( "sprites/glow_test02.vmt", 3 );
				pGluon1->SetEndWidth( pGluon1->GetWidth() );
				pGluon2->SetEndWidth( pGluon2->GetWidth() );
				switch (m_hAntiquark->m_iColour)
				{
					case CHARGE_ANTIBLUE:
						pGluon1->SetColor(255,255,0);
						break;
					case CHARGE_ANTIGREEN:
						pGluon1->SetColor(255,0,255);
						break;
					case CHARGE_ANTIRED:
						pGluon1->SetColor(0,255,255);
						break;
					case CHARGE_RED:
						pGluon1->SetColor(255,0,0);
						break;
					case CHARGE_GREEN:
						pGluon1->SetColor(0,255,0);
						break;
					case CHARGE_BLUE:
						pGluon1->SetColor(0,0,255);
						break;
				}
				switch (m_hQuarkTarget->m_iColour)
				{
					case CHARGE_ANTIBLUE:
						pGluon2->SetColor(255,255,0);
						break;
					case CHARGE_ANTIGREEN:
						pGluon2->SetColor(255,0,255);
						break;
					case CHARGE_ANTIRED:
						pGluon2->SetColor(0,255,255);
						break;
					case CHARGE_RED:
						pGluon2->SetColor(255,0,0);
						break;
					case CHARGE_GREEN:
						pGluon2->SetColor(0,255,0);
						break;
					case CHARGE_BLUE:
						pGluon2->SetColor(0,0,255);
						break;
				}
				//EntsInit seems unsuitable, perhaps because the annihilating quarks are removed in the next tick
				pGluon1->PointsInit(m_hAntiquark->GetAbsOrigin(), pChangeColour->GetAbsOrigin());
				pGluon2->PointsInit(m_hQuarkTarget->GetAbsOrigin(), pChangeColour->GetAbsOrigin());
				pGluon1->LiveForTime(1);
				pGluon2->LiveForTime(1);
			}
			pChangeColour->m_iColour = m_hQuarkTarget->m_iColour;
			
			for(int i=0;i<m_vhQuarks.Count();i++) {
				if(m_vhQuarks[i]->m_bAntiParticle == m_hTarget->m_bAntiParticle) {
					//This time select the one that's the same "antiness" as the target - it needs to be added into the hadron
					m_hMyQuark = m_vhQuarks[i];
				}
			}
			if(!m_hMyQuark) {
				if(debug_nucleus.GetBool()) {Msg("WARNING: Meson %i can't find its non-antiquark!!!\n", entindex());}
				return;
			}
			int iQuarkIndex = m_hTarget->m_vhQuarks.Find(m_hQuarkTarget);
			m_hTarget->m_vhQuarks[iQuarkIndex] = m_hMyQuark;
			
			m_vecQuarkTargetOrigin = m_hQuarkTarget->GetLocalOrigin();
			
			//Annihilation
			UTIL_Remove(m_hQuarkTarget);
			UTIL_Remove(m_hAntiquark);
			
			//Turn quark interactions back on in the target.
			m_hTarget->m_flNextQuarkIntTime = gpGlobals->curtime + 2.0f;
			m_hTarget->m_bQuarkIntDisabled = false;
			
			//Finally, we have to pop the quark into the position previously filled by the now-annihilated quark
			m_vhQuarks.RemoveAll();
			m_hMyQuark->SetLocalVelocity(m_vecQuarkTargetOrigin - m_hQuarkTarget->GetLocalOrigin());
		}
	}
}

void CMeson::ParticleUpdateSprite() {
	BaseClass::ParticleUpdateSprite();

	if( m_hSprite && !m_bColourUpToDate )
	{
		if( !m_bHitTarget ) {
			m_hSprite->SetColor(64,64,64); //Same as W boson
		} else {
			m_hSprite->SetBrightness(0, 0.15f); //Fade out once we've hit our target
		}
		m_bColourUpToDate = true;
	}
}

static bool g_bBetaPermitted = false; //Becomes true when the command is given, then becomes false as soon as some particle decays.
static void PermitBeta(const CCommand &args) {
	g_bBetaPermitted = true;
}
static ConCommand beta_decay("beta_decay", PermitBeta, "Causes a particle to beta-decay if it is energetically favourable to do so");

class CNucleon : public CBaseHadron { //A base class for nucleons that puts it on the list of particles in the nucleus
	DECLARE_CLASS(CNucleon, CBaseHadron);
	
public:
	DECLARE_DATADESC();
	
	bool m_bProton; //True if proton, false if neutron
	CNucleon();
	~CNucleon();
	
	void BetaDecay(float flEnergy);
	void Spawn();
	void UpdateOnRemove();
	
	virtual void ApplyParticlePhysics();
	virtual void ParticleUpdateSprite();
	virtual bool IsProton(){return m_bProton;}
};

LINK_ENTITY_TO_CLASS(particle_proton, CNucleon);
LINK_ENTITY_TO_CLASS(particle_neutron, CNucleon);

BEGIN_DATADESC(CNucleon)
	DEFINE_FIELD(m_bProton, FIELD_BOOLEAN),
END_DATADESC()

static CUtlVector<CNucleon *> g_vpNucleons;

ConVar use_wapstra_binding_energy( "use_wapstra_binding_energy", "0" );

static double NuclearBindingEnergy(int iProtons = -1) {
	int Z = iProtons;
	if(Z < 0) { //Use the actual number of protons
		Z = 0;
		for(int i=0;i<g_vpNucleons.Count();i++) {
			if(g_vpNucleons[i]->IsProton()) {
				Z++;
			}
		}
	}
	int A = g_vpNucleons.Count();
	int N = A - Z;
	
	//Liquid-drop model
	//Constants
	float aV, aS, aC, aA, delta;
	if(use_wapstra_binding_energy.GetBool()) {
		//Wapstra constants
		aV = 14.1;
		aS = 13;
		aC = 0.595;
		aA = 19;
		delta = 33.5;
	} else {
		//Rohlf constants
		aV = 15.75;
		aS = 17.8;
		aC = 0.711;
		aA = 23.7;
		delta = 11.18;
	}
	
	//Terms
	float flVolume = aV * A;
	float flSurface = aS * pow(A, 2.0f/3.0f);
	float flCoulomb = aC * Z*(Z-1) / pow(A, 1.0f/3.0f);
	float flAsymmetry = aA * Square<int>(A - 2*Z) / A;
	float flPairing = (A%2) ? 0 : ((Z%2) ? delta : -delta); //%2 checks if odd - if A is odd, it's zero, otherwise if Z is odd, it's +delta, otherwise it's -delta
	
	return (flVolume - flSurface - flCoulomb - flAsymmetry - flPairing) * ELEMENTARY_CHARGE * 1e6; //This was all in MeV, so multiply by proton/electron charge x 10^6 to get it in joules
}

double TotalNuclearEnergy(int iProtons) {
	int Z = iProtons;
	if(Z < 0) { //Use the actual number of protons
		Z = 0;
		for(int i=0;i<g_vpNucleons.Count();i++) {
			if(g_vpNucleons[i]->IsProton()) {
				Z++;
			}
		}
	}
	int A = g_vpNucleons.Count();
	int N = A - Z;
	//if(debug_nucleus.GetBool()){Msg("TotalNuclearEnergy: A=%i, N=%i, Z=%i\n", A, N, Z);}
	return ( (Z * (PROTON_MASS/MASS_SCALE)) + (N * (NEUTRON_MASS/MASS_SCALE)) ) * C_SQ - NuclearBindingEnergy(Z);
}

static void RecalculateNuclearOrigin(void) {
	g_flNextNuclearOriginTime = gpGlobals->curtime + 1.0f; //Every second
	
	g_vecNuclearOrigin = vec3_origin;
	for(int i=0;i<g_vpNucleons.Count();i++) {
		g_vecNuclearOrigin += g_vpNucleons[i]->GetAbsOrigin();
	}
	g_vecNuclearOrigin /= g_vpNucleons.Count();
	
	if(debug_nucleus.GetBool()){Msg("Nuclear Origin is (%.2f, %.2f, %.2f)\n", g_vecNuclearOrigin.x, g_vecNuclearOrigin.y, g_vecNuclearOrigin.z);}
	if(draw_bohr_radius.GetBool()) {
		NDebugOverlay::Line(g_vecNuclearOrigin + Vector(A0,0,0), g_vecNuclearOrigin - Vector(A0,0,0), 255, 0, 64, 0, 1.0f);
		NDebugOverlay::Line(g_vecNuclearOrigin + Vector(0,A0,0), g_vecNuclearOrigin - Vector(0,A0,0), 255, 0, 64, 0, 1.0f);
		NDebugOverlay::Line(g_vecNuclearOrigin + Vector(0,0,A0), g_vecNuclearOrigin - Vector(0,0,A0), 255, 0, 64, 0, 1.0f);
	}
}

CNucleon::CNucleon() {
	g_vpNucleons.AddToHead(this);
	m_flAppropriateMass = PROTON_MASS;
	m_flCharge = ELEMENTARY_CHARGE;
}

CNucleon::~CNucleon() {
	g_vpNucleons.FindAndRemove(this);
}

void CNucleon::Spawn() {
	m_bProton = FClassnameIs(this, "particle_proton");
	if( !m_bProton ) {
		m_flCharge = 0.0f;
		m_flAppropriateMass = NEUTRON_MASS;
	}
	
	//Create quarks:
	for (int i=1;i<=3;i++) {
		CQuark *pQuark = (CQuark *)CreateEntityByName("particle_quark");
		if(debug_nucleus.GetBool()){Msg("Nucleon %i created quark %i\n", entindex(), pQuark->entindex());}
		pQuark->SetOwnerEntity(this);
		pQuark->SetParent(this);
		if(debug_nucleus.GetBool()){Msg("Quark %i's owner and parent set\n", pQuark->entindex());}
		pQuark->m_iColour = i; //1=red, 2=green, 3=blue
		if(debug_nucleus.GetBool()){Msg("Quark %i has colour %i\n", pQuark->entindex(), i);}
		Vector vecOrigin;
		AngleVectors(QAngle(120*i, 0, 0), &vecOrigin);
		vecOrigin *= 40; //Keep it far enough from the centre so they are disinct
		if(debug_nucleus.GetBool()){Msg("Quark %i's origin determined...\n", pQuark->entindex());}
		pQuark->SetLocalOrigin(vecOrigin);
		if(debug_nucleus.GetBool()){Msg("Quark %i's origin set\n", pQuark->entindex());}
		pQuark->m_flAppropriateMass = 17.83f; //Scaled rest mass of up/down quark
		if(debug_nucleus.GetBool()){Msg("Quark %i's mass set\n", pQuark->entindex());}
		if((m_bProton && i==3) || (!m_bProton && i>=2)) { //Third in proton, both second and third in neutron
			pQuark->m_flIsospin = -0.5f; //Down quark
		} else {
			pQuark->m_flIsospin = 0.5f; //Up quark
		}
		if(debug_nucleus.GetBool()){Msg("Quark %i's isospin set\n", pQuark->entindex());}
		m_vhQuarks.AddToHead(pQuark);
		if(debug_nucleus.GetBool()){Msg("Quark %i placed in vector...\n", pQuark->entindex());}
		pQuark->Spawn();
		if(debug_nucleus.GetBool()){Msg("Quark %i spawned!\n", pQuark->entindex());}
	}
	
	BaseClass::Spawn();
	
	if(VPhysicsGetObject()) {
		VPhysicsGetObject()->EnableCollisions(false); //Nucleons are very close together, so they *will* collide
	}
}

void CNucleon::UpdateOnRemove() {
	g_vpNucleons.FindAndRemove(this);
	BaseClass::UpdateOnRemove();
}

void CNucleon::BetaDecay(float flEnergy) {
	g_bBetaPermitted = false; //No more decay until the command is issued again.
	
	//Pick a quark to change flavour:
	CQuark *pQuark = NULL;
	if(m_bProton) {
		for(int i=0;i<m_vhQuarks.Count();i++) {
			pQuark = m_vhQuarks[i];
			if(pQuark->m_flIsospin > 0) { //Up quark
				break; //Leave for loop - we have an appropriate quark
			}
		}
	} else {
		for(int i=0;i<m_vhQuarks.Count();i++) {
			pQuark = m_vhQuarks[i];
			if(pQuark->m_flIsospin < 0) { //Down quark
				break; //Leave for loop - we have an appropriate quark
			}
		}
	}
	
	m_bProton = !m_bProton;
	pQuark->m_flIsospin *= -1;
	m_flAppropriateMass = m_bProton ? PROTON_MASS : NEUTRON_MASS;
	m_flCharge = m_bProton * ELEMENTARY_CHARGE;
	
	//Physical part
	if( VPhysicsGetObject() ) {
		VPhysicsGetObject()->SetMass(m_flAppropriateMass);
	}
	
	Vector vecBosonOrigin = pQuark->GetAbsOrigin();
	QAngle angBosonDirection = vec3_angle;
	if( vecBosonOrigin != g_vecNuclearOrigin ) {
		VectorAngles(vecBosonOrigin - g_vecNuclearOrigin, angBosonDirection);
	}
	CWeakBoson *pWeak = CWeakBoson::CreateWeak(flEnergy, m_bProton ? -1 : 1, vecBosonOrigin, angBosonDirection, this);
	
	m_bColourUpToDate = m_bScaleUpToDate = false;
}

void CNucleon::ParticleUpdateSprite() {
	BaseClass::ParticleUpdateSprite();

	if( m_hSprite && !m_bColourUpToDate )
	{
		if( m_bProton ) {
			m_hSprite->SetColor(255,128,128); //Pale red, so it's not mixed up with red colour-charge
		} else {
			m_hSprite->SetColor(128,128,255); //Pale blue, so it's not mixed up with blue colour-charge
		}
		m_bColourUpToDate = true;
	}
}

void CNucleon::ApplyParticlePhysics() {
	//Decide whether or not to do a pion interaction
	if(!m_bTargettedByMeson && g_vpNucleons.Count() > 1 && gpGlobals->curtime > m_flNextQuarkIntTime && random->RandomFloat() > 0.75f) {
		bool bDecided = false;
		int iTries = 0;
		while(!bDecided && iTries < 5) {
			int iNucleonIndex = random->RandomInt(1, g_vpNucleons.Count());
			CNucleon *pNucleon = g_vpNucleons[iNucleonIndex - 1];
			if(pNucleon != this && !pNucleon->m_bTargettedByMeson) {
				if(debug_nucleus.GetBool()) {Msg("Nucleon %i has picked a target for a pion interaction!\n", entindex());}
				CQuark *pFirstQuark;
				CQuark *pSecondQuark;
				//Randomly select which quark will *not* be taking part in the interaction
				switch(random->RandomInt(0, 2)) {
					case 0:
						pFirstQuark = m_vhQuarks[1];
						pSecondQuark = m_vhQuarks[2];
						break;
					case 1:
						pFirstQuark = m_vhQuarks[0];
						pSecondQuark = m_vhQuarks[2];
						break;
					case 2:
						pFirstQuark = m_vhQuarks[0];
						pSecondQuark = m_vhQuarks[1];
				}
				
				Vector vecQuarkDest = pFirstQuark->GetLocalOrigin(); //This comes in later when we're moving the new quark in to take the old one's place
				
				int iGluonColour1 = -pFirstQuark->m_iColour; //We'll let the first one be the anticolour
				int iGluonColour2 = pSecondQuark->m_iColour;
				
				CQuark *pAntiQuark = (CQuark *)CreateEntityByName("particle_quark");
				if(debug_nucleus.GetBool()){Msg("Nucleon %i created quark %i\n", entindex(), pAntiQuark->entindex());}
				pAntiQuark->SetOwnerEntity(this);
				pAntiQuark->SetParent(this);
				if(debug_nucleus.GetBool()){Msg("Quark %i's owner and parent set\n", pAntiQuark->entindex());}
				pAntiQuark->m_iColour = iGluonColour1;
				if(debug_nucleus.GetBool()){Msg("Quark %i has colour %i\n", pAntiQuark->entindex(), iGluonColour1);}
				pAntiQuark->SetLocalOrigin(3.0f*pFirstQuark->GetLocalOrigin());
				if(debug_nucleus.GetBool()){Msg("Quark %i's origin set\n", pAntiQuark->entindex());}
				pAntiQuark->m_flAppropriateMass = 17.83f; //Scaled rest mass of up/down quark
				if(debug_nucleus.GetBool()){Msg("Quark %i's mass set\n", pAntiQuark->entindex());}
				pAntiQuark->m_flIsospin = pFirstQuark->m_flIsospin;
				pAntiQuark->ApplyParticlePhysics(); //Make it an antiquark straight away.
				if(debug_nucleus.GetBool()){Msg("Quark %i's isospin set\n", pAntiQuark->entindex());}
				pAntiQuark->Spawn();
				if(debug_nucleus.GetBool()){Msg("Quark %i spawned!\n", pAntiQuark->entindex());}
				
				CQuark *pNewQuark = (CQuark *)CreateEntityByName("particle_quark");
				if(debug_nucleus.GetBool()){Msg("Nucleon %i created quark %i\n", entindex(), pNewQuark->entindex());}
				pNewQuark->SetOwnerEntity(this);
				pNewQuark->SetParent(this);
				if(debug_nucleus.GetBool()){Msg("Quark %i's owner and parent set\n", pNewQuark->entindex());}
				pNewQuark->m_iColour = iGluonColour2;
				if(debug_nucleus.GetBool()){Msg("Quark %i has colour %i\n", pNewQuark->entindex(), iGluonColour2);}
				pNewQuark->SetLocalOrigin(1.5f*(pFirstQuark->GetLocalOrigin() + pSecondQuark->GetLocalOrigin()));
				if(debug_nucleus.GetBool()){Msg("Quark %i's origin set\n", pNewQuark->entindex());}
				pNewQuark->m_flAppropriateMass = 17.83f; //Scaled rest mass of up/down quark
				if(debug_nucleus.GetBool()){Msg("Quark %i's mass set\n", pNewQuark->entindex());}
				pNewQuark->m_flIsospin = pFirstQuark->m_flIsospin;
				if(debug_nucleus.GetBool()){Msg("Quark %i's isospin set\n", pNewQuark->entindex());}
				pNewQuark->Spawn();
				if(debug_nucleus.GetBool()){Msg("Quark %i spawned!\n", pNewQuark->entindex());}
				
				if(hadron_mode.GetInt() >= 1) {
					CBeam *pGluon1 = CBeam::BeamCreate( "sprites/glow_test02.vmt", 3 );
					CBeam *pGluon2 = CBeam::BeamCreate( "sprites/glow_test02.vmt", 3 );
					pGluon1->SetEndWidth( pGluon1->GetWidth() );
					pGluon2->SetEndWidth( pGluon2->GetWidth() );
					switch (iGluonColour1)
					{
						case CHARGE_ANTIBLUE:
							pGluon1->SetColor(255,255,0);
							break;
						case CHARGE_ANTIGREEN:
							pGluon1->SetColor(255,0,255);
							break;
						case CHARGE_ANTIRED:
							pGluon1->SetColor(0,255,255);
							break;
						case CHARGE_RED:
							pGluon1->SetColor(255,0,0);
							break;
						case CHARGE_GREEN:
							pGluon1->SetColor(0,255,0);
							break;
						case CHARGE_BLUE:
							pGluon1->SetColor(0,0,255);
							break;
					}
					switch (iGluonColour2)
					{
						case CHARGE_ANTIBLUE:
							pGluon2->SetColor(255,255,0);
							break;
						case CHARGE_ANTIGREEN:
							pGluon2->SetColor(255,0,255);
							break;
						case CHARGE_ANTIRED:
							pGluon2->SetColor(0,255,255);
							break;
						case CHARGE_RED:
							pGluon2->SetColor(255,0,0);
							break;
						case CHARGE_GREEN:
							pGluon2->SetColor(0,255,0);
							break;
						case CHARGE_BLUE:
							pGluon2->SetColor(0,0,255);
							break;
					}
					pGluon1->EntsInit(pSecondQuark, pAntiQuark);
					pGluon2->EntsInit(pSecondQuark, pNewQuark);
					pGluon1->LiveForTime(1);
					pGluon2->LiveForTime(1);
				}
				
				pSecondQuark->m_iColour = pFirstQuark->m_iColour;
				int iQuarkIndex = m_vhQuarks.Find(pFirstQuark);
				m_vhQuarks[iQuarkIndex] = pNewQuark;
				
				CMeson::CreateMeson(pFirstQuark, pAntiQuark, pNucleon, GetAbsOrigin() + 2.0f*pFirstQuark->GetLocalOrigin());
				
				m_flNextQuarkIntTime = gpGlobals->curtime + 2.0f;
				
				if(hadron_mode.GetInt() >= 1) {
					//Move the new quark in to take the old quark's place:
					pNewQuark->SetMoveType(MOVETYPE_NOCLIP);
					pNewQuark->SetLocalVelocity((vecQuarkDest - pNewQuark->GetLocalOrigin()) / 2.0f); //It'll be there in two seconds' time
					m_iQuarkToStop = iQuarkIndex;
				} else {
					pNewQuark->SetLocalOrigin(vecQuarkDest); //Just put it there
				}
				
				if(debug_nucleus.GetBool()){Msg("Nucleon %i has decided to do a pion interaction!\n", entindex());}
				bDecided = true;
			} else {
				iTries++;
			}
		}
	}
	
	BaseClass::ApplyParticlePhysics(); //This handles quark interaction scheduling
	
	if(g_bBetaPermitted) {
		//Figure out if it's energetically favourable to beta-decay:
		int Z = 0;
		for(int i=0;i<g_vpNucleons.Count();i++) {
			if(g_vpNucleons[i]->IsProton()) {
				Z++;
			}
		}
		int A = g_vpNucleons.Count();
		int N = A - Z;
		
		//If I'm a proton, decay will subtract from Z; if I'm a neutron, decay will add to Z.
		double dEnergyDifference = m_bProton ? (TotalNuclearEnergy(Z) - TotalNuclearEnergy(Z-1)) : (TotalNuclearEnergy(Z) - TotalNuclearEnergy(Z+1));
		if( dEnergyDifference > 0 ) {
			//Beta-decay will decrease the nuclear energy
			if(m_bProton) {
				if(debug_nucleus.GetBool()){Msg("Beta+ decay of %i is energetically favourable! (%.30f)\n", entindex(), dEnergyDifference);}
			} else {
				if(debug_nucleus.GetBool()){Msg("Beta- decay of %i is energetically favourable! (%.30f)\n", entindex(), dEnergyDifference);}
			}
			
			BetaDecay(dEnergyDifference);
		} else {
			if(debug_nucleus.GetBool()){Msg("Beta decay of %i is NOT energetically favourable! (%.30f)\n", entindex(), dEnergyDifference);}
		}
	}
}
#endif
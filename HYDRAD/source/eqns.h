// ****
// *
// * Class definition of the time-dependent hydrodynamic
// * equations, inherited by the adaptive mesh class
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 08/26/2021
// *
// ****


#include "cell.h"

// The radiation time-scale may not be adequate to control stability in the case
// of a large rate of energy change. The radiation time-scale can be grossly
// overestimated. Enabling this directive prevents catastrophic cooling to 
// negative energies / temperatures by setting d/dt = 0 when necessary
#define ENFORCE_POSITIVE_ELECTRON_ENERGY


// **** EQUATIONS CLASS ****

// Define the equations class
class CEquations {

    private:

    // Piece-wise polynomial fit to the field-aligned gravitational acceleration
    PPIECEWISEFIT pGravity;

#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
    // Piece-wise polynomial fit to the field-aligned magnetic field strength
    PPIECEWISEFIT pMagneticField;
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

#ifdef USE_JB
    double Tc, Te_max, old_Tc = 0.0;
#endif // USE_JB
    double lower_radiation_temperature_boundary;

    void Initialise( void );
    void FreeAll( void );

    // Function for finding the gravitational acceleration from the polynomial fit
    double CalculateGravity( double x );

    // Function for finding the smallest time-scale
    void GetSmallestTimeScale( double *delta_t, int iFirstStep );

	// The maximum current refinement level
    int iMaxRL;

#ifdef USE_KINETIC_MODEL
    // Functions for the Spitzer-Harm part of the solution
    // Tabulated values from tables I and II (for Z = 1) in Spitzer & Harm, 1953, Phys. Rev., 89, 977
    double SH_Table[51][3];
    void Get_SH_Table( void );
#endif // USE_KINETIC_MODEL

    public:

    // User specifiable and loop parameters
    PARAMETERS Params;

    // Pointer to the heating model
    PHEAT pHeat;

    // Pointers to the radiation models
    PRADIATION pRadiation, pRadiation2;

#if defined(OPTICALLY_THICK_RADIATION) || defined(BEAM_HEATING)
#ifdef OPTICALLY_THICK_RADIATION
    // Pointers to the ions for which optically-thick radiative emission
    // will be calculated
    POPTICALLYTHICKION pHI, pMgII, pCaII;
#ifdef NLTE_CHROMOSPHERE
    // Pointer to the radiative rates for the 6-level hydrogen atom
    PRADIATIVERATES pRadiativeRates;
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
#ifdef BEAM_HEATING
	// Pointers to the particular radiation objects which contain the models for hydrogen and helium
	PRADIATION pRadiation_H, pRadiation_He;
#endif // BEAM_HEATING
    // Pointer to the centre of the row at the current time (approx. the loop apex)
    PCELL pCentreOfCurrentRow;
#endif // OPTICALLY_THICK_RADIATION || BEAM_HEATING

    // Pointer to the left-most cell at the previous time (the start of the previous row)
    PCELL pStartOfPreviousRow;

    // Pointer to the left-most and right-most cells at the current time (the start and end of the current row)
    PCELL pStartOfCurrentRow, pEndOfCurrentRow;

    // Pointer to the active cell
    PCELL pActiveCell;

    // Constructor
    CEquations( void );

    // Destructor
    ~CEquations( void );

#if defined (OPENMP) || defined(USE_KINETIC_MODEL)
    // Pointer to an indexed list of cells
    PCELL *ppCellList;
    void CreateIndexedCellList( void );
#endif // OPENMP || USE_KINETIC_MODEL

#ifdef USE_KINETIC_MODEL
    void CalculateKineticModel( int iFirstStep );
    void CalculateNonMaxDFN( void );
#endif // USE_KINETIC_MODEL

    // Function for calculating physical quantities
    void CalculatePhysicalQuantities( void );

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
    void CalculateInitialPhysicalQuantities( void );
    void InitialiseRadiativeRates( void );
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

    // Function for evaluating the terms of the equations
    void EvaluateTerms( double current_time, double *delta_t, int iFirstStep );

	// Functions declared public: because needed in CMesh to calculate upper boundary
	// conditions in ghost cells for open flux tubes
#ifdef USE_POLY_FIT_TO_MAGNETIC_FIELD
	// Function for finding the magnetic field from the polynomial fit
    double CalculateMagneticField( double x );
    // Function for finding the cross-section from the polynomial fit
    double CalculateCrossSection( double x );
#endif // USE_POLY_FIT_TO_MAGNETIC_FIELD

	// Functions to provide a 2nd order accurate numerical integration of
    // the system of equations
    void Half_Time_Step( PCELLPROPERTIES pNewCellProperties, double delta_t );
    void Full_Time_Step( PCELLPROPERTIES CellProperties, double delta_t );

	int GetMaxRL( void ) { return iMaxRL; }

};

// ****
// *
// * Element Class Definition for Radiative Emission Model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 08/26/2021
// *
// ****


#include "config.h"

#ifdef TIME_VARIABLE_ABUNDANCES
    #define CORONAL_ABUNDANCE_FACTOR 4.0   // The abundance factor (enhancement of low FIP elements) that is "coronal"
    #define PHOTOSPHERIC_ABUNDANCE_FACTOR 1.0  // The abundance factor for "photospheric"
    #define LOW_FIP_THRESHOLD 10.0  // The threshold in eV below which elements are considered "low FIP"
#endif // TIME_VARIABLE_ABUNDANCES


class CElement {

    private:

    // The atomic number of the element
    int Z;

    // The abundance of the element relative to hydrogen
    double fAbund;

    #ifdef TIME_VARIABLE_ABUNDANCES
    // The first ionization potential (FIP) of the element
    double fFIP;
    #endif // TIME_VARIABLE_ABUNDANCES

    // The number of ions and their spectroscopic numbers
    int NumIons, *pSpecNum;
	
    // The temperature and density values in log_10 form
    int NumTemp, NumDen;
    double *pTemp, *pDen;
		
    // Pointer to an array of pointers, each pointing to the emissivity
    // data for an individual ion held in a NumTemp * NumDen size array
    double **ppEmiss;
	
    // Pointers to two arrays of pointers.  One points to the total ionisation rate
    // and the other to the total recombination rate of each ion at a specified
    // temperature
    double **ppIonRate, **ppRecRate;
	
    // Pointer to an array of pointers, each pointing to the fractional
    // population of an individual ion at a specified temperature
    double **ppIonFrac;
	
    // Pointers to the factor phi( n, T ) for each ion and the total phi( n, T )
    // for the element
    double **ppPhi, *pTotalPhi;
	
    // Function to open and read the ranges data file
    void OpenRangesFile( char *szRangesFilename );

    // Function to open and read the abundances data file
    void OpenAbundanceFile( char *szAbundFilename );

    // Function to open and read the emissivity data file
    void OpenEmissivityFile( char *szEmissFilename );

    // Function to open and read the total ionisation and recombination rates file
    void OpenRatesFile( char *szRatesFilename );

    // Function to open and read the ionisation balance file
    void OpenIonFracFile( char *szIonFracFilename );
    
    #ifdef TIME_VARIABLE_ABUNDANCES
    // Function to open and read the first ionization potential file
    void OpenFIPFile( char *szFIPFilename );
    #endif // TIME_VARIABLE_ABUNDANCES

    // Functions to calculate the factor phi( n, T ), which is multiplied by n^2 to calculate the radiated energy
    void CalculatePhi( void );
    void CalculateTotalPhi( void );

    // Function to free all allocated memory
    void FreeAll( void );

    // Function to return the required emissivity values
    double GetIonEmissivity( int iIon, double flog_10T, double flog_10n );

    public:
	
    // Constructor
    CElement( int iZ, char *szRangesFilename, char *szAbundFilename, char *szEmissFilename, char *szRatesFilename, char *szIonFracFilename );
	
    // Destructor
    ~CElement( void );

    // Function to initialise the element object
    void Initialise( int iZ, char *szRangesFilename, char *szAbundFilename, char *szEmissFilename, char *szRatesFilename, char *szIonFracFilename );

    // Function to return the element abundance
    double GetAbundance( void );

    // Function to return the total ionisation and total recombination rate
    // for a particular ion at a specified temperature and density
    void GetRates( int iIon, double flog_10T, double *pfIonRate, double *pfRecRate );
    void GetRates( int iIon, double flog_10T, double flog_10n, double *pfIonRate, double *pfRecRate );

    // Function to return the fractional population of a particular ion at a
    // specified temperature and density in equilibrium
    double GetEquilIonFrac( int iIon, double flog_10T );
    double GetEquilIonFrac( int iIon, double flog_10T, double flog_10n );

    // Functions to calculate the emissivity in equilibrium (this number includes multiplication by the ion fraction)
    // Multiply by the number density squared to obtain the energy radiatied in units of erg cm^-3 s^-1
    double GetEmissivity( int iIon, double flog_10T, double flog_10n );
    double GetEmissivity( double flog_10T, double flog_10n );

    // Function to calculate the rate of change with respect to time of the fractional population of the ions and the characteristic time-scale
	void Getdnibydt( double flog_10T, double flog_10n, double *pni0, double *pni1, double *pni2, double *pni3, double *pni4, double *s, double *s_pos, double *pv, double delta_s, double *pdnibydt, double *pTimeScale );

    // Functions to calculate the emissivity away from equilibrium (this number includes multiplication by the ion fraction)
    // Multiply by the number density squared to obtain the energy radiatied in units of erg cm^-3 s^-1
    double GetEmissivity( int iIon, double flog_10T, double flog_10n, double ni );
    double GetEmissivity( double flog_10T, double flog_10n, double *pni );
    
    #ifdef TIME_VARIABLE_ABUNDANCES
    // Function to return the element FIP
    double GetFIP( void );
    #endif // TIME_VARIABLE_ABUNDANCES

};

typedef CElement* PELEMENT;
typedef CElement** PPELEMENT;

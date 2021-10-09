// ****
// *
// * Class definition of the adaptive mesh
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 10/08/2021
// *
// ****

#include <time.h>

#include "eqns.h"


// **** ADAPTIVE MESH CLASS ****

// Define the adaptive mesh class
class CAdaptiveMesh : private CEquations {

    private:

	// The mesh (model) time elapsed and the current time-step
    double mesh_time, mesh_delta_t;

    // The file number for the next set of output profiles
    int iFileNumber;

    void CreateInitialMesh( void );

    void ZeroCellProperties( PCELLPROPERTIES pCellProperties );

    void Adapt( void );
	void EnforceBoundaryConditions( void );

    void FreeCurrentRow( void );
    void FreePreviousRow( void );

    void Solve( void );
#ifdef ADAPT
    void Integrate( bool bAdapt );
#else // ADAPT
    void Integrate( void );
#endif // ADAPT
#ifdef OPENMP
	void ShowProgress( double *ptimer );
#else // OPENMP
	void ShowProgress( clock_t *ptimer );
#endif // OPENMP
    void WriteToFile( void );
#ifdef UPDATE_HYDRAD_CONFIG
	void UpdateHYDRADConfig( void );
#endif // UPDATE_HYDRAD_CONFIG

    public:

    // Constructor
    CAdaptiveMesh( void );
	
    // Destructor
    ~CAdaptiveMesh( void );
	
};

// Define a type for the adaptive mesh class
typedef CAdaptiveMesh* PMESH;

// ****
// *
// * A heating code to simulate different forms of heat deposition
// *
// * Class function bodies
// *
// * (c) Dr. Stephen J. Bradshaw
// *     
// * Date last modified: 02/21/2020
// *
// ****


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "heat.h"
#include "../../Resources/source/gammabeta.h"
#include "../../Resources/source/constants.h"
#include "../../Resources/source/file.h"
#include "../../Resources/source/fitpoly.h"

#ifdef ALFVEN_WAVE_HEATING
#include <algorithm>
#endif // ALFVEN_WAVE_HEATING

CHeat::CHeat( void )
{
Initialise();
}

CHeat::~CHeat( void )
{
FreeAll();
}

void CHeat::Initialise( void )
{
GetHeatingData();
#ifdef BEAM_HEATING
	GetBeamHeatingData();
#endif // BEAM_HEATING
#ifdef OPTICALLY_THICK_RADIATION
	GetVALHeatingData();
#endif // OPTICALLY_THICK_RADIATION

#ifdef ALFVEN_WAVE_HEATING
    GetWaveData();
#endif // ALFVEN_WAVE_HEATING
}

void CHeat::FreeAll( void )
{
#ifdef OPTICALLY_THICK_RADIATION
int i;

for( i=0; i<2; i++ )
    free( ppVALHeating[i] );
free( ppVALHeating );
#endif // OPTICALLY_THICK_RADIATION

#ifdef BEAM_HEATING
if( iBeamHeatingDP )
{
	free( pfBeamSpectralIndex );
	free( pfBeamCutOff );
	free( pfBeamEnergyFlux );
	free( pfBeamTime );
}
#ifdef OPTICALLY_THICK_RADIATION
	#ifdef NLTE_CHROMOSPHERE
		if( pfQbeam )
			delete[] pfQbeam;
	#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
#endif // BEAM_HEATING

if( !NumActivatedEvents )
    return;

// Location, scale-length and maximum energy
free( s0episodic );
free( sHepisodic );
free( E0episodic );

// Start and end times of rise phase
free( tsRepisodic );
free( teRepisodic );

// Start and end times of decay phase
free( tsDepisodic );
free( teDepisodic );

#ifdef ALFVEN_WAVE_HEATING
delete fStartTime;
delete fPulseDuration;
delete fStartingPosition;
delete fTotalPoyntingFlux;
delete fPulseFrequency;
delete fPulseWaveNumber;
delete fPulseSign;
delete fEndTime;
#endif // ALFVEN_WAVE_HEATING

}

void CHeat::GetHeatingData( void )
{
FILE *pConfigFile;
int i;

// Open and read the configuration file
pConfigFile = fopen( "Heating_Model/config/heating_model.cfg", "r" );

// Get the quiescent heating parameter values
ReadDouble( pConfigFile, &s0quiescent );
ReadDouble( pConfigFile, &sHquiescent );
ReadDouble( pConfigFile, &E0quiescent );

// Get the episodic heating event parameter values
fscanf( pConfigFile, "%i", &NumActivatedEvents );

if( !NumActivatedEvents )
{
    // Close the configuration file
    fclose( pConfigFile );
    return;
}

// Allocate sufficient memory to store the positioning and timing information for each event

// Location, scale-length and maximum energy
s0episodic = (double*)malloc( sizeof(double) * NumActivatedEvents );
sHepisodic = (double*)malloc( sizeof(double) * NumActivatedEvents );
E0episodic = (double*)malloc( sizeof(double) * NumActivatedEvents );

// Start and end times of rise phase
tsRepisodic = (double*)malloc( sizeof(double) * NumActivatedEvents );
teRepisodic = (double*)malloc( sizeof(double) * NumActivatedEvents );

// Start and end times of decay phase
tsDepisodic = (double*)malloc( sizeof(double) * NumActivatedEvents );
teDepisodic = (double*)malloc( sizeof(double) * NumActivatedEvents );

for( i=0; i<NumActivatedEvents; i++ )
{
    // Location, scale-length and maximum energy
    ReadDouble( pConfigFile, &(s0episodic[i]) );
    ReadDouble( pConfigFile, &(sHepisodic[i]) );
    ReadDouble( pConfigFile, &(E0episodic[i]) );

    // Start and end times of rise phase
    ReadDouble( pConfigFile, &(tsRepisodic[i]) );
    ReadDouble( pConfigFile, &(teRepisodic[i]) );

    // Start and end times of decay phase
    ReadDouble( pConfigFile, &(tsDepisodic[i]) );
    ReadDouble( pConfigFile, &(teDepisodic[i]) );
}

// Close the configuration file
fclose( pConfigFile );
}

#ifdef BEAM_HEATING
void CHeat::GetBeamHeatingData( void )
{
FILE *pConfigFile;
char szBuffer[256];
int i;

// Open and read the configuration file
pConfigFile = fopen( "Heating_Model/config/beam_heating_model.cfg", "r" );

// Get the number of tabulated values for the beam heating parameters
fscanf( pConfigFile, "%i", &iBeamHeatingDP );
if( !iBeamHeatingDP ) return;	// There is no beam heating in the current run

// Get the header information
fscanf( pConfigFile, "%s", szBuffer );
fscanf( pConfigFile, "%s", szBuffer );
fscanf( pConfigFile, "%s", szBuffer );
fscanf( pConfigFile, "%s", szBuffer );
fscanf( pConfigFile, "%s", szBuffer );
fscanf( pConfigFile, "%s", szBuffer );
fscanf( pConfigFile, "%s", szBuffer );

pfBeamTime = (double*)malloc( sizeof(double) * iBeamHeatingDP );
pfBeamEnergyFlux = (double*)malloc( sizeof(double) * iBeamHeatingDP );
pfBeamCutOff = (double*)malloc( sizeof(double) * iBeamHeatingDP );
pfBeamSpectralIndex = (double*)malloc( sizeof(double) * iBeamHeatingDP );

if( iBeamHeatingDP == 1 )
{
	// The beam heating parameters are time-independent
	// The first time value in the beam configuration file will be the beam duration in this case
	ReadDouble( pConfigFile, &(pfBeamTime[0]) );
	// Get the remaining beam parameter values
	ReadDouble( pConfigFile, &(pfBeamEnergyFlux[0]) );
	ReadDouble( pConfigFile, &(pfBeamCutOff[0]) );
	ReadDouble( pConfigFile, &(pfBeamSpectralIndex[0]) );
} else {
	// The beam heating parameters are time-dependent
	// Get the tabulated beam parameter values
	for( i=0; i<iBeamHeatingDP; i++ )
	{
		ReadDouble( pConfigFile, &(pfBeamTime[i]) );
		// Get the remaining beam parameter values
		ReadDouble( pConfigFile, &(pfBeamEnergyFlux[i]) );
		ReadDouble( pConfigFile, &(pfBeamCutOff[i]) );
		ReadDouble( pConfigFile, &(pfBeamSpectralIndex[i]) );
	}
}

fclose( pConfigFile );

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
	pfQbeam = NULL;
	fAverageElectronEnergy = 0.0;
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
}
#endif // BEAM_HEATING

#ifdef OPTICALLY_THICK_RADIATION
void CHeat::GetVALHeatingData( void )
{
FILE *pFile;
int i;

// Open and read the configuration file
pFile = fopen( "Radiation_Model/atomic_data/OpticallyThick/VAL_atmospheres/VAL.heat", "r" );

// Get the number of data points in the file
fscanf( pFile, "%i", &iVALHeatingDP );

// Allocate sufficient memory to hold the heating data
ppVALHeating = (double**)malloc( sizeof(double*) * 2 );
for( i=0; i<2; i++ )
    ppVALHeating[i] = (double*)malloc( sizeof(double) * iVALHeatingDP );

for( i=0; i<iVALHeatingDP; i++ )
{
    // Array index [0][i] contain the column mass density and [1][i] contain the volumetric heating rate
    ReadDouble( pFile, &(ppVALHeating[0][i]) );
    ReadDouble( pFile, &(ppVALHeating[1][i]) );
}

fclose( pFile );
}
#endif // OPTICALLY_THICK_RADIATION

double CHeat::CalculateQuiescentHeating( double s )
{
double fQuiescent = 0.0, term1;

// Left footpoint
if( E0quiescent )
{
    term1 = ( s - s0quiescent );
    term1 *= term1;

    fQuiescent = E0quiescent * exp( - term1 / (2.0*sHquiescent*sHquiescent) );
}

return fQuiescent;
}

double CHeat::CalculateEpisodicHeating( double s, double t )
{
double fEpisodic = 0.0, term1, term2;
int i;

for( i=0; i<NumActivatedEvents; i++ )
{
    if( t >= tsRepisodic[i] && t <= teDepisodic[i] )
    {
        if( t >= tsRepisodic[i] && t <= teRepisodic[i] )
        {
            term1 = ( s - s0episodic[i] );
            term1 *= term1;

            term2 = ( t - tsRepisodic[i] ) / ( teRepisodic[i] - tsRepisodic[i] );

            fEpisodic += E0episodic[i] * exp( - term1 / (2.0*sHepisodic[i]*sHepisodic[i]) ) * term2;
        }
        else if( t > teRepisodic[i] && t < tsDepisodic[i] )
        {
            term1 = ( s - s0episodic[i] );
        	term1 *= term1;

            fEpisodic += E0episodic[i] * exp( - term1 / (2.0*sHepisodic[i]*sHepisodic[i]) );
        }
        else if( t >= tsDepisodic[i] && t <= teDepisodic[i] )
        {
            term1 = ( s - s0episodic[i] );
            term1 *= term1;

            term2 = 1.0 - ( ( t - tsDepisodic[i] ) / ( teDepisodic[i] - tsDepisodic[i] ) );

            fEpisodic += E0episodic[i] * exp( - term1 / (2.0*sHepisodic[i]*sHepisodic[i]) ) * term2;
        }
    }
}

return fEpisodic;
}

double CHeat::CalculateHeating( double s, double t )
{
double fQuiescentHeating = 0.0, fEpisodicHeating = 0.0;

if( E0quiescent )
    fQuiescentHeating = CalculateQuiescentHeating( s );

if( NumActivatedEvents )
    fEpisodicHeating = CalculateEpisodicHeating( s, t );

return( fQuiescentHeating + fEpisodicHeating );
}

#ifdef BEAM_HEATING
void CHeat::CalculateBeamParameters( double t, double *pBeamParams )
{
double x[3], y[3];
int i;

if( !iBeamHeatingDP ) return;
else if( iBeamHeatingDP == 1 ) {
	if( t > pfBeamTime[0] ) return;
	pBeamParams[0] = pfBeamEnergyFlux[0]; 
	pBeamParams[1] = pfBeamCutOff[0];
	pBeamParams[2] = pfBeamSpectralIndex[0];
} else {
	if( t < pfBeamTime[0] || t > pfBeamTime[iBeamHeatingDP-1] ) return;

	// Find the time interval for the interpolation of the parameter values
	for( i=0; i<iBeamHeatingDP-1; i++ )
		if( t < pfBeamTime[i+1] ) break;
	x[1] = pfBeamTime[i];
	x[2] = pfBeamTime[i+1];

	// Energy flux
	y[1] = pfBeamEnergyFlux[i];
	y[2] = pfBeamEnergyFlux[i+1];
	LinearFit( x, y, t, &(pBeamParams[0]) );

	// Cut-off
	y[1] = pfBeamCutOff[i];
	y[2] = pfBeamCutOff[i+1];
	LinearFit( x, y, t, &(pBeamParams[1]) );

	// Spectral index
	y[1] = pfBeamSpectralIndex[i];
	y[2] = pfBeamSpectralIndex[i+1];
	LinearFit( x, y, t, &(pBeamParams[2]) );
}

return;
}

double CHeat::CalculateBeamHeating( double t, double *pBeamParams, double nds, double Nstar, double n_e, double n_H, double x )
{
double energy_flux, cutoff_energy, delta;
double mu0 = 1.0;	// The cosine of the pitch angle is hard-wired for the moment
double fBeamHeating;
double term1, term2;

	// Trap the special cases
	if( ( !iBeamHeatingDP ) || ( iBeamHeatingDP == 1 && t > pfBeamTime[0] ) || ( iBeamHeatingDP > 1 && ( t < pfBeamTime[0] || t > pfBeamTime[iBeamHeatingDP-1] ) ) ) return 0.0;

	// The total energy flux in the beam at the injection site (erg cm^-2 s^-1)	
	energy_flux = pBeamParams[0];

	// Cut-off energy of the beam (erg)
	// 1.602e-9 erg / keV
	cutoff_energy = (1.602e-9) * pBeamParams[1];

	// The power-law index of the beam's energy distribution
	// delta must strictly be greater than 2
	delta = pBeamParams[2];

	//Average electron energy -> Derived from weighted average
	double avg_energy = ( ( 1.0 - delta ) / ( 2.0 - delta ) ) * cutoff_energy;

	// The three Coulomb logarithms used:
	double Lambda1 = 66.0 + ( 1.5 * log(avg_energy) ) - ( 0.5 * log(n_e) );
	double Lambda2 = 25.1 + log(avg_energy);

	// The special values of gamma and beta as defined in H&F 1994 or Emslie 1978 naming them g and b to avoid confusion with the mathematical Gamma and Beta functions
	double g = ( x * Lambda1 ) + ( ( 1.0 - x ) * Lambda2 );
	// double b = 2.0; // As prescribed by H&F
	
	// Nc and Ncstar, as defined in H&F 1994
	// double Nc = ( mu0 * cutoff_energy * cutoff_energy ) / ( g * (2.0 + (b/2.0)) * 2.0 * _PI_ * pow( ELECTRON_CHARGE, 4.0 ) );
	// Simplified to:
	double Nc = ( mu0 * cutoff_energy * cutoff_energy ) / ( g * (1.0006128424679109518248146922857e-36) );

	// double Ncstar = ( mu0 * cutoff_energy * cutoff_energy ) / ( Lambda1 * (2.0 + (b/2.0)) * 2.0 * _PI_ * pow( ELECTRON_CHARGE, 4.0 ) );
	// Simplified to:	
	double Ncstar = Nc * ( g / Lambda1 );

	// The lower bound energy of the heating integral
	// double integral_lowenergy = sqrt( (2.0+b/2.0) * g * 2.0 * _PI_ * pow(ELECTRON_CHARGE,4.0) * nds / mu0 );
	// Simplified to:
	double integral_lowenergy = sqrt( ( g * nds * (1.0006128424679109518248146922857e-36) ) / mu0 );
	if( isnan( integral_lowenergy ) ) return 0.0;
	
	term1 = (3.335376141559703172749382307619e-37) * n_H * g * ( delta - 2.0 ) * energy_flux * pow( (Nstar/Ncstar), (-delta/2.0) ) * Beta( (delta/2.0), 0.33333333 );
	term2 = 2.0 * cutoff_energy * cutoff_energy;

	if ( cutoff_energy > integral_lowenergy )
		term1 *= incompleteBeta( (delta/2.0), 0.33333333, (nds/Nc) );

	fBeamHeating = term1 / term2;
	if( isnan(fBeamHeating) || fBeamHeating < 0.0 ) fBeamHeating = 0.0;

	return( fBeamHeating );
}

#ifdef OPTICALLY_THICK_RADIATION
#ifdef NLTE_CHROMOSPHERE
void CHeat::InitQbeam( double fAvgEE, int iNumCells )
{
	// Set the average electron energy
	fAverageElectronEnergy = fAvgEE;

	// If memory has already been allocated to a beam at a previous time-step then free that memory ready for the beam at the current time-step
	if( pfQbeam )
		delete[] pfQbeam;

	// Reset the index to the beam quantities
	iQbeamIndex = 0;
	 
	iQbeamIndex_max = iNumCells * 2;
	pfQbeam = new double[iQbeamIndex_max];
}

void CHeat::SetQbeam( double s, double Qbeam )
{
	if( !pfQbeam ) return;

	pfQbeam[iQbeamIndex] = s;
	pfQbeam[iQbeamIndex+1] = Qbeam;

	iQbeamIndex += 2;
	if( iQbeamIndex == iQbeamIndex_max )
		iQbeamIndex = 0;
}

double CHeat::GetQbeam( double s )
{
	if( !pfQbeam ) return 0.0;

	double x[3], y[3], Qbeam;
	int i;

#ifdef OPENMP
	for( i=0; i<=iQbeamIndex_max-2; i+=2 )
#else // OPENMP
	for( i=iQbeamIndex; i<=iQbeamIndex_max-2; i+=2 )
#endif // OPENMP
		if( s <= pfQbeam[i] ) break;

	if( !i ) i+=2;
	else if( i > iQbeamIndex_max-2) i = iQbeamIndex_max-2;

	x[1] = pfQbeam[i-2];
	x[2] = pfQbeam[i];
	y[1] = pfQbeam[i-1];
	y[2] = pfQbeam[i+1];
	LinearFit( x, y, s, &Qbeam );
	if( Qbeam < 0.0 ) Qbeam = 0.0;

#ifdef OPENMP
#else // OPENMP
	iQbeamIndex = i;
#endif // OPENMP

	return Qbeam;
}

double CHeat::GetAvgEE( void )
{
	return fAverageElectronEnergy;	
}

void CHeat::WriteQbeam( void )
{
	FILE *pFile;
	int i;
	
	pFile = fopen( "Heating_Model/config/Qbeam.dat", "w" );
		fprintf( pFile, "%.16e\n", fAverageElectronEnergy );
		fprintf( pFile, "%i\n", iQbeamIndex_max>>1 );
		for( i=0; i<iQbeamIndex_max; i+=2 )
			fprintf( pFile, "%.16e\t%.16e\n", pfQbeam[i], pfQbeam[i+1] );
	fclose( pFile );
}

void CHeat::ReadQbeam( void )
{
	FILE *pFile;
	double fAvgEE, fs, fQbeam;
	int iNumCells, i;
	
	pFile = fopen( "Heating_Model/config/Qbeam.dat", "r" );
		ReadDouble( pFile, &fAvgEE );
		fscanf( pFile, "%i", &iNumCells );
		InitQbeam( fAvgEE, iNumCells );
		for( i=0; i<iNumCells; i++ )
		{
			ReadDouble( pFile, &fs );
			ReadDouble( pFile, &fQbeam );
			SetQbeam( fs, fQbeam );
		}
	fclose( pFile );
}
#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

#endif // BEAM_HEATING

#ifdef OPTICALLY_THICK_RADIATION
double CHeat::CalculateVALHeating( double flog10_rho_c )
{
int i;
double x[3], y[3], fVALHeating;

if( !iVALHeatingDP ) return 0.0;

// Trap limits
if( flog10_rho_c < ppVALHeating[0][0] ) return pow( 10.0, ppVALHeating[1][0] );
else if( flog10_rho_c > ppVALHeating[0][iVALHeatingDP-1] ) return pow( 10.0, ppVALHeating[1][iVALHeatingDP-1] );

for( i=0; i<iVALHeatingDP; i++ )
{
    if( flog10_rho_c <= ppVALHeating[0][i] ) break;
}

// Trap the special cases
if( i == 0 ) i = 1;
else if( i == iVALHeatingDP ) i = iVALHeatingDP - 1;

x[1] = ppVALHeating[0][i-1];
x[2] = ppVALHeating[0][i];
y[1] = ppVALHeating[1][i-1];
y[2] = ppVALHeating[1][i];

LinearFit( x, y, flog10_rho_c, &fVALHeating );

return pow( 10.0, fVALHeating );
}
#endif // OPTICALLY_THICK_RADIATION

#ifdef ALFVEN_WAVE_HEATING

// Reads in a file with the wave pulse data and sets up the pulses
void CHeat::GetWaveData( void )
{
    FILE *pConfigFile;
    int i, j, minimum_index = 0;
    double minimum_time = 1.0e100, temp_time;
    std::vector<double> row;
    
    // Open and read the configuration file
    pConfigFile = fopen( "Heating_Model/config/wave_pulses.cfg", "r" );
    
    // Read in the number of pulses
    if( fscanf( pConfigFile, "%i", &iNumPulses) != 1 )
    {
        fprintf( stderr, "Unable to open wave pulses data... aborting \n");
        exit(1);
    }
    
    // If zero, close the file and return nothing -- no wave pulses!
    if( !iNumPulses)
    {
        fclose( pConfigFile );
        return;
    }
    
    // Allocate sufficient memory to store the information for each pulse
    fStartTime = new double[iNumPulses];
    fPulseDuration = new double[iNumPulses];
    fStartingPosition = new double[iNumPulses];
    fTotalPoyntingFlux = new double[iNumPulses];
    fPulseFrequency = new double[iNumPulses];
    fPulseWaveNumber = new double[iNumPulses];
    fPulseSign = new double[iNumPulses];
    fEndTime = new double[iNumPulses];
    
    for( i=0; i<iNumPulses; i++ )
    {
        // Duration, injection location, total Poynting flux, and direction
        ReadDouble( pConfigFile, &(fStartTime[i]) );
        ReadDouble( pConfigFile, &(fPulseDuration[i]) );
        ReadDouble( pConfigFile, &(fStartingPosition[i]) );
        ReadDouble( pConfigFile, &(fTotalPoyntingFlux[i]) );
        ReadDouble( pConfigFile, &(fPulseFrequency[i]) );
        ReadDouble( pConfigFile, &(fPulseWaveNumber[i]) );
        ReadDouble( pConfigFile, &(fPulseSign[i]) );
        
        fEndTime[i] = fStartTime[i] + fPulseDuration[i];
        
        // Find the start time and index of the first pulse
        if( fStartTime[i] < minimum_time )
        {
            minimum_time = fStartTime[i];
            minimum_index = i;
        }
        
        fRayPositions.push_back(row);
        fRayPoyntingFlux.push_back(row);
        fRayFrequency.push_back(row);
        fRayWaveNumber.push_back(row);
        fRaySign.push_back(row);
    }
    
    // Close the configuration file
    fclose( pConfigFile );
    
    iRayFileNumber = 1;
    
    // Save the start time of the first ray
    NewRayTime = minimum_time;
    NewPulseIndex = minimum_index;
    
    // Calculate the start time of the second ray
    NextRayTime = NewRayTime + 1.0 / fRaysPerSecond;
    NextPulseIndex = minimum_index;
    
    // Calculate the starting times of all rays for all pulses -- saves a headache later on
    for( i=0; i<iNumPulses; i++ )
    {
        if( i != minimum_index )
        {
            if (fStartTime[i] < NextRayTime)
                NextPulseIndex = i;
            NextRayTime = std::min( NextRayTime, fStartTime[i] );
        }
        
        temp_time = fStartTime[i];
        
        fRayStartTimes.push_back( fStartTime[i]);
        iPulseIndices.push_back( i );
        for( j=1; j < (int) (fRaysPerSecond * fPulseDuration[i]); j++ )
        {
            // Store the pulse index of each ray
            iPulseIndices.push_back( i );
            
            temp_time += 1.0/fRaysPerSecond;
            // Store the time when the ray is initialized:
            fRayStartTimes.push_back( temp_time );
            
        }
        
        boolPulseInitialized.push_back( 0 );
    }
    
}

// Calculate the Coulomb logarithm
double CHeat::CalculateCoulombLog( double n_e, double n_H, double T_e, double T_i, double x )
{
    // 5.446170207e-4 = m_e / m_p
    if( T_i * 5.446170207e-4 < T_e && T_e < 116049.0 )
        return 23. - log( sqrt(n_e) * pow(T_e,-1.5) );
    
    else if( T_i * 5.446170207e-4 < 116049. && T_e > 116049.0 )
        return 24. - log( sqrt(n_e) / T_e);
    
    else if( T_e < T_i * 5.446170207e-4 )
        return 30. - log( sqrt(n_H * x) * pow(T_i,-1.5) );
    
    else
    {   // Maybe add an actual error that crashes the code.
        // This should only occur with unphysical temperatures, that is, never.
        printf("Error in Coulomb logarithm routine -> Coulomb logarithm outside of normal bounds\n");
        return 20.;
    }
    
}

// Calculate the magnetic field strength at a position s
double CHeat::CalculateMagneticField( double s, double fLoopLength )
{
    // A simple 5th order polynomial fit to the pressure gradient in the chromosphere, and
    // assume constant B in the corona.
    //
    // Assumes ~ 1000 G at the photosphere, ~ 100 G in the corona, and a power-law
    // scaling in between taken from Russell & Fletcher 2013.
    // This is completely arbitrary!
    
    // This should be updated to work properly with HYDRAD's field scaling, using a given
    // photospheric field strength.
    if( s < 3.e8 )
        return 1031.61 - 1.33014e-5 * s + 8.48378e-14 * s * s - 3.03276e-22 * s * s * s + 5.88384e-31 * s * s * s * s - 4.7167e-40 * s * s * s * s * s;
    else if( s > fLoopLength - 3.e8)
    {
        double temp_s = fLoopLength - s;
        return 1031.61 - 1.33014e-5 * temp_s + 8.48378e-14 * temp_s * temp_s - 3.03276e-22 * temp_s * temp_s * temp_s + 5.88384e-31 * temp_s * temp_s * temp_s * temp_s - 4.7167e-40 * temp_s * temp_s * temp_s * temp_s * temp_s;
    }
    else
        return 107.9;
}

// Calculate the loss in Poynting flux of a ray due to damping
void CHeat::UpdateRayPoyntingFlux( int pulse_index, int ray_index, double poynting_flux, double dz, double inverseL_DT )
{
    fRayPoyntingFlux[pulse_index][ray_index] = poynting_flux * exp( -1. * dz * inverseL_DT );
}

// Return the current number of pulses
int CHeat::GetNumberOfPulses( void )
{
    return iNumPulses;
}

// Initialize a new ray
void CHeat::InitializeRay( int pulse_index )
{
    // If this is the first ray in a pulse to be initialized:
    if( !(boolPulseInitialized[pulse_index]) )
    {
        boolPulseInitialized[pulse_index] = 1;
        iNumRaysInPulse.push_back( 1 );
    }
    else
        iNumRaysInPulse[pulse_index]++;
    
    // Initialize the parameters of the ray.  The frequency, wavenumber, and sign
    // could probably be stored in the pulse data, rather than individually for each ray,
    // but this allows it to be easily generalized if the need arises.
    // For example, reflection could be done by initializing a new pulse, or
    // simply changing the ray's parameters.
    fRayPositions[pulse_index].push_back( fStartingPosition[pulse_index]+1.e2 );
    fRayPoyntingFlux[pulse_index].push_back( fTotalPoyntingFlux[pulse_index]  );
    fRayFrequency[pulse_index].push_back( fPulseFrequency[pulse_index] );
    fRayWaveNumber[pulse_index].push_back( fPulseWaveNumber[pulse_index] );
    fRaySign[pulse_index].push_back( fPulseSign[pulse_index] );
    
    nCurrentRays++;
}

void CHeat::DeleteRay( int pulse_index, int ray_index )
{
    fRayPositions[pulse_index].erase( fRayPositions[pulse_index].begin() + ray_index );
    fRayPoyntingFlux[pulse_index].erase( fRayPoyntingFlux[pulse_index].begin() + ray_index );
    fRayFrequency[pulse_index].erase( fRayFrequency[pulse_index].begin() + ray_index );
    fRayWaveNumber[pulse_index].erase( fRayWaveNumber[pulse_index].begin() + ray_index );
    fRaySign[pulse_index].erase( fRaySign[pulse_index].begin() + ray_index );
    
    nCurrentRays--;
    iNumRaysInPulse[pulse_index]--;
    
    // If all the rays have been deleted, in a pulse, delete the pulse
    if (iNumRaysInPulse[pulse_index] == 0)
        boolPulseInitialized[pulse_index] = 0;
}

// Return the current position of ray with index ray_index
double CHeat::GetRayPosition( int pulse_index, int ray_index )
{
    return (fRayPositions[pulse_index][ray_index]);
}

// Increase the ray position by distance
void CHeat::SetRayPosition( int pulse_index, int ray_index, double distance )
{
    fRayPositions[pulse_index][ray_index] += distance;
}

// Return the frequency of ray with index ray_index
double CHeat::GetRayFrequency( int pulse_index, int ray_index )
{
    return (fRayFrequency[pulse_index][ray_index]);
}

// Return the wave number of ray with index ray_index
double CHeat::GetRayWaveNumber( int pulse_index, int ray_index )
{
    return (fRayWaveNumber[pulse_index][ray_index]);
}

// Return the Poynting flux of ray with index ray_index
double CHeat::GetRayPoyntingFlux( int pulse_index, int ray_index )
{
    return (fRayPoyntingFlux[pulse_index][ray_index]);
}

// Return the direction of ray with index ray_index
double CHeat::GetRaySign( int pulse_index, int ray_index )
{
    return (fRaySign[pulse_index][ray_index]);
}

// Gets the next time for a ray to be initialized
double CHeat::GetNewRayTime( void )
{
    return NewRayTime;
}

// Sets the next time for rays to be initialized
//  and initializes the ray(s) at the current time
void CHeat::SetNewRayTime( double current_time )
{
    std::vector<double> RayIndices;
    double max_time, temp_time;
    
    long unsigned int j, counter;  // since we compare to a vector size,
    // this throws a warning if just declared as "int"
    
    max_time = *std::max_element(fRayStartTimes.begin(), fRayStartTimes.end());
    if( current_time > max_time )
    {
        NewRayTime = 1.0e100;
        wave_tracing_delta_t = 1.0e100;
    }
    else
    {
        // NewRayTime is time of the very next ray to be initialized.
        // NextRayTime is the ray after that.  Better name for variables?
        NewRayTime = NextRayTime;
        NewPulseIndex = NextPulseIndex;
        
        counter = 0;
        for( j=0; j<iPulseIndices.size(); j++ )  // First, find the number of rays starting at this time
        {                                         // Probably just 1, but this is more general and safer.
            if( fRayStartTimes[j] == NextRayTime )
            {
                counter++;
                RayIndices.push_back( iPulseIndices[j] );
            }
        }
        
        for( j=0; j<counter; j++ )
        {
            InitializeRay( RayIndices[j] );  // I've left InitializeRay as its own function
            // so that it can be used in a more general sense when
            // the code is eventually modified to include reflection.
            // Its block of code could be copied+pasted here otherwise
            // which would save a bit of runtime.
            nCurrentRays++;
        }
        
        // Find the next time a ray starts:
        temp_time = 1.0e100;
        for( j=0; j<iPulseIndices.size(); j++ )
        {
            
            if( fRayStartTimes[j] > NewRayTime && fRayStartTimes[j] < temp_time )
            {
                temp_time = fRayStartTimes[j];
                NextPulseIndex = iPulseIndices[j];
            }
        }
        
        NextRayTime = temp_time;
    }
}

double CHeat::GetRayTimeScale( void )
{
    return wave_tracing_delta_t;
}

// Returns a boolean of whether a pulse has been initialized in the simulation yet
int CHeat::IsPulseInitialized( int pulse_index )
{
    return boolPulseInitialized[pulse_index];
}

// Returns the number of rays that have been initialized in a given pulse
int CHeat::GetNumRaysInPulse( int pulse_index )
{
    return iNumRaysInPulse[pulse_index];
}

double CHeat::GetNextOutputTime( void )
{
    return fNextOutputTime;
}

void CHeat::SetNextOutputTime( double fTime )
{
    fNextOutputTime = fTime;
}

void CHeat::IncrementRayFileNumber( void )
{
    iRayFileNumber++;
}

int CHeat::GetRayFileNumber( void )
{
    return iRayFileNumber;
}
#endif // ALFVEN_WAVE_HEATING

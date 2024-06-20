// ****
// *
// * A heating code to simulate different forms of heat deposition
// *
// * Class function definitions
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Date last modified: 02/21/2020
// *
// ****

#include "config.h"

#ifdef ALFVEN_WAVE_HEATING
#include <vector>
#endif // ALFVEN_WAVE_HEATING

#ifdef OPTICALLY_THICK_RADIATION
	#ifdef NLTE_CHROMOSPHERE
		#ifdef BEAM_HEATING
			#define Qbeam_CUT_OFF		1E-10
		#endif // BEAM_HEATING
	#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION

class CHeat {

    private:

    // Quiescent heating parameters
	// Location, scale-length and maximum energy
    	double s0quiescent, sHquiescent, E0quiescent;
	
    // Episodic heating parameters
	// The number of episodic heating events to be activated
    	int NumActivatedEvents;

        // Location, scale-length and maximum energy
        double *s0episodic, *sHepisodic, *E0episodic;

        // Start and end times of rise phase
        double *tsRepisodic, *teRepisodic;
        // Start and end times of decay phase
        double *tsDepisodic, *teDepisodic;

#ifdef BEAM_HEATING
    // Beam heating parameters
        int iBeamHeatingDP;
        double *pfBeamTime, *pfBeamEnergyFlux, *pfBeamCutOff, *pfBeamSpectralIndex;
#ifdef OPTICALLY_THICK_RADIATION
	#ifdef NLTE_CHROMOSPHERE
		int iQbeamIndex, iQbeamIndex_max;
		double *pfQbeam, fAverageElectronEnergy;
	#endif // NLTE_CHROMOSPHERE
#endif // OPTICALLY_THICK_RADIATION
#endif // BEAM_HEATING

#ifdef OPTICALLY_THICK_RADIATION
    // Background volumetric heating required to maintain the VAL atmosphere in equilibrium
        int iVALHeatingDP;
        double **ppVALHeating;
#endif // OPTICALLY_THICK_RADIATION

    void Initialise( void );
    void FreeAll( void );

    // Time independent / dependent (thermal) heating functionality
    void GetHeatingData ( void );

#ifdef BEAM_HEATING
    // Beam heating functionality
    void GetBeamHeatingData( void );
#endif // BEAM_HEATING

#ifdef OPTICALLY_THICK_RADIATION
    // Heating to maintain the lower atmosphere
    void GetVALHeatingData( void );
#endif // OPTICALLY_THICK_RADIATION

// Wave heating parameters
#ifdef ALFVEN_WAVE_HEATING
    // Pulse properties, to be read from file:
    int iNumPulses;             // The number of pulses of waves
    double *fStartTime;         // The time at which the pulse begins (s)
    double *fPulseDuration;     // The duration of each pulse (s)
    double *fStartingPosition;  // The injection location of each pulse (cm)
    double *fTotalPoyntingFlux; // The total initial Poynting flux of each pulse (erg/s/cm^2)
    double *fPulseFrequency;    // The frequency of the pulse (Hz)
    
    double *fPulseWaveNumber;   // The perpendicular wave-number of the pulse at the PHOTOSPHERE (cm^-1)
    // NB: Reep et al lists the APEX wave-number.  Conversion is simple:
    // k_p = k_a * B_p / B_a
    // -- This scales in position with B-field.
    // I prefer this to be programmed in as a photospheric value because it's easier to specify
    // a photospheric B-strength since we don't have to locate the loop apex.
    // Functionally, they're equivalent.
    
    double *fPulseSign;         // The initial direction of the pulse
    double *fEndTime;           // The time at which the pulse ends (s)
    
    std::vector<double> fRayStartTimes;     // The starting times of all the rays, to be sorted in order
                                            // to keep track of when to start a new ray
    
    std::vector<int> iPulseIndices;      // A list of which ray belongs to which pulse
    std::vector<int> boolPulseInitialized; // Has a wave pulse been initialized yet?
    std::vector<int> iNumRaysInPulse;       // The number of rays initialized in the pulse
    
    // Ray properties:
    int nCurrentRays = 0;       // The number of rays currently being traced
    double NewRayTime;          // The next time to initialize a new ray
    double NextRayTime;         // The next time after NewRayTime to initialize a ray
    int NextPulseIndex, NewPulseIndex;
    
    std::vector<std::vector<double> > fRayPositions;
    std::vector<std::vector<double> > fRayPoyntingFlux;
    std::vector<std::vector<double> > fRayFrequency;
    std::vector<std::vector<double> > fRayWaveNumber;
    std::vector<std::vector<double> > fRaySign;
    
    // The number of rays per second of injection duration, pre-defined as 50
    const double fRaysPerSecond = 50.0;
    
    // The time-scale necessary to resolve the initialization of rays
    double wave_tracing_delta_t = 0.5 / fRaysPerSecond;
    
    double fNextOutputTime;
    int iRayFileNumber;
    
    void GetWaveData ( void );
    
#endif // ALFVEN_WAVE_HEATING

    public:

    CHeat( void );
    ~CHeat( void );

    // Time independent / dependent (thermal) heating functionality
    double CalculateQuiescentHeating( double s );
    double CalculateEpisodicHeating( double s, double t );
    double CalculateHeating( double s, double t );

#ifdef BEAM_HEATING
    // Beam heating functionality
    void CalculateBeamParameters( double t, double *pBeamParams );
    double CalculateBeamHeating( double t, double *pBeamParams, double nds, double Nstar, double n_e, double n_H, double x );
		#ifdef OPTICALLY_THICK_RADIATION
			#ifdef NLTE_CHROMOSPHERE
				// Functions that will provide access to the beam energy input (at the previous time-step) for the non-thermal collision frequency
				void InitQbeam( double fAvgEE, int iNumCells );
				void SetQbeam( double s, double Qbeam );
				double GetQbeam( double s );
				double GetAvgEE( void );
				void WriteQbeam( void );
				void ReadQbeam( void );
			#endif // NLTE_CHROMOSPHERE
		#endif // OPTICALLY_THICK_RADIATION
#endif // BEAM_HEATING

#ifdef OPTICALLY_THICK_RADIATION
    // Heating to maintain the lower atmosphere
    double CalculateVALHeating( double flog10_rho_c );
#endif // OPTICALLY_THICK_RADIATION

#ifdef ALFVEN_WAVE_HEATING
    // Physical calculations:
    double CalculateCoulombLog( double n_e, double n_H, double T_e, double T_i, double x );
    double CalculateMagneticField( double s, double fLoopLength );
    void UpdateRayPoyntingFlux( int pulse_index, int num, double poynting_flux, double dz, double inverseL_DT );
    
    // Pulse and ray manipulations
    int GetNumberOfPulses( void );
    void InitializeRay( int pulse_index );
    void DeleteRay( int pulse_index, int ray_index );
    double GetRayPosition( int pulse_index, int ray_index );
    void SetRayPosition( int pulse_index, int ray_index, double distance );
    double GetRayFrequency( int pulse_index, int ray_index );
    double GetRayWaveNumber( int pulse_index, int ray_index );
    double GetRayPoyntingFlux( int pulse_index, int ray_index );
    double GetRaySign( int pulse_index, int ray_index );
    double GetNewRayTime( void );
    void SetNewRayTime( double current_time );
    double GetRayTimeScale( void );
    int IsPulseInitialized( int pulse_index );
    int GetNumRaysInPulse( int pulse_index );
    
    // ----  These four can be removed when there is a proper file printing function
    // ----  added to mesh.cpp.  Used to write file with ray values.
    double GetNextOutputTime( void );
    void SetNextOutputTime( double fTime );
    void IncrementRayFileNumber( void );
    int GetRayFileNumber( void );
    // -----
    
#endif //ALFVEN_WAVE_HEATING

};

typedef CHeat* PHEAT;

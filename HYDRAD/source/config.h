// ****
// *
// * #defines for configuring the hydrodynamic model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by HYDRAD_GUI on 12-09-2018 14:23:26
// *
// ****


// **** Output ****
#define OUTPUT_EVERY_N_TIME_STEPS 1000
#define WRITE_FILE_PHYSICAL
// **** End of Output ****

// **** Physics ****
// Heating //
#include "../../Heating_Model/source/config.h"
// End of Heating //
// Radiation //
#include "../../Radiation_Model/source/config.h"
// End of Radiation //
// Thermal Conduction //
#define HEAT_FLUX_LIMITING_COEFFICIENT 0.167
#define TIME_STEP_LIMIT 1E-10
// End of Thermal Conduction //
// Collisions //
#include "collisions.h"
// End of Collisions //
// Flux Tube //
// End of Flux Tube //
// **** End of Physics ****

// **** Solver ****
#define SAFETY_RADIATION 0.1
#define SAFETY_CONDUCTION 1.0
#define SAFETY_ADVECTION 1.0
#define SAFETY_VISCOSITY 1.0
#define TIME_STEP_INCREASE_LIMIT 1.05
#define MINIMUM_RADIATION_TEMPERATURE 2E4
#define ZERO_OVER_TEMPERATURE_INTERVAL 5E2
#define MINIMUM_TEMPERATURE 1E4
// **** End of Solver ****

// **** Grid ****
#define ADAPT
#define MAX_REFINEMENT_LEVEL 12
#define REFINE_ON_DENSITY
#define REFINE_ON_ELECTRON_ENERGY
#define ADAPT_EVERY_N_TIME_STEPS 10
#define MIN_FRAC_DIFF 0.05
#define MAX_FRAC_DIFF 0.1
#define ENFORCE_CONSERVATION
// **** End of Grid ****

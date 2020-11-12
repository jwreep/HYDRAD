// ****
// *
// * #defines for configuring the radiation model
// *
// * (c) Dr. Stephen J. Bradshaw
// *
// * Source code generated by HYDRAD_GUI on 11-11-2020 22:03:25
// *
// ****

// **** Physics ****
// Radiation //
#define USE_POWER_LAW_RADIATIVE_LOSSES
// End of Radiation //
// Collisions //
#include "../../HYDRAD/source/collisions.h"
// End of Collisions //
// Flux Tube //
// End of Flux Tube //
// **** End of Physics ****

// **** Solver ****
#define MAX_OPTICALLY_THIN_DENSITY 1.0E12
#define SAFETY_ATOMIC 1.0
#define CUTOFF_ION_FRACTION 1E-15
#define EPSILON_D 0.1
#define EPSILON_R 1.8649415311920072
// **** End of Solver ****

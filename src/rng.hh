#ifndef RNG_HH
#define RNG_HH

/* #include <gsl/gsl_rng.h> */
/* #include <gsl/gsl_randist.h> */

#include <cstdlib>
#include <cmath>

inline void seed(long s) { srand(s); }
inline double frand() {
    return (double)rand() / RAND_MAX;
}
inline double uniform(const double min, const double max) { 
    return frand() * (max - min) + min; 
}
inline double normal(const double sigma) { 
    double u = frand();
    double v = frand();
    return sqrt(-2.0 * log(u)) * sin(2.0 * M_PI * v) * sigma;
}

#endif
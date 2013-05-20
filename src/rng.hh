#ifndef RNG_HH
#define RNG_HH

#include <cstdlib> /* rand */
#include <cmath>   /* M_PI */

inline void seed(long s) { srand(s); }

/*! Returns a floating-point number between 0 and 1
 * @return A single sample from U[0, 1) */
inline double frand() {
    return (double)rand() / RAND_MAX;
}

/*! Sample from a Uniform distribution between min and max
 * @param minimum the minimum value (inclusive)
 * @param maximum the maximum value (exclusive)
 * @return A single sample from U[min, max) */
inline double uniform(const double min, const double max) { 
    return frand() * (max - min) + min; 
}

/*! Sample from Normal distribution using Box-Muller transform 
 * @param sigma The standard deviation
 * @return A single sample from N(0, sigma) */
inline double normal(const double sigma) { 
    double u = frand();
    double v = frand();
    return sqrt(-2.0 * log(u)) * sin(2.0 * M_PI * v) * sigma;
}

#endif
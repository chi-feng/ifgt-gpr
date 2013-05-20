#ifndef IFGT_HH
#define IFGT_HH

#include <iostream>
#include <fstream> 
#include <cmath>
#include <climits>
#include <vector>
#include "linalg.hh"

class IFGT {

private:
    /*! dimension of each source */
    size_t dim;                             
    /*! source locations, each location being a Vector */
    std::vector<Vector> sources;
    /*! vector of source weights */            
    std::vector<double> weights;
    /*! vector of distance to nearest center for each source */            
    std::vector<double> distances;          
    /*! get the index of the source marked as a center by id */
    std::vector<int> centerToSource;        
    /*! number of centers, aka K */
    size_t numCenters;
    /*! get the index of the center corresponding to a source id */
    std::vector<int> sourceToCenter;    
    /*! contains C_alpha^k of all K centers */   
    std::vector<Vector> centerCoefficients;
    /*! degree of maximum taylor expansion */
    size_t degree;
    /*! number of expansion monomials when expanded in
        graded lexicographic ordering */          
    size_t numExpansionTerms;
    /*! desired precision */
    double epsilon;                   
    /*! desired bandwidth */
    double bandwidth;                 
    /*! desired cluster radius */
    double clusterRadius;
    /*! desired cutoff radius */
    double cutoffRadius;
    /*! current radius */ 
    double radius;
    /*! highest number of neighbors around any point */ 
    size_t maxNeighbors;
    /*! verbosity */
    int verbose;
    
    inline double norm(size_t x, size_t y);
    void computeClusters();
    void computeCoefficients(); 
    int farthestSource();
    
public:
    IFGT(const std::vector<Vector>& sources, const Vector& weights, const double bandwidth, const size_t degree, const double clusterRadius, const double cutoffRadius, const int verbose);
    static double errorBound(const size_t p, const double radius, const double cutoff);
    void evaluate(const std::vector<Vector>& targets, std::vector<double>& result);
    void directEvaluate(const std::vector<Vector>& targets, std::vector<double>& result);
    void writeClustersToFile(const char* filename);
    static inline size_t binomialCoefficient(size_t n, size_t k) {
        assert(k <= n);
        if (k > n - k) { // (n, k) = (n, n-k), so use smaller n-k
            k = n - k;
        }
        int c = 1;
        for (size_t i = 1; i < k + 1; i++) {
            c *= n - (k - i);
            c /= i;
        }
        return c;
    };
};

#endif
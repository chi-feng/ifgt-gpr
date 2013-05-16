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
    
    size_t dim;                       /*! dimension */
    std::vector<Vector> sources;           /*! source locations */
    std::vector<double> weights;           /*! source locations */
    std::vector<double> distances;         /*! distance to nearest center */
    std::vector<int> centerToSource;       /*! center -> source mapping */
    size_t numCenters;
    std::vector<int> sourceToCenter;       /*! source -> center mapping */
    std::vector<Vector> centerCoefficients;
    
    size_t degree;                    /*! degree of taylor expansion */
    size_t numExpansionTerms;
    double epsilon;                   /*! desired precision */
    double bandwidth;                 /*! desired bandwidth */
    double clusterRadius;
    double cutoffRadius;
    double radius;
    size_t maxNeighbors;
    
    inline double norm(size_t x, size_t y);
    void computeClusters();
    void computeCoefficients(); 
    int farthestSource();
    
public:
    IFGT(const std::vector<Vector>& sources, const Vector& weights, const double bandwidth, const size_t degree, const double clusterRadius, const double cutoffRadius);
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
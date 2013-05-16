
#include "IFGT.hh"
#include "rng.hh"
#include "execution_timer.h"

/*! Constructor for the Improved Fast Gauss Transform
 \param sources A vector of N sources (dim-dimensional vectors)
 \param weights A vector of weights q_i length N
 \param bandwidth The bandwidth of the Gauss Transform
 \param degree The maximum degree (p)
 \param clusterRadius The desired cluster radius in units of bandwidth  
        (centers will be added until radius falls below this value)
 \param cutoffRadius The desired cutoff radius in units of bandwidth */
IFGT::IFGT(
    const Matrix& sources, 
    const Vector& weights,
    const double bandwidth, 
    const size_t degree,
    const double clusterRadius, 
    const double cutoffRadius
    const int verbose = 0) {
    
    this->sources = sources;
    this->weights = weights;
    this->bandwidth = bandwidth;
    this->degree = degree;
    this->clusterRadius = clusterRadius;
    this->cutoffRadius = cutoffRadius;
    
    dim = sources[0].size();
    sourceToCenter.resize(sources.size(), 0);
    distances.resize(sources.size(), 0.0);
    
    computeClusters();
    numExpansionTerms = binomialCoefficient(degree - 1 + dim, dim);
    computeCoefficients();
}

/*! Compute the error bound of the IFGT
 * \param degree The degree of the Taylor expansion
 * \param radius Target cluster radius (a multiple of bandwidth)
 * \param cutoff Cutoff distance (a multiple of bandwidth)
 * \return Error bound of the IFGT */
double IFGT::errorBound(const size_t degree, const double radius, const double cutoff) {
    // compute the leading term 2^p/p! in the truncation error
    double leading = 2.0;
    for (size_t i = 2; i <= degree; i++) {
        leading *= 2.0 / i;
    }
    double truncationError = leading * pow(radius, degree) * pow(cutoff, degree);
    double cutoffError = exp(-cutoff * cutoff);
    // return Q(truncationError + cutoffError)
    return 0.5 * erfc((truncationError + cutoffError) / sqrt(2));
}

/*! Add centers until cluster radius falls below threshold */
void IFGT::computeClusters() {
    // initially use first source as a center
    centerToSource.push_back(0);
    for (size_t i = 0; i < sources.size(); ++i) {
        distances[i] = norm(i, 0);
        sourceToCenter[i] = 0;
    }
    // add the farthest source incrementally until radius < clusterRadius
    int farthest = farthestSource();
    radius = sqrt(distances[farthest]);
    while (radius > clusterRadius * bandwidth) {
        centerToSource.push_back(farthest);
        for (size_t i = 0; i < sources.size(); i++) {
            double distance = norm(i, centerToSource.back());
            if (distance < distances[i]) {
                distances[i] = distance;
                sourceToCenter[i] = centerToSource.size() - 1;
            }
        }
        farthest = farthestSource();
        radius = sqrt(distances[farthest]);
    }
    numCenters = centerToSource.size();
}

/*! Computes the squared norm between two sources
 * \param x Index of first source
 * \param y Index of second source
 * \return |x-y|^2 */
inline double IFGT::norm(const size_t x, const size_t y) {
  double delta, total = 0;
  for (size_t i = 0; i < dim; i++) {
      delta = sources[x][i] - sources[y][i];
      total += delta * delta;
  }
  return total;
}

/*! Find the source that is farthest from any center */
int IFGT::farthestSource() {
    int farthestSource = 0;
    double farthestDistance = distances[farthestSource];
    for (size_t i = 1; i < sources.size(); i++) {
        if (distances[i] > farthestDistance) {
            farthestDistance = distances[i];
            farthestSource = i;
        }
    }
    return farthestSource;
}

/*! Compute the coefficients in (3.14) */
void IFGT::computeCoefficients() { 
    centerCoefficients.resize(numCenters, Vector(numExpansionTerms, 0.0)); 
    
    Vector sourceCenterMonomials(numExpansionTerms);   
    Vector sourceCenterDifferences(dim);
     
    Vector constantTerms(numExpansionTerms);

    Vector heads(dim, 0.0);
    Vector factorialTerms(numExpansionTerms);
    
    // precompute constants 
    double bandwidthSquared = bandwidth * bandwidth; 
    
    for (size_t i = 0; i < sources.size(); i++) {
        size_t k = sourceToCenter[i];
        for (size_t d = 0; d < dim; d++) {
            sourceCenterDifferences[d] = (sources[i][d] - sources[k][d]) / bandwidth;
            heads[d] = 0.0;
        }
        // compute source-center monomials
        // loop over graded lexicographic order
        sourceCenterMonomials[0] = 1.0;
    	for (size_t p = 1, t = 1, tail = 1; p < degree; p++, tail = t) {
    		for (size_t d = 0; d < dim; d++){
    			size_t head = heads[d];
    			heads[d] = t;
    			for (size_t j = head; j < tail; j++, t++) {
    			    sourceCenterMonomials[t] = sourceCenterDifferences[d] * sourceCenterMonomials[j];
    			}
    		}
    	}
        
        double gaussianTerm = weights[i] * exp(-distances[i] / bandwidthSquared); 
        
        for (size_t alpha = 0; alpha < numExpansionTerms; alpha++) {
            centerCoefficients[k][alpha] += gaussianTerm * sourceCenterMonomials[alpha];
        }
    }
    
    // compute the 2^|alpha|/alpha! constant at the front of C_alpha^k
    // loop over graded lexicographic order
    heads[0] = 0;
    factorialTerms[0] = 0;
	constantTerms[0] = 1.0;
    for (size_t p = 1, t = 1, tail = 1; p < degree; p++, tail = t) {
        for (size_t d = 0; d < dim; d++) {
            size_t head = heads[d];
            heads[d] = t;
            for (size_t j = head; j < tail; j++, t++) {
                if (d + 1 >= dim) // reset factorial terms 
                    factorialTerms[t] = 1; 
                else // increment factorial terms
                    factorialTerms[t] = (j < heads[d+1]) ? factorialTerms[j] + 1 : 1;
                constantTerms[t] = 2.0 * constantTerms[j] / factorialTerms[t];
            }
        }
    }
    
    // multiply sum inside C_alpha^k by the constant term
    for (size_t k = 0; k < numCenters; k++) {
        for (size_t alpha = 0; alpha < numExpansionTerms; alpha++) {
            centerCoefficients[k][alpha] *= constantTerms[alpha];
        }
    }
}

/*! Evaluate the sum of Gaussians using (3.10) 
 * \param targets The points y_j that the Gauss transform will be evaluated at
 * \param gaussTransform The transformed values will be stored in gaussTransform */
void IFGT::evaluate(const std::vector<Vector>& targets, Vector& gaussTransform) {
    std::cout << "IFGT: Evaluating Gauss transform" << std::endl;
    size_t numTargets = targets.size();
    gaussTransform.resize(numTargets, 0);
    
    Vector targetCenterMonomials(numExpansionTerms);
    Vector targetCenterDifferences(dim);
    
    Vector heads(dim);
    
    // precompute constants 
    double cutoffRadiusSquared = cutoffRadius * cutoffRadius * bandwidth * bandwidth;
    double bandwidthSquared = bandwidth * bandwidth; 
    maxNeighbors = 0;
    size_t neighbors;
    // loop through each target y_j (3.11)
    for (size_t j = 0; j < numTargets; j++) {
        neighbors = 0;
        // find clusters centers c_k within cutoff radius 
        for (size_t k = 0; k < numCenters; k++) {
            int center = centerToSource[k];
            double targetCenterDistance = 0.0;
            for (size_t d = 0; d < dim; d++) {
                targetCenterDifferences[d] = targets[j][d] - sources[center][d];
                targetCenterDistance += targetCenterDifferences[d] * targetCenterDifferences[d];
                // exit if running sum of square distance is already over threshold
                if (targetCenterDistance > cutoffRadiusSquared)
                    break;
            }
            if (targetCenterDistance < cutoffRadiusSquared) {
                neighbors++;
                // compute target-center monomials
            	targetCenterMonomials[0] = 1.0;
            	for (size_t p = 1, t = 1, tail = 1; p < degree; p++, tail=t){
            		for (size_t d = 0; d < dim; d++){
            			size_t head = heads[d];
            			heads[d] = t;
            			for (size_t j = head; j < tail; j++, t++)
            				targetCenterMonomials[t] = targetCenterDifferences[d] * targetCenterMonomials[j];
            		}
            	}
                // putting it all together
                double gaussianTerm = exp(-targetCenterDistance / bandwidthSquared);
                for (size_t alpha = 0; alpha < numExpansionTerms; alpha++) {
                    gaussTransform[j] += centerCoefficients[k][alpha] * gaussianTerm * targetCenterMonomials[alpha];
                }
            }
        }
        if (neighbors > maxNeighbors) {
            maxNeighbors = neighbors;
        }
    }
    std::cout << "IFGT: Evaluating Gauss transform...Done (maxNeighbors = " << maxNeighbors << ")" << std::endl;
    long direct = sources.size() * numTargets * dim; 
    long ifgt = sources.size() * numCenters 
        + numCenters * numExpansionTerms 
        + numTargets * maxNeighbors * numExpansionTerms;
    std::cout << "IFGT: Estimated speedup " << direct / ifgt << "x" << std::endl;
}

void IFGT::writeClustersToFile(const char* filename) {
    std::ofstream file;
    file.open(filename);
    file << centerToSource.size() << ' ' << sources.size() << std::endl; 
    for (size_t i = 0; i < centerToSource.size(); i++)
        file << i << ' ' << centerToSource[i] << std::endl;
    for (size_t i = 0; i < sources.size(); i++) {
        file << i << ' ' << sourceToCenter[i];
        for (size_t j = 0; j < dim; j++)
            file << ' ' << sources[i][j];
        file << std::endl;
    }
    file.close();
}

/*! Evaluate the discrete gauss transform directly by summing N*M square exponentials 
 * \param targets The points y_j to evaluate (length M)
 * \param gaussTransform The output vector that the results will be saved into */
void IFGT::directEvaluate(const std::vector<Vector>& targets, Vector& gaussTransform) {
    gaussTransform.resize(targets.size(), 0.0);
    double bandwidthSquared = bandwidth * bandwidth;
    for (size_t j = 0; j < targets.size(); j++) {
        for (size_t i = 0; i < sources.size(); i++) {
            double differenceSquared = 0.0;
            for (size_t d = 0; d < dim; d++) {
                double difference = sources[i][d] - targets[j][d];
                differenceSquared += difference * difference;
            }
            gaussTransform[j] += weights[i] * exp(-differenceSquared / bandwidthSquared);
        }
    }
}

/*! Evaluate the 1D square exponential covariance kernel for the Gaussian Process Regression example 
 * \param x1 The x-coordinate of the first point
 * \param x2 The x-coordinate of the second point
 * \param sigma_f The maximum covariance between the two points
 * \param length The length-scale of the Gaussian process
 * \return k(x1, x2) where k is a square exponential covariance kernel */
inline double kernel1D(const double x1, const double x2, const double sigma_f, const double length) {
    double delta = x1 - x2;
    return sigma_f * sigma_f * exp(-delta * delta / (2.0 * length * length));
}

// Solve matrix equation using conjugate-gradient
// and be really clever and use IFGT when evaluating A * p
size_t invertMatrixCGIFGT(const Matrix& sources, const Vector& b, Vector& x, const double length, const double epsilon = 1.0e-16, const size_t MAX_ITER = 1000) {
    Vector r(b.size()); // residue
    Vector p(b.size()); // std::vector along the search direction
    Vector q(b.size()); // A.p
    double e1, e2 = 0;  // epsilons
    double alpha, beta; // minimize x_{k+1} = x_k + alpha p_k

    // EVALUATE A * x
    Vector Ax;
    IFGT Ax_ifgt(sources, x, sqrt(2) * length, 10, 0.2, 4);
    Ax_ifgt.evaluate(sources, Ax);

    r = b - Ax;


    p = r;
    e1 = r * r;
    size_t iter = 0;
    for (iter = 0; iter < MAX_ITER; iter++) {
        

        // EVALUATE A * p
        Vector Ap;
        IFGT Ap_ifgt(sources, p, sqrt(2) * length, 10, 0.2, 4);
        Ap_ifgt.evaluate(sources, Ap);
        
        // q = A * p;
        q = Ap;
        
        alpha = e1 / (p * q);
        x += alpha * p;
        r -= alpha * q;
        e2 = r * r;
        std::cout << "epsilon_2 = " << e2 << std::endl;
        if (fabs(e2) < epsilon) break;
        beta = e2 / e1;
        p = r + beta * p;
        e1 = e2;
    }
    return iter;
}

void gpr_tests() {
    // parameters
    double sigma_f = 0.5; // max covariance
    double sigma_n = 0.1; // noise
    double length = 0.04;

    size_t N = 25;
    size_t dim = 1;
    size_t M = 201;
    
    std::ofstream file;
    file.open("out/gpr");

    // generate some noisy data 
    
    // print to file
    file << N << std::endl;
    Matrix x(N, Vector(dim, 0.0));
    Vector y(N);
    for(size_t i = 0; i < N; ++i) {
        x[i][0] = uniform(0, 1);
        double noisefree = sin(3.0 * M_PI * x[i][0]);
        y[i] = noisefree + normal(sigma_n);
        file << x[i][0] << " " << y[i] << " " << noisefree << std::endl;
    }
    
    std::cout << "constructing covariance matrix" << std::endl;
    
    // Construct NxN covariance matrix
    Matrix K;
    K.resize(N, Vector(N));
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            if (i == j) {
                K[i][j] = sigma_f * sigma_f + sigma_n * sigma_n;
            } else {
                K[i][j] = kernel1D(x[i][0], x[j][0], sigma_f, length);
            }
        }
    }
    
    std::cout << "inverting covariance matrix" << std::endl;
    
    Matrix Kinv; 
    invertMatrix(K, Kinv); 
    
    std::cout << "evaluating interpolants" << std::endl;
    
    file << M << std::endl;
    Vector meanValues(M);
    Vector variances(M);
    Vector noiseFree(M);
    double xstar; 
    for (size_t i = 0; i < M; i++) {
        xstar = 1.0 / (M-1) * i;
        Vector Kstar(N);
        for (size_t j = 0; j < N; j++) {
            Kstar[j] = kernel1D(x[j][0], xstar, sigma_f, length);
        }
        noiseFree[i] = sin(3.0 * M_PI * xstar);
        double mean = Kstar * Kinv * y;
        meanValues[i] = mean;
        variances[i] = K[0][0] - Kstar * Kinv * Kstar;
    }
    
    // Now we do it with CG
    // first compute quantity (Kinv y) by solving Kx=y for x
    // check first with CG without using IFGT 
    Vector Kinvy(N, 1.0);
    int iter = invertMatrixCG(K, y, Kinvy, 1e-16);
    std::cout << "invertMatrixCG used " << iter << " iterations." << std::endl;

    // check against solution mean
    Vector meanValuesCG(M);
    for (size_t i = 0; i < M; i++) {
        xstar = 1.0 / (M-1) * i;
        Vector Kstar(N);
        for (size_t j = 0; j < N; j++) {
            Kstar[j] = kernel1D(x[j][0], xstar, sigma_f, length);
        }
        double mean = Kstar * Kinvy;
        meanValuesCG[i] = mean;

        Vector KinvKstar(N, 1.0);
        iter = invertMatrixCG(K, Kstar, KinvKstar, 1e-15);
        std::cout << "invertMatrixCG used " << iter << " iterations." << std::endl;
        
        variances[i] = K[0][0] - Kstar * KinvKstar;
    }
    
    Vector difference = meanValuesCG - meanValues;
    difference /= meanValues;
    for (size_t i = 0; i < M; i++) {
        xstar = 1.0 / (M-1) * i;
        file << xstar << " " << meanValues[i] << " " << variances[i] << " " 
             << noiseFree[i] << " " << meanValuesCG[i] << std::endl;
    }
    
    
    // Now we do it with CG AND IFGT
    // first compute quantity (Kinv y) by solving Kx=y for x
    Vector Ktest;
    Vector test(N, 1.0 * sigma_f * sigma_f);
    Vector test_unweighted(N, 1.0);
    IFGT Ktest_ifgt(x, test, sqrt(2) * length, 10, 0.2, 4);
    Ktest_ifgt.directEvaluate(x, Ktest);
    printVector(Ktest);
    printVector(K * test_unweighted);
    
    Kinvy.resize(N, 1.0 * sigma_f * sigma_f);
    iter = invertMatrixCGIFGT(x, y, Kinvy, length);
    std::cout << "invertMatrixCGIFGT used " << iter << " iterations." << std::endl;
    
    
}

void ifgt_tests() {
    
    size_t N = 20000;
    size_t dim = 2;
    
    Matrix sources(N, Vector(dim, 0.0));
    Vector weights(N, 1.0);
    
    for (size_t i = 0; i < N; i++) {
        if (uniform(0, 1) < 0.5) {
            for (size_t j = 0; j < dim; j++)
                sources[i][j] = 0.8 + normal(0.1);
        } else {
            for (size_t j = 0; j < dim; j++)
                sources[i][j] = 0.2 + normal(0.1);
        }
        // for (size_t j = 0; j < dim; j++)
        //    sources[i][j] = uniform(0,1);
    }
    

    Vector gaussTransform(N, 0.0);
    Vector directGaussTransform(N, 0.0);
    Vector error;
    
    double ifgt_time, direct_time;
    
    timer_start();
    IFGT ifgt(sources, weights, 0.1, 5,  0.8, 2);
    ifgt.evaluate(sources, gaussTransform);
    timer_stop();
    ifgt_time = timer_value();
    
    std::cout << "ifgt took " << ifgt_time << " sec" << std::endl;
    
    timer_start();
    ifgt.directEvaluate(sources, directGaussTransform);
    timer_stop();
    direct_time = timer_value();
    std::cout << "direct took " << direct_time << " sec" << std::endl;
    
    std::cout << "actual speedup: " << direct_time / ifgt_time << std::endl;
    
    error = gaussTransform - directGaussTransform;
    error = (1.0 / maxValue(directGaussTransform)) * error;
    
    ifgt.writeClustersToFile("out/clusters");
    
    std::ofstream file;
    file.open("out/2d");
    for (size_t i = 0; i < sources.size(); i++) {
        file << sources[i][0] << " " << sources[i][1] << " " << gaussTransform[i] << " " << error[i] << std::endl;
    }
    file.close();
}

void matrix_tests() {
    /* finite difference test for matrix methods */
    size_t N = 120;
    Matrix A; 
    A.resize(N, Vector(N, 0.0));
    for (size_t i = 0; i < N; i++) {
        A[i][i] = 2;
        if (i > 0) A[i-1][i] = -1;
        if (i+1 < N) A[i+1][i] = -1;
    }
    Vector b(N, 0.0);
    b[N/2] = 1.0;
    Vector x(b.size(), 1.0);
    int iters = invertMatrixCG(A, b, x, 1e-13);
    std::cout << "CG took " << iters << " steps." << std::endl;
    printVector(x);
}

int main(int argc, char **argv) {

    srand(time(NULL));

    ifgt_tests();
    // gpr_tests();
    
    return 0;
}


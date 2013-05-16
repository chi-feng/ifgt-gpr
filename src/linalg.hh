#ifndef LINALG_HH
#define LINALG_HH

#include <iostream>
#include <fstream> 
#include <vector>
#include <math.h>
#include <assert.h>

/* poor-man's linear algebra */
typedef std::vector<double> Vector;
typedef std::vector<Vector> Matrix;

/*! Matrix-vector multiplication 
 * \param A The matrix
 * \param x The vector
 * \return The matrix-vector product between A and x
 */
inline Vector operator*(const Matrix& A, const Vector& x) {
    assert(A.size() > 0);
    assert(A[0].size() == x.size());
    Vector product(A.size());
    for (size_t i = 0; i < A.size(); i++) {
        product[i] = 0.0;
        for (size_t j = 0; j < x.size(); j++) {
            product[i] += A[i][j] * x[j];
        }
    }
    return product;
} 

/*! Vector-matrix multiplication 
 * \param x The vector
 * \param A The Matirx
 * \return The matrix-vector product between x and A
 */
inline Vector operator*(const Vector& x, const Matrix& A) {
    assert(A.size() > 0);
    assert(A.size() == x.size());
    Vector product(x.size());
    for (size_t i = 0; i < A.size(); i++) {
        product[i] = 0.0;
        for (size_t j = 0; j < x.size(); j++) {
            product[i] += A[j][i] * x[i];
        }
    }
    return product;
} 

/*! Vector addition 
 * \param a The first vector
 * \param b The second vector
 * \return The sum of the first and second vector */
inline Vector operator+(const Vector& a, const Vector& b) {
    assert(a.size() > 0); 
    assert(b.size() > 0);
    assert(a.size() == b.size());
    Vector sum(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        sum[i] = a[i] + b[i];
    }
    return sum;
}


/*! Vector subtraction 
 * \param a The first vector
 * \param b The second vector
 * \return The difference of the first and second vector */
inline Vector operator-(const Vector& a, const Vector& b) {
    assert(a.size() > 0); 
    assert(b.size() > 0);
    assert(a.size() == b.size());
    Vector sum(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        sum[i] = a[i] - b[i];
    }
    return sum;
}


/*! Scalar multiplication of a vector 
 * \param a The scalar
 * \param x The vector
 * \return The vector x scaled by scalar a */
inline Vector operator*(const double a, const Vector& x) {
    assert(x.size() > 0); 
    assert(a == a); // for NaN 
    Vector product(x.size());
    for (size_t i = 0; i < x.size(); i++) {
        product[i] = x[i] * a;
    }
    return product;
}


/*! Increment a vector by another vector
 * \param a The vector
 * \param b The vector to be added
 * \return The reference to the original vector, incremented */
inline Vector& operator+=(Vector& a, const Vector& b) {
    assert(a.size() > 0); 
    assert(b.size() > 0);
    assert(a.size() == b.size());
    for (size_t i = 0; i < a.size(); i++) {
        a[i] += b[i];
    }
    return a;
}

/*! Decrement a vector by another vector
 * \param a The vector
 * \param b The vector to be subtracted
 * \return The reference to the original vector, decremented */
inline Vector& operator-=(Vector& a, const Vector& b) {
    assert(a.size() > 0); 
    assert(b.size() > 0);
    assert(a.size() == b.size());
    for (size_t i = 0; i < a.size(); i++) {
        a[i] -= b[i];
    }
    return a;
}

/*! Pointwise multiplication of two vectors
 * \param a The first vector
 * \param b The second vector
 * \return The reference to the original vector, 
 * multiplied pointwise by hte second vector */
inline Vector& operator*=(Vector& a, const Vector& b) {
    assert(a.size() > 0); 
    assert(b.size() > 0);
    assert(a.size() == b.size());
    for (size_t i = 0; i < a.size(); i++) {
        a[i] *= b[i];
    }
    return a;
}

/*! Pointwise division of two vectors
 * \param a The first vector
 * \param b The second vector
 * \return The reference to the original vector, 
 * divided pointwise by hte second vector */
inline Vector& operator/=(Vector& a, const Vector& b) {
    assert(a.size() > 0); 
    assert(b.size() > 0);
    assert(a.size() == b.size());
    for (size_t i = 0; i < a.size(); i++) {
        a[i] /= b[i];
    }
    return a;
}


/*! Inner (dot) product between two vectors
 * \param a The first vector
 * \param b The second vector
 * \return The dot product between the first and second vector */
inline double operator*(const Vector& a, const Vector& b) {
    assert(a.size() > 0); 
    assert(b.size() > 0);
    assert(a.size() == b.size());
    double product = 0;
    for (size_t i = 0; i < a.size(); i++) {
        product += a[i] * b[i];
    }
    return product;
}

/*! Prints a vector to the screen 
 * \param x The vector to be printed */
void printVector(const Vector& x) {
    std::cout.precision(3);
    std::cout << std::fixed;
    std::cout << x[0];
    if (x.size() > 1) {
        for (size_t i = 1; i < x.size(); i++) {
            std::cout << ", " << x[i];
        }
    }
    std::cout << std::endl;
}

/*! Prints a matrix to the screen 
 * \param A The matrix to be printed */
void printMatrix(const Matrix& A) {
    std::cout.precision(3);
    std::cout << std::fixed;
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < A[0].size(); j++) {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

/*! Prints a matrix to the screen 
 * \param A The matrix to be printed */
double maxValue(const Vector& x) {
    double max = x[0];
    for (size_t i = 1; i < x.size(); i++) {
        if (x[i] > max) 
            max = x[i];
    }
    return max;
}

// Invert matrix using Gaussian Eliminiation 
void invertMatrix(const Matrix& A, Matrix& Ainv) {
    Matrix Acopy(A);
    Ainv.resize(A.size(), Vector(A.size(), 0.0));

    for (size_t i = 0; i < A.size(); i++) {
         for (size_t j = 0; j < A.size(); j++) {
             if (i == j) 
                 Ainv[i][i] = 1.0; 
             else
                 Ainv[i][j] = 0.0;
         }
     }
     
     for (size_t i = 0; i < A.size(); i++) {
         for (size_t j = 0; j < i; j++) {
             double scale = Acopy[i][j] / Acopy[j][j];
             for (size_t k = 0; k < A.size(); k++) {
                 Acopy[i][k] -= scale * Acopy[j][k];
                 Ainv[i][k] -= scale * Ainv[j][k];
             }
         }
     }

     for (size_t i = 0; i < A.size(); i++) {
         for (size_t j = i + 1; j < A.size(); j++) {
             double scale = Acopy[i][j] / Acopy[j][j];
             for (size_t k = 0; k < A.size(); k++) {
                 Acopy[i][k] -= scale * Acopy[j][k];
                 Ainv[i][k] -= scale * Ainv[j][k];
             }
         }

     }

     for (size_t i = 0; i < A.size(); i++) {
         double scale = 1.0 / Acopy[i][i];
         for (size_t k = 0; k < A.size(); k++) {
             Acopy[i][k] *= scale;
             Ainv[i][k] *= scale;
         }
     }
     
}

/*! Solve matrix equation Ax=b using conjugate-gradient method
 * \param A The matrix
 * \param b The rhs
 * \param x The initial guess for x
 * \param epsilon The threshold for epsilon_2 to exit inner loop
 * \param MAX_ITER the maximum number of iterations
 * \return The number of iterations */
size_t invertMatrixCG(const Matrix& A, const Vector& b, Vector& x, const double epsilon = 1.0e-16, const size_t MAX_ITER = 1000) {
    Vector r(b.size()); // residue
    Vector p(b.size()); // std::vector along the search direction
    Vector q(b.size()); // A.p
    double e1, e2 = 0;  // epsilons
    double alpha, beta; // minimize x_{k+1} = x_k + alpha p_k
    r = b - A * x;
    p = r;
    e1 = r * r;
    size_t iter = 0;
    for (iter = 0; iter < MAX_ITER; iter++) {
        q = A * p;
        alpha = e1 / (p * q);
        x += alpha * p;
        r -= alpha * q;
        e2 = r * r;
        if (fabs(e2) < epsilon) break;
        beta = e2 / e1;
        p = r + beta * p;
        e1 = e2;
    }
    return iter;
}

#endif
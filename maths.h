#ifndef MATHS_H_
#define MATHS_H_

#include <vector>
using matrix = std::vector<std::vector<double>>;
using vectord = std::vector<double>;
using vectort = std::vector<size_t>;
constexpr double epsilon = 1e-10;
 
struct LU
{
    matrix L;
    matrix U;
    vectort P;
};

class Maths
{
public:
    // operations on matrices
    static double multVectors(const vectord &a, const vectord &b);
    static matrix multMatrix(const matrix &x, const matrix &y);
    static vectord multVectorMatrix(const matrix &x, const vectord &y);
    static matrix transpose(const matrix &x);

    // the core of algorithm
    static LU LUdecomposition(const matrix &x);
    static vectord forwardSubstitution(const matrix &L, const vectord &b);
    static vectord backSubstitution(const matrix &U, const vectord &z);

    static double abs(const double &x);
};

#endif
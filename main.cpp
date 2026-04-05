#include <iostream>
#include <vector>
#include <chrono>
#include "maths.h"

class LinearRegression
{
private:
    const matrix& x;
    const vectord& targets;

public:
    LinearRegression(const matrix& x, const vectord& targets) : x(x), targets(targets) {}

    vectord LR()
    {

        matrix tx = Maths::transpose(x);
        matrix A = Maths::multMatrix(tx, x);
        vectord b = Maths::multVectorMatrix(tx, targets);
        vectord bp(b.size());

        LU decomposition = Maths::LUdecomposition(A);

        for (size_t i = 0; i < b.size(); i++)
        {
            bp[i] = b[decomposition.P[i]];
        }

        vectord z = Maths::forwardSubstitution(decomposition.L, bp);
        vectord w = Maths::backSubstitution(decomposition.U, z);

        return w;
    }
};

int main()
{

    auto start = std::chrono::high_resolution_clock::now();

    matrix x = {{3, 3, 1, 1}, {3, 5, 2, 1}, {2, 3, 4, 1}, {2, 4, 1, 1}};
    vectord targets = {30, 70, 35, 90};

    std::cout << "matrix x: " << '\n';
    for (size_t i = 0; i < x.size(); i++)
    {
        for (size_t j = 0; j < x[0].size(); j++)
            std::cout << x[i][j] << " ";
        std::cout << "\n";
    }

    LinearRegression linearRegression(x, targets);
    vectord weights = linearRegression.LR();

    std::cout << "\n";

    std::cout << "weights : " << '\n';
    for (size_t i = 0; i < weights.size(); i++)
    {
        std::cout << weights[i] << " ";
    }
    std::cout << "\n";

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Training took: " << duration.count() << " microseconds" << "\n";
}
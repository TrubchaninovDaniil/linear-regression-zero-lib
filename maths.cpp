#include "maths.h"

double Maths::multVectors(const vectord &a, const vectord &b)
{
    double c = 0;
    for (size_t i = 0; i < a.size(); i++)
        c += a[i] * b[i];

    return c;
}

matrix Maths::multMatrix(const matrix &x, const matrix &y)
{

    size_t M = x.size();
    size_t K = x[0].size();
    size_t N = y[0].size();

    matrix c(M, vectord(N, 0));

    for (size_t i = 0; i < M; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            double sum = 0;
            for (size_t k = 0; k < K; k++)
            {
                sum += x[i][k] * y[k][j];
            }
            c[i][j] = sum;
        }
    }
    return c;
}

vectord Maths::multVectorMatrix(const matrix &x, const vectord &y)
{
    size_t M = x.size();
    size_t cols = x[0].size();

    vectord c(M, 0);

    for (size_t i = 0; i < M; i++)
    {
        double sum = 0;
        for (size_t k = 0; k < cols; k++)
        {
            sum += x[i][k] * y[k];
        }
        c[i] = sum;
    }
    return c;
}

matrix Maths::transpose(const matrix &x)
{
    matrix xt(x[0].size(), vectord(x.size(), 0));
    for (size_t i = 0; i < x.size(); i++)
    {
        for (size_t j = 0; j < x[0].size(); j++)
        {
            xt[j][i] = x[i][j];
        }
    }
    return xt;
}

LU Maths::LUdecomposition(const matrix &x)
{
    matrix U = x;
    matrix L(x.size(), vectord(x.size(), 0));
    vectort P(x.size());

    for (size_t i = 0; i < L.size(); i++)
    {
        L[i][i] = 1;
        P[i] = i;
    }

    for (size_t row = 0; row < U.size() - 1; row++)
    {
        size_t max_row = row;
        for (size_t i = row + 1; i < U.size(); i++)
        {
            if (Maths::abs(U[i][row]) > Maths::abs(U[max_row][row]))
            {
                max_row = i;
            }
        }

        if (max_row != row)
        {
            std::swap(U[row], U[max_row]);
            std::swap(P[row], P[max_row]);

            for (size_t j = 0; j < row; j++)
            {
                std::swap(L[row][j], L[max_row][j]);
            }
        }

        for (size_t row2 = row + 1; row2 < U.size(); row2++)
        {

            if (Maths::abs(U[row][row]) > epsilon)
            {
                double k = U[row2][row] / U[row][row];
                L[row2][row] = k;
                for (size_t col = row; col < U.size(); col++)
                {
                    U[row2][col] -= U[row][col] * k;
                }
            }
        }
    }
    return {L, U, P};
}

vectord Maths::forwardSubstitution(const matrix &L, const vectord &b)
{
    double sum = 0;
    vectord z(b.size(), 0);
    for (size_t i = 0; i < b.size(); i++)
    {
        sum = 0;
        for (size_t j = 0; j < i; j++)
        {
            sum += L[i][j] * z[j];
        }
        z[i] = b[i] - sum;
    }
    return z;
}

vectord Maths::backSubstitution(const matrix &U, const vectord &z)
{
    double sum = 0;

    int n = z.size();
    vectord w(z.size(), 0);
    for (int i = n - 1; i >= 0; i--)
    {
        sum = 0;
        for (int j = i + 1; j < n; j++)
        {
            sum += U[i][j] * w[j];
        }
        if (Maths::abs(U[i][i]) > epsilon)
        {
            w[i] = ((z[i] - sum) / U[i][i]);
        } else {
            w[i] = 0;
        }
    }
    return w;
}

double Maths::abs(const double &x)
{
    return x >= 0 ? x : -x;
}

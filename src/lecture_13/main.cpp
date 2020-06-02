#include <iostream>
#include <fstream>
#include <nr3.h>
#include "ludcmp.h"
#include "banded.h"
#include <fstream>

std::ostream& operator<<(std::ostream& os, const MatDoub& x)
{
    os << "Matrix:" << "\n";
    for (size_t i = 0; i < x.nrows(); i++)
    {
        for (size_t j = 0; j < x.ncols(); j++)
        {
            os << x[i][j] << ", ";
        }
        os << "\n";
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const VecDoub& x)
{
    os << "Vector:" << "\n";

    for (size_t i = 0; i < x.size(); i++)
    {
        os << x[i] << ", ";
    }

    return os << "\n";
}


Int main(Int argc, Char ** argv)
{
    // Start and End of x and y (is a square...)

    std::fstream file("lecture13.csv", ios::out);

    VecDoub x(2);
    x[0] = 0;
    x[1] = 1;

    VecDoub y(2);
    y[0] = 0;
    y[1] = 1;

    auto lambda1 = [&](int N, Doub lambda)
    {
        Doub h = 1.0/N;
        Int n = N-1;
        MatDoub A(std::pow(n,2),(2*n+1), 0.0);

        Int row = 0;

        // This is to fill the banded matrix
        for (Int i = 0; i < n; i++)
        {
            for (Int j = 0; j < n; j++, row++)
            {

                // Diagonal
                A[row][n] = 4 + lambda*std::pow(h,2);

                // Check Left
                if ( j > 0 )
                {
                    A[row][n-1] = -1;
                }
                // Check Right
                if ( j < n-1 )
                {
                    A[row][n+1] = -1;
                }

                if ( i > 0 )
                {
                    A[row][0] = -1;
                }

                if ( i < n-1 )
                {
                    A[row][2*n] = -1;
                }

            }
        }

        return A;
    };

    auto f = [](Doub x, Doub y)
    {
        return 1+x+y;
    };

    auto lambda2 = [&](int N)
    {
        Doub h = 1.0/N;
        Doub xstep = (x[1] - x[0])/N;
        Doub ystep = (y[1] - y[0])/N;

        VecDoub z(std::pow(N-1,2));

        Int row = 0;

        // This is to fill the RHS vector
        for (Int i = 0; i < N-1; i++)
        {
            for (Int j = 0; j < N-1; j++, row++)
            {
                z[row] = std::pow(h,2)*f(x[0] + (1+i)*xstep, y[0] + (1+j)*xstep);
            }
        }
        return z;
    };
    
    static std::vector<double> richardson;
    Int N = 4;

    for (Int i = 0; i < 7; i++, N*=2)
    {
        MatDoub A = lambda1(N, 0);
        VecDoub b = lambda2(N);
        VecDoub z(b.size());

        Bandec ban(A, N-1, N-1);
        ban.solve(b, z);

        Int index_0505 = std::pow(N-1,2)/2;

        richardson.push_back(z[index_0505]);

    }

    std::cout.precision(5);
    N = 4;

    for (Int i = 0; i < (int)richardson.size(); i++, N*=2)
    {
        Doub error      = (i == 0)   ? NaN : richardson[i-1] - richardson[i];
        Doub error_2    = (i < 2) ? NaN : log2(( richardson[i-2] - richardson[i-1] ) / error);
        std::cout       << "N : " << std::setw(6) << N;
        std::cout       << " - " <<"Estimate : " << std::setw(10) << richardson[i] << " - " << " A[i] - A[i-1] " << std::setw(12) <<  error;
        std::cout       << " - " << " Estimate of Order: " << std::setw(10) << error_2;
        std::cout       << " - " << " The error: " << std::setw(10) << (richardson[i-1] - richardson[i]) / (4.0 - 1.0) << std::endl;
    }

    }


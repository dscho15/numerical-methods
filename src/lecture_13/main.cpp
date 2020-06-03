#include <iostream>
#include <fstream>
#include <nr3.h>
#include "ludcmp.h"
#include "banded.h"
#include <fstream>

struct pde_struct
{
    Doub x0, x1, y0, y1;
    Int N;
    Doub h;
    Doub d_x;
    Doub d_y;
    Doub lambda;

    Doub f(Doub x, Doub y)
    {
        return 1+x+y;
    }

    pde_struct(Doub x0_, Doub x1_, Doub y0_, Doub y1_, Int N_, Doub lambda_) : x0(x0_), x1(x1_), y0(y0_), y1(y1_), N(N_), lambda(lambda_)
    {
        h = 1.0/N;
        d_x = (x1-x0)/N;
        d_y = (y1-y0)/N;
    }
};

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

    auto initalize_A_matrix = [](pde_struct & pde)
    {
        Int n = pde.N-1;
        MatDoub A(std::pow(n,2),(2*n+1), 0.0);
        // This is to fill the banded matrix
        Int row = 0;
        for (Int i = 0; i < n; i++)
        {
            for (Int j = 0; j < n; j++, row++)
            {
                // Diagonal
                A[row][n] = 4 + pde.lambda*std::pow(pde.h,2);
                if ( j > 0 ) // Check left
                {
                    A[row][n-1] = -1;
                }
                if ( j < n-1 ) // Check Right
                {
                    A[row][n+1] = -1;
                }
                if ( i > 0 ) //Check lower
                {
                    A[row][0]  =  -1;
                }
                if ( i < n-1 ) //Check right
                {
                    A[row][2*n] = -1;
                }
            }
        }
        return A;
    };

    auto initalize_y_vector = [](pde_struct & pde)
    {
        Int n = pde.N-1;
        Doub h = 1.0/pde.N;

        VecDoub z(std::pow(n,2));

        Int row = 0;

        // This is to fill the RHS vector
        for (Int i = 0; i < n; i++)
        {
            for (Int j = 0; j < n; j++, row++)
            {
                z[row] = std::pow(pde.h,2)*pde.f(pde.x0 + (1+i)*pde.d_x, pde.y0 + (1+j)*pde.d_y);
                /*
                if ( j > 0 ) // Check left
                {
                    z[row] = -1;
                }
                if ( j < n-1 ) // Check Right
                {
                    z[row] = -1;
                }
                if ( i > 0 ) //Check lower
                {
                    z[row]  =  -1;
                }
                if ( i < n-1 ) //Check right
                {
                    z[row] = -1;
                }
                */
            }
        }
        return z;
    };
    
    static std::vector<double> richardson;
    Doub x0 = 0, x1 = 1, y0 = 0, y1 = 1, lambda = 0;
    Int N = 4;

    for (Int i = 0; i < 5; i++, N*=2)
    {
        pde_struct pde(x0, x1, y0, y1, N, lambda);
        MatDoub A = initalize_A_matrix(pde);
        VecDoub b = initalize_y_vector(pde);
        VecDoub z(b.size());

        Bandec ban(A, N-1, N-1);
        ban.solve(b, z);

        Int index_0505 = std::pow(N-1,2)/2;

        richardson.push_back(z[index_0505]);

    }

    std::cout.precision(5);
    N = 4;
    
    file << "N" << ", " << "A[N]" << ", " << "A[N-1] - A[N]" << ", " << "Estimate of Order" << ", " << "Error " << std::endl;

    for (Int i = 0; i < (int)richardson.size(); i++, N*=2)
    {
        Doub error      = (i == 0)   ? NaN : richardson[i-1] - richardson[i];
        Doub error_2    = (i < 2) ? NaN : log2(( richardson[i-2] - richardson[i-1] ) / error);
        std::cout       << "N : " << std::setw(6) << N;
        std::cout       << " - " <<"Estimate : " << std::setw(10) << richardson[i] << " - " << " A[i] - A[i-1] " << std::setw(12) <<  error;
        std::cout       << " - " << " Estimate of Order: " << std::setw(10) << error_2;
        std::cout       << " - " << " The error: " << std::setw(10) << (richardson[i-1] - richardson[i]) / (4.0 - 1.0) << std::endl;

        file << N << ", " << richardson[i] << ", " << error << ", " << error_2 << ", " << (richardson[i-1] - richardson[i]) / (4.0 - 1.0) << std::endl;
    }

    }


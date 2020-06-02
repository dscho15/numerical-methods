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
    // Goal is move from somewhere to..

    std::fstream file("lecture14.csv", ios::out);

    VecDoub x(2);
    x[0] = 0;
    x[1] = 1;

    VecDoub t(2);
    t[0] = 0;
    t[1] = 20;

    auto a = [](Doub x, Doub t)
    {
        return 0;
    };

    auto b = [](Doub x, Doub t)
    {
        return 1;
    };

    auto f = [](Doub x, Doub t)
    {
        return 10 * x * ( 1 - x ) * exp(-t/10);
    };

    auto u_init = [](Int N, VecDoub x)
    {
        VecDoub z(N-1);
        Doub x_d = (x[1] - x[0])/N;

        auto g = [](Doub x)
        {
            return std::pow(x,4);
        };
        
        for (size_t i = 0; i < z.size(); i++)
        {
            Doub x_i = x[0] + x_d * (i+1);
            
            // The function that is used here is g(x)
            z[i] = g(x_i);
        }

        return z;
    };

    auto construct_matrix = [](Int N, VecDoub t, VecDoub x, Doub alpha)
    {

        MatDoub A(N-1, N-1);
        Doub delta_x = (x[1] - x[0])/N;
        Doub delta_t = delta_x;
        Doub r = alpha * delta_t / (std::pow(delta_x,2));

        for (size_t i = 0; i < N-1; i++)
        {
            A[i][i] = 1+r;

            if(i > 0)
            {
                A[i][i-1] = -.5 * r;
            }

            if(i < N-2)
            {
                A[i][i+1] = -.5 * r;
            }
        }

        std::cout << A << std::endl;
        return A;
    };

    auto construct_vector = [&](Int N, VecDoub t, VecDoub x, VecDoub z, Doub t_i, Doub alpha)
    {

        VecDoub y(N-1);
        Doub delta_x = (x[1] - x[0])/N;
        Doub delta_t = delta_x;

        std::cout << delta_x << std::endl;
        std::cout << delta_t << std::endl;
        Doub r = alpha * delta_t / (std::pow(delta_x,2));

        for (size_t i = 0; i < N-1; i++)
        {
            if(i == 0)
            {
                y[i] = .5 * r * z[i+1] + (1 - r) * z[i] + .5 * r * ( a(x[0], t_i) + a(x[0], t_i + delta_t) ) + .5 * delta_t * ( f(x[0], t_i) + f(x[0], t_i + delta_t ) );
            }
            else if(i == N-2)
            {
                y[i] = .5 * r * z[i-1] + (1 - r) * z[i] + .5 * r * ( b(x[1], t_i) + b(x[1], t_i + delta_t) ) + .5 * delta_t * ( f(x[0], t_i) + f(x[0], t_i + delta_t ) );
            }
            else
            {
                Doub x_i = x[0] + (i+1)*delta_x;    
                y[i] = .5 * r * z[i+1] + (1 - r) * z[i] + .5 * r * z[i-1] + .5 * delta_t * ( f(x_i, delta_t) + f(x_i, delta_t) );
            }
        }

        std::cout << y << std::endl;
        return y;
    };

    VecDoub z   = u_init(4, x);
    MatDoub tri = construct_matrix(4, t, x, 1.0);
    VecDoub y   = construct_vector(4, t, x, z, 1, 1.0);

    return 0;
}


#include <iostream>
#include <fstream>
#include <nr3.h>
#include "ludcmp.h"
#include "banded.h"
#include "tridag.h"
#include <functional>

std::ostream &operator<<(std::ostream &os, const MatDoub &x)
{
    os << "Matrix:"
       << "\n";
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

std::ostream &operator<<(std::ostream &os, const VecDoub &x)
{
    os << "Vector:"
       << "\n";

    for (size_t i = 0; i < x.size(); i++)
    {
        os << x[i] << ", ";
    }

    return os << "\n";
}

struct pde_info
{

    Doub t0, t1, x0, x1;
    Int N;
    Doub d_x;
    Doub d_t;
    Doub alpha;

    Doub a(Doub t)
    { 
        return 0.0;
    };

    Doub b(Doub t)
    {
        return 1.0;
    }

    Doub g(Doub x)
    {
        return std::pow(x,4);
    }

    Doub f(Doub x, Doub t)
    {
        return 10 * x * (1 - x) * exp(-t / 10);
    }

    pde_info( Doub t0_, Doub t1_, Doub x0_, Doub x1_, Int N_, Doub alpha_) : t0(t0_), t1(t1_), x0(x0_), x1(x1_), N(N_), alpha(alpha_)
    {
        // c = 1
        d_x = (x1-x0)/N;
        d_t = d_x;
    }

};


Int main(Int argc, Char ** argv)
{
    // Start and End of x and y (is a square...)
    // Goal is move from somewhere to..

    std::fstream file("lecture14.csv", ios::out);

    auto u_init = [](pde_info & pde)
    {
        VecDoub z(pde.N+1);
        
        z[0]        = pde.a(0); 
        z[pde.N]    = pde.b(0);

        for (size_t i = 1   ; i < z.size()-1   ; i++)
        {
            // The function that is used here is g(x)
            Doub x_i     = pde.x0 + pde.d_x * (i);
            z[i]         = pde.g(x_i);
        }

        return z;
    };

    auto construct_tri = [](pde_info & pde, VecDoub & a, VecDoub & b, VecDoub & c)
    {
        // Resize...
        a.resize(pde.N+1);
        b.resize(pde.N+1);
        c.resize(pde.N+1);

        // Calculate consts
        Doub r = pde.alpha * pde.d_t / (std::pow(pde.d_x,2));
        
        // Defined aforehand, hence this is the upper and lower line of the "pool"
        b[0] = b[pde.N] = 1.0;
        a[0] = c[0] = a[pde.N] = c[pde.N] = 0;

        for (size_t i = 1; i < b.size()-1; i++)
        {
            b[i] = 1+r;
            a[i] = c[i] = -.5 * r;
        }
    };

    auto construct_vec = [](pde_info &pde, Doub t_i, const VecDoub & u,  VecDoub &y) 
    {
        // Resize...
        y.resize(pde.N + 1);

        // Calculate consts
        Doub r          = pde.alpha * pde.d_t / (std::pow(pde.d_x, 2));

        // Definitions, that is, the upper and lower limit of the "pool"
        y[0]            = pde.a(t_i);
        y[pde.N]        = pde.b(t_i);

        // Dynamics definitions, that is, the within values of the "pool"
        for (size_t i = 1; i < y.size() - 1; i++)
        {
            Doub x_i    = pde.x0 + pde.d_x * (i);
            y[i]        = .5*r * ( u[i-1] + u[i+1] ) + (1-r) * u[i]+ .5 * pde.d_t * ( pde.f(x_i, t_i) + pde.f(x_i, t_i - pde.d_t) ) ;
        }
    };

    VecDoub a, b, c, y;

    Doub x0 = 0.0, x1 = 1.0, t0 = 0.0, t1 = 20.0, alpha = 1.0, atol = 1e-4;
    Int N = 2;
    Doub A1 = NaN, A2 = NaN, A3 = NaN, error = NaN;

    file << "N" << ", " << "A[N]" << ", " << "A[N-1] - A[N]" << ", " << "Estimate of Order" << ", " << "Error " << std::endl;

    while( error > atol || N <= 4 )
    {
        pde_info pde(t0, t1, x0, x1, N, alpha);

        // Time start
        Doub t = t0 + pde.d_t;

        // Construct the initial u(xj,0) vector
        // and then the corresponding triangular matrix
        auto u = u_init(pde);
        construct_tri(pde, a, b, c);

        // The update loop to do LUD
        while (t <= t1)
        {
            construct_vec(pde, t, u, y);
            tridag(a, b, c, y, u);
            t += pde.d_t;
        }

        A3 = A2;
        A2 = A1;
        A1 = u[u.size()/2];

        error = abs(A1 - A2);

        Doub error_1 = (N == 2) ? NaN : A2 - A1;
        
        Doub error_2 = (N < 8) ? NaN : log2((A3 - A2) / error_1 );

        std::cout << "N : " << std::setw(6) << N;
        std::cout << " - "
                  << "Estimate : " << std::setw(10) << A1 << " - "
                  << " A[i-1] - A[i] " << std::setw(12) << error;
        std::cout << " - "
                  << " Estimate of Order: " << std::setw(10) << error_2;
        std::cout << " - "
                  << " The error: " << std::setw(10) << (error_1) / (4.0 - 1.0) << std::endl;

        error = (error_1) / (4.0 - 1.0);

        file << N << ", " << A1 << ", " << error << ", " << error_2 << ", " << (error_1) / (4.0 - 1.0) << std::endl;

        N *= 2;
    }

    return 0;
}


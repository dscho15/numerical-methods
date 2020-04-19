#include <iostream>
#include <functional>
#include <cmath>
#include <vector>
#include <nr3.h>
#include <quadrature.h>
#include <derule.h>
#include <svd.h>
#include <fstream>

const Doub T1      = 1000;
const Doub T2      = 500;
const Doub eps1    = 0.80;
const Doub eps2    = 0.60;
const Doub sigma   = 0.1712*pow(10,-9);
const Doub d       = 1.0;
const Doub w       = 1.0;
const Doub a       = -w/2;
const Doub b       = w/2;
const Doub eps1_b  = (1-eps1);
const Doub eps2_b  = (1-eps2);
const Doub K1      = eps1*sigma*std::pow(T1,4);
const Doub K2      = eps2*sigma*std::pow(T2,4);

Doub trapezoidal(int k)
{
    Int N = std::pow(2,k);
    Doub h = (b-a)/N;

    NRvector<Doub> xy(2*(N+1));

    for (Int i = 0; i < N+1; i++)
    {
        xy[i] = a+i*(b-a)/N;
        xy[i+N+1] = xy[i];
    }

    NRmatrix<Doub> A(2*(N+1), 2*(N+1));
    NRvector<Doub> x(2*(N+1));
    NRvector<Doub> y(2*(N+1));

    NRmatrix<Doub> U(2*(N+1), N+1);
       
    NRvector<Doub> I(2*(N+1));
    NRvector<Doub> x_h(2*(N+1));
    NRvector<Doub> Q(2);
    Q[0] = 0; Q[1] = 0;

    auto F = [&](Doub x, Doub y)
    {
        return 0.5 / pow(d*d + pow(x-y, 2), 1.5);
    };

    for (Int i = 0; i < 2*(N+1); i++)
    {
        if(i < N+1)
        {
            y[i] = -K1;
            for (Int j = 0; j < 2*(N+1); j++)
            {
                if(j < N+1)
                    if(j == 0 || j == N)
                    {
                        A[i][j] = 0.5*h*eps1_b*F(xy[i],xy[j]);
                        U[i][j] = 0.5*h*F(xy[i],xy[j]);
                        U[i+N+1][j] = 0.5*h*F(xy[i],xy[j]);
                    }
                    else
                    {
                        A[i][j] = h*eps1_b*F(xy[i],xy[j]);
                        U[i][j] = h*F(xy[i],xy[j]);
                        U[i+N+1][j] = h*F(xy[i],xy[j]);
                    }
                else
                    if(i == j-(N+1)) 
                        A[i][j] = -1.0;
                    else
                        A[i][j] = 0.0;
            }  
        }
        else
        {
            y[i] = -K2;
            for (Int j = 0; j < 2*(N+1); j++)
            {
                if(j < N+1)
                    if(i-(N+1) == j) 
                        A[i][j] = -1.0;
                    else
                        A[i][j] = 0.0;
                else
                    if(j == N+1 || j == 2*N+1 )
                        A[i][j] = h*0.5*eps2_b*F(xy[i],xy[j]);
                    else
                        A[i][j] = h*eps2_b*F(xy[i],xy[j]);
            }
        }
    }

    SVD svd(A);
    svd.solve(y, x);

    for (size_t i = 0; i < 2*(N+1); i++)
    {
        double sum = 0;
        for (size_t j = 0; j < N+1; j++)
        {
            if(i < N+1)
                sum += U[i][j] * x[j];
            else
                sum += U[i][j] * x[j+N+1];
        }

        I[i] = h*sum;
        x_h[i] = h*x[i];

        if(i == 0 || i == N || i == N+1 || i == 2*N+1 )
        {
            I[i] *= 0.5;
            x_h[i] *= 0.5;
        }

        if(i < N+1)
            Q[0] += x_h[i] - I[i];
        else
            Q[1] += x_h[i] - I[i];
    }

    std::cout << "Q:= " << Q[0] << " \t " << Q[1] << std::endl;

}



int main()
{
    std::cout.precision(5);
    
    trapezoidal(1);
    trapezoidal(2);
    trapezoidal(3);
    trapezoidal(4);
    trapezoidal(5);
    trapezoidal(6);
    trapezoidal(7);
    trapezoidal(8);
    trapezoidal(9);

    return 1;
}
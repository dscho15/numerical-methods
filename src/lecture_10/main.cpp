#include <iostream>
#include <fstream>
#include <functional>
#include <nr3.h>
#include <new>
#include "ludcmp.h"
#include "qrdcmp.h"
#include <roots_multidim.h>

Doub function_calls = 0;

VecDoub func(VecDoub & x)
{
    
    VecDoub f(x.size());

    if( x.size() == 0 )
        return VecDoub(2);
    
    function_calls += 1;

    f[0] = x[0]*x[1];
    f[1] = -x[0]*x[0];

    return f;
}

VecDoub func(VecDoub_I & x)
{
    
    VecDoub f(x.size());

    if( x.size() == 0 )
        return VecDoub(2);
    
    function_calls += 1;

    f[0] = x[0]*x[1];
    f[1] = -x[0]*x[0];

    return f;
}


VecDoub func(VecDoub && x)
{
    VecDoub f(x.size());

    if( x.size() == 0 )
        return VecDoub(2);

    f[0] = x[0]*x[1];
    f[1] = -x[0]*x[0];

    function_calls += 1;

    return f;
}

VecDoub operator/(const VecDoub & a, const Doub && b)
{
    VecDoub y(a.size());
    for (size_t i = 0; i < a.size(); i++)
        y[i] = a[i] / b;
    return y;
}

VecDoub operator*(const VecDoub & a, const VecDoub & b)
{
    VecDoub y(a.size());
    for (size_t i = 0; i < a.size(); i++)
        y[i] = a[i]*b[i];
    return y;
}

VecDoub operator*(const VecDoub & a, const VecDoub && b)
{
    VecDoub y(a.size());
    for (size_t i = 0; i < a.size(); i++)
        y[i] = a[i]*b[i];
    return y;
}

VecDoub operator*(const VecDoub & a, const Doub & b)
{
    VecDoub y(a.size());
    for (size_t i = 0; i < a.size(); i++)
        y[i] = a[i]*b;
    return y;
}

VecDoub operator*(const VecDoub & a, const Doub && b)
{
    VecDoub y(a.size());
    for (size_t i = 0; i < a.size(); i++)
        y[i] = a[i]*b;
    return y;
}

VecDoub operator+(const VecDoub & a, const VecDoub & b)
{
    VecDoub y(a.size());
    for (size_t i = 0; i < a.size(); i++)
        y[i] = a[i] + b[i];
    return y;
}

VecDoub operator+(const VecDoub & a, const VecDoub && b)
{
    VecDoub y(a.size());
    for (size_t i = 0; i < a.size(); i++)
        y[i] = a[i] + b[i];
    return y;
}

VecDoub operator-(const VecDoub & a, const VecDoub & b)
{
    VecDoub y(a.size());
    for (size_t i = 0; i < a.size(); i++)
        y[i] = a[i] - b[i];
    return y;
}

VecDoub operator-(const VecDoub & a, const VecDoub && b)
{
    VecDoub y(a.size());
    for (size_t i = 0; i < a.size(); i++)
        y[i] = a[i] - b[i];
    return y;
}

void richardson_extrapolation(Doub N, Doub A_h3, Doub A_h2, Doub A_h1)
{
    std::cout << "N: " << setw(15) << N << " - A: " << setw(15) << A_h3 << " - A(last)-A(new): " << setw(15) 
    << A_h2 - A_h3 << " - Rich-alpha^k: " << setw(15)<< log2((A_h1 - A_h2) / (A_h2 - A_h3)) << " - f: " << setw(15)
    << function_calls << std::endl;
    
    function_calls = 0;
};

void euler(const VecDoub y, Doub max)
{
    Doub N = 4;

    std::ofstream file("lecture10_euler.csv");
    VecDoub x_prev[2];
    VecDoub x = y;
    x_prev[0].resize(2);
    x_prev[1].resize(2);

    x_prev[0][0] = NAN; x_prev[0][1] = NAN;
    x_prev[1][0] = NAN; x_prev[1][1] = NAN;

    std::cout << "euler ---------------------------" << std::endl;

    for (size_t i = 0; i < 18; i++)
    {
        Doub h = max/N;
        for (size_t j = 0; j < N; j++)
        {
            Doub t = h*j;
            x = x + func(x)*h;
        }

        richardson_extrapolation(N, x[0], x_prev[0][0], x_prev[1][0]);
        x_prev[1] = x_prev[0];
        x_prev[0] = x;
        x[0] = 1.0;
        x[1] = 1.0;
        N *= 2;
    }
}

void midpoint(const VecDoub & y, Doub max)
{
    Doub N = 4;

    std::ofstream file("lecture10_midpoint.csv");
    VecDoub k1(y.size());
    VecDoub k2(y.size());
    VecDoub x_prev[2];
    VecDoub x = y;
    x_prev[0].resize(2);
    x_prev[1].resize(2);

    x_prev[0][0] = NAN; x_prev[0][1] = NAN;
    x_prev[1][0] = NAN; x_prev[1][1] = NAN;

    std::cout << "midpoint ---------------------------" << std::endl;

    for (size_t i = 0; i < 13; i++)
    {
        Doub h = max/N;
        for (size_t j = 0; j < N; j++)
        {
            Doub t  = h*j;
            k1 = func(x);
            k2 = func(x+k1*h/2.0);
            x = x + k2*h;
            //file << t << ", " << x[0] << ", " << x[1] << std::endl;
        }

        richardson_extrapolation(N, x[0], x_prev[0][0], x_prev[1][0]);
        x_prev[1] = x_prev[0];
        x_prev[0] = x;
        x[0] = 1.0;
        x[1] = 1.0;
        N *= 2;
    }
}

void runge_kutta_4th(VecDoub & y, Doub max)
{
    Doub N = 4;

    std::ofstream file("lecture10_runge_kutta.csv");

    // x_half 
    VecDoub k1(y.size());
    VecDoub k2(y.size());
    VecDoub k3(y.size());
    VecDoub k4(y.size());
    VecDoub x_prev[2];
    VecDoub x = y;
    x_prev[0].resize(2);
    x_prev[1].resize(2);

    x_prev[0][0] = NAN; x_prev[0][1] = NAN;
    x_prev[1][0] = NAN; x_prev[1][1] = NAN;

    std::cout << "runge_kutta_4th ---------------------------" << std::endl;

    for (size_t i = 0; i < 14; i++)
    {
        Doub h = max/N;
        for (size_t j = 0; j < N; j++)
        {
            // This is the second parameter!!!
            Doub t  = h*j;
            // Apply Runge Kutta Formulas to find 
            // next value of y 
            k1 = func(x);
            k2 = func(x+k1*h*.5);
            k3 = func(x+k2*h*.5);
            k4 = func(x+k3*h);
            x = x + (k1+k2*2+k3*2+k4)*h*1/6.0;
            //file << t << ", " << x[0] << ", " << x[1] << std::endl;
        }

        richardson_extrapolation(N, x[0], x_prev[0][0], x_prev[1][0]);
        x_prev[1] = x_prev[0];
        x_prev[0] = x;
        x[0] = 1.0;
        x[1] = 1.0;
        N *= 2;
    }
}

void leap_frog(VecDoub & x, Doub max)
{
    Doub t = 0;
    Doub h = 0.001;
    Doub it = max/h;

    std::ofstream file("lecture10_leap_frog.csv");

    // y_n-1 y_n
    VecDoub y[2];

    y[0] = x;
    y[1] = x + func(x)*h;


    for (size_t i = 1; i <= it/2; i++)
    {
        // This is the second parameter!!!
        t = 2*h*i;
        
        // Apply Runge Kutta Formulas
        x = y[0] + func(y[1])*2.0*h;

        y[0] = y[1];
        y[1] = x;

        file << t << ", " << x[0] << ", " << x[1] << std::endl;
    }
}

void trapezoidal(VecDoub & y, Doub max)
{
    Doub N = 4;
    VecDoub x_prev[2];
    VecDoub x = y;
    x_prev[0].resize(2);
    x_prev[1].resize(2);

    x_prev[0][0] = NAN; x_prev[0][1] = NAN;
    x_prev[1][0] = NAN; x_prev[1][1] = NAN;
    std::cout << "trapezoidal ---------------------------" << std::endl;

    for (size_t i = 0; i < 14; i++)
    {

        Doub h = max/N;

        for (size_t j = 0; j < N; j++)
        {
            VecDoub dx = func(x);
            VecDoub x_nxt = x+dx*h;

            function<VecDoub(VecDoub_I &)> phi = [&](VecDoub_I &z) -> VecDoub 
            {
                return z - x - ( dx + func(z) ) * h * 0.5;
            };

            bool newt_failed;
            newt(x_nxt, newt_failed, phi);
            x = x_nxt;
            Doub t = j * h;
        }

        richardson_extrapolation(N, x[0], x_prev[0][0], x_prev[1][0]);
        x_prev[1] = x_prev[0];
        x_prev[0] = x;
        x[0] = 1.0;
        x[1] = 1.0;
        N *= 2;
    }
}


Int main(Int argc, Char ** argv)
{

    VecDoub x(2);
    
    x[0] = 1.0;
    x[1] = 1.0;

    euler(x, 20.0);
    
    x[0] = 1.0;
    x[1] = 1.0;

    midpoint(x, 20.0);

    x[0] = 1.0;
    x[1] = 1.0;
    
    runge_kutta_4th(x, 20.0);

    x[0] = 1.0;
    x[1] = 1.0;
    leap_frog(x, 20.0);

    x[0] = 1.0;
    x[1] = 1.0;

    trapezoidal(x, 20);

}


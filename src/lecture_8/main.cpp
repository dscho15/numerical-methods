#include <iostream>
#include <functional>
#include <cmath>
#include <vector>

typedef double (*exercise)(double x);
typedef double (*integration)(const double x_0, const double x_n, const int steps, exercise func);

double exercise1(double x)
{
    return std::cos(std::pow(x,2))*std::exp(-x);
}

double exercise2(double x)
{
    return std::sqrt(x)*std::cos(std::pow(x,2))*std::exp(-x);
}

double exercise3(double x)
{
    return 1/std::sqrt(x)*std::cos(std::pow(x,2))*std::exp(-x);
}

double exercise4(double x)
{
    return 1000*std::exp(-1/x)*std::exp(1/(1-x));
}




template<class T>
double trapezoidal_integration(const double x_0, const double x_n, const int steps, T func )
{

    double step_size = (x_n-x_0)/steps;
    double sum = 0;

    for (size_t i = 0; i < steps; i++)
    {
        sum += step_size*(0.5*func(x_0 + step_size*i)+0.5*func(x_0+step_size*(i+1)));
        //std::cout << x_0 + step_size*i << " " << x_0+step_size*(i+1) << std::endl;
    }

    return sum;
}

template<class T>
double rectangular_integration(const double x_0, const double x_n, const int steps, T func )
{
    double step_size = (x_n-x_0)/steps;
    double step_pos = 0.5*step_size + x_0;
    double sum = step_size*func(step_pos);

    for (size_t i = 1; i < steps; i++)
    {
        step_pos += step_size;
        sum += step_size*func(step_pos);
        //std::cout << step_pos << std::endl;
    }
    
    return sum;
}

// Richardson Extrapolation

template<class T>
double simpsons_integration(const double x_0, const double x_n, const int steps, T func )
{
    double step_size = (x_n-x_0)/steps;
    double sum = 0;

    for (size_t i = 0; i < steps; i++)
    {
        sum += step_size / 6.0 * ( func( x_0 + step_size * i ) + 4.0 * func( x_0 + step_size * ( i + 0.5 ) ) + func( x_0 + step_size * ( i + 1 ) ) );
       //std::cout << x_0 + step_size * i << " " << x_0 + step_size * ( i + 0.5 ) << " " << x_0 + step_size * ( i + 1 ) << std::endl;
    }   

    return sum;
}

template<class T>
double error_est(const double x_0, const double x_n, const double alpha, T func, exercise ptr)
{
    std::vector<double> A;
    A.push_back(func(x_0, x_n, 1, ptr));
    A.push_back(func(x_0, x_n, 2, ptr));
    int inc = 2;
    int j = 0;
    double k = 0;

    std::cout << "Estimating error now :=" << std::endl;

    for (size_t i = 1; i < 10; i++)
    {
        inc += std::pow(2,i);
        A.push_back(func(x_0, x_n, inc, ptr));
        k = round(log2((A[j]-A[j+1])/(A[j+1]-A[j+2])));
        std::cout << "k := " << k << std::endl;
        std::cout << "epsilon := " << (A[j+1]-A[j])/(std::pow(alpha,k)-1) << std::endl;
        j++;
    }

    return (A[j+1]-A[j])/(std::pow(alpha,k)-1);
}


int main()
{
    std::cout.precision(15);

    integration inte = &trapezoidal_integration;
    double ex1_trap = error_est(0.0,1.0, 2.0, inte, exercise1);

    inte = &rectangular_integration;
    double ex1_rect = error_est(0.0, 1.0, 2.0, inte, exercise1);

    inte = &simpsons_integration;
    double ex1_simp = error_est(0.0, 1.0, 2.0, inte, exercise1);

    //
    inte = &trapezoidal_integration;
    double ex2_trap = error_est(0.0,1.0, 2.0, inte, exercise2);

    inte = &rectangular_integration;
    double ex2_rect = error_est(0.0, 1.0, 2.0, inte, exercise2);

    inte = &simpsons_integration;
    double ex2_simp = error_est(0.0, 1.0, 2.0, inte, exercise2);

    //
    inte = &trapezoidal_integration;
    double ex3_trap = error_est(0.001,1.0, 2.0, inte, exercise3);

    inte = &rectangular_integration;
    double ex3_rect = error_est(0.0, 1.0, 2.0, inte, exercise3);

    inte = &simpsons_integration;
    double ex3_simp = error_est(0.001, 1.0, 2.0, inte, exercise3);

    //
    inte = &trapezoidal_integration;
    double ex4_trap = error_est(0.001, 1.0, 2.0, inte, exercise4);

    inte = &rectangular_integration;
    double ex4_rect = error_est(0.001, 1.0, 2.0, inte, exercise4);

    inte = &simpsons_integration;
    double ex4_simp = error_est(0.001, 1.0, 2.0, inte, exercise4);


    return 1;
}
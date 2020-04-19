#include <iostream>
#include <functional>
#include <cmath>
#include <vector>
#include <nr3.h>
#include <quadrature.h>
#include <derule.h>

Ldoub exercise1(Ldoub x, Ldoub delta)
{
    return cos(x*x)*exp(-x);
}

Ldoub exercise2(Ldoub x, Ldoub delta)
{
    return sqrt(x)*cos(x*x)*exp(-x);
}

Ldoub exercise3(Ldoub x, Ldoub delta)
{
    return 1000*exp(-1/x)*exp(-1/(1-x));
}

Ldoub exercise4(Ldoub x, Ldoub delta)
{
    if(abs(x) < 0.001)
        return 1/sqrt(delta)*cos(delta*delta)*exp(-delta);
    return 1/sqrt(x)*cos(x*x)*exp(-x);
}

void error_est(DErule<Ldoub (Ldoub x, Ldoub Delta)> & exe)
{
    auto epsilon_trim_lambda = [&]()
    {
        auto q = [](Ldoub t)
        {
            return exp(-2.0*sinh(t));
        };

        auto del = [&](Ldoub t)
        {
            return (exe.b-exe.a)*q(t)/(1.0+q(t));
        };

        auto x = [&](Ldoub t)
        {
            return 1/2*(exe.b-exe.a) + 1/2*(exe.b-exe.a)*tanh(sinh(t));
        };

        return (exe.func(x(-exe.hmax+del(-exe.hmax)), del(-exe.hmax)) + exe.func(x(exe.hmax+del(exe.hmax)), del(exe.hmax))) * exp(-exp(exe.hmax));
    };

    auto epsilon_disc_lambda = [&]()
    {
        return ( exp(-M_PI*M_PI*exe.n/exe.hmax) );
    };

    for (size_t i = 0; i < 5; i++)
    {
        Ldoub A_last = exe.s;
        exe.next();
        std::cout << "A: " << exe.s << " \t |A-A_i|: " << abs(exe.s-A_last) << " \t epsilon_trim " << epsilon_trim_lambda() << " \t epsilon_disc " << epsilon_disc_lambda() << std::endl;
    }    
}



int main()
{
    std::cout.precision(15);

    DErule exe1(exercise1, 0, 1, 4.3);
    error_est(exe1);

    std::cout << "<---------------------->" << std::endl;

    DErule exe2(exercise2, 0, 1, 4.3);
    error_est(exe2);

    std::cout << "<---------------------->" << std::endl;

    DErule exe3(exercise3, 0, 1, 4.3);
    error_est(exe3);

    std::cout << "<---------------------->" << std::endl;

    DErule exe4(exercise4, 0.0, 1.0, 4.3);
    error_est(exe4);
    std::cout << "<---------------------->" << std::endl;
    
    return 1;
}
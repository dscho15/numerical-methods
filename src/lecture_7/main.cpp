#include <iostream>
#include <curses.h>
#include <cmath>
#include <nr3.h>
#include <ludcmp.h>
#include <qrdcmp.h>
#include <roots_multidim.h>
#include <utilities.h>

// Fixed values by default
constexpr double v      = 120;
constexpr double k      = 2.5;
constexpr double w      = 4;
constexpr double alpha  = 2*10e-7;


// Fixed Values by usr
double d      = 0;
double n      = 0;


VecDoub F(const VecDoub & state)
{
    VecDoub out(8);

    Doub p   = state[0],    L     = state[1],
         x   = state[2],    a     = state[3],
         phi = state[4],    theta = state[5],
         L0  = state[6],    H     = state[7];

    out[0] = p                  - a * ( std::cosh( x / a ) - 1 );
    out[1] = L                  - 2 * a * std::sinh( x / a );
    out[2] = d                  - ( 2 * x + 2 * k * std::cos(theta) );
    out[3] = n                  - ( p + k * std::sin( theta ) );
    out[4] = std::tan(phi)      - std::sinh(x/a);
    out[5] = std::tan(theta)    - (1 + v / ( w * L0 ) * std::tan( phi ) );
    out[6] = L                  - L0 * ( 1 + alpha * H );
    out[7] = H                  - w * L0 / ( 2 * std::sin( phi ) );

    return out;
}

void fill_vec8(VecDoub & q_init, const Doub & p, const Doub & L, const Doub & x, const Doub & a, const Doub & phi, const Doub & theta, const Doub & L0, const Doub & H)
{
    q_init.resize(8);
    q_init[0] = p;
    q_init[1] = L;
    q_init[2] = x;
    q_init[3] = a;
    q_init[4] = phi;
    q_init[5] = theta;
    q_init[6] = L0;
    q_init[7] = H;
}

void print_state(const VecDoub & state)
{
    printf("| p = %.5f | L = %.5f | x = %.5f | a = %.5f | \u03C6 = %.5f | \u03B8 = %.5f | L0 = %.5f | H = %.5f | \n", state[0], state[1], state[2], state[3], state[4], state[5], state[6], state[7]);
}


int main()
{
    VecDoub q_init(8);

    Doub p, L, x, a, phi, theta, L0, H;

    p = 0.5; L = 28; x = 15; a = 40; H = 100;
    phi = M_PI_4; theta = M_PI_2; L0 = 30; H = 100;

    fill_vec8(q_init, p, L, x, a, phi, theta, L0, H);

    for ( double n_t : {5.0, 2.0, 1.0, 0.2, 0.1} )
    {
        d   = 30;
        n   = n_t;

        // Use the previous defined q_init rather than filling the vector once agin
        Bool check = true;
        newt(q_init, check, F);

        print_state(q_init);
        util::print(F(q_init), "Error plugging state into F(q):");
        std::cout << "---------------------------------------------------------------------------------------------------------------------------" << std::endl;

    }

    // The results seems to increase in H, when n is decreasing
    // this is reasonable, hence p

    // What happens if your function is very expensive to compute?
    // If it was very expensive, then 


    // What happens if your function is stochastic?
    
    // Assume you’re able to compute the Jacobian.  
    // Does this help in the above two cases?

    return 0;
}


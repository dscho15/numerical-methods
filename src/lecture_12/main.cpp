#include <iostream>
#include <fstream>
#include "nr3.h"
#include "tridag.h"
#include "ludcmp.h"

using namespace std;

struct two_point_boundary_problem
{
    Doub F(Doub yd, Doub y, Doub x)
    {
        return sin(yd) - cos(y) + 2*x ;
    };

    Doub F_y(Doub yd, Doub y, Doub x)
    {
        return sin(y);
    };

    Doub F_yd(Doub yd, Doub y, Doub x)
    {
        return cos(yd);
    };

    Doub jacobian(VecDoub & y)
    {
        // Defines
        Int N = y.size();
        a.resize(N);
        b.resize(N);
        c.resize(N);

        Doub h = (x1 - x0)/((Doub)(N+1.0));
        Doub dy = (y[1] - alpha) / (2*h);

        b[0] = 2 + std::pow(h,2) * F_y(dy, y[0], x0 + h);
        c[0] = -1 + (.5*h) * F_yd(dy, y[0], x0 + h);
        
        dy = (beta - y[N-2])/(2*h);

        b[N-1] = 2 + std::pow(h,2) * F_y(dy, y[N-1], x1 - h);
        a[N-1] = -1 - (.5*h) * F_yd(dy, y[N-1], x1 - h);

        for (Int i = 1; i < N-1; i++)
        {
            Doub x = x0 + h*(i+1);
            dy     = (y[i+1] - y[i-1])/(2*h);
            a[i]   = -1.0 - .5*h * F_yd(dy, y[i], x);
            b[i]   = 2.0 + std::pow(h, 2.0) * F_y(dy, y[i],  x);  
            c[i]   = -1.0 + .5*h * F_yd(dy, y[i], x);
        }
    };

    VecDoub func(VecDoub_I & y)
    {
        // Okay, so basically y[0] = alpha, y[1] = x, y[2] = beta,
        // but y is not incorporating this, so basically, y a size smaller (N=2 -> N=1)
        // so y[0] = x 

        Int N = y.size();
        VecDoub phi(N);
        Doub h = (x1 - x0)/((Doub)(N+1.0));

        if(N == 1)
        {
            Doub dy = (beta - alpha)/(2*h);
            phi[0]  = -alpha + y[0] - beta + std::pow(h,2) * F(dy, y[0], x0 + h); 
        }
        else
        {
            Doub dy   = (y[1]-alpha)/(2*h);
            phi[0]    = 2*y[0] - y[1] + std::pow(h,2) * F(dy, y[0], x0 + h) - alpha;

            dy        = (beta - y[N-2])/(2*h);
            phi[N-1]  = 2*y[N-1] - y[N-2] + std::pow(h,2) * F(dy, y[N-1], x1 - h) - beta;

            for (size_t i = 1; i < N-1; i++)
            {
                Doub x    = x0 + (i+1)*h;
                dy        = (y[i+1] - y[i-1])/(2*h);
                phi[i]    = 2*y[i] - y[i-1] - y[i+1] + std::pow(h,2) * F(dy, y[i], x);
            }
        }

        return phi;
    };

    two_point_boundary_problem(Doub x0_, Doub x1_, Doub alpha_, Doub beta_) : x0(x0_), x1(x1_), alpha(alpha_), beta(beta_) 
    {

    };

    // Used for tridiag
    VecDoub a, b, c;

    Doub x0, x1, alpha, beta;
};

template <class T>
void lnsrch(VecDoub_I &xold, const Doub fold, VecDoub_I &g, VecDoub_IO &p,
VecDoub_O &x, Doub &f, const Doub stpmax, Bool &check, T &func) {
	const Doub ALF=1.0e-4, TOLX=numeric_limits<Doub>::epsilon();
	Doub a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
	Doub rhs1,rhs2,slope=0.0,sum=0.0,temp,test,tmplam;
	Int i,n=xold.size();
	check=false;
	for (i=0;i<n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=0;i<n;i++)
			p[i] *= stpmax/sum;
	for (i=0;i<n;i++)
		slope += g[i]*p[i];
	if (slope >= 0.0) throw("Roundoff problem in lnsrch.");
	test=0.0;
	for (i=0;i<n;i++) {
		temp=abs(p[i])/MAX(abs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;
	for (;;) {
		for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
		f=func(x);
		if (alam < alamin) {
			for (i=0;i<n;i++) x[i]=xold[i];
			check=true;
			return;
		} else if (f <= fold+ALF*alam*slope) return;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(f-fold-slope));
			else {
				rhs1=f-fold-alam*slope;
				rhs2=f2-fold-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc < 0.0) tmplam=0.5*alam;
					else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
					else tmplam=-slope/(b+sqrt(disc));
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = f;
		alam=MAX(tmplam,0.1*alam);
	}
}


struct NRfmin {
	VecDoub fvec;
    two_point_boundary_problem * ode2;
	Int n;
	NRfmin(two_point_boundary_problem * ode2_) : ode2(ode2_)
    {}
	Doub operator() (VecDoub_I &x) {
		n=x.size();
		Doub sum=0;
		fvec=ode2->func(x);
		for (Int i=0;i<n;i++) sum += SQR(fvec[i]);
		return 0.5*sum;
	}
};

void newt(VecDoub_IO &x, Bool &check, two_point_boundary_problem * ode2) 
{
	const Int MAXITS=200;
	const Doub TOLF=1.0e-8,TOLMIN=1.0e-12,STPMX=100.0;
	const Doub TOLX=numeric_limits<Doub>::epsilon();
	Int i,j,its,n=x.size();
	Doub den,f,fold,stpmax,sum,temp,test;
	VecDoub g(n),p(n),xold(n);
	MatDoub fjac(n,n);
	NRfmin fmin(ode2);
	VecDoub &fvec=fmin.fvec;
	f=fmin(x);
	test=0.0;
	for (i=0;i<n;i++)
		if (abs(fvec[i]) > test) test=abs(fvec[i]);
	if (test < 0.01*TOLF) {
		check=false;
		return;
	}
	sum=0.0;
	for (i=0;i<n;i++) sum += SQR(x[i]);
	stpmax=STPMX*MAX(sqrt(sum),Doub(n));
	for (its=0;its<MAXITS;its++) 
    {
        ode2->jacobian(x);
        for (i=0;i<n;i++) 
		{
			sum=0.0;
			if(i == 0)
				sum += ode2->b[i]*fvec[i] + ode2->c[i]*fvec[i+1];
			else if( i == n-1)
				sum += ode2->b[i]*fvec[i] + ode2->a[i]*fvec[i-1];
			else
				sum += ode2->a[i]*fvec[i-1] + ode2->b[i]*fvec[i] + ode2->c[i]*fvec[i+1];
			g[i]=sum;
		}
		for (i=0;i<n;i++) xold[i]=x[i];
		fold=f;
		for (i=0;i<n;i++) p[i] = -fvec[i];
        tridag(ode2->a, ode2->b, ode2->c, p, p);
		lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin);
		test=0.0;
		for (i=0;i<n;i++)
			if (abs(fvec[i]) > test) test=abs(fvec[i]);
		if (test < TOLF) {
			check=false;
			return;
		}
		if (check) {
			test=0.0;
			den=MAX(f,0.5*n);
			for (i=0;i<n;i++) {
				temp=abs(g[i])*MAX(abs(x[i]),1.0)/den;
				if (temp > test) test=temp;
			}
			check=(test < TOLMIN);
			return;
		}
		test=0.0;
		for (i=0;i<n;i++) {
			temp=(abs(x[i]-xold[i]))/MAX(abs(x[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX)
			return;
	}
	throw("MAXITS exceeded in newt");
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

    Int N = 2;
    Doub x0 = 0, x1 = 2, alpha = 0, beta = 1;
    Doub A3 = NaN, A2 = NaN, A1 = NaN, error = NaN, error_1 = NaN, error_2 = NaN, atol = 1e-4;
    std::fstream file("lecture12.csv", ios::out);
    file << "N" << ", " << "A[N]" << ", " << "A[N-1] - A[N]" << ", " << "Estimate of Order" << ", " << "Error " << std::endl;

    VecDoub x(N-1);

    // Initial guess....
    for (size_t i = 0; i < x.size(); i++)
    {
        x[i] = alpha + i*(beta-alpha)/(x.size());
    }

    while( error > atol || N <= 8 )
    {
        
        two_point_boundary_problem ode2(x0, x1, alpha, beta);
        bool check;

        for (size_t i = 0; i < x.size(); i++)
        {
            x[i] = alpha + i*(beta-alpha)/(x.size());
        }

        newt(x, check, &ode2);

        A3 = A2;
        A2 = A1;
        A1 = x[x.size()/2.0];

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
        error = abs((error_1) / (4.0 - 1.0));

        file << N << ", " << A1 << ", " << error << ", " << error_2 << ", " << (error_1) / (4.0 - 1.0) << std::endl;

        VecDoub x_cpy = x;
        N *= 2;
        x.resize(N-1);
        
        Doub itp = (x_cpy[0] - alpha)/2.0;
        x[0] = itp;
        itp = (x_cpy[x_cpy.size()-1] - beta)/2.0;
        x[N-2] = itp;

        for (size_t i = 1; i < N-2; i++)
        {
            if ( (i % 2) == 1)
            {
                x[i] = x_cpy[i/2.0];
            }
            else
            {
                itp = x_cpy[(i+1)/2.0] - x[i-1];
                x[i] = itp;
            }
        }
    }
}


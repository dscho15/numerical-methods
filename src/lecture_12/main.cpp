#include <iostream>
#include <fstream>
#include <functional>
#include <nr3.h>
#include "ludcmp.h"
#include "tridag.h"
#include "qrdcmp.h"
#include <fstream>

constexpr Doub alpha_= 0.0;
constexpr Doub beta_ = 1.0;
constexpr Doub a_    = 0.0;
constexpr Doub b_    = 2.0;

VecDoub func(VecDoub_I &y) 
{
    VecDoub phi(y.size());
    Doub h = (b_-a_)/y.size();

    auto lambda = [](Doub x, Doub dy, Doub y)
    {
        return 2*x + sin( dy ) - cos(y);
    };

    for (size_t i = 0; i < y.size(); i++)
    {
        Doub x = (b_-a_)/(y.size()+1) * (i+1); 
        Doub dy = 0;
        if( i == 0)
        {
            dy      = (y[i+1] - alpha_)/(2*h);
            phi[i]  = y[i+1] - 2 * y[i] - std::pow(h,2) * lambda(x, dy, y[i]) + alpha_;
        }
        else if ( i == y.size() )
        {
            dy      = (beta_ - y[i-1])/(2*h);
            phi[i]  = -2*y[i] + y[i-1] - std::pow(h,2) * lambda(x, dy, y[i]) + beta_;
        }
        else
        {
            dy      = (y[i+1] - y[i-1])/(2*h);
            phi[i]  = y[i-1] - 2 * y[i] + y[i+1] - std::pow(h,2) * lambda(x, dy, y[i]);
        }

    }
    return phi;
}

template <class T>
void lnsrch(VecDoub_I &xold, const Doub fold, VecDoub_I &g, VecDoub_IO &p,
VecDoub_O &x, Doub &f, const Doub stpmax, Bool &check, T &func) 
{
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

template <class T>
struct NRfdjac {
	const Doub EPS;
	T &func;
	VecDoub a;
	VecDoub b;
	VecDoub c;
	NRfdjac(T &funcc) : EPS(1.0e-8), func(funcc) 
	{}
	void jacobian (VecDoub_I &y, VecDoub_I &fvec) 
	{
		Int n 		= y.size();
        Doub zero 	= 0.0;
		a.resize(n);
		b.resize(n);
		c.resize(n);
		
		//MatDoub df(n, n, zero);

        auto ydd_yd = [](Doub y_d, Doub y, Doub x)
        {
            return cos(y_d);
        };

        auto ydd_y = [](Doub y_d, Doub y, Doub x)
        {
            return sin(y);
        };

		for (Int i = 0; i < n; i++)
        {
            Doub h         = (b_-a_)/(n+1);
            Doub x         = a_+h*(i+1);
            Doub y_d       = 0;

			a[i] = 0; b[i] = 0; c[i] = 0;

            if(i == 0)
            {
                y_d         = (y[i+1]-alpha_)/(2.0*h);
				b[i]    = -(2.0 + std::pow(h, 2.0) * ydd_y(y_d, y[i], x));
                c[i]    = -(-1.0 + h / 2.0 * ydd_yd(y_d, y[i], x));
            }
            else if(i == n-1)
            {
                y_d         = (beta_ - y[i-1])/( 2.0 * h  );
				a[i]    = -(-1.0 - h / 2.0 * ydd_yd(y_d, y[i], x));
                b[i]    = -(2.0 + std::pow(h, 2.0) * ydd_y(y_d, y[i], x));
            }
            else
            {
                y_d         = (y[i+1] - y[i-1]) / (2.0 * h);
				a[i]   = -(-1.0 - h / 2.0 * ydd_yd(y_d, y[i], x));
                b[i]   = -(2.0 + std::pow(h, 2.0) * ydd_y(y_d, y[i],  x));  
                c[i]   = -(-1.0 + h / 2.0 * ydd_yd(y_d, y[i], x));
            }
        }
	}
};

template <class T>
struct NRfmin {
	VecDoub fvec;
	T &func;
	Int n;
	NRfmin(T &funcc) : func(funcc){}
	Doub operator() (VecDoub_I &x) {
		n=x.size();
		Doub sum=0;
		fvec=func(x);
		for (Int i=0;i<n;i++) sum += SQR(fvec[i]);
		return 0.5*sum;
	}
};

template <class T>
void newt(VecDoub_IO &x, Bool &check, T &vecfunc) 
{
	const Int MAXITS=200;
	const Doub TOLF=1.0e-8,TOLMIN=1.0e-12,STPMX=100.0;
	const Doub TOLX=numeric_limits<Doub>::epsilon();
	Int i,j,its,n=x.size();
	Doub den,f,fold,stpmax,sum,temp,test;
	VecDoub g(n),p(n),xold(n);
	NRfmin<T> fmin(vecfunc);
	NRfdjac<T> fdjac(vecfunc);
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
		fdjac.jacobian(x,fvec);

		for (i=0;i<n;i++) 
		{
			sum=0.0;

			if(i == 0)
				sum += fdjac.b[i]*fvec[i] + fdjac.c[i]*fvec[i+1];
			else if( i == n-1)
				sum += fdjac.b[i]*fvec[i] + fdjac.a[i]*fvec[i-1];
			else
				sum += fdjac.a[i]*fvec[i-1] + fdjac.b[i]*fvec[i] + fdjac.c[i]*fvec[i+1];

			g[i]=sum;
		}

		for (i=0;i<n;i++) xold[i]=x[i];
		fold=f;
		for (i=0;i<n;i++) p[i] = -fvec[i];

		tridag(fdjac.a, fdjac.b, fdjac.c, p, p);
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

Doub richardson_extrapolation(Doub N, Doub A_h3, Doub A_h2, Doub A_h1)
{
	auto order = round(log2((A_h1 - A_h2) / (A_h2 - A_h3)));
    std::cout << "N: " << setw(10) << N << " - A: " << setw(15) << A_h3 << " - A(last)-A(new): " << setw(15) 
    << A_h2 - A_h3 << " - Rich-alpha^k: " << setw(15)<< order << " - Error Estimate: " << setw(15) << ( A_h2 - A_h1 ) / (std::pow(2,order) - 1) << std::endl;
	return ( A_h2 - A_h1 ) / (std::pow(2,order) - 1);
};

void relax()
{
	// Init 1
    VecDoub_IO y(3);

    for (size_t i = 0; i < y.size(); i++) y[i] = alpha_ + (beta_ - alpha_)/(y.size()+1) * (i+1);

	// Init 2
	int n = 2;
	VecDoub y_save_1 = y;
	VecDoub y_save_2 = y;
	ofstream log;
	log.open("lec12_file.csv");
	Int least_nb_of_runs = 0;

	while(true)
	{
		// Run Newtons Method, should compare with y_new_1 after
		Bool chk = false;
    	newt(y, chk, func);

		// Compare with y_new_1
		auto k = richardson_extrapolation(y.size(), y[(y.size()-1)/2], y_save_1[(y_save_1.size()-1)/2], y_save_2[(y_save_2.size()-1)/2]);

		Doub h = (b_-a_)/(y.size()+1);
		for (size_t i = 0; i < y.size(); i++)
		{
			Doub x = a_ + h*(i+1);
			log << x << ", " << y[i] << "\n";
		}

		/*
		* M = readmatrix("lec12_file.csv")
	    * plot(M(:,1),M(:,2),'*')
		*/

		least_nb_of_runs = least_nb_of_runs + 1;

		// Be careful with tuning this error constant
		const Doub epsilon = 1e-4;
		if ( least_nb_of_runs++ > 2 && k < epsilon) break;

		y_save_2 = y_save_1;
		y_save_1 = y;

		// Interpolate for next Run
		VecDoub temp;
		
		// Setup for interpolation
		n *= 2;
		temp.resize( n + y.size() );

		for (size_t i = 0; i < y.size(); i++)
		{
			if ( i == 0)
			{
				temp[i		 ]  = 1.0 / 2.0 * ( y[i] + alpha_);
				temp[i+1	 ]  = y[i];
			}
			else if ( i == ( y.size() - 1 ) )
			{
				temp[2*i    ] 	= 1.0 / 2.0 * ( y[i] + beta_);
			}
			else
			{
				temp[2*i 	 ]  = 1.0 / 2.0 *( y[i-1] + y[i] );
				temp[2*i + 1]	= y[i];
			}
		}

		// Overwrite last value
		y = temp;
	}

	log.close();
}

Int main(Int argc, Char ** argv)
{

    VecDoub x(2);

    VecDoub y0(2);
    y0[0] = 0;
    y0[1] = 1;

    VecDoub x0(2);
    x0[0] = 0;
    x0[1] = 2;

    relax();

}


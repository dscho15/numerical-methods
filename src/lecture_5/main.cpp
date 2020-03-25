#include <iostream>
#include <cmath>
#include <nr3.h>
#include <roots.h>

template <class T>
double bisection(T &func, const double x1, const double x2, const double xacc) {
	const int JMAX=50;
	double dx,xmid,rtb;
	double f=func(x1);
	double fmid=func(x2);
	if (f*fmid >= 0.0) 
        throw("Root must be bracketed for bisection in bisection");
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (int j=0;j<JMAX;j++) 
	{
		fmid=func(xmid=rtb+(dx *= 0.5));
		if (fmid <= 0.0) 
            rtb=xmid;
		if (abs(dx) < xacc || fmid == 0.0) 
            return rtb;

		printf("%f \t - \t %f \t - \t %d \n", rtb, dx, j);
	}

	throw("Too many bisections in bisection");
}

template <class T>
double secant(T &func, const double x1, const double x2, const double xacc) 
{
	const int MAXIT=30;
	double xl,rts;
	double fl=func(x1);
	double f=func(x2);
	if (std::abs(fl) < std::abs(f)) {
		rts=x1;
		xl=x2;
		std::swap(fl,f);
	} 
    else 
    {
		xl=x1;
		rts=x2;
	}
	for (int j=0;j<MAXIT;j++) {
		double dx=(xl-rts)*f/(f-fl);
		xl=rts;
		fl=f;
		rts += dx;
		f=func(rts);

        printf("%f  - %f    - %d \n", xl, dx, j);

		if ( std::abs(dx) < xacc || f == 0.0) 
            return rts;
	}
	throw("Maximum number of iterations exceeded in secant");
}

template <class T>
double true_false(T &func, const double x1, const double x2, const double xacc) {
	const int MAXIT=30;
	double xl, xh, del, fl=func(x1), fh=func(x2);
	if (fl*fh > 0.0) 
        throw("Root must be bracketed in true_false");
	if (fl < 0.0) 
    {
		xl=x1;
		xh=x2;
	} 
    else 
    {
		xl=x2;
		xh=x1;
		std::swap(fl,fh);
	}
	double dx=xh-xl;
	for (int j=0;j<MAXIT;j++) 
    {
		double rtf=xl+dx*fl/(fl-fh);
		double f=func(rtf);
		if (f < 0.0) {
			del=xl-rtf;
			xl=rtf;
			fl=f;
		} 
        else 
        {
			del=xh-rtf;
			xh=rtf;
			fh=f;
		}
		dx=xh-xl;
        printf("%f  - %f    - %d \n", xl, dx, j);
		if (abs(del) < xacc || f == 0.0) 
            return rtf;
	}
	throw("Maximum number of iterations exceeded in true_false");
}

template <class T>
double ridder(T &func, const double x1, const double x2, const double xacc) 
{
	const int MAXIT=60;
	double fl=func(x1);
	double fh=func(x2);
	if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
		double xl=x1;
		double xh=x2;
		double ans=-9.99e99;
		for (int j=0;j<MAXIT;j++) {
			double xm=0.5*(xl+xh);
			double fm=func(xm);
			double s=sqrt(fm*fm-fl*fh);
            double dx = xh-xl;
            printf("%f  - %f    - %d \n", xl, dx, j);
			if (s == 0.0) return ans;
			double xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
			if (abs(xnew-ans) <= xacc) return ans;
			ans=xnew;
			double fnew=func(ans);
			if (fnew == 0.0) return ans;
			if (SIGN(fm,fnew) != fm) {
				xl=xm;
				fl=fm;
				xh=ans;
				fh=fnew;
			} else if (SIGN(fl,fnew) != fl) {
				xh=ans;
				fh=fnew;
			} else if (SIGN(fh,fnew) != fh) {
				xl=ans;
				fl=fnew;
			} else throw("never get here.");
			if (abs(xh-xl) <= xacc) return ans;
		}
		throw("ridder exceed maximum iterations");
	}
	else {
		if (fl == 0.0) return x1;
		if (fh == 0.0) return x2;
		throw("root must be bracketed in ridder.");
	}
}


template <class T, class T_dot>
double newton_raphson(T &func, T_dot &func_dot, const double xi, const double xacc)
{
    constexpr int max_it = 30;
    double x1, dx, x0 = xi;

    for (int i = 0; i < max_it; i++)
    {
        x1 = x0 - func(x0)/func_dot(x0);
        dx = x1-x0;
        printf("%f  - %f    - %d \n", x0, dx, i);
        if(abs(dx) < xacc) return x1;
        x0 = x1;
    }
    throw("newton_raphson exceed maximum it");
}

double func(double x)
{
    return x-cos(x);
}

double func_dot(double x)
{
    return 1+sin(x);
}

int main()
{
    // Question 1
    double xl = 0, xh = M_PI_2, xacc = 1e-8;

    printf("Start of Bisection := \n");
    double bisection_q1 = bisection(func, xl, xh, xacc);
    std::cout << bisection_q1 << "\n" << std::endl;

    // Question 2
                                xacc = 1e-16;
    printf("Start of secant := \n");
    double secant_q1 = secant(func, xl, xh, xacc);
    std::cout << secant_q1 << "\n" << std::endl;

    // Question 3
    printf("Start of true_false := \n");
    double true_false_q1 = true_false(func, xl, xh, xacc);
    std::cout << true_false_q1 << "\n" << std::endl;

    // Question 4
    printf("Start of ridder := \n");
    double ridder_q1 = ridder(func, xl, xh, xacc);
    std::cout << ridder_q1 << "\n" << std::endl;

    // Question 5
    printf("Newton Rawphson := \n");
    double newton_q1 = newton_raphson(func, func_dot, xl, xacc);
    std::cout << newton_q1 << "\n" << std::endl;




    return 1;
}
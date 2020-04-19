template<class T>
struct DErule : Quadrature {
	Ldoub a,b,hmax,s;
	T &func;

	DErule(T &funcc, const Ldoub aa, const Ldoub bb, const Ldoub hmaxx=3.7)
		: func(funcc), a(aa), b(bb), hmax(hmaxx) {n=0;}

	Ldoub next() {
		Ldoub del,fact,q,sum,t,twoh;
		Int it,j;
		n++;
		if (n == 1) {
			fact=0.25;
			return s=hmax*2.0*(b-a)*fact*func(0.5*(b+a),0.5*(b-a));
		} else {
			for (it=1,j=1;j<n-1;j++) it <<= 1;
			twoh=hmax/it;
			t=0.5*twoh;
			for (sum=0.0,j=0;j<it;j++) {
				q=exp(-2.0*sinh(t));
				del=(b-a)*q/(1.0+q);
				fact=q/SQR(1.0+q)*cosh(t);
				sum += fact*(func(a+del,del)+func(b-del,del));
				t += twoh;
			}
			return s=0.5*s+(b-a)*twoh*sum;
		}
	}
};

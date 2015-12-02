#include "nr_headers/nr3.h"
#include "nr_headers/erf.h"
#include "nr_headers/ludcmp.h"
#include "nr_headers/qrdcmp.h"
#include "nr_headers/roots_multidim.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>

#include "headers/write_matrix.h"

using namespace::std;

struct Funcvec{
	double jeo,jee,jie,jio,jei,jii;
	double theta_e,theta_i;
	int K;
	double mo;
	double sqK;
	double E,I,Je,Ji;
	Funcvec(const double moo,const double jeoo, const double jeee, const double jiee,
			const double jioo, const double jeii, const double jiii,
			const double theta_ee,const double theta_ii,int KK) :
			mo(moo),jeo(jeoo),jee(jeee),jie(jiee),jio(jioo),jei(jeii),
			jii(jiii), theta_e(theta_ee), theta_i(theta_ii),K(KK),
			sqK(sqrt((double)K)), E(jeo),I(jio),Je(-1*jei),Ji(-1*jii) {}
	Erf erfunction;
	VecDoub operator()(const VecDoub_I& x) {
		VecDoub f(2);
		f[0] = E*mo+x[0]-Je*x[1]-(theta_e-sqrt(x[0]+Je*Je*x[1])*erfunction.inverfc(x[0]))/sqK;
		f[1] = I*mo+x[0]-Ji*x[1]-(theta_i-sqrt(x[0]+Ji*Ji*x[1])*erfunction.inverfc(x[1]))/sqK;
		return f;
	}
};
//VecDoub funcv(const VecDoub_I& x)
//{
//
//	Erf erfunc;
//	VecDoub f(2);
//	double me = x[0];
//	double mi = x[1];
//	f[0] = E*mo+me-Je*mi-(theta_e-sqrt(me+Je*Je*mi)*erfunc.inverfc(me))/sqK;
//	f[1] = I*mo+me-Ji*mi-(theta_i-sqrt(me+Ji*Ji*mi)*erfunc.inverfc(mi))/sqK;
//	return f;
//}


int main()
{
	double jeo = 1.0;
	double jee = 1.0;
	double jie = 1.0;
	double jio = 0.8;
	double jei = -2.0;
	double jii = -1.8;

	double theta_e = 1.0;
	double theta_i = 0.7;

	double K = 1000;

	int N=2;
	bool check;
	VecDoub_IO x(N);
	VecDoub f(N);

	int n = 1000;
	vector<double> m(n);
	for(int i=0;i<n;i++) m[i] = 0.1+0.27*i/(double)n;
	vector<double> me(n), mi(n);
	for(int i=0;i<n;i++) {

		Funcvec	funcv(m[i],jeo,jee,jie,jio,jei,jii,theta_e,theta_i,K);

		x[0] = 1.5*m[i];
		x[1] = 1.5*m[i];
		broydn(x,check,funcv);
//		newt(x,check,funcv);
		f = funcv(x);
		if(check) cout << " shit's fucked up at i="<<i << endl;
	
		me[i] = x[0];
		mi[i] = x[1];
	}

	write_matrix(me,n,"me.csv");
	write_matrix(mi,n,"mi.csv");
	write_matrix(m,n,"mo.csv");


	return 0;
}






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
	Funcvec(double moo,double jeoo,double jeee,double jiee,
			double jioo,double jeii,double jiii, double theta_ee,
			double theta_ii,int KK) :
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

struct Funcvec0 {
	double jeo, jee,jie,jio,jei,jii;
	double mo;
	double E,I,Je,Ji;
	Funcvec0(double moo, double jeoo, double jeee,double jiee,
			double  jioo, double jeii, double jiii) :
			mo(moo), jeo(jeoo), jee(jeee),jie(jiee),
			jio(jioo),jei(jeii),jii(jiii), E(jeo),
			I(jio), Je(-1.*jei), Ji(-1.*jii) {}
	VecDoub operator()(const VecDoub_I& x) {
		VecDoub f(2);
		f[0] = E*mo+x[0]-Je*x[1];
		f[1] = I*mo+x[0]-Ji*x[1];
		return f;
	}
};


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
	for(int i=0;i<n;i++) m[i] = 0.04+0.26*i/(double)n;
	vector<double> me(n), mi(n);
	for(int i=0;i<n;i++) {

		Funcvec	funcv(m[i],jeo,jee,jie,jio,jei,jii,theta_e,theta_i,K);

		x[0] = 1.5*m[i];
		x[1] = 1.5*m[i];
		broydn(x,check,funcv);
		f = funcv(x);
		if(check) cout << " shit's fucked up at i="<<i << endl;
	
		me[i] = x[0];
		mi[i] = x[1];
	}

	write_matrix(me,n,"me.csv");
	write_matrix(mi,n,"mi.csv");
	write_matrix(m,n,"mo.csv");

	VecDoub_IO xx(N);
	VecDoub ff(N);
	vector<double> mm(n);
	for(int i=0;i<n;i++) mm[i] = 0.3*i/(double)n;
	vector<double> mme(n), mmi(n);
	for(int i=0;i<n;i++) {
		Funcvec0 funcv0(mm[i],jeo,jee,jie,jio,jei,jii);
		xx[0] = 1.5*mm[i];
		xx[1] = xx[0];
		broydn(xx,check,funcv0);
		if(check) cout << " error at i= " << i << endl;
		mme[i] = xx[0];
		mmi[i] = xx[1];
	}

	write_matrix(mme,n,"mme.csv");
	write_matrix(mmi,n,"mmi.csv");
	write_matrix(mm,n,"mmo.csv");


	return 0;
}






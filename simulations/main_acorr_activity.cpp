/*
	calculates average autocorr 
	of the activity of units

*/



// nr_headers (nr_headers/)
#include "nr_headers/nr3.h"
#include "nr_headers/ran.h"
#include "nr_headers/fourier.h"
#include "nr_headers/correl.h"

// std headers
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <time.h>

//other headers (headers/)
#include "headers/generate_matrix.h"
#include "headers/write_matrix.h"
#include "headers/other.h"
#include "headers/check_empty.h"

using namespace::std;

typedef std::vector<std::vector<int> > Matint;
int main(int argc,char *argv[])
{
	int seed;
	seed = 11111111;


	int N = 10000;
	int Ne = N;
	int Ni = N;
	int No = N;

	int K = 1000;

	double Jeo = 1.0;
	double Jee = 2.0;
	double Jie = 1.0;
	double Jio = 0.8;
	double Jei = -2.0;
	double Jii = -1.8;

	double theta_e = 1.0;
	double theta_i = 0.8;
	double D = 0.3;
	double mo = 0.08;
	double tau = .9;

	double me0 = 0;
	double mi0 = 0;

	int tmax = 100*N;

	Ranq1 r(seed);

	Matint EE = connectivity_matrix(Ne,Ne,K,r);
	Matint EI = connectivity_matrix(Ne,Ni,K,r);
	Matint IE = connectivity_matrix(Ni,Ne,K,r);
	Matint II = connectivity_matrix(Ni,Ni,K,r);
	Matint EO = connectivity_matrix(Ne,No,K,r);
	Matint IO = connectivity_matrix(Ni,No,K,r);

	vector<double> the(Ne,0.0);
	vector<double> thi(Ni,0.0);
	double sqK = sqrt((double)K);
	for(int i=0;i<Ne;i++) the[i] = (theta_e + r.doub()*D)*sqK;
	for(int i=0;i<Ni;i++) thi[i] = (theta_i + r.doub()*D)*sqK;


	vector<int> nwe(Ne,0);
	vector<int> nwi(Ni,0);
	vector<int> nwo(No,0);
//	for(int i=0;i<Ne;i++) if(me0>r.doub()) nwe[i] = 1;
//	for(int i=0;i<Ni;i++) if(mi0>r.doub()) nwi[i] = 1;
	for(int i=0;i<No;i++) if(mo>r.doub()) nwo[i] = 1;

//	vector<double> node_e_in(tmax,0.0);
//	vector<double> node_i_in(tmax,0.0);
//	vector<double> node_spike(tmax,0.0);
	
	int tmax2 = pow2(tmax);
	VecDoub state(tmax2,0.0);
	VecDoub acorrnode(tmax2,0.0);
	vector<double> acorr(tmax2/2,0.0);
	
	// start simulation
	int ie, ii, nwe_ie,nwi_ii;
	double change_nwe, change_nwi;
	double current,currentEE, currentEI, currentIE,
		currentII, currentEO, currentIO;

	int nfull =0;
	int naverage = 50;
	vector<int> empty_nwe(Ne,0.0);
//	vector<int> empty_nwi(Ni,0.0);
	for(int inode=0;inode<naverage;inode++) {
			nwe = empty_nwe;
			nwi = empty_nwe;						
			for(int t=1;t<10*N;t++) {

			//update exitatory population
			ie = r.int64() % (Ne-1);
			currentEE = Jee*dotproduct(EE[ie],nwe,Ne);
			currentEI = Jei*dotproduct(EI[ie],nwi,Ni);
	//		currentEO = Jeo*mo*K;
			currentEO = Jeo*dotproduct(EO[ie],nwo,No);
			current = currentEE+currentEO+currentEI;
			if(current>the[ie]){
				nwe_ie = 1;
				if(nwe[ie] == 1){
					 change_nwe = 0;
				} else {
					change_nwe = 1;
					nwe[ie] = 1;
//					if(ie == inode) node_spike[t] = 1;
				}
			} else{
				nwe_ie = 0;
				if(nwe[ie] == 0) {
					change_nwe = 0;
				} else {
					change_nwe = -1;
					nwe[ie] = 0;
				}
			}
			
			// update inhibitory population
			ii = r.int64() % (Ni-1);
			currentIE = Jie*dotproduct(IE[ii],nwe,Ne);
			currentII = Jii*dotproduct(II[ii],nwi,Ni);
	//		currentIO = Jio*mo*K;
			currentIO = Jio*dotproduct(IO[ii],nwo,No);
			current = currentIE+currentII+currentIO;
			if(current>thi[ii]){
				nwi_ii = 1;
				if(nwi[ii] == 1) {
					change_nwi =0;
				} else {
					change_nwi = 1;
					nwi[ii] = 1;
				}	
			} else {
				nwi_ii = 0;
				if(nwi[ii] == 0) {
					change_nwi = 0;
				} else {
					change_nwi = -1;
					nwi[ii] = 0;
				}
			}

			// extra inhibitory update
			if( (1/tau - 1)> r.doub()) {
				ii = r.int64() % (Ni-1);
				currentIE = Jie*dotproduct(IE[ii],nwe,Ne);
				currentII = Jii*dotproduct(II[ii],nwi,Ni);
	//			currentIO = Jio*mo*K;
				currentIO = Jio*dotproduct(IO[ii],nwi,No);
				current = currentIE+currentII+currentIO;
				if(current>thi[ii]){
					nwi_ii = 1;
					if(nwi[ii] == 1) {
						change_nwi =0;
					} else {
						change_nwi = 1;
						nwi[ii] = 1;
					}	
				} else {
					nwi_ii = 0;
					if(nwi[ii] == 0) {
						change_nwi = 0;
					} else {
						change_nwi = -1;
						nwi[ii] = 0;
					}
				}
			}
		}
		for(int t=1;t<tmax;t++) {

			//update exitatory population
			ie = r.int64() % (Ne-1);
			currentEE = Jee*dotproduct(EE[ie],nwe,Ne);
			currentEI = Jei*dotproduct(EI[ie],nwi,Ni);
	//		currentEO = Jeo*mo*K;
			currentEO = Jeo*dotproduct(EO[ie],nwo,No);
			current = currentEE+currentEO+currentEI;
			if(current>the[ie]){
				nwe_ie = 1;
				if(nwe[ie] == 1){
					 change_nwe = 0;
				} else {
					change_nwe = 1;
					nwe[ie] = 1;
//					if(ie == inode) node_spike[t] = 1;
				}
			} else{
				nwe_ie = 0;
				if(nwe[ie] == 0) {
					change_nwe = 0;
				} else {
					change_nwe = -1;
					nwe[ie] = 0;
				}
			}
			
			// update inhibitory population
			ii = r.int64() % (Ni-1);
			currentIE = Jie*dotproduct(IE[ii],nwe,Ne);
			currentII = Jii*dotproduct(II[ii],nwi,Ni);
	//		currentIO = Jio*mo*K;
			currentIO = Jio*dotproduct(IO[ii],nwo,No);
			current = currentIE+currentII+currentIO;
			if(current>thi[ii]){
				nwi_ii = 1;
				if(nwi[ii] == 1) {
					change_nwi =0;
				} else {
					change_nwi = 1;
					nwi[ii] = 1;
				}	
			} else {
				nwi_ii = 0;
				if(nwi[ii] == 0) {
					change_nwi = 0;
				} else {
					change_nwi = -1;
					nwi[ii] = 0;
				}
			}

			// extra inhibitory update
			if( (1/tau - 1)> r.doub()) {
				ii = r.int64() % (Ni-1);
				currentIE = Jie*dotproduct(IE[ii],nwe,Ne);
				currentII = Jii*dotproduct(II[ii],nwi,Ni);
	//			currentIO = Jio*mo*K;
				currentIO = Jio*dotproduct(IO[ii],nwi,No);
				current = currentIE+currentII+currentIO;
				if(current>thi[ii]){
					nwi_ii = 1;
					if(nwi[ii] == 1) {
						change_nwi =0;
					} else {
						change_nwi = 1;
						nwi[ii] = 1;
					}	
				} else {
					nwi_ii = 0;
					if(nwi[ii] == 0) {
						change_nwi = 0;
					} else {
						change_nwi = -1;
						nwi[ii] = 0;
					}
				}
			}

			state[t] = nwe[inode];
		}
		bool empty = true;
		empty = check_empty(state,tmax2);
		if(!empty) {
			nfull++;
			correl(state,state,acorrnode);
			double acorrnode0 = acorrnode[0];
			for(int i=0;i<tmax2/2;i++) {
				acorr[i] += acorrnode[i]/acorrnode0;
				acorr[i] += acorrnode[i];
			}
		}

	}
	cout << nfull << endl;
	for(int i=0;i<tmax2/2;i++) acorr[i]/=(double)nfull;
	write_matrix(acorr,tmax2/2,"acorr.csv");
	
	return 0;
}


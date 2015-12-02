// nr_headers (nr_headers/)
#include "nr_headers/nr3.h"
#include "nr_headers/ran.h"


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

using namespace::std;

typedef std::vector<std::vector<int> > Matint;
int main()
{
	time_t tstart = time(NULL);

	int inode = 0;
	ofstream inode_out("inode.txt");
	inode_out << inode;

	int N = 5000;
	int Ne = N;
	int Ni = N;
	int No = N;

	int K = 500;

	double Jeo = 1.0;
	double Jee = 1.0;
	double Jie = 1.0;
	double Jio = 0.8;
	double Jei = -2.5;
	double Jii = -1.8;

	double theta_e = 1.0;
	double theta_i = 0.8;
	double D = 0.;
	double mo = 0.08;
	double tau = .9;

	double me0 = 0;
	double mi0 = 0;

	int tmax = 200*N;

	int seed = 123456789;
	Ranq1 r(seed);

	ofstream log("log.txt");
	log << endl;
	time_t t1,t2;
	t1 = time(NULL);

	Matint EE = connectivity_matrix(Ne,Ne,K,r);
	Matint EI = connectivity_matrix(Ne,Ni,K,r);
	Matint IE = connectivity_matrix(Ni,Ne,K,r);
	Matint II = connectivity_matrix(Ni,Ni,K,r);
//	Matint EO = connectivity_matrix(Ne,No,K,r);
//	Matint IO = connectivity_matrix(Ni,No,K,r);

	t2 = time(NULL);
	log << "generating connectivity matrices took: " << difftime(t2,t1) << " sec" << endl;	

	vector<double> the(Ne,0.0);
	vector<double> thi(Ni,0.0);
	double sqK = sqrt((double)K);
	for(int i=0;i<Ne;i++) the[i] = (theta_e + r.doub()*D)*sqK;
	for(int i=0;i<Ni;i++) thi[i] = (theta_i + r.doub()*D)*sqK;


	vector<int> nwe(Ne,0);
	vector<int> nwi(Ni,0);
	vector<int> nwo(No,0);
	for(int i=0;i<Ne;i++) if(me0>r.doub()) nwe[i] = 1;
	for(int i=0;i<Ni;i++) if(mi0>r.doub()) nwi[i] = 1;
	for(int i=0;i<No;i++) if(mo>r.doub()) nwo[i] = 1;

	vector<double> nwe_activity(tmax,0.0);
	vector<double> nwi_activity(tmax,0.0);
	vector<double> node_e_in(tmax,0.0);
	vector<double> node_i_in(tmax,0.0);
	vector<double> node_spike(tmax,0.0);
	
	t1 = time(NULL);	
	// start simulation
	int ie, ii, nwe_ie,nwi_ii;
	double change_nwe, change_nwi;
	double current,currentEE, currentEI, currentIE,
		currentII, currentEO, currentIO;

	nwe_activity[0] = average_vec(nwe,Ne);
	nwi_activity[0] = average_vec(nwi,Ni);
	for(int t=1;t<tmax;t++) {

		//update exitatory population
		ie = r.int64() % (Ne-1);
		currentEE = Jee*dotproduct(EE[ie],nwe,Ne);
		currentEI = Jei*dotproduct(EI[ie],nwi,Ni);
		currentEO = Jeo*mo*K;
		current = currentEE+currentEO+currentEI;
		if(current>the[ie]){
			nwe_ie = 1;
			if(nwe[ie] == 1){
				 change_nwe = 0;
			} else {
				change_nwe = 1;
				nwe[ie] = 1;
				if(ie == inode) node_spike[t] = 1;
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
		nwe_activity[t] = nwe_activity[t-1]+change_nwe/Ne;
		
		// update inhibitory population
		ii = r.int64() % (Ni-1);
		currentIE = Jie*dotproduct(IE[ii],nwe,Ne);
	
		if(t==1000) cout << currentIE << endl;
		currentII = Jii*dotproduct(II[ii],nwi,Ni);
		currentIO = Jio*mo*K;
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
		nwi_activity[t] = nwi_activity[t-1]+change_nwi/Ni;

		// extra inhibitory update
		if( (1/tau - 1)> r.doub()) {
			ii = r.int64() % (Ni-1);
			currentIE = Jie*dotproduct(IE[ii],nwe,Ne);
			currentII = Jii*dotproduct(II[ii],nwi,Ni);
			currentIO = Jio*mo*K;
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
			nwi_activity[t] = nwi_activity[t]+change_nwi/Ni;
		}


		// calculate inode currents
		currentEE = Jee*dotproduct(EE[inode],nwe,Ne);
		currentEI = Jei*dotproduct(EI[inode],nwi,Ni);
		currentEO = Jeo*mo*K;
		current = currentEE+currentEO+currentEI;
		node_e_in[t] = currentEE+currentEO;
		node_i_in[t] = currentEI;
	}

	t2 = time(NULL);
	log<< "simulation took: " << difftime(t2,t1) << " sec" << endl;

	t1 = time(NULL);
	write_matrix(nwe_activity,tmax,"Eactivity.csv");
	write_matrix(nwi_activity,tmax,"Iactivity.csv");
	write_matrix(node_e_in,tmax,"Ein.csv");
	write_matrix(node_i_in,tmax,"Iin.csv");
	write_matrix(node_spike,tmax,"spike.csv");

//	write_matrix(EE,Ne,Ne,"EE.csv");
//	write_matrix(EI,Ne,Ni,"EI.csv");
//	write_matrix(IE,Ni,Ne,"IE.csv");
//	write_matrix(II,Ni,Ni,"II.csv");
//	write_matrix(EO,Ne,No,"EO.csv");
//	write_matrix(IO,Ni,No,"IO.csv");
	write_matrix(the,Ne,"the.csv");
	write_matrix(thi,Ni,"thi.csv");

	log << "writing matrices took: " << difftime(time(NULL),t1) << " sec" << endl;
	log << endl << "total time: " << difftime(time(NULL),tstart) << " sec" << endl;
	log << endl;
	return 0;
}


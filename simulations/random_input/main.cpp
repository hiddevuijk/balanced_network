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
#include "headers/network_struct.h"
#include "headers/advance_network.h"
#include "headers/generate_matrix.h"
#include "headers/write_matrix.h"
#include "headers/read_input.h"
#include "headers/other.h"

using namespace::std;

int main(int argc,char *argv[])
{
	int N,K,tt,seed;
	double theta_e,theta_i,D,mo,mf;
	Network NW;
	
	int ti,tf;
	read_input(NW,N,K,theta_e,theta_i,D,mo,tt,seed,"input.txt");

	// if seed is input from command line use that:
	if(argc ==2) {
		seed = stoi(argv[1]);
	}
	string outname;
	if(argc == 2) {
		outname = argv[1];
		outname += "h.csv";
	} else {
		outname  = "nwe_act_acorr.csv";
	}

	int tmax = 	tt*N;
	int tmax2 = 2*tmax;
	tmax2 = pow2(tmax2);
	int t_to_stable = 10*N;

	Ranq1 r(seed);

	// setup network
	NW.Ne = N;
	NW.Ni = N;
	NW.No = N;

	NW.EE = connectivity_matrix(NW.Ne,NW.Ne,K,r);
	NW.EI = connectivity_matrix(NW.Ne,NW.Ni,K,r);
	NW.IE = connectivity_matrix(NW.Ni,NW.Ne,K,r);
	NW.II = connectivity_matrix(NW.Ni,NW.Ni,K,r);
	NW.EO = connectivity_matrix(NW.Ne,NW.No,K,r);
	NW.IO = connectivity_matrix(NW.Ni,NW.No,K,r);

	NW.the = vector<double> (NW.Ne,0.0);
	NW.thi = vector<double> (NW.Ni,0.0);
	double sqK = sqrt((double)K);
	for(int i=0;i<NW.Ne;i++) NW.the[i] = (theta_e + r.doub()*D)*sqK;
	for(int i=0;i<NW.Ni;i++) NW.thi[i] = (theta_i + r.doub()*D)*sqK;

	int inode =0;

	// set up starting state of the network
	State ST;
	ST.nwe = vector<int> (NW.Ne,0);
	ST.nwi = vector<int> (NW.Ni,0);
	ST.nwo = vector<int> (NW.No,0);
	for(int i=0;i<NW.No;++i){
		if(r.doub()<mo) ST.nwo[i]=1;
	}
	ST.node_spike = vector<int> (tmax,0);
	ST.nwe_activity = vector<double> (tmax,0.0);
	ST.nwi_activity = vector<double> (tmax,0.0);

	ST.inode = 0;
	
	// start simulation
	//advance until t_to_stable
	for(int i=0;i<t_to_stable*10;++i) {
		advanceNW_stable_state(NW,ST,r);
		for(int i=0;i<NW.Ne;++i){
			if(ST.nwe[i] ==1) NW.the[i] +=0.01;
			if(ST.nwe[i]==0 && NW.the[i] > 0.1) NW.the[i] -=0.1;
		}	
	}

	// activities at t=0
	ST.nwe_activity[0] = average_vec(ST.nwe,NW.Ne);
	ST.nwi_activity[0] = average_vec(ST.nwi,NW.Ni);
	// evolve and save values
	for(int t=1;t<tmax;t++) {
		// advance network 1 time step
		advanceNW(NW,ST,t,r);	

	}


	VecDoub nwe_act(tmax2,0.0);
	VecDoub nwe_act_acorr(tmax2,0.0);
	double m = mean(ST.nwe_activity,tmax);
	for(int i=0;i<tmax;++i) {
		nwe_act[i] = ST.nwe_activity[i] - m;
	}
	correl(nwe_act,nwe_act,nwe_act_acorr);
	
	// write results
	write_matrix(nwe_act_acorr,tmax2/2,outname);
	write_matrix(ST.nwe_activity,tmax,"nwe_act.csv");
	return 0;
}


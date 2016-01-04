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
#include "headers/network_struct.h"
#include "headers/advance_network.h"
#include "headers/generate_matrix.h"
#include "headers/write_matrix.h"
#include "headers/read_input.h"

using namespace::std;

int main(int argc,char *argv[])
{
	int N,K,tt,seed;
	double theta_e,theta_i,D,mo,mf;
	Network NW;
	int sti_max;
	int ti,tf;
	read_input(NW,N,K,theta_e,theta_i,D,mo,tt,seed,sti_max,"input.txt");

	// if seed is input from command line use that:
	if(argc ==2) {
		seed = stoi(argv[1]);
	}

	int tmax = 	tt*N;
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

	// set up starting state of the network
	State ST;
	ST.inode = 0;
	ST.current_sti = 0;
	ST.total_sti = 0;
	ST.nwe = vector<int> (NW.Ne,0);
	ST.nwi = vector<int> (NW.Ni,0);
	ST.nwo = vector<int> (NW.No,0);

	for(int i=0;i<NW.No;++i){
		if(r.doub()<mo) ST.nwo[i]=1;
	}

	// start simulation
	//advance until t_to_stable
	for(int i=0;i<t_to_stable;++i) {
		advanceNW_stable_state(NW,ST,r);
	}

	// evolve and save values
	int t=0;
	bool abort = false;
	while(ST.total_sti<sti_max){
		// advance network 1 time step
		advanceNW(NW,ST,t,r);
		++t;
		if(t/N == 100 and ST.total_sti <3) {
			cerr << seed << " aborted, les than 3 spikes after t=100 \n";
			abort = true;
			break;
		}
	}

	// write results
	if(!abort) {
		string name = "sti";
		if(argc>1) name = name + to_string(seed);
		name = name + ".csv";
		ofstream outstream(name);
		outstream << setprecision(16);
		for(int i=0;i<ST.total_sti;++i)
			outstream << ST.sti[i]/(double)N << ';';
		cout << name <<'\n';
	}
	return 0;
}


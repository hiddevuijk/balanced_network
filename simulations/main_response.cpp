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
#include "headers/read_input_perturbation.h"

using namespace::std;

int main(int argc,char *argv[])
{
	int N,K,tt,seed;
	double theta_e,theta_i,D,mo,mf;
	Network NW;
	
	int ti,tf;
	read_input(NW,N,K,theta_e,theta_i,D,mo,mf,ti,tf,tt,seed,"input.txt");
	double a = mf-mo;
	ti = ti*N;
	tf = tf*N;

	// if seed is input from command line use that:
	if(argc == 1) {
		seed = 123456789;
	} else if(argc ==2) {
		seed = stoi(argv[1]);
	} else cerr << argv[0] << " does not know what to do with input";


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

	int inode =0;

	// set up starting state of the network
	State ST;
	ST.nwe = vector<int> (NW.Ne,0);
	ST.nwi = vector<int> (NW.Ni,0);
	ST.nwo = vector<int> (NW.No,0);

	ST.node_spike = vector<int> (tmax,0);
	ST.nwe_activity = vector<double> (tmax,0.0);
	ST.nwi_activity = vector<double> (tmax,0.0);

	ST.inode = 0;
	
	bool node_1 = true;
	int i;
	// start simulation

	//advance until t_to_stable
	for(int i=0;i<t_to_stable;++i) {
		advanceNW_stable_state(NW,ST,r);
	}

	// activities at t=0
	ST.nwe_activity[0] = average_vec(ST.nwe,NW.Ne);
	ST.nwi_activity[0] = average_vec(ST.nwi,NW.Ni);

	// evolve and save values
	for(int t=1;t<tmax;t++) {

	
		// advance network 1 time step
		advanceNW(NW,ST,t,r);	
		
		if(t>ti && t<tf && r.doub()<a) {
			int i;
			do{
				i = r.int64() % (NW.No-1);
				if(ST.nwo[i] == 1) {
					node_1 = true;
				} else {
					node_1 = false;
				}
			} while(node_1);
			ST.nwo[i] =1;
		}
	}



	// write results
	write_matrix(ST.nwe_activity,tmax,"Eactivity.csv");
	write_matrix(ST.nwi_activity,tmax,"Iactivity.csv");

	return 0;
}


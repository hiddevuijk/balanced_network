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

using namespace::std;

int main(int argc,char *argv[])
{
	int seed;
	if(argc == 1) {
		seed = 123456789;
	} else if(argc ==2) {
		seed = stoi(argv[1]);
	} else cerr << argv[0] << " does not know what to do with input";


	int inode = 0;
	ofstream inode_out("inode.txt");
	inode_out << inode;

	int N = 1000;
	int K = 10;
	double theta_e = 1.0;
	double theta_i = 0.8;
	double D = 0.3;
	double mo = 0.08;

	int tmax = 10*N;
	int t_to_stable = 10*N;

	Ranq1 r(seed);
	// setup network
	Network NW;
	NW.Ne = N;
	NW.Ni = N;
	NW.No = N;

	NW.Jeo = 1.0;
	NW.Jee = 2.0;
	NW.Jie = 2.0;
	NW.Jio = 0.8;
	NW.Jei = -2.0;
	NW.Jii = -1.8;

	NW.tau = .9;

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
	ST.nwe = vector<int> (NW.Ne,0);
	ST.nwi = vector<int> (NW.Ni,0);
	ST.nwo = vector<int> (NW.No,0);

	ST.node_spike = vector<int> (tmax,0);
	ST.nwe_activity = vector<double> (tmax,0.0);
	ST.nwi_activity = vector<double> (tmax,0.0);

	ST.inode = 0;

	vector<double> node_e_in(tmax,0.0);
	vector<double> node_i_in(tmax,0.0);
	


	// start simulation
	for(int i=0;i<t_to_stable;++i) {
		advanceNW_stable_state(NW,ST,r);
	}


	ST.nwe_activity[0] = average_vec(ST.nwe,NW.Ne);
	ST.nwi_activity[0] = average_vec(ST.nwi,NW.Ni);

	for(int t=1;t<tmax;t++) {

		// advance network 1 time step
		advanceNW(NW,ST,t,r);	

		// calculate inode currents
		node_e_in[t] = NW.Jee*dotproduct(NW.EE[ST.inode],ST.nwe,NW.Ne);
		node_e_in[t] +=  NW.Jeo*dotproduct(NW.EO[ST.inode],ST.nwo,NW.No);
		node_i_in[t] = NW.Jei*dotproduct(NW.EI[ST.inode],ST.nwi,NW.Ni);
	}


	write_matrix(ST.nwe_activity,tmax,"Eactivity.csv");
	write_matrix(ST.nwi_activity,tmax,"Iactivity.csv");
	write_matrix(node_e_in,tmax,"Ein.csv");
	write_matrix(node_i_in,tmax,"Iin.csv");
	write_matrix(ST.node_spike,tmax,"spike.csv");

	write_matrix(NW.the,NW.Ne,"the.csv");
	write_matrix(NW.thi,NW.Ni,"thi.csv");

	return 0;
}


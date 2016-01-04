#ifndef GUARD_network_struct_h
#define GUARD_network_struct_h

#include <vector>

struct Network {
	std::vector<std::vector<int> > EE, EI, IE, II,EO,IO;
	double Jee, Jei, Jie, Jii, Jeo, Jio;
	int Ne,Ni,No;
	std::vector<double> the, thi;
	double tau;

};

struct State {
	std::vector<int>  nwe,nwi,nwo;
	std::vector<double> nwe_activity, nwi_activity;
	int inode;
	std::vector<int> node_spike;
};

#endif


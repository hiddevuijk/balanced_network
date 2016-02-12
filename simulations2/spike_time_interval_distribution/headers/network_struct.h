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
	int inode;
	std::vector<int> nwe,nwi,nwo;
	int current_sti;
	int total_sti;
	std::vector<int> sti;
};

#endif


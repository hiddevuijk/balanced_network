#ifndef GUARD_advance_network_h
#define GUARD_advance_network_h

#include "network_struct.h"
#include "generate_matrix.h"

#include <vector>

template<class R>
int find(const std::vector<int>& v,int N, int a, R r)
{
	int max = 1000;
	int m=0;
	int ii;
	do{
		ii = r.int64() %(N-1);
	} while(v[ii] != a && m++<maxx);
	return ii;
}



// ONLY WORKS FOR Ne=Ni !!!!!
template<class R>
void advanceNW(const Network& NW, State& ST, int t, R& r)
{

	
	int i, nwe_ie,nwi_ii;
	double change_nwe, change_nwi;
	double current,currentEE, currentEI, currentIE,
		currentII, currentEO, currentIO;

	// update excitatory node ie
	i = r.int64() % (NW.Ne-1);
	currentEE = NW.Jee*dotproduct(NW.EE[i],ST.nwe,NW.Ne);
	currentEI = NW.Jei*dotproduct(NW.EI[i],ST.nwi,NW.Ni);
	currentEO = NW.Jeo*dotproduct(NW.EO[i],ST.nwo,NW.No);
	current = currentEE+currentEO+currentEI;
	if(current>NW.the[i]){
		nwe_ie = 1;
		if(ST.nwe[i] == 1){
			 change_nwe = 0;
		} else {
			change_nwe = 1;
			ST.nwe[i] = 1;
			if(i == ST.inode) ST.node_spike[t] = 1;
		}
	} else{
		nwe_ie = 0;
		if(ST.nwe[i] == 0) {
			change_nwe = 0;
		} else {
			change_nwe = -1;
			ST.nwe[i] = 0;
		}
	}
	ST.nwe_activity[t] = ST.nwe_activity[t-1]+change_nwe/NW.Ne;
		
	if(NW.tau<1.) {
		// update inhibitory population
		i = r.int64() % (NW.Ni-1);
		currentIE = NW.Jie*dotproduct(NW.IE[i],ST.nwe,NW.Ne);
		currentII = NW.Jii*dotproduct(NW.II[i],ST.nwi,NW.Ni);
		currentIO = NW.Jio*dotproduct(NW.IO[i],ST.nwo,NW.No);
		current = currentIE+currentII+currentIO;
		if(current>NW.thi[i]){
			nwi_ii = 1;
			if(ST.nwi[i] == 1) {
				change_nwi =0;
			} else {
				change_nwi = 1;
				ST.nwi[i] = 1;
			}	
		} else {
			nwi_ii = 0;
			if(ST.nwi[i] == 0) {
				change_nwi = 0;
			} else {
				change_nwi = -1;
				ST.nwi[i] = 0;
			}
		}
		ST.nwi_activity[t] = ST.nwi_activity[t-1]+change_nwi/NW.Ni;

		// extra inhibitory update
		if( (1./NW.tau - 1)> r.doub()) {
			i = r.int64() % (NW.Ni-1);
			currentIE = NW.Jie*dotproduct(NW.IE[i],ST.nwe,NW.Ne);
			currentII = NW.Jii*dotproduct(NW.II[i],ST.nwi,NW.Ni);
			currentIO = NW.Jio*dotproduct(NW.IO[i],ST.nwi,NW.No);
			current = currentIE+currentII+currentIO;
			if(current>NW.thi[i]){
				nwi_ii = 1;
				if(ST.nwi[i] == 1) {
					change_nwi =0;
				} else {
					change_nwi = 1;
					ST.nwi[i] = 1;
				}	
			} else {
				nwi_ii = 0;
				if(ST.nwi[i] == 0) {
					change_nwi = 0;
				} else {
					change_nwi = -1;
					ST.nwi[i] = 0;
				}
			}
			ST.nwi_activity[t] = ST.nwi_activity[t]+change_nwi/NW.Ni;
		}

	}
	if(NW.tau>1. && (1./NW.tau) > r.doub() ) {
		// update inhibitory population
		i = r.int64() % (NW.Ni-1);
		currentIE = NW.Jie*dotproduct(NW.IE[i],ST.nwe,NW.Ne);
		currentII = NW.Jii*dotproduct(NW.II[i],ST.nwi,NW.Ni);
		currentIO = NW.Jio*dotproduct(NW.IO[i],ST.nwo,NW.No);
		current = currentIE+currentII+currentIO;
		if(current>NW.thi[i]){
			nwi_ii = 1;
			if(ST.nwi[i] == 1) {
				change_nwi =0;
			} else {
				change_nwi = 1;
				ST.nwi[i] = 1;
			}	
		} else {
			nwi_ii = 0;
			if(ST.nwi[i] == 0) {
				change_nwi = 0;
			} else {
				change_nwi = -1;
				ST.nwi[i] = 0;
			}
		}
		ST.nwi_activity[t] = ST.nwi_activity[t-1]+change_nwi/NW.Ni;

	}
	if(NW.tau==1.) {
		// update inhibitory population
		i = r.int64() % (NW.Ni-1);
		currentIE = NW.Jie*dotproduct(NW.IE[i],ST.nwe,NW.Ne);
		currentII = NW.Jii*dotproduct(NW.II[i],ST.nwi,NW.Ni);
		currentIO = NW.Jio*dotproduct(NW.IO[i],ST.nwo,NW.No);
		current = currentIE+currentII+currentIO;
		if(current>NW.thi[i]){
			nwi_ii = 1;
			if(ST.nwi[i] == 1) {
				change_nwi =0;
			} else {
				change_nwi = 1;
				ST.nwi[i] = 1;
			}	
		} else {
			nwi_ii = 0;
			if(ST.nwi[i] == 0) {
				change_nwi = 0;
			} else {
				change_nwi = -1;
				ST.nwi[i] = 0;
			}
		}
		ST.nwi_activity[t] = ST.nwi_activity[t-1]+change_nwi/NW.Ni;

	}

}


template<class R>
void advanceNW_stable_state(const Network& NW, State& ST, R& r)
{

	
	int i, nwe_ie,nwi_ii;
	double change_nwe, change_nwi;
	double current,currentEE, currentEI, currentIE,
		currentII, currentEO, currentIO;

	// update excitatory node ie
	i = r.int64() % (NW.Ne-1);
	currentEE = NW.Jee*dotproduct(NW.EE[i],ST.nwe,NW.Ne);
	currentEI = NW.Jei*dotproduct(NW.EI[i],ST.nwi,NW.Ni);
	currentEO = NW.Jeo*dotproduct(NW.EO[i],ST.nwo,NW.No);
	current = currentEE+currentEO+currentEI;
	if(current>NW.the[i]){
		ST.nwe[i] = 1;
	} else {
		ST.nwe[i] = 0;
	}
		
	if(NW.tau<1.) {
		// update inhibitory population
		i = r.int64() % (NW.Ni-1);
		currentIE = NW.Jie*dotproduct(NW.IE[i],ST.nwe,NW.Ne);
		currentII = NW.Jii*dotproduct(NW.II[i],ST.nwi,NW.Ni);
		currentIO = NW.Jio*dotproduct(NW.IO[i],ST.nwo,NW.No);
		current = currentIE+currentII+currentIO;
		if(current>NW.thi[i]){
			ST.nwi[i] = 1;
		} else {
			ST.nwi[i] = 0;
		}

		// extra inhibitory update
		if( (1./NW.tau - 1)> r.doub()) {
			i = r.int64() % (NW.Ni-1);
			currentIE = NW.Jie*dotproduct(NW.IE[i],ST.nwe,NW.Ne);
			currentII = NW.Jii*dotproduct(NW.II[i],ST.nwi,NW.Ni);
			currentIO = NW.Jio*dotproduct(NW.IO[i],ST.nwi,NW.No);
			current = currentIE+currentII+currentIO;
			if(current>NW.thi[i]){
				ST.nwi[i] = 1;
			} else {
				ST.nwi[i] = 0;
			}
		}

	}
	if(NW.tau>1. && (1./NW.tau) > r.doub() ) {
		// update inhibitory population
		i = r.int64() % (NW.Ni-1);
		currentIE = NW.Jie*dotproduct(NW.IE[i],ST.nwe,NW.Ne);
		currentII = NW.Jii*dotproduct(NW.II[i],ST.nwi,NW.Ni);
		currentIO = NW.Jio*dotproduct(NW.IO[i],ST.nwo,NW.No);
		current = currentIE+currentII+currentIO;
		if(current>NW.thi[i]){
			ST.nwi[i] = 1;
		} else {
			ST.nwi[i] = 0;
		}		
	}
	if(NW.tau==1.) {
		// update inhibitory population
		i = r.int64() % (NW.Ni-1);
		currentIE = NW.Jie*dotproduct(NW.IE[i],ST.nwe,NW.Ne);
		currentII = NW.Jii*dotproduct(NW.II[i],ST.nwi,NW.Ni);
		currentIO = NW.Jio*dotproduct(NW.IO[i],ST.nwo,NW.No);
		current = currentIE+currentII+currentIO;
		if(current>NW.thi[i]){
			ST.nwi[i] = 1;
		} else {
			ST.nwi[i] = 0;
		}
	}

}

// advance and learn
template<class R>
void advance_learnNW(Network& NW, State& ST, double rate, int t, R& r)
{

	
	int i, nwe_ie,nwi_ii, ii;
	double change_nwe, change_nwi;
	double current,currentEE, currentEI, currentIE,
		currentII, currentEO, currentIO;

	// update excitatory node ie
	i = r.int64() % (NW.Ne-1);
	currentEE = NW.Jee*dotproduct(NW.EE[i],ST.nwe,NW.Ne);
	currentEI = NW.Jei*dotproduct(NW.EI[i],ST.nwi,NW.Ni);
	currentEO = NW.Jeo*dotproduct(NW.EO[i],ST.nwo,NW.No);
	current = currentEE+currentEO+currentEI;
	if(current>NW.the[i]){
		// break an excitatory connection
		// add inhibitory connection
		ii = find(NW.EE[i],NW.Ne,1,r);
		if(rate>r.doub()) NW.EE[i][ii] = 0;
		ii = find(NW.EI[i],NW.Ni,0,r);
		if(rate>r.doub()) NW.EI[i][ii] = 1;

		nwe_ie = 1;
		if(ST.nwe[i] == 1){
			 change_nwe = 0;
		} else {
			change_nwe = 1;
			ST.nwe[i] = 1;
			if(i == ST.inode) ST.node_spike[t] = 1;
		}
	} else{
		// add an excitatory connection
		// break inhibitory connection
		ii = find(NW.EE[i],NW.Ne,0,r);
		if(rate>r.doub()) NW.EE[i][ii] = 1;
		ii = find(NW.EI[i],NW.Ni,1,r);
		if(rate>r.doub()) NW.EI[i][ii] = 0;

		nwe_ie = 0;
		if(ST.nwe[i] == 0) {
			change_nwe = 0;
		} else {
			change_nwe = -1;
			ST.nwe[i] = 0;
		}
	}
	ST.nwe_activity[t] = ST.nwe_activity[t-1]+change_nwe/NW.Ne;
		
	if(NW.tau<1.) {
		// update inhibitory population
		i = r.int64() % (NW.Ni-1);
		currentIE = NW.Jie*dotproduct(NW.IE[i],ST.nwe,NW.Ne);
		currentII = NW.Jii*dotproduct(NW.II[i],ST.nwi,NW.Ni);
		currentIO = NW.Jio*dotproduct(NW.IO[i],ST.nwo,NW.No);
		current = currentIE+currentII+currentIO;
		if(current>NW.thi[i]){
			// break an excitatory connection
			// add inhibitory connection
			ii = find(NW.IE[i],NW.Ne,1,r);
			if(rate>r.doub()) NW.IE[i][ii] = 0;
			ii = find(NW.II[i],NW.Ni,0,r);
			if(rate>r.doub()) NW.II[i][ii] = 1;

			nwi_ii = 1;
			if(ST.nwi[i] == 1) {
				change_nwi =0;
			} else {
				change_nwi = 1;
				ST.nwi[i] = 1;
			}	
		} else {
			// add an excitatory connection
			// break inhibitory connection
			ii = find(NW.IE[i],NW.Ne,0,r);
			if(rate>r.doub()) NW.IE[i][ii] = 1;
			ii = find(NW.II[i],NW.Ni,1,r);
			if(rate>r.doub()) NW.II[i][ii] = 0;

			nwi_ii = 0;
			if(ST.nwi[i] == 0) {
				change_nwi = 0;
			} else {
				change_nwi = -1;
				ST.nwi[i] = 0;
			}
		}
		ST.nwi_activity[t] = ST.nwi_activity[t-1]+change_nwi/NW.Ni;

		// extra inhibitory update
		if( (1./NW.tau - 1)> r.doub()) {
			i = r.int64() % (NW.Ni-1);
			currentIE = NW.Jie*dotproduct(NW.IE[i],ST.nwe,NW.Ne);
			currentII = NW.Jii*dotproduct(NW.II[i],ST.nwi,NW.Ni);
			currentIO = NW.Jio*dotproduct(NW.IO[i],ST.nwi,NW.No);
			current = currentIE+currentII+currentIO;
			if(current>NW.thi[i]){
				// break an excitatory connection
				// add inhibitory connection
				ii = find(NW.IE[i],NW.Ne,1,r);
				if(rate>r.doub()) NW.IE[i][ii] = 0;
				ii = find(NW.II[i],NW.Ni,0,r);
				if(rate>r.doub()) NW.II[i][ii] = 1;


				nwi_ii = 1;
				if(ST.nwi[i] == 1) {
					change_nwi =0;
				} else {
					change_nwi = 1;
					ST.nwi[i] = 1;
				}	
			} else {
				// add an excitatory connection
				// break inhibitory connection
				ii = find(NW.IE[i],NW.Ne,0,r);
				if(rate>r.doub()) NW.EE[i][ii] = 1;
				ii = find(NW.II[i],NW.Ni,1,r);
				if(rate>r.doub()) NW.EI[i][ii] = 0;

		
				nwi_ii = 0;
				if(ST.nwi[i] == 0) {
					change_nwi = 0;
				} else {
					change_nwi = -1;
					ST.nwi[i] = 0;
				}
			}
			ST.nwi_activity[t] = ST.nwi_activity[t]+change_nwi/NW.Ni;
		}

	}
	if(NW.tau>1. && (1./NW.tau) > r.doub() ) {
		// update inhibitory population
		i = r.int64() % (NW.Ni-1);
		currentIE = NW.Jie*dotproduct(NW.IE[i],ST.nwe,NW.Ne);
		currentII = NW.Jii*dotproduct(NW.II[i],ST.nwi,NW.Ni);
		currentIO = NW.Jio*dotproduct(NW.IO[i],ST.nwo,NW.No);
		current = currentIE+currentII+currentIO;
		if(current>NW.thi[i]){
			// break an excitatory connection
			// add inhibitory connection
			ii = find(NW.IE[i],NW.Ne,1,r);
			if(rate>r.doub()) NW.IE[i][ii] = 0;
			ii = find(NW.II[i],NW.Ni,0,r);
			if(rate>r.doub()) NW.II[i][ii] = 1;


			nwi_ii = 1;
			if(ST.nwi[i] == 1) {
				change_nwi =0;
			} else {
				change_nwi = 1;
				ST.nwi[i] = 1;
			}	
		} else {
			// add an excitatory connection
			// break inhibitory connection
			ii = find(NW.IE[i],NW.Ne,0,r);
			if(rate>r.doub()) NW.IE[i][ii] = 1;
			ii = find(NW.II[i],NW.Ni,1,r);
			if(rate>r.doub()) NW.II[i][ii] = 0;


			nwi_ii = 0;
			if(ST.nwi[i] == 0) {
				change_nwi = 0;
			} else {
				change_nwi = -1;
				ST.nwi[i] = 0;
			}
		}
		ST.nwi_activity[t] = ST.nwi_activity[t-1]+change_nwi/NW.Ni;

	}
	if(NW.tau==1.) {
		// update inhibitory population
		i = r.int64() % (NW.Ni-1);
		currentIE = NW.Jie*dotproduct(NW.IE[i],ST.nwe,NW.Ne);
		currentII = NW.Jii*dotproduct(NW.II[i],ST.nwi,NW.Ni);
		currentIO = NW.Jio*dotproduct(NW.IO[i],ST.nwo,NW.No);
		current = currentIE+currentII+currentIO;
		if(current>NW.thi[i]){
			// break an excitatory connection
			// add inhibitory connection
			ii = find(NW.IE[i],NW.Ne,1,r);
			if(rate>r.doub()) NW.IE[i][ii] = 0;
			ii = find(NW.II[i],NW.Ni,0,r);
			if(rate>r.doub()) NW.II[i][ii] = 1;


			nwi_ii = 1;
			if(ST.nwi[i] == 1) {
				change_nwi =0;
			} else {
				change_nwi = 1;
				ST.nwi[i] = 1;
			}	
		} else {
			// add an excitatory connection
			// break inhibitory connection
			ii = find(NW.IE[i],NW.Ne,0,r);
			if(rate>r.doub()) NW.IE[i][ii] = 1;
			ii = find(NW.II[i],NW.Ni,1,r);
			if(rate>r.doub()) NW.II[i][ii] = 0;


			nwi_ii = 0;
			if(ST.nwi[i] == 0) {
				change_nwi = 0;
			} else {
				change_nwi = -1;
				ST.nwi[i] = 0;
			}
		}
		ST.nwi_activity[t] = ST.nwi_activity[t-1]+change_nwi/NW.Ni;

	}

}
#endif



#ifndef GUARD_advance_network_h
#define GUARD_advance_network_h

#include "network_struct.h"
#include "generate_matrix.h"

#include <vector>



// ONLY WORKS FOR Ne=Ni !!!!!
template<class R>
void advanceNW(const Network& NW, State& ST, int t, R& r)
{	
	int i, nwe_ie,nwi_ii;
	double current,currentEE, currentEI, currentIE,
		currentII, currentEO, currentIO;

	// update excitatory node ie
	i = r.int64() % (NW.Ne-1);
	currentEE = NW.Jee*dotproduct(NW.EE[i],ST.nwe,NW.Ne);
	currentEI = NW.Jei*dotproduct(NW.EI[i],ST.nwi,NW.Ni);
	currentEO = NW.Jeo*dotproduct(NW.EO[i],ST.nwo,NW.No);
	current = currentEE+currentEO+currentEI;
	if(i != ST.inode) {
		++ST.current_sti;
		if(current>NW.the[i]){
			ST.nwe[i] = 1;
		} else {
			ST.nwe[i] = 0;
		}
	} else {
		// updateing inode, check spike
		if(current>NW.the[i]) {
			if(ST.nwe[i] == 0) {
				// spike
				ST.sti.push_back(ST.current_sti);
				ST.current_sti = 0;
				++ST.total_sti;
			} else {
				//no spike
				++ST.current_sti;
			}
			ST.nwe[i] = 1;
		} else {
			++ST.current_sti;
			ST.nwe[i] = 0;
		}
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

#endif


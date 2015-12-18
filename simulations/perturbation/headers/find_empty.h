#ifndef GUARD_find_empty_h
#define GUARD_find_empty_h

#include <vector>

template<class R>
void find_empty(int& i, const std::vector<int>& v, int N, R& r, bool& empty)
{
	int MAX = 1000*N;
	empty = false;
	while(!empty && MAX>0) {
		--MAX;
		i = r.int64() % (N-1);
		if(v[i] == 0) empty = true;
	}
}



#endif

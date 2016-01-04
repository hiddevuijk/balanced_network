#ifndef GUARD_check_empty_h
#define GUARD_check_empty_h

template<class V>
bool check_empty(const V& vec, const int& N)
{
	bool empty = true;
	for(int i=0;i<N;i++) {
		if(vec[i] !=0){
			empty = false;
			break;
		}
	}
	return empty;
}




#endif


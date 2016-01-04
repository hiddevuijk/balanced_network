#ifndef GUARD_generate_matrix_h
#define GUARD_generate_matrix_h

template<class R>
std::vector<std::vector<int> > connectivity_matrix(int N1, int N2, int K, R& r)
{
	std::vector<int> row_vec;
	std::vector<std::vector<int> > C(N1,row_vec);

	double p = (double)K/N2;

	for(int i=0;i<N1;i++) {
		std::vector<int> row(N2,0);
		for(int j=0;j<N2;j++) {
			if(p>r.doub()) row[j] = 1; 
		}
		C[i] = row; 
	}

	return C;
}

int dotproduct(const std::vector<int>& v1, const std::vector<int>& v2,const int& N)
{
	int sum = 0;
	for(int i=0;i<N;i++)
		sum += v1[i]*v2[i];
	return sum;
}

double average_vec(const std::vector<int>& vec, int N) {
	double sum = 0.0;
	for(int i=0;i<N;i++) sum += vec[i];
	return sum/(double)N;
}

#endif


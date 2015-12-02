#ifndef GUARD_generate_matrix_h
#define GUARD_generate_matrix_h




template<class R>
std::vector<Eigen::VectorXi> connectivity_matrix(int N1, int N2, int K, R& r)
{
	Eigen::VectorXi row_vec(N2);
	std::vector<Eigen::VectorXi> C(N1,row_vec);

	double p = (double)K/N2;

	for(int i=0;i<N1;i++) {
		Eigen::VectorXi row(N2);
		for(int j=0;j<N2;j++) {
			if(p>r.doub()) row[j] = 1; 
		}
		C[i] = row; 
	}

	return C;
}

int dotproduct(const Eigen::VectorXi& v1, const std::vector<int>& v2,const int& N)
{
	int sum = 0;
	for(int i=0;i<N;i++)
		sum += v1[i]*v2[i];
	return sum;
}

double average_vec(const Eigen::VectorXi& vec, int N) {
	double sum = 0.0;
	for(int i=0;i<N;i++) sum += vec[i];
	return sum/(double)N;
}

#endif


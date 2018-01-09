#include <iostream>
#include <string>
#include <vector>
#include "matutil.cpp"
#include <omp.h>
#include <math.h>

using namespace std;

static long int TESTSIZE = 8;
static long int NUMBER_THREADS = 4;

template<typename S>
void forward_seq_lcs(vector<vector<long int> > &A, const S& X,const S& Y) {
	for(long int i = 1; i < A.size(); ++i)
		for(long int j = 1; j < A[0].size(); ++j)
			if(X[i-1] == Y[j-1]) {
				A[i][j] = A[i-1][j-1] + 1;
			}
			else {
				A[i][j] = max(A[i][j-1], A[i-1][j]);
			}
}

string backward_seq_lcs(vector<vector<long int> > & A, const string & s1,
		const string & s2)
{
	long int i = A.size() - 1;
	long int j = A[i].size() - 1;
	// don't know size unless seq
	string result = "";
	while (i > 0 && j > 0) {
		if (s1[i-1] == s2[j-1]) {
			result.insert(0,1,s1[i-1]);
			--j; --i;
		}
		else if (A[i-1][j] > A[i][j-1]) {
			--i;
		}
		else {
			--j;
		}
	}
	return result;
}

// seems a bit off
void print_diag(vector<vector<long int> > & A,int ysize) {
	bool printed = true;
	int stage = 0;
	int row = 0;
	while(printed) {
		printed = false;
		stage = 0;
		// prints upper triangle
		while(stage < A.size()) {
			if (stage <= ysize + 1) {
				if (row < A[stage].size()) {
					cout << A[stage][row] << " ";
					printed = true;
				}
			} else {
				
			}
			++stage;
		}
		cout << endl;
		++row;
	}
}
string backward_seq_lcs_diag(vector<vector<long int> > & A, const string & s1,
		const string & s2)
{
	//print_diag(A,s2.size());
	long int stage = A.size() - 1;
	long int elem = 0;
	string result = ""; // actually don't know size beforehand unless seq
	while (stage > s2.size() + 1 || (elem > 0 && stage-elem-1 > 0)) {
	//while (stage > s2.size() + 1 || (elem - 1 >= 0 && stage-elem-2 >= 0)) {
	//while (stage > 2) {
		//cout << "\nstage,elem: " << stage << "," << elem << "\n";
		if (stage <= s2.size() + 1) {
			if (s1[elem-1] == s2[stage-elem-1-1]) {
				result = s1[elem-1] + result;
				stage-=2; --elem;
			}
			else if (A[stage-1][elem] > A[stage-1][elem-1]) {
				--stage;
			}
			else {
				--stage; --elem;
			}
		} else {
			if (s1[stage-1+elem-s2.size()-1] == s2[s2.size()-elem-1]) {
				result = s2[s2.size()-elem-1] + result;
				stage-=2; ++elem;
			}
			else if (A[stage-1][elem] > A[stage-1][elem+1]) {
				--stage;
			}
			else {
				--stage; ++elem;
			}
		}
	}
	return result;
}

void rand_row(vector<long int> & row) {
	for (long int i = 1; i < row.size(); ++i) {
		row[i] = rand() % 100 + 1;
	}
}
static long int iterations = 0;
void forward_LTDP_lcs_row (vector<vector<long int> > & A, const string & X, const
		string & Y)
{
	// could just do have for ltdp and half for wavefront when it comes to it
	// for a simple benchmark
	long int p = NUMBER_THREADS;
	long int n = A.size() - 1; // thread 0 starts calculating from stage 1 not 0
	long int tid, first_stage, last_stage, i,j;
	for (tid = 0; tid < p; ++tid) {
		last_stage =  (tid+1)*n/p;
		rand_row(A[last_stage]);
	}
	// diff, and temp are used to check parallel
	long int diff,temp;
	bool all_converged = false;
	bool converged;
	diff = -1000;
	// it should converge in at most threads number of iterations
	while (!all_converged) {
		++iterations;
		all_converged = true;
#pragma omp parallel private(tid,first_stage,last_stage,i,diff,temp,j,converged) reduction(&&:all_converged)
		{
			tid = omp_get_thread_num();
			first_stage = 1 + tid*n/p;
			last_stage =  (tid+1)*n/p;
			for (i = first_stage; i <= last_stage; ++i) {
				converged = true;
				for(j = 1; j < A[0].size(); ++j) {
					if(X[i-1] == Y[j-1]) {
						temp = A[i-1][j-1] + 1;
					}
					else {
						temp = max(A[i][j-1], A[i-1][j]);
					}
					if (j > 1 && diff != temp - A[i][j])
						converged = false;
					diff = temp - A[i][j];
					A[i][j] = temp;
				}
				if (iterations < 2) {
					converged = false;
				}
				if (converged) {
					break;
				}
			}
			all_converged = all_converged && converged;
		}
	}
}

void rand_diag(vector<long int> & row,long stage,int ysize) {
	if (stage <= ysize+1) {
		for (long int i = 1; i < row.size()-1; ++i) {
			row[i] = rand() % 100 + 1;
		}
	} else {
		for (long int i = 0; i < row.size() > 0; ++i) {
			row[i] = rand() % 100 + 1;
		}
	}
}

// maybe should make these inline
// ignore the conv variable for wavefront only
void lcs_diag_elem_upper(vector<vector<long> > & A, vector<vector<long> > & B, const string & X,
		const string & Y, long stage, long elem, bool & conv)
{
	long temp, diff;
	if(X[elem-1] == Y[stage-elem-2]) {
		temp = A[stage-2][elem-1] + 1;
	}
	else {
		temp = max(A[stage-1][elem-1], A[stage-1][elem]);
	}
	if (elem > 1 && diff != temp - A[stage][elem])
		conv = false;
	diff = temp - A[stage][elem];
	B[stage][elem] = temp;
}
long lcs_diag_elem_lower(vector<vector<long> > & A, vector<vector<long> > & B, const string & X,
		const string & Y, long stage, long elem, bool & conv)
{
	long temp, diff;
	if(X[stage-1+elem-Y.size()-1] == Y[Y.size()-elem-1]) {
		temp = A[stage-2][elem+1] + 1;
	}
	else {
		temp = max(A[stage-1][elem+1], A[stage-1][elem]);
	}
	if (elem && diff != temp - A[stage][elem])
		conv = false;
	diff = temp - A[stage][elem];
	B[stage][elem] = temp;
}

// pass in A twice, if read and write to same array
// make this inline perhaps
bool lcs_diag(vector<vector<long> > & A, vector<vector<long> > & B, const string & X,
		const string & Y, long stage)
{
	bool conv = true;
	long elem;
	// could do register allocation on all these indices
	if (stage <= (Y.size() + 1)) {
		for(elem = 1; elem < A[stage].size() - 1; ++elem) {
			lcs_diag_elem_upper(A,B,X,Y,stage,elem,conv);
		}
	} else {
		for(elem = 0; elem < A[stage].size(); ++elem) {
			lcs_diag_elem_lower(A,B,X,Y,stage,elem,conv);
		}
	}
	return conv;
}
// inline this .. probs
void memoize_diag(vector<vector<long int> > & A, const string & X,
		const string & Y)
{
	A.resize(X.size()+Y.size()+2,vector<long int>());
	for (long int stage = 0; stage < A.size(); ++stage) {
		if (stage <= Y.size()+1) {
			long int new_size = min(stage,(long int)X.size()+1);
			A.at(stage).resize(new_size,0);
		} else {
			long int ind = 1 + Y.size() - (stage - (Y.size() + 1));
			long int new_size = min(ind,(long int)(X.size())+1);
			A.at(stage).resize(new_size,0);
		}
	}
}
// i think the lesson is that if statements are WORSE than extra variables
// probably should have just went with the index translator function ...
void forward_LTDP_lcs_diag (vector<vector<long int> > & A, const string & X,
		const string & Y)
{
	memoize_diag(A,X,Y);
	// initialize random stages
	long int p = NUMBER_THREADS;
	// there are better ways of distributing threads ofc..
	long int n = A.size() - 3; // just going to distribute stages to threads evenly
	long int tid, first_stage, last_stage, stage, elem;
	for (tid = 0; tid < p; ++tid) {
		last_stage =  2 + (tid+1)*n/p;
		// now need to init random the last two stages
		rand_diag(A[last_stage-1],last_stage-1,Y.size());
		rand_diag(A[last_stage],last_stage,Y.size());
	}
	long int diff,temp;
	bool all_converged = false;
	bool converged;
	while (!all_converged) {
		++iterations;
		all_converged = true;
#pragma omp parallel private(tid,first_stage,last_stage,stage,diff,temp,elem,converged) reduction(&&:all_converged)
		{
			tid = omp_get_thread_num();
			first_stage = 3 + tid*n/p;
			last_stage =  2 + (tid+1)*n/p;
			// doing the triangles separately, removes the if statement check
			// every iteration
			for (stage = first_stage; stage <= last_stage; ++stage) {
				converged = lcs_diag(A,A,X,Y,stage);
				if (iterations < 2) {
					converged = false;
				}
				if (converged) {
					break;
				}
			}
			all_converged = all_converged && converged;
		}
	}
}
void forward_lcs_wavefront (vector<vector<long int> > & A, const string & X,
		const string & Y)
{
	memoize_diag(A,X,Y);
	bool conv = false; // this is ignored
	long elem;
	for (long stage = 3; stage < A.size(); ++stage) {
		if (NUMBER_THREADS < 2*A[stage].size()) {
			// sequential
			lcs_diag(A,A,X,Y,stage);
		}
		else {
#pragma omp parallel for private(elem,conv)
			for (elem = 0; elem < A[stage].size(); ++elem) {
				if (stage <= (Y.size() + 1)) {
					for(elem = 1; elem < A[stage].size() - 1; ++elem) {
						lcs_diag_elem_upper(A,A,X,Y,stage,elem,conv);
					}
				} else {
					for(elem = 0; elem < A[stage].size(); ++elem) {
						lcs_diag_elem_lower(A,A,X,Y,stage,elem,conv);
					}
				}
			}
		}
	}
}
// uses nested parallelism
// keep numthreads the same
// take the sqrt, for now just assume threads can be sqr rooted
// then set it back to NUMBER_THREADS
void forward_LTDP_lcs_diag_wavefront (vector<vector<long int> > & A, const string & X,
		const string & Y)
{
	memoize_diag(A,X,Y);
	// initialize random stages
	long int p;
	p = (int)sqrt(0.5+NUMBER_THREADS);
	omp_set_num_threads(p);
	//p = NUMBER_THREADS;

	// there are better ways of distributing threads ofc..
	long int n = A.size() - 3; // just going to distribute stages to threads evenly
	long int tid, first_stage, last_stage, stage, elem;
	for (tid = 0; tid < p; ++tid) {
		last_stage =  2 + (tid+1)*n/p;
		// now need to init random the last two stages
		rand_diag(A[last_stage-1],last_stage-1,Y.size());
		rand_diag(A[last_stage],last_stage,Y.size());
	}
	long int diff,temp;
	bool all_converged = false;
	bool converged;
	while (!all_converged) {
		++iterations;
		all_converged = true;
#pragma omp parallel private(tid,first_stage,last_stage,stage,diff,temp,elem,converged) reduction(&&:all_converged)
		{
			tid = omp_get_thread_num();
			first_stage = 3 + tid*n/p;
			last_stage =  2 + (tid+1)*n/p;

			// doing the triangles separately, removes the if statement check
			// every iteration
			for (stage = first_stage; stage <= last_stage; ++stage) {
				converged = true;
				// should make inline function or just use an extra variable prev_diff
				// could do register allocation on all these indices

				long int first_element = -2;
				long int last_element = -2;
				long int tid2;
				long int n2;
				bool conv;
				// try to fix the bad distribution of threads from ltdp_diag
				if (omp_get_num_threads() < 8*A[stage].size()) {
					converged = lcs_diag(A,A,X,Y,stage);
					if (iterations < 2) {
						converged = false;
					}
					if (converged) {
						break;
					}
				} else {
#pragma omp parallel private(first_element,last_element,tid2,elem,temp,conv) reduction(&&:converged)
					{
						tid2 = omp_get_thread_num();
						conv = true;
						if (stage <= (Y.size() + 1)) {
							n2 = A[stage].size() - 2;
							first_element = 1 + tid2*n2/p;
							last_element = 0 + (tid2+1)*n2/p;

							for(elem = first_element; elem <= last_element; ++elem) {
								lcs_diag_elem_upper(A,A,X,Y,stage,elem,conv);
							}
						} else {
							n2 = A[stage].size();
							first_element = 0 + tid2*n2/p;
							last_element = -1 + (tid2+1)*n2/p;
							for(elem = first_element; elem <= last_element; ++elem) {
								lcs_diag_elem_lower(A,A,X,Y,stage,elem,conv);
							}
						}
						converged = converged && conv;
					}
				}
				if (iterations < 2) {
					converged = false;
				}
				if (converged) {
					break;
				}
			}
			all_converged = all_converged && converged;
		}
	}
	omp_set_num_threads(NUMBER_THREADS);
}

// return value is whether row is parallel
// make inline?
// if A and B are the same then well it's the same o__o, you can put it the
// same, having two of them just makes it easy heh
bool lcs_row(vector<vector<long> > & A, vector<vector<long> > & B, const string & X,
		const string & Y, long stage)
{
	bool conv = true;
	long diff;
	long temp;
	for(int j = 1; j < A[stage].size(); ++j) {
		if(X[stage-1] == Y[j-1]) {
			temp = A[stage-1][j-1] + 1;
		}
		else {
			temp = max(A[stage][j-1], A[stage-1][j]);
		}
		if (j > 1 && diff != temp - A[stage][j])
			conv = false;
		diff = temp - A[stage][j];
		B[stage][j] = temp;
	}
	return conv;
}
// this will work if you chunk it into blocks of 10
// and then each block will have a converged flag
// need an array of converged
vector<vector< long int> > forward_LTDP_lcs_row_dynamic (const string & X, const string & Y)
{
	// dynamic needs two arrays A and B
	vector<vector< long int> > A(X.size()+1, vector<long>(Y.size()+1,0));
	vector<vector< long int> > B(X.size()+1, vector<long>(Y.size()+1,0));
#pragma omp parallel for
	for (int stage = 2; stage < A.size(); ++stage) {
		rand_row(A[stage]);
	}
	long int diff,temp,i,j,m,n;
	m = A.size();
	n = A[0].size();
	bool all_converged = false;
	bool converged;
	bool write_to_B = true;
	while (!all_converged) {
		++iterations;
		all_converged = true;
		// schedule(guided) vs dynamic?
#pragma omp parallel for private(diff,converged) reduction(&&:all_converged) schedule(static)
		for (i = 1; i < m; ++i) {
			if (i == 1) {
				//cout << "there are " << omp_get_num_threads() << " threads running\n";
			}
			converged = write_to_B ? lcs_row(A,B,X,Y,i) : lcs_row(B,A,X,Y,i);
			if (iterations < 2) {
				converged = false;
			}
			all_converged = all_converged && converged;
		}
		write_to_B = write_to_B ? false : true;
		//cout << write_to_B << endl;
	}
	if (write_to_B) {
		return A;
	}
	return B;
}

void rand_sequence(string &s) {
	for(long int i = 0; i < s.size(); ++i)
		s[i] = (rand() % 4 + 'a'); 
}

// long int is fine up to 2 million so...
// omp_set_dynamic(0);     // Explicitly disable dynamic teams
// omp_set_num_threads(4); // Use 4 threads for all consecutive parallel regions
// omp_get_num_threads() // use in parallelized sections to verify num threads
// http://stackoverflow.com/questions/11095309/openmp-set-num-threads-is-not-working


// https://docs.oracle.com/cd/E19205-01/819-5270/aewbc/index.html
// may need nested parallelism to combine LTDP and wavefront
// only times forward at the moment
// make sure to be measuring complete time of calculation
// including backward phase
string timed_run(string s,const string & s1, const string & s2) {
	iterations = 0;
	if (s == "seq" || s == "ltdp_row") {
		vector<vector<long int> > A(s1.size()+1,vector<long int>(s2.size()+1,0));
		double time = -omp_get_wtime();
		if (s == "seq") {
			forward_seq_lcs(A,s1,s2);
			s = s + " ";
		} else if (s == "ltdp_row") {
			omp_set_dynamic(0);
			omp_set_num_threads(NUMBER_THREADS);
			forward_LTDP_lcs_row(A,s1,s2);
			cout << s << " converged in " << iterations << " iterations\n";
		}
		time+= omp_get_wtime();
		cout << s << " timed: " << time << endl;
		return backward_seq_lcs(A,s1,s2);
	} else if (s == "ltdp_diag") {
		vector<vector<long int> > B;
		omp_set_dynamic(0);
		omp_set_num_threads(NUMBER_THREADS);

		double time = -omp_get_wtime();
		forward_LTDP_lcs_diag(B,s1,s2);
		time+= omp_get_wtime();

		cout << s << " converged in " << iterations << " iterations\n";
		cout << s << " timed: " << time << endl;
		return backward_seq_lcs_diag(B,s1,s2);
	} else if (s == "ltdp_wavefront") {
		vector<vector<long int> > B;
		omp_set_dynamic(0);
		omp_set_num_threads(NUMBER_THREADS);
		omp_set_nested(1);

		double time = -omp_get_wtime();
		forward_LTDP_lcs_diag_wavefront(B,s1,s2);
		time+= omp_get_wtime();

		cout << s << " converged in " << iterations << " iterations\n";
		cout << s << " timed: " << time << endl;
		omp_set_nested(0);
		return backward_seq_lcs_diag(B,s1,s2);
	} else if (s == "dyn_row") {
		omp_set_num_threads(NUMBER_THREADS);
		omp_set_dynamic(1);
		double time = -omp_get_wtime();
		vector<vector<long int> > A;
		A = forward_LTDP_lcs_row_dynamic(s1,s2);
		time+= omp_get_wtime();
		cout << s << " converged in " << iterations << " iterations\n";
		cout << s << " timed: " << time << endl;
		return backward_seq_lcs(A,s1,s2);
	} else if (s == "wavefront") {
		omp_set_num_threads(NUMBER_THREADS);
		omp_set_dynamic(0);
		double time = -omp_get_wtime();
		vector<vector<long int> > A;
		forward_lcs_wavefront(A,s1,s2);
		time+= omp_get_wtime();
		cout << s << " timed: " << time << endl;
		return backward_seq_lcs_diag(A,s1,s2);
	}
}
void compare_seq_ltdp_row(string s1, string s2) {
	string res2 = timed_run("seq",s1,s2);

	string res1 = timed_run("ltdp_row",s1,s2);
	if (TESTSIZE < 20) {
		cout << s1 << endl;
		cout << s2 << endl;
		cout << res1 << endl;
		cout << res2 << endl;
	}
	if (res1 == res2) {
		cout << "results match\n";
		cout << NUMBER_THREADS << " threads were used\n";
		cout << "word1 was size: " << s1.size() << endl;
		cout << "word2 was size: " << s2.size() << endl;
	}
}
void compare_all(string & s1, string & s2) {
	int res1,res2,res3,res4,res5,res6;
	int result[32];
	int run = 0;
	result[run++]=timed_run("seq",s1,s2).size();
	result[run++]=timed_run("ltdp_row",s2,s1).size();
	result[run++]=timed_run("ltdp_diag",s1,s2).size();
	result[run++]=timed_run("ltdp_wavefront",s1,s2).size();
	result[run++]=timed_run("wavefront",s1,s2).size();
	//result[run++]=timed_run("dyn_row",s1,s2).size();
	int i;
	for (i = 1; i < run; ++i) {
		if (result[i] == result[i-1]) {
			cout << "run " << i-1 << " matches run " << i << endl;
		} else {
			cout << result[i-1] << endl;
			cout << result[i] << endl;
			break;
		}
	}
	if (i == run) {
		cout << "all results match!\n";
	}
}

int main()
{
	TESTSIZE = 20000;
	NUMBER_THREADS = 4;
	srand(13);
	// why does values of 10 and 11 for width make diag and row different?
	string s1(TESTSIZE/100,'*'); rand_sequence(s1);  // make sure this string is smaller size
	string s2(TESTSIZE,'*'); rand_sequence(s2); 
	if (TESTSIZE < 20) {
		cout << s1 << endl;
		cout << s2 << endl;
	}
	//compare_seq_ltdp_row(s1,s2);
	compare_all(s1,s2);

	return 0;
}

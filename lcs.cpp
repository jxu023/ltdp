#include <vector>
#include <ctime>
#include <cstdlib>
#include <string>
#include <iostream>
#include <thread>
#include <cassert>

#include "ltdp.cpp"
#include "ltdp-stages.cpp"
#include "naive.cpp"
#include "backtrack.cpp"

#define TESTSIZE 2

using namespace std;


void setsequence(vector<char> &s) {
  for(int i = 0; i < s.size(); ++i)
    s[i] = (rand()%10 + '0'); 
}

int main(int argc, char * argv[]) {
  //timing taken from Quentin Perez http://stackoverflow.com/questions/15092504/how-to-time-a-function-in-milliseconds-without-boosttimer
  clock_t start;
  srand(13);
  //vector<char> s1 = "abcdef";
  //vector<char> s2 = "acebdf";
  vector<char> s1(TESTSIZE); setsequence(s1); 
  vector<char> s2(TESTSIZE); setsequence(s2); 

  if(TESTSIZE < 12) 
    cout << "s1:" << string(s1.data()) <<"\ns2:"<< string(s2.data()) <<endl;
  vector<vector<int> > Cnaive(s1.size()+1, vector<int>(s2.size()+1));
  vector<vector<int> > C(s1.size()+1, vector<int>(s2.size()+1));
  
  matutil::initdiag(Cnaive,2);
  matutil::initdiag(Cnaive,3);
  matutil::printmat(Cnaive);

  start = clock();
  vector<vector<int> > Ccoal = matutil::converttodiag(Cnaive);
  cout << "Time to convert for memory " << (clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << "ms" << endl;
  matutil::printmat(Ccoal);

  /*
  start = clock();
  naive::lcs(Cnaive, s1,s2);
  cout << "naive solution took " << (clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << "ms" << endl;
  matutil::printmat(Cnaive);

  start = clock();
  ltdp::gensol(C, s1,s2, ltdp::stage::lcs);
  cout << "staged solution took " << (clock() - start) / (double)(CLOCKS_PER_SEC / 1000) <<"ms"<< endl;
  assert(matutil::equal(Cnaive, C));
  
  start = clock();
  ltdp::gensol(C, s1,s2, ltdp::stage::wavefront::lcs);
  cout << "wavefront(openmp) staged solution took " << (clock() - start) / (double)(CLOCKS_PER_SEC / 1000) 
      <<"ms NOTE: this is CPU time, Num cores: "<< thread::hardware_concurrency()  <<endl;
  assert(matutil::equal(Cnaive,C));
  */
  cout << "successfully finished!" << endl;
  //if(TESTSIZE < 12) printmat(lens);
  /*
  set<vector<char>> allstrings = backtrackAll<string>(lens, s1,s2,s1.length(), s2.length());
  for(const auto& str : allstrings)
    cout << "backtrackall found: " << str << endl;
  */
}
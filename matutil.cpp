#include<vector>
#include<iostream>
#include <cstdlib>
#include <algorithm>    // std::adjacent_find
#include <functional>   // std::not_equal_to

namespace matutil {

    /*  need to run on well 
     */
  template<typename t>
    void printmat(std::vector<std::vector<t> > mat) {
      for (int i = 0; i < mat.size(); ++i) {
        for (int j = 0; j < mat[i].size(); ++j) {
          std::cout << mat[i][j] << "\t";
        }
        std::cout << std::endl; 
      }
    }


  template<typename t>
    bool equal(std::vector<std::vector<t> > X, std::vector<std::vector<t> > Y) {
      for(int r =0; r< X.size(); ++r)
        for(int c=0;c<X[r].size();++c)
          if(X[r][c] != Y[r][c])
            return false;
      return true;
    }

  void initdiag(std::vector<std::vector<int> > &C, int si){
    for(int i = 1; i < si; ++i) {
      if(si-i < C.size() && i < C[0].size()) {
        C[si-i][i] = std::rand(); //is this the best? maybe we want to start at something that is already almost parallel...
      }
    }
  }


  //this grabbing only works from 
  //we make a copy to check for parallelness 
  std::vector<int> grabdiag(const std::vector<std::vector<int> > &C, int diagi) {
    unsigned int s= std::max(0, 1+diagi-(int)(C.size()+C[0].size())/2);
    unsigned int t= std::min(diagi, (int)C[0].size()-1); 

    std::vector<int> diag(t-s+1);
    unsigned int d = 0;
    for(unsigned int i = s; i < diagi ; ++i) {
      diag[d] = C[diagi-i][i];
      ++d;
    }
    return diag;
  }


  //convert our matrix to represent depencies in a way thats a little better for memory access
  std::vector<std::vector<int> > converttodiag(const std::vector<std::vector<int> > &C) {
    std::vector<std::vector<int> > Cnew(C.size()+C[0].size()-1);
    for(unsigned int stage = 0;stage<Cnew.size();++stage)
      Cnew.at(stage) = matutil::grabdiag(C, stage);
    return Cnew;
  }
  //there is probably a better way to do this
  //Found it: try this http://stackoverflow.com/questions/7572640/how-do-i-know-if-two-vectors-are-near-parallel
  bool checkifdiagisparallel(std::vector<std::vector<int> > &C, int diagi, std::vector<int> checkdiag) {
    unsigned int s= std::max(1, 1+diagi-(int)(C.size()+C[0].size())/2);
    unsigned int t= std::min(diagi, (int)C.at(0).size()+1); 
    if(checkdiag.size() != s-t)
      return false;
    int diff = C[diagi-s-1][s-1] - checkdiag[0];
    unsigned int di = 1;
    int prevdiff;
    for(unsigned int i = s+1; s < t; ++i) {
      prevdiff = diff;
      diff = C[diagi-i-1][i-1] - checkdiag[di];
      if(diff != prevdiff) //early abondonment 
        return false;
      ++di;
    }
    return true;
    //the line below can be used if you calculate the difference in parallel and store them in a vector named "difference"
    //Vlad from Moscow http://stackoverflow.com/questions/20287095/checking-if-all-elements-of-a-vector-are-equal-in-c
    //return std::adjacent_find( difference.begin(), difference.end(), std::not_equal_to<int>() ) == difference.end();
  }


}

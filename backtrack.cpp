
#include <vector>
#include <string>
#include <algorithm>
#include <set> 

namespace backtrack {

    template<typename S>
    std::string lcs(const std::vector<std::vector<int> >& C, const S& X,const S& Y, int i, int j) {
        if(i == 0 || j == 0)
            return "";
        else if(X[i-1] == Y[j-1])
            return lcs(C, X, Y, i-1, j-1) + X[i-1];
        else
            if(C[i][j-1] > C[i-1][j])
                return lcs(C,X,Y,i,j-1);
            else
                return lcs(C,X,Y,i-1,j);
    }

    /* bleh bleh need to replace the '{""}' with 'new std::set()' or something to run on well
    template<typename S>
    std::set<std::string> lcsAll(const std::vector<std::vector<int> >& C, const S& X, const S& Y, int i, int j) {
        if(i == 0 || j == 0)
            return {""};
        else if(X[i-1] == Y[j-1]) {
            std::set<S> R;
            std::set<S> tmp = backtrackAll(C,X,Y,i-1,j-1);
            for(S s : tmp)
            R.insert(s + X[i-1]);
            return R;  
        } else {
            std::set<S> R;
            if(C[i][j-1] >= C[i-1][j]) {
            std::set<S> tmp = backtrackAll(C,X,Y,i,j-1);
            R.insert(tmp.begin(), tmp.end());
            }
            if(C[i-1][j] >= C[i][j-1]) {
            std::set<S> tmp = backtrackAll(C,X,Y,i-1,j);
            R.insert(tmp.begin(), tmp.end());
            }
            return R;
        }
    }
    */
}
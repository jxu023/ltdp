
#include <vector>

namespace naive {
    template<typename S>
    void lcs(std::vector<std::vector<int> > &C, const S& X,const S& Y) {
        for(int i =0; i < C.size(); ++i) C[i][0] = 0; 
        for(int j =0; j < C[0].size(); ++j) C[0][j] = 0;

        for(int i = 1; i <= X.size(); ++i)
            for(int j = 1; j <= Y.size(); ++j)
            if(X[i-1] == Y[j-1])
                C[i][j] = C[i-1][j-1] +1;
            else
                C[i][j] = std::max(C[i][j-1], C[i-1][j]);
    }
}
#include <vector>
#include <iostream>
#include <omp.h>

namespace ltdp {
    namespace stage{
        //note the if inside the for, this is bad 
        //this is for parallel versions of the stages
        // there is probably a better way of doing that too (maybe add another parameter)...
        // better yet we just check if its compiled with open mp or not... http://stackoverflow.com/questions/1300180/ignore-openmp-on-machine-that-doesnt-have-it
        // it's easier to keep it like this for now for testing purposes
        namespace wavefront{
            template<typename sequence>
            void lcs(std::vector<std::vector<int> > &C, const sequence& x, const sequence& y, int si) {
                unsigned int s= std::max(1, 1+si-(int)(C.size()+C[0].size())/2);
                unsigned int t= std::min(si, (int)C.at(0).size()+1); 
                #pragma omp parallel for
                for(unsigned int i = s; i < t; ++i) {
                    int xi = (si-i-1); int yi = (i-1);
                    if(x[xi] == y[yi])
                        C[si-i][i] = C[xi][yi] + 1; //match so we increment
                    else
                        C[si-i][i] = std::max(C[si-i][yi], C[xi][i]); //no match take the max
                }
            }
        }

        template<typename sequence>
        bool lcsmulti(const sequence& prev1,  const sequence& prev0, sequence& workstage,  const sequence& x, const sequence& y) {
            
            sequence copy = workstage;
            unsigned int si = workstage.size()-1;

            //if si > halfway through then start with i =0 and go to workstage.size()-1!
            if(workstage.size() > prev0.size()) {
                for(unsigned int i = 1; i < workstage.size()-1; ++i) {
                    unsigned int xi = si-i-1; //This are wrong
                    unsigned int yi = i-1;  //This is wrong
                    //here check that xi and yi are valid 
                    if(x[xi] == y[yi]) 
                        workstage[i] = prev1[i-1] + 1; 
                    else 
                        workstage[i] = std::max(prev0[i-1], prev0[i]); // this line could give you ima
                }
            } else {
                for(unsigned int i = 0; i < workstage.size(); ++i) {
                    unsigned int xi = si-i-1; //this is wrong
                    unsigned int yi = i-1; // this is wrong
                    if(x[xi] == y[yi]) 
                        workstage[i] = prev1[i] + 1;
                    else                            
                        workstage[i] = std::max(prev0[i],prev0[i+1]);
                }
            }
            return (isparallel(copy,workstage));
        }

		template<typename sequence>
        void lcs(std::vector<std::vector<int> > &C, const sequence& x, const sequence& y, int si) {
            //std::cout << "stage:"<<si<<std::endl;
            for(unsigned int i = std::max(1, 1+si-(int)(C.size()+C[0].size())/2); i < si && i-1 < y.size() ; ++i) {//start at 1 + 1 for every row past the middle
                unsigned int xi = (si-i-1);unsigned int yi = (i-1);
                //std::cout << X[xi] << "[" << xi << "]"  << "|" << Y[yi] << "[" << yi << "]" << "\t";
                if(x[xi] == y[yi])
                    C[si-i][i] = C[xi][yi] + 1; //match so we increment
                else
                    C[si-i][i] = std::max(C[si-i][yi], C[xi][i]); //no match take the max
                
            }            
            //std::cout << std::endl;
        }
    }
}


#include <vector>
#include <algorithm>  //std::max std::min
#include <functional> //std::all_of
#include <omp.h>

#include "matutil.cpp"

#define STARTSTAGE 2

namespace ltdp {
    
    //gensol and gensol par both are for the unmodiffied matrix C
    // The third argument is the signature of a stage iterating function
    template<typename sequence>
    void gensol(std::vector<std::vector<int> > & C, const sequence& x, const sequence& y, void (*stagefun) (std::vector<std::vector<int> >&, const sequence&, const sequence&, int)) {
        unsigned int numstages = x.size() + y.size();
        for(unsigned int i =0; i < C.size(); ++i) C[i][0] = 0;  
        for(unsigned int j =1; j < C[0].size(); ++j) C[0][j] = 0; //do you see we start at 1
        for(unsigned int i =2; i <= numstages; ++i) stagefun(C,x,y,i);
    }

    template<typename sequence>
    void gensolpar(std::vector<std::vector<int> > & C, const sequence& x, const sequence& y, 
                    void (*stagefunmem) (const sequence&, const sequence&, sequence &, const sequence&, const sequence&), 
                    int nps=2) {
        unsigned int numstages = x.size() + y.size()-1;
        unsigned int stagesleft = numstages;
        std::vector<bool> converged(numstages);
        
        //this is probably not the best way to do this
        for(unsigned int i = 0; i < numstages; ++i) {
            if(i <= numstages/2)
            for(unsigned int j = 0;)
        }

        for(unsigned int i = 2; i < numstages; ++i) {
            stagefunmem(C[i-2],C[i-1],C[i],x,y);
        }

        while(stagesleft > 0) {
            stagesleft = std::count(converged.begin(),converged.end(),false); //probably a better way to do this
            for(unsigned int i = 0; i < stagesleft; ++i) {
                
            }
        }
    }


    /*
    // nps number of process
    template<typename sequence>
    void gensolpar(std::vector<std::vector<int> > & C, const sequence& x, const sequence& y, void (*stagefun) (std::vector<std::vector<int> >&, const sequence&, const sequence&, int), int nps=2) {
        unsigned int numstages = x.size() + y.size();
        std::vector<bool> converged(numstages);
        std::vector<int> prevsol;
        
        unsigned int incr = (numstages - STARTSTAGE)/nps;

        // https://www.dartmouth.edu/~rc/classes/intro_openmp/schedule_loops.html
        omp_set_num_threads(nps);//we want the number of threads to be the number of parallel stages
        #pragma omp parallel private(prevsol)
        {
            //we still init our matrix in here because we already made the threads might as well exploit them
            #pragma omp for
            for(unsigned int i =0; i < C.size(); ++i) C[i][0] = 0;  
            #pragma omp for
            for(unsigned int j =1; j < C[0].size(); ++j) C[0][j] = 0; 
            #pragma omp for
            for(unsigned int stage = STARTSTAGE-1 + incr; stage < numstages; stage += incr) { //we don't have to init the first row and we are initing TWO previous rows... dependencies for the diagonal... this might need to be changed for other dependency types
                matutil::initstage(C, stage-1);
                matutil::initstage(C, stage-2);
            }

            //run each one up until just before the next thread started
            #pragma omp for 
            for(unsigned int stage = STARTSTAGE-1 + incr; stage < numstages; stage++) { //we don't have to init the first row and we are initing TWO previous rows... dependencies for the diagonal... this might need to be changed for other dependency types
                stagefun(C,x,y,stage);
            }
            //everything from the stage0 to stage0+incr-1 is converged, lets mark them
            #pragma omp for
            for(unsigned int run = 0; run < incr; ++run) converged[run] = true;
            /*
            //continue running, this time save the old values until we converged
            bool conv = false;
            do {


                conv = std::all_of(converged.begin(), converged.end(), true);
            } while(conv);

            //we have to make a 2d vector which holds the stage for each of our wavefronts before we over write it 
            //does it have to be 2d or can our threads each be responsible for their own
            //we will resize our vector... keep in mind we are getting the thread infront of us's solution...
            //should replace this with a for loop...
            unsigned int tid = omp_get_thread_num();
            unsigned int targetstage = ((tid+1)*incr);
            if(targetstage < numstages) {
                prevsol = matutil::grabdiag(C, targetstage);
                stagefun(C,x,y,targetstage);
                if(matutil::checkifdiagisparallel(C, targetstage, prevsol)) {
                    //it is parallel! mark the stage as complete and all of the stages up to this
                    //if its not the first stage we also mark every stage since it started as converged
                    //mark all prevoius stages up to stage-incr as marked
                    for(int i = 0; i< incr; ++i) converged[targetstage-i] = true;
                }
            }
            

            //then we will have another loop to check for parallelness...
            //mark all the finished stages (up to the next frontier)
            
            //continue to iterate until all stages finish
            
        }
    }
    */


}


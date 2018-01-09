#include<iostream>
#include<omp.h>

using namespace std;

int main() {
    omp_set_num_threads(2);//we want the number of threads to be the number of parallel stages
    #pragma omp parallel
    {
        #pragma omp for schedule(static)
        for(int i = 0; i < 10; ++i)
            cout << i << endl;
    }
}
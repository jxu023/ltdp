#a: lcs.cpp ltdp.cpp ltdp-stages.cpp matutil.cpp backtrack.cpp naive.cpp
	#g++ -std=c++11 -fopenmp lcs.cpp
a: main.cpp matutil.cpp
	g++ -fopenmp main.cpp -o main -g
run:
	./main

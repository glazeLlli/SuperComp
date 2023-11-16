all:
	g++ -O3 -fopenmp ompParallel.cpp -std=c++11 -o para40
	./para40

seqGen:
	g++ -O3 sequential.cpp -std=c++11 -o seq

submitSeq:
	./seqSubmit.sh

ompGen:
	g++  -O3 -fopenmp ompParallelSub.cpp -std=c++11 -o para40Sub

ompGen80:
	g++  -O3 -fopenmp ompParallelSub.cpp -std=c++11 -o para80Sub

ompGen160:
	g++  -O3 -fopenmp ompParallelSub.cpp -std=c++11 -o para160Sub

submitOmp:
	./ompSubmit.sh

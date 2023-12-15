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

MPIGen:
	g++ -O3 -std=c++11 -fopenmp mpiPure.cpp -o mpiPure -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm -DUSE_MPI

pureSub:
	./MPIPureSub.sh

HybridGen40:
	g++ -O3 -std=c++11 -fopenmp mpiHybrid.cpp -o Hybrid40 -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm -DUSE_MPI

HybridGen80:
	g++ -O3 -std=c++11 -fopenmp mpiHybrid.cpp -o Hybrid80 -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm -DUSE_MPI

HybridGen160:
	g++ -O3 -std=c++11 -fopenmp mpiHybrid.cpp -o Hybrid160 -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm -DUSE_MPI

HybridSub1:
	./HybridSub1.sh

HybridSub2:
	./HybridSub2.sh

HybridSub3:
	./HybridSub3.sh
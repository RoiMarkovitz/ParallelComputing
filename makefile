build:
	mpicxx -fopenmp -c main.c -o main.o
	mpicxx -fopenmp -c cFunctions.c -o cFunctions.o
	mpicxx -fopenmp -c mpiUtil.c -o mpiUtil.o
	nvcc -I./inc -c cudaFunctions.cu -o cudaFunctions.o
	mpicxx -fopenmp -o alignment_score main.o cFunctions.o mpiUtil.o cudaFunctions.o  /usr/local/cuda-11.0/lib64/libcudart_static.a -ldl -lrt
	 
clean:
	rm -f *.o ./alignment_score

run:
	mpiexec -np 2 ./alignment_score < input.txt 

runOn2:	
	mpiexec -np 2 -machinefile  mf  -map-by  node  ./alignment_score



# ParallelComputing
 
A project for the course "Parallel and distributed computing" at the academy.
The project is written using MPI, OpenMP and Cuda.

The project deals with the comparison of series of letters in the English ABC.
When comparing 2 series of strings, we calculate a "score" of the comparison. The score is called "alignment score". The higher the score, the more similar the series are considered.
There are four signs for similarity of letters. If both letters are the same they are marked as "$". If both letters are not the same but belong to the first group of strings, then they are marked as "%". If they do not belong to the first group of strings but belong to the second group of strings then they are marked as "#". If they dont belong to any group then they are marked as space.
For a series of sequence letters we will define the Mutant Sequence to be marked with MS (k) as the series of letters obtained by adding hyphen after the k-place in a seuqence when k = 1, 2 ...

For a given strings Seq1, Seq2 where Seq2 is shorter, we have to find the offset and mutation level for which the alignment score will be maximum.
Parallel version should work faster than serial version on **one computer.** (if it was for two computers, then implementation was different. we consider here one GPU and one CPU as our resources).

### The program

The program creates a number of processes as defined for it in the makefile. Process number 0 reads the data from a text file and sends relevant data to the rest of the processes. The division of the work between the processes is done statically, and if there is a remainder in the division it will be divided equally between the first processes whose number is smaller than the remainder of the division. To reduce the frequency of the use of MPI_Send and MPI_Recieve communication, I decided to use the MPI_Pack functions for process number 0 and MPI_Unpack for the rest of the processes. After sending and receiving the data, process number 0 transfers his works and other data to the GPU and the GPU performs a calculation for all the results. Once the results are obtained, the process opens threads to find the maximum result, offset and mutation level for each string. Later on, process 0 collects the results obtained from the rest of the processes and prints them. In parrallel to process number 0 work, the rest of the processes will find the maximum result, offset and mutation level for each of the works (strings) assigned to them by using threads.

### The algorithm and data structures

- To avoid many unnecessary calculations, I have written down a separate program that creates a 26 by 26 size matrix of integers in the range of 0 to 3 that represent the signs.   I defined this matrix in the main program as a matrix of constants. From the matrix its possible to retrieve the correct sign without using many unnecessary calculations. 
  Also to keep things ordered, each result is kept in a structure that includes three fields: score, offset, mutant_lvl. 
- The algorithm compares the first string in the program with a second type string. In the outer loop the level of the mutation is determined, starting from 1.   
  In the inner loop at the first level the offset is determined starting from 0. In the inner loop at the second level a comparison of each pair of characters is made and         their value is added to the current result. At the end of the run of the second level internal loop, we get the result for a particular mutation level and offset. Then a check   is performed for whether the result obtained is greater than the current maximum, if so - the current maximum result is updated and also the level of mutation and the offset     that yields the maximum result.

### Parralization

- The rest of the processes: To the algorithm described above I added another external loop which uses threads. Each thread works on a different string depending on the index     assigned to it from the system. This way the work will be divided and carried out by several threads instead of one single thread as happens in the serial version. 
- Process 0: First turns to the GPU. For each string a grid is opened in two-dimensional topology where the x-axis represents the offset number and the y-axis represents the       mutant_lvl + 1. In each block a result of offset and mutant_lvl is calculated according to its identification in the grid. The number of threads in each block is based on the   number of character comparisons to be performed. If the number of comparisons to be made is above 1024 (maximum number of threads within a block), some threads will perform a   calculation or an additional number of calculations. When the GPU finishes performing the calculations, it sends the results to the processor. A function that works in           parallel recieves the results and obtains for each string the maximum result, the offset and the corresponding mutation level.




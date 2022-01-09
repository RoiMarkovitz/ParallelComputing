#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "cFunctions.h"
#include "mpiUtil.h"

// Submitted and developed by Roi Markovitz


int main(int argc, char *argv[])
{        
    char** work_arr;  // sequence2 strings. the work_arr of every process
    AlignmentScore* best_scores; // each cell represents the best offset, mutant_lvl and score of a sequence 2 string. 
    int num_of_sequences, work_arr_size; // num_of_sequences is the number of sequence 2 strings. work_arr_size is the number of sequence 2 strings of a process.
    int my_rank, num_procs; // my_rank is the number of a process. num_procs is the overall number of processes.
    int weights[NUM_SIGNS] = {0}; // weights to calculate alignment score.
    char* seq1; // sequence1 string. all sequence2 strings will be compared to it
    double t; // time variable to calculate serial and parallel times
    
    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    MPI_Datatype alignment_score_type;
    create_alignment_score_type(&alignment_score_type); 
    
    if(my_rank == ROOT) // its process 0
    // reads information from standard input 
    // gives work to other processes
    // activates CUDA for calculations 
    // writes results to standard output
    { 
      AlignmentScore* score_array; // scores of all seq2 strings
      int* work_array_sizes, *score_array_sizes; 
      int score_array_total_size = 0;
      if(readInputFromFile(&work_arr, &best_scores, weights, &seq1, &num_of_sequences, &work_array_sizes, &work_arr_size, num_procs, my_rank))
        MPI_Abort(MPI_COMM_WORLD, 1);      

      rootSendGeneralDataToAllOtherProcesses(num_of_sequences, num_procs, seq1, weights);
        
      rootSendSeq2ToAllOtherProcesses(work_arr_size, num_procs, work_array_sizes, work_arr);

      if(initAlignmentScoreArray(&score_array, &score_array_sizes, work_arr, seq1, work_arr_size, &score_array_total_size))
        MPI_Abort(MPI_COMM_WORLD, 1); 
     
      t = MPI_Wtime(); // Parallel code time start
    
     // start cuda  
      if(startCudaAlignmentScoreCalculation(work_arr, work_arr_size, seq1, weights, score_array, score_array_sizes, score_array_total_size))
              MPI_Abort(MPI_COMM_WORLD, 1);

      findBestScoresParallel(best_scores, score_array, score_array_sizes, work_arr_size); 
                                                                
      // recieve best scores from processes      
      int offset = work_arr_size;
      for (int worker_id = 1; worker_id < num_procs; worker_id++)   
      {
           MPI_Recv(best_scores + offset, work_array_sizes[worker_id], alignment_score_type, worker_id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
           offset += work_array_sizes[worker_id];
      }     
                  
      printf("Parallel run time is %lf\n", MPI_Wtime()-t);
      printResults(best_scores, num_of_sequences); 
               
      // zero best_scores before serial run
      zeroAlignmentScores(best_scores, num_of_sequences);
               
      t = MPI_Wtime(); // Serial code time start            
      for(int i=0; i<num_of_sequences; i++)
         serialAlignmentScoreAlgorithm(seq1, work_arr[i], &best_scores[i], weights); 
      printf("Serial run time is %lf\n", MPI_Wtime()-t);
      printResults(best_scores, num_of_sequences); 

      freeRootMemory(work_arr, num_of_sequences, seq1, best_scores, work_array_sizes, score_array, score_array_sizes);  

    } // end if process 0  
    else // the other processes
    {    
      
      if(processRecieveGeneralDataFromRoot(&num_of_sequences, &seq1, weights))
        MPI_Abort(MPI_COMM_WORLD, 1);
   
      calculateNumWorks(&work_arr_size, num_of_sequences, num_procs, my_rank);
         
      if(initWorkerProcessesArrays(&work_arr, work_arr_size, &best_scores))
        MPI_Abort(MPI_COMM_WORLD, 1);
      
      if(processRecieveSeq2DataFromRoot(work_arr_size, work_arr))
        MPI_Abort(MPI_COMM_WORLD, 1);
        
      parallelThreadsAlignmentScoreAlgorithm(seq1, work_arr, best_scores, weights, work_arr_size);
       
      MPI_Send(best_scores, work_arr_size, alignment_score_type, ROOT, 0, MPI_COMM_WORLD); 
  
      freeWorkerMemory(work_arr, work_arr_size, seq1, best_scores);
    } // end if other processes
                   
    MPI_Finalize();

    return 0;
}

 
   
   
   



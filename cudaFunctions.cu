#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <iostream>
#include "cFunctions.h"

                             
__global__ void calculateAlignmentScores(char* dev_seq1, char* dev_seq2, int seq2_len, int* dev_weights, AlignmentScore* dev_score_array, int current_location_in_score_array, int* dev_signs_match_matrix, int max_threads_per_block)
{    
    // each block calculates a score for a pair of offset and mutant_lvl 
    int score_index = blockIdx.y * gridDim.x + blockIdx.x + current_location_in_score_array; 
    int thread_local_id = threadIdx.x;
    int skip_letter = 0;
    int cycle = 0; // number of cycles of max_threads_per_block
    int i,j;
     
    // max num of compares is seq2_len 
    if (thread_local_id < seq2_len) 
    { 	
        do
        {
            // in the grid -> blockIdx.y + 1 = mutant_lvl = k
		    if (thread_local_id >= blockIdx.y + 1) // ignore "-" comparison                
			    skip_letter = 1;

            i = dev_seq1[thread_local_id+blockIdx.x+skip_letter] - 'A'; 
            j = dev_seq2[thread_local_id] - 'A';
            atomicAdd(&dev_score_array[score_index].score, dev_weights[dev_signs_match_matrix[i*NUM_LETTERS + j]]);
            cycle++;
            thread_local_id = threadIdx.x + cycle*max_threads_per_block;
        } while (thread_local_id < seq2_len); // check if thread should work again     
    }
             
    if(threadIdx.x == 0) 
    {
        dev_score_array[score_index].offset = blockIdx.x; 
        dev_score_array[score_index].mutant_lvl = blockIdx.y + 1; // blockIdx.y + 1 = mutant_lvl of the block in the grid. 
    }

}

int checkStatus(cudaError_t cudaStatus, char* dev_seq1, char* dev_seq2, AlignmentScore* dev_score_array, int* dev_weights, int* dev_signs_match_matrix, std::string err)
{
    if(cudaStatus != cudaSuccess)
    {
        std::cout << err <<std::endl;

        if (dev_seq1 != NULL)
            free(dev_seq1);
        if (dev_seq2 != NULL)      
            free(dev_seq2);
        if (dev_score_array != NULL) 
            free(dev_score_array);
        free(dev_weights);
        free(dev_signs_match_matrix);
           
        return 1;
    }
    return 0; // all is fine
}

 
int startCudaAlignmentScoreCalculation(char** work_arr, int work_arr_size, char* seq1, int* weights, AlignmentScore* score_array, int* score_array_sizes, int score_array_total_size)
{
    char* dev_seq1 = NULL;  
    char* dev_seq2 = NULL; 
    AlignmentScore* dev_score_array = NULL;
    int seq1_len, seq2_len;
    int* dev_weights = 0;
    int* dev_signs_match_matrix = 0;
    int size, offsets, mutants;
    int max_threads_per_block;
    int current_location_in_score_array = 0;
  
    cudaError_t cudaStatus;

    cudaDeviceProp prop;
    cudaStatus = cudaGetDeviceProperties(&prop,0);

    // SIGNS_MATCH_MATRIX allocation and copy from host to device
    size = NUM_LETTERS * NUM_LETTERS * sizeof(int);
    cudaStatus = cudaMalloc((void**)&dev_signs_match_matrix, size);
    if(checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda malloc for dev_signs_match_matrix failed!"))
        return 1;
    cudaStatus = cudaMemcpy(dev_signs_match_matrix, SIGNS_MATCH_MATRIX, size, cudaMemcpyHostToDevice);
    if(checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda memcpy SIGNS_MATCH_MATRIX to device failed!"))
        return 1; 
    
    // weights memory allocation and copy from host to device
    size = sizeof(int)*NUM_SIGNS;
    cudaStatus = cudaMalloc((void**)&dev_weights, size);
    if(checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda malloc for dev_weights failed!"))
        return 1;
    cudaStatus = cudaMemcpy(dev_weights, weights, size, cudaMemcpyHostToDevice);
    if(checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda memcpy weights to device failed!"))
        return 1;

    // seq1 memory allocation and copy from host to device
    seq1_len = strlen(seq1);
    size = sizeof(char)* (seq1_len + 1);
    cudaStatus = cudaMalloc((void**)&dev_seq1, size);
    if(checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda malloc for dev_seq1 failed!"))
        return 1;
    cudaStatus = cudaMemcpy(dev_seq1, seq1, size, cudaMemcpyHostToDevice);
    if(checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda memcpy seq1 to device failed!"))
        return 1;

    // dev_score_array memory allocation and copy from host to device
    size = sizeof(AlignmentScore) * score_array_total_size;
    cudaStatus = cudaMalloc((void**)&dev_score_array, size);
    if(checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda malloc for dev_score_array failed!"))
        return 1;
    cudaStatus = cudaMemcpy(dev_score_array, score_array, size, cudaMemcpyHostToDevice);
    if(checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda memcpy score_array to device failed!"))
        return 1;
    
    // loops over seq2 strings
    for (int i = 0; i < work_arr_size; i++)
    {
        seq2_len = strlen(work_arr[i]); 
        size = sizeof(char)* (seq2_len + 1);   
        cudaStatus = cudaMalloc((void**)&dev_seq2, size);
        if(checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda malloc for dev_seq2 failed!"))
            return 1;  
        cudaStatus = cudaMemcpy(dev_seq2, work_arr[i], size, cudaMemcpyHostToDevice);
        if(checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda memcpy work_arr[i] to device failed!"))
            return 1;
    
        offsets = seq1_len - seq2_len;
        mutants = seq2_len; 
        dim3 numBlocks (offsets, mutants);
        max_threads_per_block = prop.maxThreadsPerBlock < mutants ? prop.maxThreadsPerBlock : mutants;
      
        // activate kernel
        calculateAlignmentScores<<<numBlocks,max_threads_per_block>>>(dev_seq1, dev_seq2, seq2_len, dev_weights, dev_score_array, current_location_in_score_array, dev_signs_match_matrix, max_threads_per_block);
        cudaStatus = cudaDeviceSynchronize();
        if(checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda Kernel for calculateAlignmentScores failed!"))
            return 1;
        
        current_location_in_score_array += score_array_sizes[i];
    }
  
    // copy data from device to host
    size = sizeof(AlignmentScore) * score_array_total_size;
    cudaStatus = cudaMemcpy(score_array, dev_score_array, size, cudaMemcpyDeviceToHost);
    if(checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda memcpy dev_score_array to host failed!"))
        return 1;

	// free cude memory   
    cudaStatus = cudaFree(dev_signs_match_matrix);
    if (checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda free dev_signs_match_matrix failed!"))
       return 1;
    cudaStatus = cudaFree(dev_weights);
    if (checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda free dev_weights failed!"))
        return 1;
    cudaStatus = cudaFree(dev_seq1);
    if (checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda free dev_seq1 failed!"))
        return 1;    
    cudaStatus = cudaFree(dev_seq2);
    if (checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda free dev_seq2 failed!"))
        return 1;   
    cudaStatus = cudaFree(dev_score_array);
    if (checkStatus(cudaStatus, dev_seq1, dev_seq2, dev_score_array, dev_weights, dev_signs_match_matrix, "Cuda free dev_score_array failed!"))
        return 1;
 
    return 0; // all is fine
}



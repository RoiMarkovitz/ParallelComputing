#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "mpiUtil.h"
#include <stddef.h>
#include "cFunctions.h"


char buff[BUFF_SIZE];
int position; // tracks buffer filled size, current position in buffer


void create_alignment_score_type(MPI_Datatype* alignment_score_type)
{
    int block_sizes[NUM_OF_ATTRIBUTES] = {1, 1, 1};
    MPI_Aint disp[NUM_OF_ATTRIBUTES];
    MPI_Datatype types[NUM_OF_ATTRIBUTES] = {MPI_INT, MPI_INT, MPI_INT};
 
    disp[0] = offsetof(AlignmentScore, score);
    disp[1] = offsetof(AlignmentScore, offset);
    disp[2] = offsetof(AlignmentScore, mutant_lvl);
   
    MPI_Type_create_struct(NUM_OF_ATTRIBUTES, block_sizes, disp ,types, alignment_score_type);
    MPI_Type_commit(alignment_score_type);
}

void rootSendGeneralDataToAllOtherProcesses(int num_of_sequences, int num_procs, char* seq1, int* weights)
{  
    int lenseq1 = strlen(seq1) + 1;
    position = 0;
    MPI_Pack(&num_of_sequences, 1, MPI_INT, buff, BUFF_SIZE, &position, MPI_COMM_WORLD);     
    MPI_Pack(&lenseq1, 1, MPI_INT, buff, BUFF_SIZE, &position, MPI_COMM_WORLD);
    MPI_Pack(seq1, lenseq1, MPI_CHAR, buff, BUFF_SIZE, &position, MPI_COMM_WORLD);
    MPI_Pack(weights, NUM_SIGNS, MPI_INT, buff, BUFF_SIZE, &position, MPI_COMM_WORLD);

    for (int worker_id = 1; worker_id < num_procs; worker_id++)
        MPI_Send(buff, position, MPI_PACKED, worker_id, 0, MPI_COMM_WORLD);
}

int processRecieveGeneralDataFromRoot(int* num_of_sequences, char** seq1, int* weights)
{
    int lenseq1;
    position = 0;
    MPI_Recv(buff, BUFF_SIZE, MPI_PACKED, ROOT, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Unpack(buff, BUFF_SIZE, &position, num_of_sequences, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buff, BUFF_SIZE, &position, &lenseq1, 1, MPI_INT, MPI_COMM_WORLD);
    *seq1 = (char*)malloc(sizeof(char)* lenseq1);
    if(!(*seq1))
    {
        fprintf(stderr, "seq1 allocation failed\n");  
        return 1;
    }
    MPI_Unpack(buff, BUFF_SIZE, &position, (*seq1), lenseq1, MPI_CHAR, MPI_COMM_WORLD);
    MPI_Unpack(buff, BUFF_SIZE, &position, weights, NUM_SIGNS, MPI_INT, MPI_COMM_WORLD); 
    return 0;
}

void rootSendSeq2ToAllOtherProcesses(int work_arr_size, int num_procs, int* work_array_sizes, char** work_arr)
{
    int k = work_arr_size; // so process 0 dedicated works wont be sent to other processes
    int lenseq2;    
    for (int worker_id = 1; worker_id < num_procs; worker_id++)
    {
        position = 0; 
        for (int j = 0; j < work_array_sizes[worker_id]; j++) 
        {
            lenseq2 = strlen(work_arr[k]) + 1;
            MPI_Pack(&lenseq2, 1, MPI_INT, buff, BUFF_SIZE, &position, MPI_COMM_WORLD);
            MPI_Pack(work_arr[k], lenseq2, MPI_CHAR, buff, BUFF_SIZE, &position, MPI_COMM_WORLD);
            k++;       
        }
        if (work_array_sizes[worker_id] > 0)
            MPI_Send(buff, position, MPI_PACKED, worker_id, 0, MPI_COMM_WORLD);                     
    }
}

int processRecieveSeq2DataFromRoot(int work_arr_size, char** work_arr)
{    
    position = 0;
    int lenseq2;   
    if(work_arr_size > 0)
    {
        MPI_Recv(buff, BUFF_SIZE, MPI_PACKED, ROOT, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < work_arr_size; i++)
        {
            MPI_Unpack(buff, BUFF_SIZE, &position, &lenseq2, 1, MPI_INT, MPI_COMM_WORLD);
            work_arr[i] = (char*)malloc(sizeof(char) * lenseq2);
            if(!work_arr[i])
            {
                fprintf(stderr, "seq2 string allocation failed\n");  
                return 1;
            }
            MPI_Unpack(buff, BUFF_SIZE, &position, work_arr[i], lenseq2, MPI_CHAR, MPI_COMM_WORLD);         
        }                     
    }
    return 0;
}








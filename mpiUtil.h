#pragma once

#define ROOT 0

void create_alignment_score_type(MPI_Datatype* alignment_score_type);
void rootSendGeneralDataToAllOtherProcesses(int num_of_sequences, int num_procs, char* seq1, int* weights);
int processRecieveGeneralDataFromRoot(int* num_of_sequences, char** seq1, int* weights);
void rootSendSeq2ToAllOtherProcesses(int work_arr_size, int num_procs, int* work_array_sizes, char** work_arr);
int processRecieveSeq2DataFromRoot(int work_arr_size, char** work_arr);



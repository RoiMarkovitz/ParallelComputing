#include "cFunctions.h"
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


int readInputFromFile(char*** work_arr, AlignmentScore** best_scores, int* weights, char** seq1, int* num_of_sequences, int** work_array_sizes, int* work_arr_size, int num_procs, int my_rank)
{
      char input_buf[BUFF_SIZE]; 
      char temp_seq[MAX_LETTERS_SEQ1];    
        
      fgets(input_buf, BUFF_SIZE, stdin); 
      int num_elements_returned = sscanf(input_buf, "%d %d %d %d", &weights[0], &weights[1], &weights[2], &weights[3]);
    
      for(int i=0 ; i<NUM_SIGNS; i++)
      {
        if (weights[i] < 0)
        {
          fprintf(stderr, "All weights must be above or equal to 0\n");  
          return 1;
        }
      }

      for(int i=1; i<NUM_SIGNS; i++)
          weights[i] *=-1;
        
      scanf("%s", temp_seq);
      
      *seq1 = (char*) malloc(sizeof(char) * (strlen(temp_seq) + 1));
      if(!(*seq1))
      {
        fprintf(stderr, "seq1 allocation failed\n");  
        return 1;
      }
      strcpy((*seq1), temp_seq);
      captalizeLetters(*seq1);

      scanf("%d", num_of_sequences);
      
      if (*num_of_sequences < 1)
      {
          fprintf(stderr, "Must have at least one sequence to compare with\n");  
          return 1;
      }

      if(initBestScores(best_scores, *num_of_sequences))
        return 1;
      
      *work_array_sizes = (int*)malloc(sizeof(int) * num_procs);
      if(!(*work_array_sizes))
      {
        fprintf(stderr, "work_array_sizes allocation failed\n");  
        return 1;
      }
      calculateArrayWorksSizes(*work_array_sizes, *num_of_sequences, num_procs);
      *work_arr_size = (*work_array_sizes)[0];
       
      *work_arr = (char**)malloc(sizeof(char*) * (*num_of_sequences));  
      if(!(*work_arr))
      {
        fprintf(stderr, "work_arr allocation failed\n");  
        return 1;
      }
        
      for(int i=0; i<*num_of_sequences; i++)
      {
          scanf("%s", temp_seq);      
          (*work_arr)[i] = (char*)malloc(sizeof(char) * (strlen(temp_seq) + 1));
           if(!(*work_arr)[i])
          {
            fprintf(stderr, "seq2 string allocation failed\n");  
            return 1;
          }
                 
          strcpy((*work_arr)[i], temp_seq);
          captalizeLetters((*work_arr)[i]);
            
          if (strlen((*work_arr)[i]) >= strlen((*seq1)))
          {
            fprintf(stderr, "Sequence to compare with must be shorter than the first sequence\n");  
            return 1;
          }        
      }

      return 0; // everything is fine
}

void captalizeLetters(char* string)
{
    int i = 0;
    
    while (string[i]) 
    { 
        string[i] = toupper(string[i]);
        i++;
    }
}

// used by every process to calculate its work_arr_size
void calculateNumWorks(int* work_arr_size, int num_of_sequences, int num_procs, int my_rank)
{
    if (num_of_sequences >= num_procs)
    {
        *work_arr_size = num_of_sequences / num_procs; 
        int remains = num_of_sequences % num_procs;
        if (remains != 0)
          if(my_rank < remains)
             *work_arr_size+=1;
    }
    else // there are less sequences than num_procs
    {
      if (my_rank < num_of_sequences)
         *work_arr_size = 1;
      else
         *work_arr_size = 0;
    }
}

// caulcautes the work array size of every process for future useage of process 0
void calculateArrayWorksSizes(int* work_array_sizes, int num_of_sequences, int num_procs)
{
    for (int my_rank = 0; my_rank < num_procs; my_rank++)
      calculateNumWorks(&work_array_sizes[my_rank], num_of_sequences, num_procs, my_rank);    
}

int initAlignmentScoreArray(AlignmentScore** score_array, int** score_array_sizes, char** work_arr, char* seq1, int work_arr_size, int* score_array_total_size)
{   
    *score_array_sizes = (int*)malloc(sizeof(int) * work_arr_size);
    if(!(*score_array_sizes))
    {
        fprintf(stderr, "score_array_sizes allocation failed\n");  
        return 1;
    }
    int num_mutants;
    int num_offsets;
    int num_scores_in_seq;
    for (int i = 0; i < work_arr_size; i++)
    {
      num_mutants = strlen(work_arr[i]); // num of mutants is the length of seq2
      num_offsets = strlen(seq1) - num_mutants; // num of offsets in each mutant is (length of seq1 - length of seq2)
      num_scores_in_seq = num_mutants * num_offsets;  // num of muntants * num of offsets is the number of scores calculated in a single seq2    
      *score_array_total_size += num_scores_in_seq;
      (*score_array_sizes)[i] = num_scores_in_seq;  
    }

    *score_array = (AlignmentScore*)malloc(sizeof(AlignmentScore) * (*score_array_total_size));
    if(!(*score_array))
    {
        fprintf(stderr, "score_array allocation failed\n");  
        return 1;
    }  
    zeroAlignmentScores(*score_array, *score_array_total_size);

    return 0; // all is fine
}

int initWorkerProcessesArrays(char*** work_arr, int work_arr_size, AlignmentScore** best_scores)
{
    *work_arr = (char**)malloc(sizeof(char*) * work_arr_size);
    if(!(*work_arr))
    {
        fprintf(stderr, "work_arr allocation failed\n");  
        return 1;
    }

    if(initBestScores(best_scores, work_arr_size))
        return 1;
      
    return 0;
}

int initBestScores(AlignmentScore** best_scores, int work_arr_size)
{
    *best_scores = (AlignmentScore*)malloc(sizeof(AlignmentScore) * work_arr_size);
     if(!(*best_scores))
     {
        fprintf(stderr, "best_scores allocation failed\n");  
        return 1;
     }
     
    zeroAlignmentScores(*best_scores, work_arr_size);
   
    return 0;
}

void zeroAlignmentScores(AlignmentScore* best_scores, int work_arr_size)
{
    for (int i = 0; i < work_arr_size; i++)
    {
        best_scores[i].offset = 0; best_scores[i].mutant_lvl = 0; best_scores[i].score = 0;
    }
}

void serialAlignmentScoreAlgorithm(char* seq1, char* seq2, AlignmentScore* result, int* weights)
{
      int max_alignment_score;   
      int skip_letter; // to ignore "-" comparison 
      int str2len = strlen(seq2);
      int str1len = strlen(seq1);
    
      for(int mutant_lvl = 1; mutant_lvl <= str2len; mutant_lvl++)  // mutations loop   
      {    
        for(int offset = 0; offset < str1len - str2len; offset++) // offsets loop
        {
            skip_letter = 0;
                        
            for(int seq1_index = 0; seq1_index < str2len; seq1_index++) // chars match loop 
            {
              if (seq1_index == mutant_lvl) // ignore "-" comparison                
                skip_letter = 1;                          

              result->score += weights[SIGNS_MATCH_MATRIX[seq1[seq1_index+offset+skip_letter] - 'A'][seq2[seq1_index] - 'A']];        
                                                                          
            } // end for char to char comparison between seq1 and seq2
                
            if(mutant_lvl == 1 && offset == 0)
            {
              max_alignment_score = result->score;
              result->offset = offset;
              result->mutant_lvl = mutant_lvl; 
            }
                    
            if (result->score > max_alignment_score)
            {
              max_alignment_score = result->score;
              result->offset = offset;
              result->mutant_lvl = mutant_lvl;        
            }

            result->score = 0;
                        
        } // end for offsets
      } // end for mutations
      result->score = max_alignment_score;      
}

void parallelThreadsAlignmentScoreAlgorithm(char* seq1, char** work_arr, AlignmentScore* results, int* weights, int work_arr_size)
{   
      #pragma omp parallel for default(none) shared(seq1, work_arr, results, weights, work_arr_size) 
      for (int i = 0; i < work_arr_size; i++) // sequence2 strings loop   
        serialAlignmentScoreAlgorithm(seq1, work_arr[i], &results[i], weights);              
}

void findBestScoresParallel(AlignmentScore* best_scores, AlignmentScore* score_array, int* score_array_sizes, int work_arr_size)
{
      int max_alignment_score;
      int start = 0; int end = 0;
       
      #pragma omp parallel for default(none) shared(best_scores, score_array, score_array_sizes, work_arr_size) private(max_alignment_score) firstprivate(start, end) 
      for (int i = 0; i < work_arr_size; i++)
      {     
        for (int j = 0; j < i; j++)
          start += score_array_sizes[j]; // location where scores of seq2[i] start
           
        max_alignment_score = score_array[start].score; // the first score in the seq2[i] is currently the max
        updateBestScores(&best_scores[i], max_alignment_score, score_array[start].offset, score_array[start].mutant_lvl);
        end = start + score_array_sizes[i];
        start++;
        for (; start < end; start++)
        {
             if (score_array[start].score > max_alignment_score)
             {        
                max_alignment_score = score_array[start].score;
                updateBestScores(&best_scores[i], max_alignment_score, score_array[start].offset, score_array[start].mutant_lvl); 
             }    
        }  
        start = 0;
        end = 0;        
      }        
}

void updateBestScores(AlignmentScore* best_scores, int score, int offset, int mutant_lvl)
{
    best_scores->score = score;
    best_scores->offset = offset;
    best_scores->mutant_lvl = mutant_lvl;
}

void printResults(AlignmentScore* all_results, int size)
{
    for(int i=0; i< size; i++)
        printf("n = %d\t k = %d\n", all_results[i].offset, all_results[i].mutant_lvl);

   // printf("n = %d\t k = %d\t score = %d\n", all_results[i].offset, all_results[i].mutant_lvl, all_results[i].score);   
}

void freeWorkerMemory(char** work_arr, int work_arr_size, char* seq1, AlignmentScore* best_scores)
{
    for (int i = 0; i < work_arr_size; i++)
        free(work_arr[i]);
   
    free(best_scores);
    free(seq1);
    free(work_arr);
}

void freeRootMemory(char** work_arr, int work_arr_size, char* seq1, AlignmentScore* best_scores, int* work_array_sizes, AlignmentScore* score_array, int* score_array_sizes)
{
    for (int i = 0; i < work_arr_size; i++)
        free(work_arr[i]);
   
    free(best_scores);
    free(seq1);
    free(work_arr);
    free(work_array_sizes);
    free(score_array);
    free(score_array_sizes);
}


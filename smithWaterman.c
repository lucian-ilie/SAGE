#include <string.h>
#include<stdio.h>
#include<stdlib.h>
#include "sage.h"

double similarityScore(char a,char b);
double findArrayMaximum(double array[],int length);
int checkAlignment(double matchScore,double mismatchPenalty, double gapPenalty, char *sequence1, char *sequence2);

int ind;
double mismatchPenaltyGlobal,matchScoreGlobal, gapPenaltyGlobal;

/*int main()
{
	char *string_1,*string_2;
	string_1="ACACCCAC";
	string_2="CCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCACCCCCCCCCCCCCCCCACCCCC";
	int score=check_alignment(2,-10,-2, string_1,string_2);
	printf("String lengths: %d, %d\n",strlen(string_1),strlen(string_2));
	printf("Score %d\n",score);
}*/

int checkAlignment(double matchScore, double mismatchPenalty, double gapPenalty, char *sequence1,char *sequence2)
{
	int return_value;
	if(strlen(sequence1)<2)
	{
		sequence1=NULL;
	}
	if(strlen(sequence2)<2)
	{
		sequence2=NULL;
	}
	if(sequence1==NULL && sequence2==NULL)
	{
		return_value=0;
		return return_value;
	}
	if(sequence1==NULL && sequence2!=NULL)
	{
		return_value=strlen(sequence2)*gapPenalty;
		return return_value;
	}
	if(sequence1!=NULL && sequence2==NULL)
	{
		return_value=strlen(sequence1)*gapPenalty;
		return return_value;
	}	
	int i,j;		
	mismatchPenaltyGlobal=mismatchPenalty;
	matchScoreGlobal=matchScore;
	gapPenaltyGlobal=gapPenalty;
	char *seq_a,*seq_b; 
	seq_a=sequence1;  
	seq_b=sequence2;
	int N_a=strlen(seq_a);
	int N_b=strlen(seq_b);
	double** H;
	int dimension1_max=N_a+1;
	int dimension2_max=N_b+1;
	H = (double**)malloc(dimension1_max*sizeof(double*));
	for (i=0;i<dimension1_max;i++) 
	{
  		H[i]=(double*)malloc(dimension2_max*sizeof(double));
	}   
	for(i=0;i<dimension1_max;i++)
	{
		for(j=0;j<dimension2_max;j++)
		{
			H[i][j]=0.;
		}
	}
	for(i=1;i<dimension1_max;i++)
	{
		H[i][0]=i*gapPenaltyGlobal;
	}
	for(j=1;j<dimension2_max;j++)
	{
		H[0][j]=j*gapPenaltyGlobal;
	}
	double temp[4];
	int **I_i,**I_j;     // Index matrices to remember the 'path' for backtracking
	I_i = (int**)malloc(dimension1_max*sizeof(int*));
	for (i=0;i<dimension1_max;i++) 
	{
  		I_i[i]=(int*)malloc(dimension2_max*sizeof(int));
	}
	I_j = (int**)malloc(dimension1_max*sizeof(int*));
	for (i=0;i<dimension1_max;i++) 
	{
  		I_j[i]=(int*)malloc(dimension2_max*sizeof(int));
	}
	for(i=0;i<dimension1_max;i++)
	{
		for(j=0;j<dimension2_max;j++)
		{
			I_i[i][j]=0;
			I_j[i][j]=0;
	    	}
	}
  	// here comes the actual algorithm
	for(i=1;i<=N_a;i++)
	{
		for(j=1;j<=N_b;j++)
		{
			temp[0] = H[i-1][j-1]+similarityScore(seq_a[i-1],seq_b[j-1]); 
			temp[1] = H[i-1][j]+gapPenaltyGlobal;
			temp[2] = H[i][j-1]+gapPenaltyGlobal;
			temp[3] = -1000000.;
			H[i][j] = findArrayMaximum(temp,4);
			switch(ind)
			{
				case 0:// score in (i,j) stems from a match/mismatch
					I_i[i][j] = i-1;
					I_j[i][j] = j-1;
					break;
				case 1:// score in (i,j) stems from a deletion in sequence A
					I_i[i][j] = i-1;
					I_j[i][j] = j;
					break;
				case 2:// score in (i,j) stems from a deletion in sequence B
					I_i[i][j] = i;
					I_j[i][j] = j-1;
					break;
				case 3:// (i,j) is the beginning of a subsequence
					I_i[i][j] = i;
					I_j[i][j] = j;	
					break;
			}
		}
	}
	// Print the matrix H to the console
/*	printf("**********************************************\n");;
	printf("The scoring matrix is given by  \n\n");
	for(i=0;i<=N_a;i++)
	{
		for(j=0;j<=N_b;j++)
		{
			printf("%2.0f ",H[i][j]);
		}
		printf("\n");
	}

	printf("\n\n");
	for(i=0;i<=N_a;i++)
	{
		for(j=0;j<=N_b;j++)
		{
			printf("(%2d,%2d) ",I_i[i][j], I_j[i][j]);
		}
		printf("\n");
	}
*/


//	int H_max=H[N_a][N_b];
	int i_max=N_a;
	int j_max=N_b;
	// Backtracking from H_max
	int current_i=i_max,current_j=j_max;
	int next_i=I_i[current_i][current_j];
	int next_j=I_j[current_i][current_j];
	int tick=0;
	//char consensus_a[N_a+N_b+2],consensus_b[N_a+N_b+2];
	char *consensus_a,*consensus_b;
	if((consensus_a = malloc((N_a+N_b+2)*sizeof(consensus_a)))==NULL)
	{
		printError(102, "consensus_a");
		exit(EXIT_FAILURE);
	}
	if((consensus_b = malloc((N_a+N_b+2)*sizeof(consensus_b)))==NULL)
	{
		printError(102, "consensus_b");
		exit(EXIT_FAILURE);
	}
	while(1)
	{
		//printf("%d, %d, %d, %d\n",current_i,current_j, next_i,next_j);
		if(next_i==current_i)
			consensus_a[tick] = '-';// deletion in A
		else
			consensus_a[tick] = seq_a[current_i-1];   // match/mismatch in A
		if(next_j==current_j)
			consensus_b[tick] = '-';                  // deletion in B
		else
			consensus_b[tick] = seq_b[current_j-1];   // match/mismatch in B
		if(current_i==0 && current_j==0)
			break;		
		else if(current_i==0 && current_j!=0)
			current_j--;
		else if(current_i!=0 && current_j==0)
			current_i--;
		else
		{
			current_i = next_i;
			current_j = next_j;
		}
		next_i = I_i[current_i][current_j];
		next_j = I_j[current_i][current_j];
		tick++;
	}
	//if((fpo_output_file=fopen(output_text_file_name,"a"))==NULL) printError(OPEN_FILE,output_text_file_name);
	//fprintf(fpo_output_file,"%s\n",seq_a);
	//fprintf(fpo_output_file,"%s\n\n",seq_b);
	//for(i=tick-1;i>=0;i--)
	//{
	//	fprintf(fpo_output_file,"%c",consensus_a[i]);
	//}
	//fprintf(fpo_output_file,"\n");
	//for(i=tick-1;i>=0;i--)
	//{
	//	if(consensus_a[i]==consensus_b[i])
	//		fprintf(fpo_output_file,"|");
	//	else
	//		fprintf(fpo_output_file," ");
	//}
	//fprintf(fpo_output_file,"\n");
	//for(i=tick-1;i>=0;i--)
	//{
	//	fprintf(fpo_output_file,"%c",consensus_b[i]);
	//}
	//fprintf(fpo_output_file,"\n");
	int match=0,mismatch=0,gap=0;
	for(i=tick-1;i>=0;i--)
	{
		if(consensus_a[i]==consensus_b[i])
			match++;
		else if(consensus_a[i]=='-' || consensus_b[i]=='-')
			gap++;	
		else if(consensus_a[i]!=consensus_b[i])
			mismatch++;		
	}
	//fprintf(fpo_output_file,"Match = %d, Mismatch = %d Gap = %d\n",match,mismatch,gap); 
	return_value=H[N_a][N_b];
	for (i=0;i<dimension1_max;i++) //free all memory
	{
  		free(H[i]);
		free(I_i[i]);
		free(I_j[i]);
	}
	free(H);
	free(I_i);
	free(I_j);
	free(consensus_a);
	free(consensus_b);
	//fclose(fpo_output_file);
	return return_value;
}


double similarityScore(char a,char b)
{
	double result;
	if(a==b)
	{
		result=matchScoreGlobal;
	}
	else
	{
		result=mismatchPenaltyGlobal;
	}
	return result;
}

double findArrayMaximum(double array[],int length)
{
	double max = array[0];            // start with max = first element
	int i;
	ind=0;
	for(i = 1; i<length; i++)
	{
		if(array[i] > max)
		{
			max = array[i];
			ind = i; 
		}
	}
	return max;                    // return highest value in array
}


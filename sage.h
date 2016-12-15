#define contigLengthThreshold 100 // minimum length of string to report as a contig
#define P 0.001 // This is the probability to delete simple edges with long delta.
#define minDelta 1000// min delta to rely on a-statistics value

unsigned int readLength;//length of the reads
unsigned int minOverlap; // minimum overlap length
unsigned int *Mean,*standardDeviation;
char outputTextFileName[10000];
FILE *fpo_output_file; 

void printError(int err_no, char *mssg);
int checkAlignment(double matchScore, double mismatchPenalty, double gapPenalty, char *sequence1,char *sequence2);
void myfprintf(char * format, ...);

/* ============================================================================================
   Error numbers;                                                                   
   ============================================================================================ */
#define OPEN_FILE 101 // Error opening file
#define MEM_ALLOC 102 // Error allocating memory
#define OTHERS 200 // Other errors

/* ============================================================================================
	sage.c version 1
	Written by: Md. Bahlul Haider, Phd. Student, the University of Western Ontario
	Last updated: 30 July 2012

	This program does the following tasks.

	=> Builds a hash table of all 28bp prefix and suffix of the reads and their reverse complements.

	=> Find overlapping reads using hash table and inset edges in the graph.

	=> Remove transitive edges using Gene Meyrs algorithm while building the overlap graph.

	=>Use CS2 to to find flow in the graph.
	Details of CS2 can be found at: "An Efficient Implementation of a Scaling Minimum-Cost Flow
	Algorithm"-by Andrew V. Goldberg, J. Algorithms 22(1): 1-29 (1997). CS2 was downloaded from
	http://www.igsystems.com/cs2/

	NOTES:
	=>The output graphs (.gdl files) of this program can be viewed by AiSee. Demo version of
	AiSee can be downloaded from http://www.aisee.com/
   ============================================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>
#include <unistd.h> //To get the host name
#include <stdint.h> // Defines the data types int8_t, int64_t etc.
#include <float.h>
#include <sys/stat.h> //Defines time_t, st_size etc.
#include <stdarg.h> //Allows functions to accept an indefinite number of arguments
#include <zlib.h> // This is needed for kseq.h
#include "kseq.h" // Header file for reading dataset in blocks.
#include "sage.h" //Header file for sage.

#define DEBUGGING 0 //Print everything in debugging mode. if set to 1 then the program will save all intermediate files in the disk
#define BUILDGRAPH 1 // 1 = build the graph and save it in a file 0 = load graph from file.
#define SAVEGRAPH 0 //1 = save the graph 0 = do not save the graph.
#define TEST 0 // TEST=1 does not print the graph in every steps.


/* ============================================================================================
   This structure stores the graph with the strings in each of the edges.
   ============================================================================================ */
struct edgeListStructure // store list of edges of a node in the overlap graph
{
	uint64_t ID;
	uint64_t from_ID;
	uint8_t typeOfEdge;
	int8_t isReducible;
	uint64_t lengthOfEdge;
	float flow;
	uint64_t *listOfReads;
	struct edgeListStructure *next;
	struct edgeListStructure *previous;
	struct edgeListStructure *twinEdge;
};
typedef struct edgeListStructure Edge;


struct readToEdgeStructure //
{
	uint64_t *locationForward;
	uint64_t *locationReverse;
	Edge *edge;
	struct readToEdgeStructure *next;
};

struct matePairStructure // Store mate pair information
{
	uint64_t ID;
	uint64_t freq;
	uint8_t flag;
	int8_t type1;
	int8_t type2;
	int8_t library;
	struct matePairStructure *next;
};

struct indexes
{
	uint64_t firstNode;
	uint64_t lastNode;
	uint64_t position;
	struct indexes *next;
};

/* ============================================================================================
   Global variables.
   ============================================================================================ */
int flag, numberOfLibrary=0,*support,*frequency, aStatisticsThreshold=3,hashStringLength, iteration,numberOfSameBasepairThreshold;
uint64_t numberOfUniqueReads, numberOfReads,count=0;
uint64_t genomeSize=0, longestContig=0,numberOfContigs=0,N_50=0,N_90=0,N_50_N=0,N_90_N=0,base_pair_covered=0;
extern unsigned int readLength, minOverlap, *Mean, *standardDeviation;
extern char outputTextFileName[10000];
extern FILE *fpo_output_file;
float genomeCovered;
Edge **supportedEdg1,**supportedEdg2,**graph,***allPairPaths;
char directoryName[10000];
unsigned int *upperBoundOfInsert,*lowerBoundOfInsert, minimumUpperBoundOfInsert,maximumUpperBoundOfInsert;
uint64_t *allPairPathLengths,totalPaths,arraySize=1000,pair_of_support=0;
uint64_t longestPathLength=0,defaultArraySize=10000000, readArrayLength, **listOfHashReadsInt,highestPaths=0;
uint8_t **readsInt=NULL,**readsReverseInt=NULL,*markedNodes,*selfOverlappingReads;
uint64_t numberOfHashMiss=0,numberOfHashMatch=0,**graphEconomy;
struct readToEdgeStructure **read_to_edge;
struct matePairStructure **mate_pair;
struct indexes **indx;
uint8_t *exploredReads;
int maxSearch=0;
uint64_t sizeOfHashTable,precomputeHash;




/* ============================================================================================
   Function declarations.
   ============================================================================================ */
double main_cs2(char *input_filename,char * output_filename);
void buildOverlapGraphEconomy(void);
int insertEdge(uint64_t vertex_u,uint64_t vertex_v, int type,uint64_t *listForward, uint64_t *listReverse,float flow,uint64_t delta);
void printGraph(char *graph_file);
int deleteEdge(Edge* edge);
uint64_t contractCompositePaths(void);
uint64_t removeCycles(void);
int findMultiplePaths(Edge *edge);
uint64_t genomeSizeEstimation(void);
void computeMinCostFlow(void); // Convex cost function with a_statistics.
void quicksortReads(uint8_t **reads,int64_t left, int64_t right);
void printError(int err_no, char *mssg);
uint64_t removeC2S1Edges(void);
int64_t getIdOfRead(char * read);
void mapReadsToEdges(void);
void deleteSimpleEdgesWithLongDelta(void);
uint64_t reduceTrees(void);
uint64_t removeTriangles(void);
void readDatasetInBytes(char *filename);
int mergeEdges(Edge *edge1,Edge *edge2, int typeOfMerge);
uint64_t reduceLoops(void);
uint64_t removeCliques(void);
void memoryInitialize(void);
char* readGenome(char* filename);
void itoa(int n, char s[]);
void reverse(char s[]);
uint64_t deleteEdgesWithSimilarStrings(void);
void updateGraph(void);
int findPathsBetweenMatepairs(uint64_t mate_pair_1,uint64_t mate_pair_2,int type1,int type2,int library);
void cleanGraph(void);
int isUniqueInEdge(Edge *edge,uint64_t read);
int matechEdgeType(Edge *edge1,Edge *edge2);
uint64_t insertPath(Edge **edges,int *support_flag, uint64_t top,uint64_t path_found);
int64_t* findDistanceOnEdge(Edge *edge, uint64_t read);
uint64_t resolvePairsOfEdges(Edge **supportedEdge1,Edge **supportedEdge2, uint64_t supported_pairs);
void markSimilarEdges(void);
void savePath(Edge **edges, int distance,uint64_t top);
uint64_t exploreGraph(uint64_t node,int level);
uint64_t findAllPairPaths(void);
void mapReadLocations(void);
uint64_t partitionPaths(int64_t left, int64_t right);
void quicksortPaths(int64_t left, int64_t right);
void insertIntoList(Edge *v);
int reverseEdgeType(int type);
Edge* whichEdgeToDelete(Edge *uv,Edge *uw,Edge *vw);
int checkFlow(Edge **edges, uint64_t top, Edge *next_edge);
uint64_t startingPosition(uint64_t firstReadOnEdge,uint64_t second_read);
void mapMatePairs(char *filename,int library);
void reduceTransitiveEdges(Edge *edge);
void freeAllocatedMemory(void);
int combinedEdgeType(Edge *edge1, Edge* edge2);
int checkIfAllReadsPresent(void);
void mergeContigs(void);
uint64_t mergeFinal(void);
uint64_t removeDeadEnds(void);
uint64_t removeBubbles(void);
void deleteDelta(void);
void meanSdEstimation(void);
void computeMeanSD(int mu,int sd,int library, int *returnMu,int *returnSD);
int deleteTemporaryFiles(void);
void getDegreeStatistics(void);
void indexPaths(uint64_t num);
void freeIndexedPaths(void);
void markTransitiveEdge(uint64_t from_ID);
void sortEdgesEconomy(int type);
void quicksortEdgesEconomy(uint64_t *listForSorting,int64_t left, int64_t right,int type);
uint8_t *charsToBytes(char *read);
uint64_t hashTableInsert(uint64_t *value,uint64_t readNumber,uint8_t type);
uint64_t hashTableSearch(uint64_t *value);
void hashPrefixesAndSuffix(void);
int insertEdgeEconomy(uint64_t vertex_u,uint64_t vertex_v, int delta, int type);
void freeHashTable(void);
int64_t findNextPrimeFromFile(int64_t value);
void convertGraph(void);
uint64_t insertReadIntoList(char *read1, char *read2,uint64_t *tableSize);
int isGoodRead(char *read);
void computeReverseComplementOfReads(void);
void readOverlapGraph(void);
void saveOverlapGraph(void);
int compareStringInReads(char *string1,char *string2,int start);
int compareStringInBytes(uint8_t *read2bits1,uint8_t *read2bits2,int start);
uint64_t get64BitInt(uint8_t *read,int start,int length);
void quicksortContigs(Edge **contigEdges,int64_t left, int64_t right);
void saveContigsToFile(char *graph_file,Edge **contig_edges,uint64_t from, uint64_t to);
uint64_t findGenomeSize(uint64_t previousEstimation);
void quicksortListOfLongEdges(Edge **edge_1,Edge **edge_2,int * pair_support, int64_t *gapDistance, int64_t left, int64_t right);
int makeDirectory(char *directoryName);
void sortEdgeListOfOverlapGraph(void);
void sortListOfNeighborsOverlapGraph(Edge *edge,uint64_t number,uint64_t from_ID);
void quicksortEdges(Edge ** listForSorting,int64_t left, int64_t right);
Edge** getListOfFeasibleEdges(Edge *edge, uint64_t *number);
Edge** getListOfCompositeEdges(uint64_t *number);
uint64_t mergeAccordingToSupport(Edge **edge_1,Edge **edge_2,int *pair_support,int64_t *gapDistance, int total_edge_pair_to_merge,int support_threshold);
int compareEdges(uint64_t a,uint64_t b,int type);
uint64_t *getListOfReads(Edge *edge1,Edge *edge2);
char *stringInEdge(Edge *edge,uint64_t *gapCount, uint64_t *nCount);
char *bytesToChars(uint8_t *readInt);
int stringCompareInBytes(uint8_t *read1,uint8_t *read2);
void stringCopyInBytes(uint8_t *read1,uint8_t *read2);
int mergeDisconnectedEdges(Edge *edge1,Edge *edge2,int gap);
uint64_t *getListOfReadsInDisconnectedEdges(Edge *edge1,Edge *edge2,int gap,uint64_t *length);
int stringOverlapSize(char *string1, char *string2);
int insertAllEdgesOfRead(uint64_t read);
int removeTransitiveEdges(uint64_t read);
char* reverseComplement(char *string);
uint64_t convertStringToNumber(char *str,int length);
uint64_t *getHashStringInInt(uint64_t read,int type);
uint64_t* get64Bit2Int(uint8_t *read,int start,int length);
uint64_t getHashValue(uint64_t *value);
void organizeReads(void);
void whereIs(uint64_t read);
int hashSupportInsert(Edge ** supportedEdge1, Edge ** supportededge2,Edge *edge1,Edge *edge2);

/* ============================================================================================
   This function will return minimum of a and b.
   ============================================================================================ */
inline int minimumValue(int a, int b)
{
    return a > b ? b : a;
}

/* ============================================================================================
   This function will return maximum of a and b.
   ============================================================================================ */

inline int maximumValue(int a, int b)
{
    return a > b ? a : b;
}

KSEQ_INIT(gzFile,gzread)

/* ============================================================================================
   Main function
   ============================================================================================ */
int main(int argc, char* argv[])
{
	int flg,a,i;
	uint64_t node_resolved;
	time_t starting_time=time (NULL);
	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo=localtime(&rawtime);
	if (argc < 4)
	{
		printf("Usage: ./SAGE <inputFile(s)> <outputDirName> <minOverlap> \n");
		exit(EXIT_FAILURE);
	}

	numberOfLibrary = argc - 3;
	minOverlap=atoi(argv[2+numberOfLibrary]);
	
	//make directory for output files
	strcat(directoryName,argv[1+numberOfLibrary]); 
	int dirNameLength = strlen(directoryName);
	if(directoryName[dirNameLength-1]!='/')
		strcat(directoryName,"/");
	makeDirectory(directoryName);
	
	//get first input file name to use for output names
	char inName[10000], outNameRev[10000];
	strcpy(inName,argv[1]);
	int inNameLength = strlen(inName), endOfName = 0, totalNameLength = 0, firstDot=0;
	while(endOfName==0 && inNameLength >= 0)
	{
		if(inName[inNameLength]!='/' && firstDot==1)
		{
			outNameRev[totalNameLength] = inName[inNameLength];
			totalNameLength++;
		}
		
		if(inName[inNameLength]=='/')
			endOfName = 1;
		if(inName[inNameLength]=='.' && firstDot==0)
			firstDot = 1;
		inNameLength--;
	}
	
	char outName[totalNameLength+1];
	int nameCounter, revNameLoc = totalNameLength-1;	
	for(nameCounter = 0; nameCounter <= totalNameLength; nameCounter++)
	{
		outName[nameCounter] = outNameRev[revNameLoc];
		revNameLoc--;
	}
	
	
	strcpy(outputTextFileName,directoryName); strcat(outputTextFileName,outName); strcat(outputTextFileName,"_SAGE_output_"); strcat(outputTextFileName,argv[2+numberOfLibrary]); strcat(outputTextFileName,".txt");// main ouput file
	if((fpo_output_file=fopen(outputTextFileName,"a"))==NULL) printError(OPEN_FILE,outputTextFileName);
	char hostName[10000]; gethostname(hostName,10000);
	myfprintf("***********************************************************************************************************\n");
	myfprintf("Name Length: %d\n",totalNameLength);
	myfprintf("\tEXECUTING PROGRAM: SAGE version 1.0\n\t         EXE FILE: %s\n",argv[0]);
	myfprintf("\t   NO. OF LIBRARY: %d\n",numberOfLibrary);
	myfprintf("\t        LIBRARIES:");
	for(i=1;i<1+numberOfLibrary;i++)
		myfprintf(" (%2d) %s",i, argv[i]);
	myfprintf("\n");
	myfprintf("\t        DIRECTORY: %s\n",argv[1+numberOfLibrary]);
	myfprintf("\t  MINIMUM OVERLAP: %d\n",minOverlap);
	myfprintf("\t             TIME: %s",asctime(timeinfo));
	myfprintf("\t         HOSTNAME: %s\n",hostName);
	if(DEBUGGING)
		myfprintf("\t        DEBUGGING: YES\n");
	else
		myfprintf("\t        DEBUGGING: NO\n");

	if(BUILDGRAPH)
		myfprintf("\t       BUILDGRAPH: YES\n");
	else
		myfprintf("\t       BUILDGRAPH: NO\n");

	if(SAVEGRAPH)
		myfprintf("\t        SAVEGRAPH: YES\n");
	else
		myfprintf("\t        SAVEGRAPH: NO\n");
	myfprintf("***********************************************************************************************************\n");
	if(minOverlap<16)
	{
		myfprintf("Too small minimum overlap length, it must be larger than 15.\n");
		exit(0);
	}
	if(BUILDGRAPH) // Build the graph from the scratch.
	{
		for(i=1;i<1+numberOfLibrary;i++)
		{
			myfprintf("\nReading Library %d\n",i);
			readDatasetInBytes(argv[i]); // uint8_t **readsInt and uint8_t **readsReverseInt from input file.
		}
		organizeReads();
		buildOverlapGraphEconomy(); //Free the read after hash Table is built.
		getDegreeStatistics(); // Mark highly dense parts of the graph.
		convertGraph(); // Convert economy graph to overlap graph.
		saveOverlapGraph(); // Save the graph to disk.
	}
	else // Load the graph from file.
	{
		for(i=1;i<1+numberOfLibrary;i++)
		{
			myfprintf("\nReading Library %d\n",i);
			readDatasetInBytes(argv[i]); // char **reads and **readsReverse from input file.
		}
		organizeReads();
		readOverlapGraph(); // Read the graph from the file.
	}
	while(1)
	{
		a=contractCompositePaths(); // Contract composite paths. u----v----w => u----w
		a+=removeDeadEnds(); // Remove dead-ends from the graph
		a+=removeBubbles(); // Remove bubbles from the graph.
		if(a==0) // if nothing is reduced is this iteration.
			break;
	}
	deleteDelta(); // Do not need it?
	numberOfReads-=checkIfAllReadsPresent(); // Subtract the missing reads from the number of reads.
	genomeSize=genomeSizeEstimation(); // Estimate the size of the genome.

	computeMinCostFlow();// Convex Cost with a_stats

	removeC2S1Edges();
	deleteSimpleEdgesWithLongDelta();
	for(i=1;i<1+numberOfLibrary;i++)
	{
		myfprintf("\nMapping Library %d\n",i);
		mapMatePairs(argv[i],i); // Loads the reads from file again.
	}
	meanSdEstimation(); // Estimate the mean and standard deviation of insert size.
	removeTriangles(); // Now the graph can have multiple edges between same endpoints.
	reduceTrees(); // In-tree and out-tree reduction.
	char singleEndFileName[10000];
	strcpy(singleEndFileName,"SingleEnd");
	printGraph(singleEndFileName); // Save the singled end assembly results.
	FILE* fp_summary_se=fopen("summarySingleEnd.txt","a");
	fprintf(fp_summary_se,"n:N50 %10"PRIu64" N50: %10"PRIu64" MAX: %8"PRIu64" SUM: %10"PRIu64" ( %5ld min.) %s\n",N_50_N,N_50,longestContig,base_pair_covered,(time(NULL)-starting_time)/60,argv[1+numberOfLibrary]);
	fclose(fp_summary_se);
	iteration=1; flg=0;
	markSimilarEdges(); // Mark similar edges (bubble like edges with similar string).
	memoryInitialize(); // Initialize some memory for the paths etc.
	while(1)
	{

		time_t second_counter=time(NULL);
		myfprintf("\n\n\n***********************************************************************************************************\n                                                ITERATION %2d\n***********************************************************************************************************",iteration);
		updateGraph(); // Simplify the graph if possible.
		node_resolved=findAllPairPaths(); // Merge pair of edges according to mate pair support/
		myfprintf("%"PRIu64" pair(s) of edges merged in iteration %2d in %4ld seconds\n",node_resolved,iteration, time(NULL)-second_counter);
		if(node_resolved==0 || iteration>=10)
		{
			if(flg==0)
			{
				flg=1; deleteEdgesWithSimilarStrings();
			}
			else
			{
				break;
			}
		}
		iteration++;
	}
	int removeCount;
	myfprintf("\n***********************************************************************************************************\n\n\n");
	char afterSupportFilename[10000];
	strcpy(afterSupportFilename,"afterSupport");
	if(TEST)
		printGraph(afterSupportFilename); // Save the singled end assembly results.
	do // Final clean up of the graph.
	{
		cleanGraph();
		removeCount=reduceLoops();
		removeCount+=removeTriangles();
		removeCount+=reduceTrees();
	}while(removeCount>0);
	mergeContigs(); // Merge pair of supported edges that does not have common endpoint.
	deleteTemporaryFiles(); // Delete any file with .tmp extention.
	char pairedEndFileName[10000];
	strcpy(pairedEndFileName,outName); strcat(pairedEndFileName,"_SAGE_assembly_"); strcat(pairedEndFileName,argv[2+numberOfLibrary]);
	printGraph(pairedEndFileName); // contig file in one directoryName
	freeAllocatedMemory();
	myfprintf("Highest number of paths: %"PRIu64"\n    Longest path Length: %"PRIu64"\n",highestPaths,longestPathLength);
	myfprintf("Total time needed %4ld sec. = %4ld Minutes.\n",time (NULL)-starting_time,(time (NULL)-starting_time)/60);
	//FILE* fp_summary_pe=fopen("summaryPairedEnd.txt","a");
	//fprintf(fp_summary_pe,"n:N50 %10"PRIu64" N50: %10"PRIu64" MAX: %8"PRIu64" SUM: %10"PRIu64" ( %5ld min.) %s\n",N_50_N,N_50,longestContig,base_pair_covered,(time(NULL)-starting_time)/60,argv[1+numberOfLibrary]);
	//fclose(fp_summary_pe);
	myfprintf("\n\n########################################### END ###########################################################\n");
	fclose(fpo_output_file);
	return 0;
}

/* ============================================================================================
   This function allocated all global variables.
   ============================================================================================ */
void memoryInitialize(void)
{
	if((supportedEdg1=(Edge**)malloc(defaultArraySize*sizeof(Edge*)))==NULL) printError(MEM_ALLOC, "supportedEdg1");
	if((supportedEdg2=(Edge**)malloc(defaultArraySize*sizeof(Edge*)))==NULL) printError(MEM_ALLOC, "supportedEdg2");
	if((support = (int*)malloc(defaultArraySize*sizeof(int)))==NULL) printError(MEM_ALLOC, "support");
}

/* ============================================================================================
   Compare two strings in bytes
   ============================================================================================ */
int stringCompareInBytes(uint8_t *read1,uint8_t *read2)
{
	uint64_t i;
	for(i=0;i<readArrayLength;i++)
	{
		if(read1[i]<read2[i])
			return -1;
		if(read1[i]>read2[i])
			return 1;
	}
	return 0;
}

/* ============================================================================================
   Copy string in bytes.
   ============================================================================================ */

void stringCopyInBytes(uint8_t *read1,uint8_t *read2)
{
	uint64_t i;
	for(i=0;i<readArrayLength;i++)
		*(read1+i)=*(read2+i);
}

/* ============================================================================================
   Memory efficient reading of the dataset. Only reads in blocks. Stores the reads in uint8_t **
   ============================================================================================ */
void readDatasetInBytes(char *filename)
{
	myfprintf("In function readDatasetBytes().\n");
	time_t seconds_s=time(NULL),seconds_temp=time(NULL);
	static uint64_t readInFile=1,tableSize;
	int l;
	char *read1=NULL,*read2=NULL;
	gzFile fp; kseq_t *seq;
	if((fp=gzopen(filename,"r"))==NULL) printError(OPEN_FILE,filename);
	seq = kseq_init(fp);
	myfprintf("Reading from file...\n");
	if(readInFile!=1)
	{
		if((read1=(char *)malloc((readLength+1)*sizeof(char)))==NULL) printError(MEM_ALLOC, "read1");
		if((read2=(char *)malloc((readLength+1)*sizeof(char)))==NULL) printError(MEM_ALLOC, "read2");
	}
	while ((l=kseq_read(seq)) >= 0)
	{
		if(readInFile==1)
		{
			readLength=strlen(seq->seq.s);
			if(readLength<=minOverlap)
			{
				myfprintf("Read Length should be > minimum overlap length.\nExiting program...\n");
				exit(0);
			}
			numberOfSameBasepairThreshold=(int)(readLength*80/100);
			readArrayLength=(int)ceil((double)readLength/4.0);
			if((read1=(char *)malloc((readLength+1)*sizeof(char)))==NULL) printError(MEM_ALLOC, "read1");
			if((read2=(char *)malloc((readLength+1)*sizeof(char)))==NULL) printError(MEM_ALLOC, "read2");
		}
		if(readInFile%2==1)
		{
			strcpy(read1,seq->seq.s);
		}
		else
		{
			strcpy(read2,seq->seq.s);
			if(isGoodRead(read1) && isGoodRead(read2)) // If both reads are good.
			{
				numberOfReads=insertReadIntoList(read1,read2,&tableSize); // Save the reads.
			}

		}
		readInFile++; // Count the total number of reads (good+bad) in the file.
		if(readInFile%1000000==0)
		{
			myfprintf("Finished Reading %10"PRIu64" to %10"PRIu64" reads in %4ld seconds\n",readInFile-1000000,readInFile, time(NULL)-seconds_temp);
			seconds_temp=time(NULL);
		}

	}
	free(read1);free(read2);
	myfprintf("Finished reading from file in %4ld sec.\n           Total reads: %10"PRIu64"\n            Good Reads: %10"PRIu64"\n             Bad Reads: %10"PRIu64"\n            Array size: %10"PRIu64"\nMinimum overlap length: %10d BP\n", time(NULL)-seconds_s,readInFile-1,numberOfReads,readInFile-numberOfReads-1,tableSize,minOverlap);
	myfprintf("Function readDatasetInBytes() in %4ld sec.\n",time(NULL)-seconds_s);
	kseq_destroy(seq);
	gzclose(fp);
}

void organizeReads(void)
{
	myfprintf("\nIn function organizeReads().\n");
	time_t seconds_s=time(NULL);
	uint64_t i,j;
	char *readReverse,*read;
	time_t qs_start=time(NULL);
	quicksortReads(readsInt,1,numberOfReads); // Sorts all the reads alphabetically.
	myfprintf("Quicksort reads finished in %4ld sec.\n",time(NULL)-qs_start);
	if((frequency=(int*)malloc((numberOfReads+1)*sizeof(int)))==NULL) printError(MEM_ALLOC, "frequency");
	for(i=0;i<numberOfReads+1;i++)
		frequency[i]=0;
	for(i=1,j=1;i<=numberOfReads;i++) // Remove duplicates from the list and store the frequency.
	{
		if(stringCompareInBytes(readsInt[j],readsInt[i])!=0)
		{
			j++;
			if(i!=j)
				stringCopyInBytes(readsInt[j],readsInt[i]);
		}
		frequency[j]++; // Store the frequency of the reads. We will need it for genome estimation and setting bounds of flow.
	}
	numberOfUniqueReads=j;// Number of unique reads in the dataset.
	for(i=numberOfUniqueReads+1;i<=numberOfReads;i++) // Only keep the list of reads upto numberOfUniqueReads.
		free(readsInt[i]);
	if((readsReverseInt=(uint8_t **)malloc((numberOfUniqueReads+1)*sizeof(uint8_t *)))==NULL) printError(MEM_ALLOC, "readsReverseInt");
	myfprintf("Computing reverse complements of reads.\n",time(NULL)-qs_start);

	for(i=1;i<=numberOfUniqueReads;i++)
	{
		read=bytesToChars(readsInt[i]);
		readReverse=reverseComplement(read);
		readsReverseInt[i]=charsToBytes(readReverse);
		free(readReverse);
		free(read);
	}

	if(DEBUGGING)
	{
		char readFile[10000];
		strcpy(readFile,directoryName);
		strcat(readFile,"reads.txt");
		FILE *fpReads=fopen(readFile,"w");
		for(i=1;i<=numberOfUniqueReads;i++)
		{
			char *readForward, *readReverse;
			readForward=bytesToChars(readsInt[i]);
			readReverse=bytesToChars(readsReverseInt[i]);
			fprintf(fpReads,"%10"PRIu64" %s %s\n",i,readForward,readReverse);
			free(readForward);
			free(readReverse);
		}
		fclose(fpReads);
	}
	myfprintf("Number of unique reads: %"PRIu64"\n           Read length: %d BP\n",numberOfUniqueReads,readLength);
	myfprintf("Function organizeReads in %4ld sec.\n",time(NULL)-seconds_s);
}

/* ============================================================================================
   This function maps pair of reads. It will load the reads from file again.
   ============================================================================================ */
void mapMatePairs(char *filename,int library)
{
	time_t seconds_s=time(NULL);
	static int allocationFlag=1;
	uint64_t flag_insert,i,a1=0,a2=0,a3=0,a4=0,readInFile=1;
	int64_t id1,id2;
	int8_t type1,type2;
	int l;
	struct matePairStructure *u,*v,*w;
	char *read1,*read2,*tempString;
	if((read1=(char *)malloc((readLength+1)*sizeof(char)))==NULL) printError(MEM_ALLOC, "read1");
	if((read2=(char *)malloc((readLength+1)*sizeof(char)))==NULL) printError(MEM_ALLOC, "read2");
	gzFile fp; kseq_t *seq;
	myfprintf("In function mapMatePairs().\n");
	if(allocationFlag==1)
	{
		if((mate_pair=(struct matePairStructure **)malloc((numberOfUniqueReads+1)*sizeof(struct matePairStructure *)))==NULL) printError(MEM_ALLOC, "mate_pair");
		for(i=0;i<=numberOfUniqueReads;i++)
			mate_pair[i]=NULL;
		allocationFlag=0;
	}
	if((fp=gzopen(filename,"r"))==NULL) printError(OPEN_FILE,filename);
	seq = kseq_init(fp);
	myfprintf("Reading from file...\n");
	while ((l=kseq_read(seq)) >= 0)
	{
		if(readInFile%2==1)
		{

			strcpy(read1,seq->seq.s);
		}
		else
		{
			strcpy(read2,seq->seq.s);
			if(isGoodRead(read1) && isGoodRead(read2))
			{
				id1=getIdOfRead(read1); id2=getIdOfRead(read2);
				if(frequency[llabs(id1)]>10000) // We do not add high frequency reads.
				{
					tempString=bytesToChars(readsInt[llabs(id1)]);
					myfprintf("Read %lld String %s appeared %d times changed to 1000\n",llabs(id1),tempString, frequency[llabs(id1)]);
					free(tempString);
					frequency[llabs(id1)]=10000;
				}
				if(frequency[llabs(id2)]>10000)
				{
					tempString=bytesToChars(readsInt[llabs(id2)]);
					myfprintf("Read %lld String %s appeared %d times changed to 1000\n",llabs(id2),tempString, frequency[llabs(id2)]);
					free(tempString);
					frequency[llabs(id2)]=10000;
				}
				if(frequency[llabs(id1)]<10000 && frequency[llabs(id1)]<10000)
				{
					if(id1>0)
					{
						type1=1; a1++;
					}
					else
					{
						type1=-1; a2++;
					}
					if(id2>0)
					{
						type2=1; a3++;
					}
					else
					{
						type2=-1; a4++;
					}
					id1=llabs(id1);id2=llabs(id2);
					flag_insert=1;
					for(w=mate_pair[id1];w!=NULL;w=w->next)  // Is this check really necessary ??? probably taking most of the time.
					{
						if(w->ID==(uint64_t)id2 && w->type1==type1 && w->type2==type2 && w->library==library)// Already present in the list.
						{
							flag_insert=0; w->freq++;
							break;
						}
					}
					if(flag_insert==1) // Insert if not already present
					{
						if((u = (struct matePairStructure *)malloc(sizeof(struct matePairStructure)))==NULL)
							printError(MEM_ALLOC, "u");
						u->ID=id2; u->freq=1;
						u->flag=1; u->type1=type1;u->type2=type2;
						u->library=library;
						u->next=mate_pair[id1]; mate_pair[id1]=u;
					}
					flag_insert=1;
					for(w=mate_pair[id2];w!=NULL;w=w->next) // Is this check really necessary ??? probably taking most of the time.
					{
						if(w->ID==(uint64_t)id1 && w->type1==type2 && w->type2==type1 && w->library==library) //Already present
						{
							flag_insert=0;
							w->freq++; break;
						}
					}
					if(flag_insert==1) // Insert if not already present.
					{
						if((v = (struct matePairStructure *)malloc(sizeof(struct matePairStructure)))==NULL)
							printError(MEM_ALLOC, "v");
						v->ID=id1; v->freq=1;
						v->flag=1; v->type1=type2; v->type2=type1;
						v->library=library;
						v->next=mate_pair[id2]; mate_pair[id2]=v;
					}
				}
			}

		}
		readInFile++;
	}
	free(read1);free(read2);
	kseq_destroy(seq);
	gzclose(fp);
	myfprintf("Matepair orientation:\n      (+) (+): %10"PRIu64" pairs.\n      (+) (-): %10"PRIu64" pairs.\n      (-) (+): %10"PRIu64" pairs.\n      (-) (-): %10"PRIu64" pairs.\nReads in file: %10"PRIu64"\nFunction mapMatePairs() in %4ld sec.\n",a1,a2,a3,a4,readInFile-1,time(NULL)-seconds_s);

}


/* ============================================================================================
   This function will estimate mean and standard deviation from matepairs on the same edge
   ============================================================================================ */
void meanSdEstimation(void)
{
	myfprintf("\nIn function meanSdEstimation().\n");
	time_t second_s=time(NULL);
	int counter,i;
	if((Mean=(unsigned int *)malloc((numberOfLibrary+1)*sizeof(unsigned int)))==NULL) printError(MEM_ALLOC, "Mean");
	if((standardDeviation=(unsigned int *)malloc((numberOfLibrary+1)*sizeof(unsigned int)))==NULL) printError(MEM_ALLOC, "Mean");
	if((upperBoundOfInsert=(unsigned int *)malloc((numberOfLibrary+1)*sizeof(unsigned int)))==NULL) printError(MEM_ALLOC, "upperBoundOfInsert");
	if((lowerBoundOfInsert=(unsigned int *)malloc((numberOfLibrary+1)*sizeof(unsigned int)))==NULL) printError(MEM_ALLOC, "lowerBoundOfInsert");
	for(i=0;i<=numberOfLibrary;i++)
	{
		Mean[i]=0;standardDeviation[i]=0;
	}
	mapReadsToEdges();
	mapReadLocations();
	for(counter=1;counter<=numberOfLibrary;counter++)
	{
		int mu=5000,sd=5000,returnMu,returnSD; // Initial mean set to 5000
		myfprintf("Computing mean and SD of library %d\n",counter);
		for(i=1;i<=10;i++) // To avoid infinite loop
		{
			computeMeanSD(mu,sd,counter,&returnMu,&returnSD);
			if(abs(mu-returnMu)<=returnMu/100  && abs(sd-returnSD)<=returnSD/100)
			{
				mu=returnMu; sd=returnSD;
				myfprintf("Mean: %5d SD: %5d (Final extimation)\n\n",mu,sd);
				break;
			}
			else
			{
				mu=returnMu; sd=returnSD;
				myfprintf("Mean: %5d SD: %5d (Estimation %2d)\n",mu,sd,i);
			}
		}
		Mean[counter]=returnMu; standardDeviation[counter]=returnSD;
	}
	minimumUpperBoundOfInsert=1000000;
	maximumUpperBoundOfInsert=0;
	for(i=1;i<=numberOfLibrary;i++)
	{
		lowerBoundOfInsert[i]=Mean[i]-3*standardDeviation[i];
		//if(lowerBoundOfInsert[i]<0)
			lowerBoundOfInsert[i]=0;
		upperBoundOfInsert[i]=Mean[i]+3*standardDeviation[i];

		if(minimumUpperBoundOfInsert>=upperBoundOfInsert[i])
			minimumUpperBoundOfInsert=upperBoundOfInsert[i];
		if(maximumUpperBoundOfInsert<=upperBoundOfInsert[i])
			maximumUpperBoundOfInsert=upperBoundOfInsert[i];

		myfprintf("Libaray %3d: Mean: %5d SD: %5d LB: %5d UB: %5d\n",i, Mean[i],standardDeviation[i],lowerBoundOfInsert[i],upperBoundOfInsert[i]);
	}
	myfprintf("Minimum upper bound: %10d\nMaximum upper bound: %10d\nFunction meanSdEstimation() in %4ld sec.\n",minimumUpperBoundOfInsert,maximumUpperBoundOfInsert,time(NULL)-second_s);
}

/* ============================================================================================
   Next estimation of mean and standard deviation from the previous estimation.
   ============================================================================================ */
void computeMeanSD(int mu,int sd,int library,int *returnMu,int *returnSD)
{
	uint64_t i,mate_pair_1,mate_pair_2,distance;
	int64_t *distance1,*distance2;
	long int sum=0,counter=0;
	long double sqrd_error=0;
	struct matePairStructure *a;
	struct readToEdgeStructure *u,*v;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		mate_pair_1=i;
		for(a=mate_pair[mate_pair_1];a!=NULL;a=a->next)
		{
			mate_pair_2=a->ID;
			if(a->flag==0 && mate_pair_1<mate_pair_2 && a->library==library)
			{
				for(u=read_to_edge[mate_pair_1];u!=NULL;u=u->next)
				{
					for(v=read_to_edge[mate_pair_2];v!=NULL;v=v->next)
					{
						if(u->edge==v->edge || u->edge==v->edge->twinEdge) // mate_pairs on same edge
						{
							distance1=findDistanceOnEdge(v->edge,mate_pair_1);
							distance2=findDistanceOnEdge(v->edge,mate_pair_2);
							distance=(uint64_t)llabs(llabs(distance1[1])-llabs(distance2[1]));
							if(distance1[0]==1 && distance2[0]==1 && distance<(uint64_t)(4*mu))
							{
								sqrd_error+=(mu-distance)*(mu-distance);
								counter++; sum+=distance;
								free(distance1); free(distance2);
								break;
							}
							free(distance1); free(distance2);
						}
					}
				}
			}
		}
	}
	myfprintf("Mate-pairs considered: %d\n",counter);
	mu=sum/counter;	sd=sqrt(sqrd_error/(counter-1));
	*returnMu=mu; *returnSD=sd;
}


/* ============================================================================================
   This function finds ID of a read. Given a read, this function finds the ID from the sorted
   list by binary search.
   ============================================================================================ */
int64_t getIdOfRead(char *read)
{
	int64_t lb=1,ub=numberOfUniqueReads, mid;
	int comp_result, flag;
	uint8_t *readInteger;
	char *read_r=reverseComplement(read);
	if(strcmp(read,read_r)<0)
	{
		readInteger=charsToBytes(read);
		flag=1;
	}
	else
	{
		readInteger=charsToBytes(read_r);
		flag=-1;
	}
	free(read_r);
	while(lb<=ub)
	{
		mid=(ub+lb)>>1;
		comp_result=stringCompareInBytes(readInteger,readsInt[mid]);
		if(comp_result==0)
		{
			free(readInteger);
			return mid*flag;
		}
		else if(comp_result>0)
			lb=mid+1;
		else
			ub=mid-1;
	}
	free(readInteger);
	printf("Returning 0\n");
	return 0;
}

/* ============================================================================================
   This function will sort the array reads[]. reads[] is a global array, so the array was not
   passed in the function.
   ============================================================================================ */
void quicksortReads(uint8_t **reads,int64_t left, int64_t right)
{
	int64_t i=left,j=right, pivot=(left+right)>>1,index;
	uint8_t *pvt,*temp;
	if((pvt=(uint8_t *)malloc((readArrayLength)*sizeof(uint8_t)))==NULL) printError(MEM_ALLOC, "pvt");
	for(index=0;index<(int64_t)readArrayLength;index++)
	{
		pvt[index]=reads[pivot][index];
	}
	while(i<j)
	{
		while(stringCompareInBytes(reads[i],pvt)<0)
			i++;
		while(stringCompareInBytes(reads[j],pvt)>0)
			j--;
		if (i<=j)
		{
			temp=*(reads+i);
			*(reads+i)=*(reads+j);
			*(reads+j)=temp;
			i++;
			j--;
		}
	}
	free(pvt);
	if (left < j)
		quicksortReads(reads,left,j);
    	if (i < right)
		quicksortReads(reads,i,right);
}

/* ============================================================================================
   This function estimates the genome size.
   ============================================================================================ */
uint64_t genomeSizeEstimation(void)
{
	myfprintf("\nIn function genomeSizeEstimation().\n");
	time_t seconds_s=time(NULL);
	uint64_t i,currentGenomeSize,previousGenomeSize=0;
	for(i=1;i<=10;i++) // To avoid infinite loop
	{
		currentGenomeSize=findGenomeSize(previousGenomeSize);
		if(previousGenomeSize==currentGenomeSize) // Should be unchanged from previous estimation.
		{
			myfprintf("Genome Size: %10"PRIu64" (Final estimation.)\n",currentGenomeSize);
			break;
		}
		else
		{
			previousGenomeSize=currentGenomeSize;
			myfprintf("Genome Size: %10"PRIu64" (Estimation %2"PRIu64")\n",currentGenomeSize,i);
		}
	}
	myfprintf("\nFunction genomeSizeEstimation() in %4ld sec.\n",time(NULL)-seconds_s);
	return currentGenomeSize;
}


/* ============================================================================================
   Find next estimation from previous estimation.
   ============================================================================================ */
uint64_t findGenomeSize(uint64_t previousEstimation)
{
	uint64_t i,delta,k,deltaSum=0,kSum=0,currentGenomeSize,index;
	float a_statistic;
	Edge *v;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		for(v=graph[i];v!=NULL;v=v->next)
		{
			if(i<v->ID)
			{
				delta=v->lengthOfEdge-1;
				k=0;
				if(v->listOfReads[0]!=0)
				{
					for(index=1;index<=(uint64_t)(v->listOfReads[0]);index++)
						k=k+frequency[v->listOfReads[index]>>24];
				}
				if(previousEstimation!=0) // If not the first estimation
				{
					a_statistic=(float)((float)delta*(float)((float)numberOfReads/(float)previousEstimation)-(float)((float)k*log((float)2)));
					if(a_statistic>=aStatisticsThreshold && delta>=minDelta)
					{
						deltaSum=deltaSum+delta;
						kSum=kSum+k;
					}
				}
				else if(v->lengthOfEdge>500)  // Only take edges longer than 500 for the first estimation
				{
					deltaSum=deltaSum+delta;
					kSum=kSum+k;
				}
			}
		}
	}
	currentGenomeSize=(uint64_t)((float)numberOfReads/(float)kSum*(float)deltaSum);
	return currentGenomeSize;
}

/* ============================================================================================
   This function removes composite path in the graph
   ============================================================================================ */
uint64_t contractCompositePaths(void)
{
	myfprintf("\nIn function contractCompositePaths().\n");
	time_t seconds_s;
	seconds_s=time (NULL);
	uint64_t i,how_many_removed=0, unable_to_remove=0;
	int flag;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graph[i]!=NULL && graph[i]->next!=NULL && graph[i]->next->next==NULL) // Nodes with exactly two neighbors
		{
			flag=0;
			Edge *u;
			for(u=graph[graph[i]->ID];u!=NULL;u=u->next) // This loop will avoid multiple edges between 2 nodes
			{
				if(u->ID==graph[i]->next->ID)
				{
					flag=1; break;
				}
			}
			if(flag==1)
				continue;
			if(graph[i]->flow!=graph[i]->next->flow)
			{
				myfprintf("In-flow and out-flow different at node %4"PRIu64". flow: %.2f vs %.2f\n",i,graph[i]->flow,graph[i]->next->flow);
				continue;
			}
			if(mergeEdges(graph[i]->twinEdge,graph[i]->next,2))
				how_many_removed++;
			else
				unable_to_remove++;
		}
	}
	myfprintf("   Nodes removed: %10"PRIu64"\nUnable to remove: %10"PRIu64"\nFunction contractCompositePaths() in %4ld sec.\n",how_many_removed,unable_to_remove,time(NULL)-seconds_s);
	return how_many_removed;
}

/* ============================================================================================
   This function deletes the edge between vertex_a and vertex_b.
   ============================================================================================ */
int deleteEdge(Edge *edge) // Delete in constant time.
{
	if(DEBUGGING)
	{
		if(edge->from_ID==edge->ID)
			myfprintf("Deleting loop in deleteEdge() (%"PRIu64",%"PRIu64") length: %"PRIu64"\n",edge->from_ID,edge->ID,edge->lengthOfEdge);
	}
	if(edge->previous==NULL && edge->next==NULL) //Only one element
	{
		graph[edge->from_ID]=NULL;
	}
	else if(edge->previous==NULL && edge->next!=NULL) //First element
	{
		graph[edge->from_ID]=edge->next;
		edge->next->previous=NULL;
	}
	else if(edge->previous!=NULL && edge->next==NULL) //Last element;
	{
		edge->previous->next=NULL;
	}
	else // Middle element
	{
		edge->previous->next=edge->next;
		edge->next->previous=edge->previous;
	}
	free(edge->listOfReads);
	free(edge); // free the edge
	return 1;
}



/* ============================================================================================
   This function prints the overlap graph in overlap_graph.gdl file. The graph can be viewed by
   aisee (free software available at http://www.aisee.com/)
   ============================================================================================ */
void printGraph(char *graph_file)
{
	myfprintf("\nIn function printGraph().\n");
	Edge **contig_edges,*v,*u;
	time_t seconds_s=time(NULL);
	uint64_t color,i,thickness,flg,lengthOfArray=1000,numberOfContigsFound=0;
	char name_2d[10000];
	//FILE *fpo_overlap_graph;

	if((contig_edges=(Edge **)malloc(lengthOfArray*sizeof(Edge *)))==NULL) printError(MEM_ALLOC, "contig_edges");
	strcpy(name_2d,directoryName);strcat(name_2d,graph_file); strcat(name_2d,".gdl");

	/*
	if((fpo_overlap_graph=fopen(name_2d,"w"))==NULL) printError(OPEN_FILE,name_2d);

	fprintf(fpo_overlap_graph,"graph: {\nlayoutalgorithm :forcedir\nfdmax:704\ntempmax:254\ntempmin:0\ntemptreshold:3\ntempscheme:3\ntempfactor:1.08\nrandomfactor:100\ngravity:0.0\nrepulsion:161\nattraction:43\nignore_singles:yes\nnode.fontname:\"helvB10\"\nnode.shape:box\nnode.width:80\nnode.height:20\nnode.borderwidth:1\nnode.bordercolor:31\n");
	for(i=1;i<=numberOfUniqueReads;i++) // Nodes
		if(graph[i]!=NULL)
			fprintf(fpo_overlap_graph,"node: { title:\"%"PRIu64"\" label: \"%"PRIu64"\"}\n",i,i);
	fprintf(fpo_overlap_graph,"edge.arrowstyle:none\n");
	*/
	
	for(i=1;i<=numberOfUniqueReads;i++) //Edges
	{
		for(v=graph[i];v!=NULL;v=v->next)
		{
			if(i>=v->ID)
			{
				if(i==v->ID) // do not double count loops (a,a)
				{
					flg=0;
					for(u=graph[i];u!=v;u=u->next)
					{
						if(u==v->twinEdge)// twin edge already considered
						{
							flg=1; break;
						}
					}
					if(flg==1)
						continue;
				}
				count++;
				if(v->listOfReads[0]!=0) thickness=5;//Thick composite edges
				else thickness=1;//Thin simple edge
				if(v->typeOfEdge==1) color=3; //color code for aisee 3=green
				else if(v->typeOfEdge==2) color=1; //1=blue
				else color=2; //2=red
				if(N_90!=0 && v->lengthOfEdge>N_90 ) color=12; //edges longer than N90 12=darkmagenta
				if(!v->isReducible) color=31; //edges that are not reducible are color in black
				uint64_t k=0;
				k=v->listOfReads[0];
				/*
				if(v->typeOfEdge==0)
					fprintf(fpo_overlap_graph,"edge: { source:\"%"PRIu64"\" target:\"%"PRIu64"\" thickness: %"PRIu64" backarrowstyle: solid color: %"PRIu64" label: \"Flow: %8.2f Delta: %"PRIu64", k: %"PRIu64"\" }\n",i,v->ID,thickness,color,v->flow,v->lengthOfEdge,k);
				else if(v->typeOfEdge==1)
					fprintf(fpo_overlap_graph,"edge: { source:\"%"PRIu64"\" target:\"%"PRIu64"\" thickness: %"PRIu64" backarrowstyle:solid arrowstyle:solid color: %"PRIu64" label: \"Flow: %8.2f Delta: %"PRIu64", k: %"PRIu64"\"}\n",i,v->ID,thickness,color,v->flow,v->lengthOfEdge,k);
				else if(v->typeOfEdge==2)
					fprintf(fpo_overlap_graph,"edge: { source:\"%"PRIu64"\" target:\"%"PRIu64"\" thickness: %"PRIu64" color: %"PRIu64" label: \"Flow: %8.2f Delta: %"PRIu64", k: %"PRIu64"\"}\n",i,v->ID,thickness,color,v->flow,v->lengthOfEdge,k);
				else if(v->typeOfEdge==3)
					fprintf(fpo_overlap_graph,"edge: { source:\"%"PRIu64"\" target:\"%"PRIu64"\" thickness: %"PRIu64" arrowstyle:solid color: %"PRIu64" label: \"Flow: %8.2f Delta: %"PRIu64", k: %"PRIu64"\"}\n",i,v->ID,thickness,color,v->flow,v->lengthOfEdge,k);
				*/
				if(v->lengthOfEdge>=contigLengthThreshold)
				{
					++numberOfContigsFound;
					contig_edges[numberOfContigsFound]=v;
					if(numberOfContigsFound>=lengthOfArray-5)
					{
						lengthOfArray+=1000; // Increase the array size as needed.
						if((contig_edges=(Edge **)realloc(contig_edges,lengthOfArray*sizeof(Edge *)))==NULL) printError(MEM_ALLOC,"reallocating contig_edges failed");
					}
				}
			}
		}
	}
	//fprintf(fpo_overlap_graph,"}}");
	quicksortContigs(contig_edges,1,numberOfContigsFound); // Quicksort the contigs according to their lengths.
	saveContigsToFile(graph_file,contig_edges,1,numberOfContigsFound); //Save contigs in a fasta file.
	//fclose(fpo_overlap_graph);
	free(contig_edges);
	myfprintf("Function printGraph() in %4ld sec.\n",time(NULL)-seconds_s);
}


/* ============================================================================================
   Sort contigs according to their lengths.
   ============================================================================================ */
void quicksortContigs(Edge **contigEdges,int64_t left, int64_t right)
{
	int64_t i=left,j=right;
	uint64_t pivot=contigEdges[(left+right)>>1]->lengthOfEdge;
	Edge *tempEdge;
	while(i<j)
	{
		while(contigEdges[i]->lengthOfEdge>pivot)
			i++;
		while(contigEdges[j]->lengthOfEdge<pivot)
			j--;
		if (i<=j)
		{
			tempEdge=*(contigEdges+i);	*(contigEdges+i)=*(contigEdges+j);	*(contigEdges+j)=tempEdge;
			i++; j--;
		}
	}
	if (left < j)
		quicksortContigs(contigEdges,left, j);
    	if (i < right)
		quicksortContigs(contigEdges,i, right);
}


/* ============================================================================================
   Sort contigs according to their lengths.
   ============================================================================================ */

void saveContigsToFile(char *graph_file,Edge **contig_edges,uint64_t from, uint64_t to)
{
	uint64_t i,contig_length,N_X=0,N50=0,N90=0,N50_count=0,N90_count=0,gapCount,nCount;
	base_pair_covered=0;
	FILE *fpo_contigs_file;
	Edge *v;
	char contigFileName[10000],*contig_string,*stringPointer;
	strcpy(contigFileName,directoryName); strcat(contigFileName,graph_file); strcat(contigFileName,".fasta");
	if((fpo_contigs_file=fopen(contigFileName,"w"))==NULL) printError(OPEN_FILE,contigFileName);
	for(i=from;i<=to;i++)
	{
		if(contig_edges[i]->lengthOfEdge>=(int)(contigLengthThreshold))
			base_pair_covered+=contig_edges[i]->lengthOfEdge;
	}
	myfprintf("\nTable: Partial list of long contigs in printGraph().\n\t--------------------------------------------------------------------\n\t| # |       EDGE        |     LENGTH  | FLOW | GAP COUNT | N COUNT |\n\t--------------------------------------------------------------------\n");
	for(i=from;i<=to;i++)
	{
		v=contig_edges[i];
		contig_string=stringInEdge(v,&gapCount,&nCount);
		stringPointer=contig_string;
		contig_length=v->lengthOfEdge;
		if(i==from) longestContig=contig_length;
		if(contig_length>=(int)(contigLengthThreshold)) //Put contigs in contig file in sorted order
		{
			fprintf(fpo_contigs_file,">contig_%"PRIu64" flow=%.2f Edge length=%"PRIu64" String length: %d e(%"PRIu64",%"PRIu64") gapCount=%"PRIu64" nCount=%"PRIu64"\n",i,v->flow,contig_length,(int)(strlen(contig_string)),v->from_ID,v->ID,gapCount,nCount);
			contig_string=contig_string+(readLength/2);
			contig_length=contig_length+(readLength/2);
			while(contig_length-(readLength/2)>=60)
			{
				fprintf(fpo_contigs_file,"%.*s\n",60,contig_string);
				contig_string=contig_string+60; contig_length=contig_length-60;
			}
			if(contig_length-(readLength/2)>0) // Print the remaining of the contig
				fprintf(fpo_contigs_file,"%.*s\n",(int)contig_length-(readLength/2),contig_string);
		}
		free(stringPointer);
		contig_length=v->lengthOfEdge;
		if(i<=10) // Print only large 100 contigs in output file. This slows down the program if there are many contigs.
			myfprintf("\t|%2"PRIu64" | %8"PRIu64",%8"PRIu64" | %8"PRIu64" BP |%5.1f | %9"PRIu64" | %7"PRIu64" |\n",i,contig_edges[i]->from_ID,contig_edges[i]->ID,contig_length,contig_edges[i]->flow,gapCount,nCount);
		N_X+=contig_length;
		if(!N50_count && N_X>=base_pair_covered/2)
		{
			 N50=contig_length; N50_count=i;
		}
		if(!N90_count && N_X>=base_pair_covered*9/10)
		{
			N90=contig_length; N90_count=i;
		}

	}
	myfprintf("\t--------------------------------------------------------------------\nNumber of contigs: %10"PRIu64"\n",to);
	numberOfContigs=to; N_50=N50; N_90=N90;N_50_N=N50_count;N_90_N=N90_count;
	if(genomeSize)
		genomeCovered=(float)((float)base_pair_covered/(float)genomeSize)*100;
	else
		genomeCovered=0;
	fclose(fpo_contigs_file);
	myfprintf("\n      Genome Size: %10"PRIu64" BP\n  Number of Reads: %10"PRIu64"\n   Genome Covered: %10"PRIu64" BP (%.4f%%)\nNumber of Contigs: %10"PRIu64"\n   Longest Contig: %10"PRIu64"\n              N50: %10"PRIu64" BP (%6"PRIu64")\n              N90: %10"PRIu64" BP (%6"PRIu64")\n",genomeSize, numberOfReads,base_pair_covered, genomeCovered, to,longestContig,N50,N50_count,N90,N90_count);
}



/* ============================================================================================
   Insert edge in the overlap graph. It does not check for duplicates and also does not consider
   transitive edges.
   ============================================================================================ */
int insertEdge(uint64_t vertex_u,uint64_t vertex_v, int type,uint64_t *listForward, uint64_t *listReverse,float flow,uint64_t delta)
{

	Edge *v,*u;
	if((v = (Edge *)malloc(sizeof(Edge)))==NULL) printError(MEM_ALLOC, "v");
	if((u = (Edge *)malloc(sizeof(Edge)))==NULL) printError(MEM_ALLOC, "u");
	v->from_ID=vertex_u; 		u->from_ID=vertex_v;
	v->ID=vertex_v; 		u->ID=vertex_u;
	v->typeOfEdge=type; 		u->typeOfEdge=reverseEdgeType(type);
	v->twinEdge=u;			u->twinEdge=v;
	v->flow=flow; 			u->flow=flow;
	v->next=NULL; 			u->next=NULL;
	v->previous=NULL;		u->previous=NULL;
	v->listOfReads=listForward;	u->listOfReads=listReverse;
	v->isReducible=1;		u->isReducible=1;
	v->lengthOfEdge=delta; 		u->lengthOfEdge=delta;
	insertIntoList(v); 		insertIntoList(u);
	return 1;

}

/* ============================================================================================
   Mark transitive edges using Gene Myers algorithm.
   Time: O(nd) , n is the number of node and d is the highest degree.
   Memory: O(1)
   TODO: Think what happens for clique.
   ============================================================================================ */
void markTransitiveEdge(uint64_t from_ID)
{
	uint64_t i,j,a,b,type1,type2;
	for(i=1;i<=(uint64_t)(graphEconomy[from_ID][0]);i++)
		markedNodes[graphEconomy[from_ID][i]>>24]=1; // Inplay
	for(i=1;i<=(uint64_t)(graphEconomy[from_ID][0]);i++)
	{
		a=graphEconomy[from_ID][i]>>24;
		if(markedNodes[a]==1) // If Inplay
		{
			for(j=1;j<=(uint64_t)(graphEconomy[a][0]);j++)
			{
				b=graphEconomy[a][j]>>24;
				if(markedNodes[b]==1)
				{
					type1=(graphEconomy[from_ID][i]&0X0000000000C00000)>>22;
					type2=(graphEconomy[a][j]&0X0000000000C00000)>>22;
					if((type1==0 || type1==2) && (type2==0 || type2==1)) // check edge orientation
						markedNodes[b]=2; //Eleminated
					else if((type1==1||type1==3) && (type2==2 || type2==3)) //check edge orientation
						markedNodes[b]=2; //Eleminated
				}
			}
		}
	}
	for(i=1;i<=(uint64_t)(graphEconomy[from_ID][0]);i++)
	{
		if(markedNodes[graphEconomy[from_ID][i]>>24]==2) //Eliminated
		{
			graphEconomy[from_ID][i]=(graphEconomy[from_ID][i]|0X0000000000100000); //marked for deleteion
		}
	}

	for(i=1;i<=(uint64_t)(graphEconomy[from_ID][0]);i++)
		markedNodes[graphEconomy[from_ID][i]>>24]=0; // Mark as vacant.
	markedNodes[from_ID]=0; // Mark as vacant.
	exploredReads[from_ID]=2;
}


int removeTransitiveEdges(uint64_t read)
{
	uint64_t index,numberRemoved;
	uint64_t *newList;
	if((newList=(uint64_t *)malloc(1*sizeof(uint64_t )))==NULL) printError(MEM_ALLOC, "newList");
	newList[0]=0;
	for(index=1;index<=(uint64_t)(graphEconomy[read][0]);index++)
	{
		if(!(graphEconomy[read][index]&0X0000000000100000))
		{
			if((newList=(uint64_t *)realloc(newList,(newList[0]+2)*sizeof(uint64_t )))==NULL) printError(MEM_ALLOC, "realloc newList");
			newList[++newList[0]]=graphEconomy[read][index];
		}
	}
	numberRemoved=graphEconomy[read][0]-newList[0];
	free(graphEconomy[read]);
	graphEconomy[read]=newList;
	return numberRemoved;
}

/* ============================================================================================
   This function returns the reverse type of an edge.
   TODO: Do we need it? think.
   ============================================================================================ */
int reverseEdgeType(int type)
{
	if(type==0)
		return 3;
	else if(type==3)
		return 0;
	return type;
}


/* ============================================================================================
   This function inserts v in the list of edges of v->from_ID. Double linked list.
   Time: O(1)
   ============================================================================================ */
void insertIntoList(Edge *v) // insert in linear time;
{
	if(graph[v->from_ID]==NULL) //first time
		graph[v->from_ID]=v;
	else //insert in the top
	{
		v->next=graph[v->from_ID];
		graph[v->from_ID]->previous=v;
		graph[v->from_ID]=v;
	}
	return;
}


/* ============================================================================================
   This function creats the overlap graph from hash table.
   Time: O(nl^2d)
   Memory: O(nd)
   ============================================================================================ */
void buildOverlapGraphEconomy(void)
{
	myfprintf("\nIn function buildOverlapGraphEconomy().\n");
	uint64_t transitiveFreeNodes=0,nodeExplored=0,i=0;
	uint64_t hashMissTemp=0,totalEdgeInserted=0,totalEdgeInsertedTemp=0,transitiveEdgesRemoved=0,transitiveEdgesRemovedTemp=0,noEdge=0,highestNeighbours=0;
	uint64_t *queue,start=0,end=0,read1=0,read2=0,read3=0,index1=0,index2=0;
	
	time_t seconds_s=time(NULL);
	hashPrefixesAndSuffix(); // Hash prefix and suffix of reads and reverese reads
	time_t seconds_ovlp=time(NULL);
	if((exploredReads=(uint8_t *)malloc((numberOfUniqueReads+1)*sizeof(uint8_t)))==NULL) printError(MEM_ALLOC, "exploredReads");
	if((selfOverlappingReads=(uint8_t *)malloc((numberOfUniqueReads+1)*sizeof(uint8_t)))==NULL) printError(MEM_ALLOC, "selfOverlappingReads");
	if((graphEconomy=(uint64_t **)malloc((numberOfUniqueReads+1)*sizeof(uint64_t *)))==NULL) printError(MEM_ALLOC, "graphEconomy");
	if((markedNodes=(uint8_t *)malloc((numberOfUniqueReads+1)*sizeof(uint8_t)))==NULL) printError(MEM_ALLOC, "markedNodes");
	if((queue=(uint64_t *)malloc((numberOfUniqueReads+1)*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC, "queue");

	for(i=0;i<=numberOfUniqueReads;i++) // Initialize all flags
	{
		graphEconomy[i]=NULL; selfOverlappingReads[i]=0; exploredReads[i]=0; markedNodes[i]=0;
	}
	myfprintf("\n---------------------------------------------------------------------------------------\n|    FROM   |     TO     |    EDGES   |   TIME |   HS MISS  |  TR EDGE   | TR FR NODE |\n---------------------------------------------------------------------------------------\n");
	
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(exploredReads[i]==0)
		{
			//connectedComponent++;
			start=0;end=0; // Initialize queue start and end
			queue[end++]=i; // Put the read in the queue
			while(start<end) // This loop will explore all connected component starting from read i
			{
				nodeExplored++;
				
				read1=queue[start++];
				if(exploredReads[read1]==0)
					totalEdgeInserted+=insertAllEdgesOfRead(read1); // Explore current node
				if(graphEconomy[read1]!=NULL) // Read has some edges (required only for the first read when a new queue starts
				{
					if(highestNeighbours<(uint64_t)(graphEconomy[read1][0]))
						highestNeighbours=graphEconomy[read1][0];
					if(exploredReads[read1]==1) // Expore neighbours first
					{
						for(index1=1;index1<=(uint64_t)(graphEconomy[read1][0]);index1++) // Explore all neighbours
						{
							read2=graphEconomy[read1][index1]>>24;
							if(exploredReads[read2]==0) // Not explored
							{
								queue[end++]=read2; // Put in the queue
								totalEdgeInserted+=insertAllEdgesOfRead(read2);
							}
						}
						markTransitiveEdge(read1); // Mark transitive edges
					}
					if(exploredReads[read1]==2)
					{
						for(index1=1;index1<=(uint64_t)(graphEconomy[read1][0]);index1++) // Then explore all neighbour's neighbours
						{
							read2=graphEconomy[read1][index1]>>24;
							if(exploredReads[read2]==1)
							{
								for(index2=1;index2<=(uint64_t)(graphEconomy[read2][0]);index2++) // Explore all neighbours neighbours
								{
									read3=graphEconomy[read2][index2]>>24;
									if(exploredReads[read3]==0) // Not explored
									{
										queue[end++]=read3; // Put in the queue
										totalEdgeInserted+=insertAllEdgesOfRead(read3);
									}
								}
								markTransitiveEdge(read2); // Mark transitive edge
							}
						}
						
						transitiveFreeNodes++;
						
						transitiveEdgesRemoved+=removeTransitiveEdges(read1); // Remove the transitive edges
					}
				}
				else
				{
					noEdge++;
				}
				if(nodeExplored%100000==0) // print after every 100,000 reads. Will write about 500 times for the C.elegans
				{
					myfprintf("|%10"PRIu64" | %10"PRIu64" | %10"PRIu64" | %6ld | %10"PRIu64" | %10"PRIu64" | %10"PRIu64" |\n",nodeExplored-100000+1,nodeExplored, totalEdgeInserted-totalEdgeInsertedTemp,time(NULL)-seconds_ovlp, numberOfHashMiss-hashMissTemp,transitiveEdgesRemoved-transitiveEdgesRemovedTemp,transitiveFreeNodes);
					//myfprintf("|%10d | %10d | %10llu | %6ld | %10llu | %10llu | %10d | i: %d CC: %d noEdge: %d HST NBR: %d HS MTCH: %llu\n",nodeExplored-100000+1,nodeExplored, totalEdgeInserted-totalEdgeInsertedTemp,time(NULL)-seconds_ovlp, numberOfHashMiss-hashMissTemp,transitiveEdgesRemoved-transitiveEdgesRemovedTemp,transitiveFreeNodes,i,connectedComponent,noEdge,highestNeighbours,numberOfHashMatch);
					seconds_ovlp=time(NULL);
					totalEdgeInsertedTemp=totalEdgeInserted;
					hashMissTemp=numberOfHashMiss;
					transitiveEdgesRemovedTemp=transitiveEdgesRemoved;
					highestNeighbours=0;
					numberOfHashMatch=0;
				}
			}
		}
	}
	myfprintf("|%10"PRIu64" | %10"PRIu64" | %10"PRIu64" | %6ld | %10"PRIu64" | %10"PRIu64" | %10"PRIu64" |\n---------------------------------------------------------------------------------------\n\n",nodeExplored-(nodeExplored%100000)+1,nodeExplored,totalEdgeInserted-totalEdgeInsertedTemp,time(NULL)-seconds_ovlp, numberOfHashMiss-hashMissTemp,transitiveEdgesRemoved-transitiveEdgesRemovedTemp,transitiveFreeNodes);
	myfprintf("Function buildOverlapGraphEconomy(): %4ld sec.\n     Total edges inserted: %10"PRIu64"\nTotal number of hash miss: %10"PRIu64"\n  Transitive edge removed: %10"PRIu64"\n",time(NULL)-seconds_s, totalEdgeInserted,numberOfHashMiss,transitiveEdgesRemoved);
	free(markedNodes); free(exploredReads); free(queue);
	freeHashTable(); // we do not need it after building the graph.
}

/* ============================================================================================
   This function inserts all the edges of "read" and mark it as explored. At the end it sorts
   the edges according to their length.
   ============================================================================================ */

int insertAllEdgesOfRead(uint64_t read)
{
	uint64_t j,k,read1=read,totalEdgeInserted=0,numberToSearch=0;
	uint64_t *uint64Window;
	int64_t index;
	if(exploredReads[read1]!=0) // already explored
		return 0;
	else
		exploredReads[read1]=1; // mark as the node is explored

	for(j=0;j<=(uint64_t)(readLength-hashStringLength);j++)
	{
		uint64Window=get64Bit2Int(readsInt[read1],j,hashStringLength);
		index=hashTableSearch(uint64Window);
		free((uint64_t *) uint64Window);
		if(index!=-1) // found in hashtable
		{
			numberToSearch+=listOfHashReadsInt[index][0];
			numberOfHashMatch+=listOfHashReadsInt[index][0];
			for(k=1;k<=listOfHashReadsInt[index][0];k++)
			{
				uint64_t read2=listOfHashReadsInt[index][k]>>24;
				uint64_t type=listOfHashReadsInt[index][k] & 0X0000000000000003;
				int ovlpLength=-1,typeOfEdgeToInsert=-1;
				if(exploredReads[read2]) continue;
				if(type==0 && j<=(uint64_t)(readLength-minOverlap)&& compareStringInBytes(readsInt[read1],readsInt[read2],j))
				{
					ovlpLength=j;typeOfEdgeToInsert=3;
				}
				else if(type==1 && j>=(uint64_t)(minOverlap-hashStringLength) && compareStringInBytes(readsReverseInt[read1],readsReverseInt[read2],readLength-j-hashStringLength))
				{
					ovlpLength=readLength-j-hashStringLength;typeOfEdgeToInsert=0;
				}
				else if(type==2 && j<=(uint64_t)(readLength-minOverlap) && compareStringInBytes(readsInt[read1],readsReverseInt[read2],j))
				{
					ovlpLength=j;typeOfEdgeToInsert=2;
				}
				else if(type==3 && j>=(uint64_t)(minOverlap-hashStringLength) && compareStringInBytes(readsReverseInt[read1],readsInt[read2],readLength-j-hashStringLength))
				{
					ovlpLength=readLength-j-hashStringLength;typeOfEdgeToInsert=1;
				}
				if(ovlpLength!=-1)
				{

					{
						totalEdgeInserted+=insertEdgeEconomy(read1,read2,ovlpLength,typeOfEdgeToInsert);
					}
				}
			}
		}
	}
	if(graphEconomy[read1]!=NULL && graphEconomy[read1][0]>1) // Read has some edges
		quicksortEdgesEconomy(graphEconomy[read1],1,graphEconomy[read1][0],1); //sort according to length
	return totalEdgeInserted*2; // counting the twin edges
}


/* ============================================================================================
   This function inserts edges in the economy graph.
   ============================================================================================ */
int insertEdgeEconomy(uint64_t vertex_u,uint64_t vertex_v, int delta, int type)
{
	if(vertex_u==vertex_v) //do not want to add edge (a,a)
	{
		selfOverlappingReads[vertex_u]=1;
		return 0;
	}
	// FFFFFFFFFFXXXXXX represents the destination node.
	// XXXXXXXXXXCXXXXX represents the type of edge.
	// XXXXXXXXXXXFFFFF represents the length of edge.
	// XXXXXXXXXX3XXXXX is used for marking the edges.
	uint64_t edgeuv=0,edgevu=0,u=vertex_u,v=vertex_v,typeForward64=type,length64=delta,typeReverse64=reverseEdgeType(type);
	edgeuv=(v<<24)|(typeForward64<<22)|(length64);
	edgevu=(u<<24)|(typeReverse64<<22)|(length64);
	if(graphEconomy[vertex_u]==NULL) // First Edge of vertex_u
	{
		if((graphEconomy[vertex_u]=(uint64_t *)malloc(11*sizeof(uint64_t )))==NULL) printError(MEM_ALLOC, "graphEconomy[vertex_u]");
		graphEconomy[vertex_u][0]=0;
	}
	else if(graphEconomy[vertex_u][0]%10==0)
	{
		if((graphEconomy[vertex_u]=(uint64_t *)realloc(graphEconomy[vertex_u],(graphEconomy[vertex_u][0]+11)*sizeof(uint64_t )))==NULL) printError(MEM_ALLOC, "graphEconomy[vertex_u]");
	}


	if(graphEconomy[vertex_v]==NULL) // First Edge of vertex_v
	{
		if((graphEconomy[vertex_v]=(uint64_t *)malloc(11*sizeof(uint64_t )))==NULL) printError(MEM_ALLOC, "graphEconomy[vertex_v]");
		graphEconomy[vertex_v][0]=0;
	}
	else if(graphEconomy[vertex_v][0]%10==0)
	{
		if((graphEconomy[vertex_v]=(uint64_t *)realloc(graphEconomy[vertex_v],(graphEconomy[vertex_v][0]+11)*sizeof(uint64_t )))==NULL) printError(MEM_ALLOC, "graphEconomy[vertex_v]");
	}
	graphEconomy[vertex_u][++graphEconomy[vertex_u][0]]=edgeuv;
	graphEconomy[vertex_v][++graphEconomy[vertex_v][0]]=edgevu;
	return 1;

}

/* ============================================================================================
   This function inserts edges in the overlap graph.
   ============================================================================================ */
void insertEdgeInGraph(uint64_t from,uint64_t to,int overlap_size,int type)
{
	static int allocateFlag=1;
	uint64_t i;
	uint64_t *listForward,*listReverse;
	if(allocateFlag==1) // only allocate memory for the first time.
	{
		if((graph=(Edge **)malloc((numberOfUniqueReads+1)*sizeof(Edge *)))==NULL) printError(MEM_ALLOC, "graph");
		for(i=0;i<=numberOfUniqueReads;i++)
			graph[i]=NULL;
		allocateFlag=0;
	}
	if((listForward=(uint64_t *)malloc((1)*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC, "listForward");
	if((listReverse=(uint64_t *)malloc((1)*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC, "listReverse");
	listForward[0]=0;listReverse[0]=0;
	insertEdge(llabs(from),llabs(to),type,listForward,listReverse,0,overlap_size);
}

/* ============================================================================================
   This function frees the hash table and bit representation of the reads.
   ============================================================================================ */
void freeHashTable(void)
{
	uint64_t i;
	for(i=0;i<(uint64_t)sizeOfHashTable;i++) // free the reads represented by uint8_t
		if(listOfHashReadsInt[i]!=NULL)
			free(listOfHashReadsInt[i]);
	free(listOfHashReadsInt);
}


/* ==================================================================================
   This function prints an error message
   ================================================================================== */
void printError(int err_no, char *mssg){
  /* In: Error number and message. */
	switch (err_no)
	{
		case OPEN_FILE:
			myfprintf("Could not open file: %s.\n",mssg);
			deleteTemporaryFiles();
			exit(EXIT_FAILURE);
			break;
		case MEM_ALLOC:
			myfprintf("Malloc failed: %s.\n",mssg);
			deleteTemporaryFiles();
			exit(EXIT_FAILURE);
			break;
		case OTHERS:
			myfprintf("%s.\n",mssg);
			deleteTemporaryFiles();
			exit(EXIT_FAILURE);
			break;
	}
}


/* ============================================================================================
   This function finds multiple path with one simple and another composite edge from a node in
   the reduced overlap graph. simple edge of C2S1 is deleted by this function.
   ============================================================================================ */
uint64_t removeC2S1Edges(void)
{
	myfprintf("\nIn function removeC2S1Edges().\n");
	uint64_t i,simple_edge_deleted_in_this_phase=0,com_out_edge_counter=0,sim_out_edge_counter=0,com_in_edge_counter=0,sim_in_edge_counter=0;
	float flow_in,flow_out,flow_sim;
	Edge *v,*v_next;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graph[i]!=NULL)
		{
			com_out_edge_counter=0; sim_out_edge_counter=0; com_in_edge_counter=0; sim_in_edge_counter=0;
			flow_in=0;flow_out=0;flow_sim=0;
			for(v=graph[i];v!=NULL;v=v->next)
			{
				if(v->typeOfEdge==2 || v->typeOfEdge==3)// if out-edge
				{
					if(v->listOfReads[0]!=0) // composite edges
					{
						com_out_edge_counter++;
						flow_out+=v->flow;
					}
					else//simple edges
					{
						sim_out_edge_counter++;;
						flow_out=v->flow;
						flow_sim+=v->flow;
					}
				}
				else if(v->typeOfEdge==0 || v->typeOfEdge==1)// if in-edge
				{
					if(v->listOfReads[0]!=0)//composite edges
					{
						com_in_edge_counter++;
						flow_in+=v->flow;
					}
					else //simple edges
					{
						sim_in_edge_counter++;
						flow_in+=v->flow;
						flow_sim+=v->flow;
					}
				}
			}

			if(com_out_edge_counter==1 && com_in_edge_counter==1 && (sim_out_edge_counter==1 || sim_in_edge_counter==1) && flow_in==flow_out && flow_sim==0 && flow_in==1)
			{
				if(graph[i]!=NULL)
				{

					for(v=graph[i];v!=NULL;v=v_next)
					{
						v_next=v->next;
						if(v->listOfReads[0]==0)// simple edge
						{

							deleteEdge(v->twinEdge);
							deleteEdge(v);
							simple_edge_deleted_in_this_phase++;
						}
					}
				}
			}
		}
	}
	myfprintf("Simple edge deleted: %10"PRIu64"\n",simple_edge_deleted_in_this_phase);
	if(simple_edge_deleted_in_this_phase!=0)
		removeC2S1Edges();
	return 0;
}


/* ===========================================================================================
   This function will remove unnecessary mate pairs (Mate pairs that are on the same edge).
   ============================================================================================ */
void mapReadsToEdges(void)
{
	myfprintf("\nIn function mapReadsToEdges().\n");
	uint64_t i,flag,index;
	static int allocateFlag=1;
	if(allocateFlag==1)
	{
		if((read_to_edge = (struct readToEdgeStructure **)malloc((numberOfUniqueReads+1)*sizeof(struct readToEdgeStructure *)))==NULL) printError(MEM_ALLOC, "Read to Edge mapping");
		for(i=0;i<=numberOfUniqueReads;i++)
			read_to_edge[i]=NULL;
		allocateFlag=0;
	}
	time_t seconds_s=time (NULL);
	Edge *v;
	struct matePairStructure *x;
	struct readToEdgeStructure *delete_a,*delete_b,*edg,*edg_temp;
	for(i=0;i<=numberOfUniqueReads;i++) // delete previous mapping of read to edge
	{
		if(read_to_edge[i]!=NULL) // Nodes
		{
			for(delete_a=read_to_edge[i];delete_a!=NULL;delete_a=delete_b)
			{
				delete_b=delete_a->next;
				free(delete_a->locationForward);
				free(delete_a->locationReverse);
				free(delete_a);

			}
		}
		read_to_edge[i]=NULL;
	}
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graph[i]!=NULL)
		{
			for(v=graph[i];v!=NULL;v=v->next)
			{
				if(i<=v->ID)
				{
					for(index=1;index<=(uint64_t)(v->listOfReads[0]);index++)
					{
						flag=1;
						for(edg_temp=read_to_edge[v->listOfReads[index]>>24];edg_temp!=NULL;edg_temp=edg_temp->next)
						{
							if(edg_temp->edge==v || edg_temp->edge==v->twinEdge)
							{
								flag=0; //already present
								break;
							}
						}
						if(flag==1) // Not present in the list already.
						{
							if((edg=(struct readToEdgeStructure *)malloc(sizeof(struct readToEdgeStructure)))==NULL) printError(MEM_ALLOC, "read to edge");
							edg->edge=v;
							edg->locationForward=NULL;
							edg->locationReverse=NULL;
							edg->next=read_to_edge[v->listOfReads[index]>>24];
							read_to_edge[v->listOfReads[index]>>24]=edg;

						}
					}
				}
			}
		}
	}
	struct readToEdgeStructure *a,*b;
	for(i=1;i<=numberOfUniqueReads;i++) // for each read i
	{
		for(x=mate_pair[i];x!=NULL;x=x->next) // for each matepair x of read i
		{
			if(x->flag==1)
			{
				for(a=read_to_edge[i];a!=NULL;a=a->next)
				{
					for(b=read_to_edge[x->ID];b!=NULL;b=b->next)
					{
						if(a->edge==b->edge || a->edge==b->edge->twinEdge) //both same edge
							x->flag=0; //we do not delete it, just mark it
					}
				}
			}
		}
	}
	myfprintf("Function mapReadsToEdges() in %4ld sec.\n", time(NULL)-seconds_s);
}



/* ===========================================================================================
   This function will delete simple edges with long delta for which the probability of not
   having any other read in the interval delta is <= P
   ============================================================================================ */
void deleteSimpleEdgesWithLongDelta(void)
{
	time_t seconds_s=time (NULL);
	myfprintf("\nIn function deleteSimpleEdgesWithLongDelta().\n");
	uint64_t minDeltaToDelete=100000,i;
	unsigned int k;
	Edge *v,*next_v;
	float probability;
	for(k=1;k<=readLength-minOverlap;k++)
	{
		probability=pow((float)((float)(genomeSize-k)/(float)genomeSize),(float)numberOfReads);
		if(probability<=P)
		{
			minDeltaToDelete=k;
			myfprintf(" Min delta to delete: %8"PRIu64"\n         Probability: %8.6f\n",minDeltaToDelete,probability);
			break;
		}
	}
	uint64_t deletedSimpleEdge=0;
	if(minDeltaToDelete<100000)
	{
		for(i=1;i<=numberOfUniqueReads;i++)
		{
			if(graph[i]!=NULL)
			{
				for(v=graph[i];v!=NULL;v=next_v)
				{
					next_v=v->next;
					if(v->listOfReads[0]==0 && v->lengthOfEdge>=minDeltaToDelete && v->flow==0 && (v->isReducible==1 || v->isReducible==-2))
					{
							deletedSimpleEdge++;
							deleteEdge(v->twinEdge);
							deleteEdge(v);
					}
				}
			}
		}
	}
	myfprintf("Simple edges deleted: %10"PRIu64"\nFunction deleteSimpleEdgesWithLongDelta() in %4ld sec.\n",deletedSimpleEdge,time(NULL)-seconds_s);
}


/* ============================================================================================
   Deletes simple edges with long length.
   ============================================================================================ */

void deleteDelta(void)
{
	time_t seconds_s=time (NULL);
	myfprintf("\nIn function deleteDelta().\n");
	uint64_t deletedSimpleEdge=0,i;
	Edge *v,*next_v;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graph[i]!=NULL)
		{
			for(v=graph[i];v!=NULL;)
			{
				next_v=v->next;
				if(v->listOfReads[0]==0 && v->lengthOfEdge>=10 && (v->isReducible==1 || v->isReducible==-2))
				{
						deletedSimpleEdge++;
						deleteEdge(v->twinEdge);
						deleteEdge(v);
				}
				v=next_v;
			}
		}
	}
	myfprintf("Total %"PRIu64" simple edges delete in function deleteDelta() in %4ld sec.\n",deletedSimpleEdge,time(NULL)-seconds_s);
}


/* ===========================================================================================
   This function will remove dead ends in the graph. These are the reads with errors.
   ============================================================================================ */
uint64_t removeDeadEnds(void)
{
	time_t seconds_s=time (NULL);
	myfprintf("\nIn function removeDeadEnds().\n");
	uint64_t i,deleted_reads_with_error=0,flag=0;
	Edge *v,*v_next;
	int in_edge,out_edge;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graph[i]!=NULL)
		{
			in_edge=0;
			out_edge=0;
			flag=0;
			for(v=graph[i];v!=NULL;v=v->next)
			{
				if(v->listOfReads[0]!=0) //composite edge
				{
					int count=0;
					if(v->listOfReads[0]!=0)
						count=v->listOfReads[0];
					if(count>5) //composite edges with only few reads in it
					{
						flag=1;
						break;
					}
				}
				if(v->from_ID==v->ID) // if there is a loop
				{
						flag=1;
						break;
				}

				if(v->typeOfEdge==0 || v->typeOfEdge==1)
					in_edge++;
				else
					out_edge++;
			}
			if(flag==0 && ((in_edge==0 && out_edge>0)||(in_edge>0 && out_edge==0))) // only one type of edge
			{

				for(v=graph[i];v!=NULL;v=v_next)
				{
					v_next=v->next;
					deleteEdge(v->twinEdge);
					deleteEdge(v);
				}
				deleted_reads_with_error++;
			}
		}
	}
	myfprintf("Simple edges deleted: %10"PRIu64"\nFunction removeDeadEnds() in %4ld sec.\n",deleted_reads_with_error,time(NULL)-seconds_s);
	return deleted_reads_with_error;
}


/* ===========================================================================================
   This function will remove bubbles in the graph. These are the reads with errors.
   ============================================================================================ */
uint64_t removeBubbles(void)
{
	time_t seconds_s=time (NULL);
	myfprintf("\nIn function removeBubbles().\n");
	uint64_t i,deleted_reads_with_error=0;
	Edge *v;
	int in_edge,out_edge;
	Edge *in_edge_p=NULL,*out_edge_p=NULL,*edge_other=NULL;
	uint64_t a,b;
	int64_t delta_1,delta_2, no_1,no_2;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graph[i]!=NULL)
		{
			in_edge=0;
			out_edge=0;
			for(v=graph[i];v!=NULL;v=v->next)
			{
				if(v->typeOfEdge==0 || v->typeOfEdge==1)
				{
					in_edge++;
					in_edge_p=v->twinEdge;

				}
				else
				{
					out_edge++;
					out_edge_p=v;
				}
			}
			if(in_edge==1 && out_edge==1)
			{
				delta_1=in_edge_p->lengthOfEdge+out_edge_p->lengthOfEdge;
				a=in_edge_p->from_ID;
				b=out_edge_p->ID;
				for(v=graph[a];v!=NULL;v=v->next)
				{
					if(v->ID==b)
					{
						edge_other=v;
						delta_2=v->lengthOfEdge;
						no_1=1;no_2=0;
						if(in_edge_p->listOfReads[0]!=0) no_1+=in_edge_p->listOfReads[0];
						if(out_edge_p->listOfReads[0]!=0) no_1+=out_edge_p->listOfReads[0];
						if(edge_other->listOfReads[0]!=0) no_2+=edge_other->listOfReads[0];
						if(abs(delta_1-delta_2)<50) //similar length
						{
							if(no_1<(no_2/2))
							{
								deleteEdge(in_edge_p->twinEdge);
								deleteEdge(in_edge_p);
								deleteEdge(out_edge_p->twinEdge);
								deleteEdge(out_edge_p);
								deleted_reads_with_error++;

							}
							if(no_2<(no_1/2))
							{
								deleteEdge(edge_other->twinEdge);
								deleteEdge(edge_other);
								deleted_reads_with_error++;
							}
						}
						break;
					}
				}
			}
		}
	}
	myfprintf("Simple edges deleted: %10"PRIu64"\nFunction removeBubbles() in %4ld sec.\n",deleted_reads_with_error,time(NULL)-seconds_s);
	return deleted_reads_with_error;
}


/* ===========================================================================================
   This function simplifies in-tree and out-tree.
   ============================================================================================ */
uint64_t reduceTrees(void)
{
	time_t seconds_s=time (NULL);
	myfprintf("\nIn function reduceTrees().\n");
	uint64_t i,j,in_edge,out_edge,node_removed=0,sim_in,sim_out;
	float flow_in,flow_out;
	Edge *u, **in_edges, **out_edges;
	if((in_edges=(Edge **)malloc(1000*sizeof(Edge *)))==NULL) printError(MEM_ALLOC, "in_edges");
	if((out_edges=(Edge **)malloc(1000*sizeof(Edge *)))==NULL) printError(MEM_ALLOC, "out_edges");
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graph[i]!=NULL)
		{
			flow_in=0; flow_out=0;
			in_edge=0; out_edge=0;
			sim_in=0; sim_out=0;
			for(u=graph[i];u!=NULL;u=u->next)
			{
				if((u->typeOfEdge==0 || u->typeOfEdge==1))
				{
					if(u->flow>=1 && u->flow ==(float)(int)u->flow && u->ID!=u->from_ID)// && u->frist_read!=NULL)//strlen(u->string)>=readLength) // all integer flow >=1
					{
						in_edges[in_edge++]=u;
						flow_in+=u->flow;
						if(u->listOfReads[0]==0)
						{
							sim_in++;
						}
					}
					else
					{
						in_edge=100;
						out_edge=100;
						break;
					}

				}
				else if((u->typeOfEdge==2 || u->typeOfEdge==3))
				{
					if(u->flow>=1  && u->flow==(float)(int)u->flow && u->ID!=u->from_ID)// && u->frist_read!=NULL)//strlen(u->string)>=readLength) // all integer flow >=1
					{
						out_edges[out_edge++]=u;
						flow_out+=u->flow;
						if(u->listOfReads[0]==0)
						{
							sim_out++;
						}
					}
					else
					{
						in_edge=100;
						out_edge=100;
						break;
					}
				}
			}
			if(sim_in+sim_out>1)
				continue;
			if(((in_edge>1 && out_edge==1) || (in_edge==1 && out_edge>1)) && flow_in==flow_out)
			{
				if(in_edge>1 && out_edge==1)
				{
					for(j=0;j<in_edge;j++)
					{
						if(mergeEdges(in_edges[j]->twinEdge,out_edges[0],2))
							node_removed++;
					}
				}
				else if(in_edge==1 && out_edge>1)
				{
					for(j=0;j<out_edge;j++)
					{
						if(mergeEdges(in_edges[0]->twinEdge,out_edges[j],2))
							node_removed++;
					}
				}
			}
		}
	}
	free(in_edges);free(out_edges);
	myfprintf("Total %"PRIu64" pair of edges merged.\nFunction reduceTrees() in %4ld sec.\n",node_removed,time(NULL)-seconds_s);
	return node_removed;
}


/* ============================================================================================
   We may have triangle in the graph after removing composite path and removing cycle. This
   function will remove triangles.
   ============================================================================================ */
uint64_t removeTriangles(void)
{
	myfprintf("\nIn function removeTriangles().\n");
	uint64_t i,total_removed=0,unable_to_remove=0;
	time_t seconds_s,seconds_e;
	seconds_s=time (NULL);
	for (i=1;i<=numberOfUniqueReads;i++)
	{
		if(graph[i]!=NULL && graph[i]->next!=NULL && graph[i]->next->next==NULL) // Only two neighbours
		{
			if(graph[i]->flow==graph[i]->next->flow) // no edges (a,a) and not (a,b) (a,b)
			{
				if(mergeEdges(graph[i]->twinEdge,graph[i]->next,2))
					total_removed++;
				else
					unable_to_remove++;
			}
		}
	}
	seconds_e=time(NULL);
	myfprintf("   Nodes Removed: %8"PRIu64"\nUnable to remove: %8"PRIu64"\nFunction removeTriangles() in %4ld sec.\n",total_removed,unable_to_remove,seconds_e-seconds_s);
	return total_removed;
}


/* ============================================================================================
   Two edges (e1, e2) and (e2,e3) is merged to (e1,e3) here.
   ============================================================================================ */
int mergeEdges(Edge *edge1,Edge *edge2, int type_of_merge)
{
	int return_value,type=combinedEdgeType(edge1,edge2);
	float flow=0.0;
	if(edge1==edge2 || edge1==edge2->twinEdge || edge1->isReducible==0 || edge2->isReducible==0 || type==-1)
		return 0;
	if(edge1->flow==0.0 && edge2->flow==0.0) // when we call this function before flow
	{
		flow=0.0;
	}
	else
	{
		if(type_of_merge==1) // only one unit of flow is used from both edges. this is used when we resolve pairs based on support
		{
			flow=1.0;
		}
		else if(type_of_merge==2) // we take minimum flow of the two edges. type_of_merge=2 when we resolve in and out tree and also when we resolve composite path
		{
			if(edge1->flow<=edge2->flow) // minimum flow in edge1
				flow=edge1->flow;
			else // minimum flow in edge2
				flow=edge2->flow;
		}
	}
	return_value=insertEdge(edge1->from_ID,edge2->ID,type,getListOfReads(edge1,edge2),getListOfReads(edge2->twinEdge,edge1->twinEdge),flow,edge1->lengthOfEdge+edge2->lengthOfEdge);
	if(edge1->flow<=flow)
	{
		/*if(edge1->listOfReads[0]==0 && type_of_merge==1) // do not delete simple edges
		{
			edge1->twinEdge->flow=0.0; edge1->flow=0.0;
		}
		else*/
		{
			deleteEdge(edge1->twinEdge); deleteEdge(edge1);
		}
	}
	else
	{
		edge1->twinEdge->flow-=flow; edge1->flow-=flow;
	}

	if(edge2->flow<=flow)
	{
		/*if(edge2->listOfReads[0]==0 && type_of_merge==1)
		{
			edge2->twinEdge->flow=0.0; edge2->flow=0.0;
		}
		else*/
		{
			deleteEdge(edge2->twinEdge); deleteEdge(edge2);
		}
	}
	else
	{
		edge2->twinEdge->flow-=flow; edge2->flow-=flow;
	}
	return return_value;
}



/* ============================================================================================
   Merges the list of reads in two edges.
   MSB 40 bits read number
   Then 1 bit for orentation
   Then 1 bit for flag
   Then 11 bits for distance from previous
   Then 11 bits for distance to next
   ============================================================================================ */

uint64_t *getListOfReads(Edge *edge1,Edge *edge2)
{
	uint64_t *returnList,node=edge1->ID,orientation=0,distFromPrev,distToNext,flag=0;
	int i;
	if(edge1->typeOfEdge==1 || edge1->typeOfEdge==3) orientation=1;

	if((returnList=(uint64_t *)malloc((edge1->listOfReads[0]+edge2->listOfReads[0]+2)*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC, "returnList");
	returnList[0]=edge1->listOfReads[0]+edge2->listOfReads[0]+1;

	if(edge1->listOfReads[0]==0)
		 distFromPrev=edge1->lengthOfEdge; // If edge1 is simple edge
	else
		distFromPrev=(edge1->listOfReads[edge1->listOfReads[0]]&0X00000000000007FF); // If edge1 is composite edge

	if(edge2->listOfReads[0]==0)
		distToNext=edge2->lengthOfEdge; // If edge2 is simple edge
	else
		distToNext=(edge2->listOfReads[1]&0X00000000003FF800)>>11; // If edge2 is composite edge

	for(i=1;i<=(int)(edge1->listOfReads[0]);i++) // Copy list from edge1
		returnList[i]=edge1->listOfReads[i];

	returnList[edge1->listOfReads[0]+1]=(node<<24) | (orientation<<23) | (flag<<22) | (distFromPrev<<11) | distToNext; // Insert current node

	for(i=1;i<=(int)(edge2->listOfReads[0]);i++) // Copy list from edge2
		returnList[edge1->listOfReads[0]+i+1]=edge2->listOfReads[i];
	return returnList; // Return the new list
}
/* ============================================================================================
   Two edges (e1, e2) and (e2,e3) is merged to (e1,e3) here.
   ============================================================================================ */
int mergeDisconnectedEdges(Edge *edge1,Edge *edge2,int gap)
{
	int returnValue,type=0;
	uint64_t length1,length2;
	if((edge1->typeOfEdge==0 || edge1->typeOfEdge==1) && (edge2->typeOfEdge==0 || edge2->typeOfEdge==2)) type=0;
	else if((edge1->typeOfEdge==0 || edge1->typeOfEdge==1) && (edge2->typeOfEdge==1 || edge2->typeOfEdge==3)) type=1;
	else if((edge1->typeOfEdge==2 || edge1->typeOfEdge==3) && (edge2->typeOfEdge==0 || edge2->typeOfEdge==2)) type=2;
	else if((edge1->typeOfEdge==2 || edge1->typeOfEdge==3) && (edge2->typeOfEdge==1 || edge2->typeOfEdge==3)) type=3;
	uint64_t *returnList1=getListOfReadsInDisconnectedEdges(edge1,edge2,gap,&length1);
	uint64_t *returnList2=getListOfReadsInDisconnectedEdges(edge2->twinEdge,edge1->twinEdge,gap,&length2);
	if(length1!=length2) myfprintf("Checkfor length: %d %d\n",length1, length2);
	returnValue=insertEdge(edge1->from_ID,edge2->ID,type,returnList1,returnList2,1.0,length1);
	if(edge1->flow<=1.0)
	{
		deleteEdge(edge1->twinEdge);deleteEdge(edge1);
	}
	else
	{
		edge1->twinEdge->flow-=1.0; edge1->flow-=1.0;
	}
	if(edge2->flow<=1.0)
	{
		deleteEdge(edge2->twinEdge);deleteEdge(edge2);
	}
	else
	{
		edge2->twinEdge->flow-=1.0; edge2->flow-=1.0;
	}
	return returnValue;
}


/* ============================================================================================
   Merges the list of reads in two edges. when the edges does not have common node.
   0XFFFFFFFF00000000 read number.
   0X00000000FF000000 orientation of the read.
   0X0000000000FF0000 distance from previous read.
   0X000000000000FF00 distance to next read.
   0X00000000000000FF unused.
   ============================================================================================ */

uint64_t *getListOfReadsInDisconnectedEdges(Edge *edge1,Edge *edge2,int gap,uint64_t *length)
{
	uint64_t *returnList,node1,orient1,distPrev1,distNext1,node2,orient2,distPrev2=0,distNext2=0;
	int i,index=1,flag; char *string1,*string2;
	node1=edge1->ID;
	if(edge1->typeOfEdge==1 || edge1->typeOfEdge==3)
	{
		orient1=1; string1=bytesToChars(readsInt[node1]);
	}
	else
	{
		orient1=0; string1=bytesToChars(readsReverseInt[node1]);
	}
	distPrev1=(edge1->listOfReads[edge1->listOfReads[0]]&0X00000000000007FF); // If edge1 is composite edge

	node2=edge2->from_ID;
	if(edge2->typeOfEdge==2 || edge2->typeOfEdge==3)
	{
		orient2=1; string2=bytesToChars(readsInt[node2]);
	}
	else
	{
		orient2=0; string2=bytesToChars(readsReverseInt[node2]);
	}
	distNext2=(edge2->listOfReads[1]&0X00000000003FF800)>>11; // If edge2 is composite edge

	if(node1==node2) // only one read will be inserted
	{
		flag=1;
		distNext1=distNext2;
		*length=edge1->lengthOfEdge+edge2->lengthOfEdge;
	}
	else // both the nodes should be inserted
	{
		flag=2;
		int len=stringOverlapSize(string1,string2);
		if(len>0) // node1 and node 2 overlaps
		{
			distNext1=len; distPrev2=len;
			*length=edge1->lengthOfEdge+edge2->lengthOfEdge+len;
		}
		else // node1 and node2 does not overlap
		{
			if(gap<0) // Concatenate the string
			{
				distNext1=readLength; distPrev2=readLength;
				*length=edge1->lengthOfEdge+edge2->lengthOfEdge+distNext1;

			}
			else // Add N's between the reads.
			{
				distNext1=gap+readLength; distPrev2=gap+readLength;
				//distNext1=readLength; distPrev2=readLength;
				*length=edge1->lengthOfEdge+edge2->lengthOfEdge+distNext1;
			}
		}
	}

	if((returnList=(uint64_t *)malloc((edge1->listOfReads[0]+edge2->listOfReads[0]+flag+1)*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC, "returnList");
	returnList[0]=edge1->listOfReads[0]+edge2->listOfReads[0]+flag;
	for(i=1;i<=(int)(edge1->listOfReads[0]);i++) // Copy list from edge1
		returnList[index++]=edge1->listOfReads[i];
	returnList[index++]=(node1<<24) | (orient1<<23) | (distPrev1<<11) | (distNext1);// | gap; // Insert current node
	if(flag==2)
		returnList[index++]=(node2<<24) | (orient2<<23) | (distPrev2<<11) | distNext2; // Insert current node

	for(i=1;i<=(int)(edge2->listOfReads[0]);i++) // Copy list from edge2
		returnList[index++]=edge2->listOfReads[i];
	free(string1);free(string2);
	return returnList; // Return the new list
}



int stringOverlapSize(char *string1, char *string2)
{
	int i,j,misMatchCount;
	for(i=0;i<(int)(readLength-10);i++)
	{
		misMatchCount=0;
		for(j=0;j<(int)(readLength-i);j++)
		{
			if(string1[i+j]!=string2[j])
			{
				misMatchCount++;
				if(misMatchCount>((int)(readLength-i)/10)+1)
				{
					break;
				}
			}
		}
		if(j==(int)(readLength-i)) // at least 10 base pairs are same
		{
			return i;
		}
	}
	return -1; // do not overlap
}




/* ============================================================================================
   This function finds a loop like (v1->v2->v2->v3) and replace it by v1->v3.
   ============================================================================================ */
uint64_t reduceLoops(void)
{
	myfprintf("\nIn function reduceLoops().\n");
	time_t seconds_s=time(NULL);
	uint64_t i,number_of_loop_resolved=0;
	Edge *v,*ab,*bb,*bc;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		ab=NULL;bb=NULL,bc=NULL;
		if(graph[i]!=NULL && graph[i]->next!=NULL && graph[i]->next->next!=NULL && graph[i]->next->next->next!=NULL && graph[i]->next->next->next->next==NULL && (graph[i]->ID==i || graph[i]->next->ID==i || graph[i]->next->next->ID==i || graph[i]->next->next->next->ID==i))
		{
			for(v=graph[i];v!=NULL;v=v->next)
			{
				if(v->ID==i)
					bb=v;
				else if(ab==NULL)
					ab=v->twinEdge;
				else if(bc==NULL)
					bc=v;
			}
			if(bb==NULL || ab==NULL || bc==NULL)
				continue;
			if(DEBUGGING)
			{
				myfprintf("loop found in node %"PRIu64"\n",i);
			}
			if(bb->typeOfEdge==1||bb->typeOfEdge==2) // we cannot simplify if b<--------->b or b>--------<b
			{
				if(DEBUGGING)
				{
					myfprintf("Cannot resolve because edge (%"PRIu64",%"PRIu64") has type %d\n",bb->from_ID,bb->ID,bb->typeOfEdge);
				}
				continue;
			}
			if((ab->typeOfEdge==1 || ab->typeOfEdge==3)&& (bb->typeOfEdge==3))
			{
				if(mergeEdges(ab,bb,1))
					number_of_loop_resolved++;
			}
			else if((ab->typeOfEdge==0 || ab->typeOfEdge==2)&&(bb->typeOfEdge==0))
			{
				if(mergeEdges(ab,bb,1))
					number_of_loop_resolved++;
			}
			else if((ab->typeOfEdge==1 || ab->typeOfEdge==3)&& (bb->twinEdge->typeOfEdge==3))
			{
				if(mergeEdges(ab,bb->twinEdge,1))
					number_of_loop_resolved++;
			}
			else if((ab->typeOfEdge==0 || ab->typeOfEdge==2)&&(bb->twinEdge->typeOfEdge==0))
			{
				if(mergeEdges(ab,bb->twinEdge,1))
					number_of_loop_resolved++;
			}
		}
	}
	myfprintf("Loops resolved: %10"PRIu64"\nFunction reduceLoops() in %4ld sec.\n",number_of_loop_resolved, time(NULL)-seconds_s);
	return number_of_loop_resolved;
}


/* ============================================================================================
   This function reads genome from "filename" and returns the string (for testing Only)
   ============================================================================================ */
char* readGenome(char* filename)
{
	myfprintf("\nIn function readGenome().\n");
	uint64_t i,j;
	char *Genome;
	FILE *fpi;
	struct stat file_info;
	if (stat(filename,&file_info) == -1) // get information about the file
	{
		myfprintf("Error opening input file: %s\n",filename); // Error opening genome file
		Genome = NULL;
  	}
	fpi=fopen(filename,"r");
	if(fpi==NULL)
	{
		myfprintf("Error opening input file: %s\n",filename);
		Genome = NULL;
	}
	Genome = (char *)calloc(file_info.st_size+1,sizeof(char));
	uint64_t count=fread(Genome,1,file_info.st_size,fpi);
	myfprintf("Genome characters count: %"PRIu64"\n", count);
	j=0;
	for(i=0;i<(uint64_t)file_info.st_size;i++)
	{
		if(Genome[i]=='A' || Genome[i]=='C' || Genome[i]=='G' || Genome[i]=='T' || Genome[i]=='a' || Genome[i]=='c' || Genome[i]=='g' || Genome[i]=='t')
		{
			Genome[j++]=Genome[i];
		}
	}
	Genome[j]='\0';
	fclose(fpi);
	return Genome;
}

uint64_t convertStringToNumber(char *str,int length)
{
	uint64_t value=0;
	int i;
	for(i=0;i<length;i++)
	{
		if(*(str+i)=='A' || *(str+i)=='a')
			value=(value<<2);
		else if(*(str+i)=='C' || *(str+i)=='c')
			value=(value<<2 )|1;
		else if(*(str+i)=='G' || *(str+i)=='g')
			value=(value<<2)|2;
		else if(*(str+i)=='T' || *(str+i)=='t')
			value=(value<<2)|3;
		else
			return 0XFFFFFFFFFFFFFFFF;

	}
	return value;
}

/* ============================================================================================
   This function convert into to string
   ============================================================================================ */
void itoa(int n, char s[])
{
	int i, sign;
	if ((sign = n) < 0)  /* record sign */
		n = -n;          /* make n positive */
	i = 0;
	do
	{       /* generate digits in reverse order */
		s[i++] = n % 10 + '0';   /* get next digit */
	} while ((n /= 10) > 0);     /* delete it */
	if (sign < 0)
		s[i++] = '-';
	s[i] = '\0';
	reverse(s);
}

/* ============================================================================================
   This function reverses a string
   ============================================================================================ */
void reverse(char s[])
{
	int i, j;
	char c;
	for (i = 0, j = strlen(s)-1; i<j; i++, j--)
	{
		c = s[i];
		s[i] = s[j];
		s[j] = c;
	}
}

/* ============================================================================================
   if string in (a,b) and (a,b) are similar and we have (b,c) and (b,d), then we merge
   (a,b)(b,c) and (a,b)(b,d)
   ============================================================================================ */
void markSimilarEdges(void)
{
	myfprintf("\nIn function markSimilarEdges().\n");
	time_t second_s=time(NULL);
	uint64_t i, score,length1,length2,max_length,gapCount,nCount;
	Edge *v,*edge1,*edge2;
	float in_flow,out_flow;
	int edges_marked=0;
	char *string1,*string2;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graph[i]!=NULL)
		{
			in_flow=0;
			out_flow=0;
			for(v=graph[i];v!=NULL;v=v->next)
			{
				if(v->typeOfEdge==1 || v->typeOfEdge==0)
					in_flow+=v->flow;
				else
					out_flow+=v->flow;
			}
			if(in_flow!=out_flow)
				continue;
			for(edge1=graph[i];edge1!=NULL;edge1=edge1->next)
			{
				if(edge1->listOfReads[0]==0) continue;
				length1=edge1->lengthOfEdge;
				for(edge2=edge1->next;edge2!=NULL;edge2=edge2->next)
				{
					if(edge2->listOfReads[0]==0) continue;
					length2=edge2->lengthOfEdge;
					max_length=maximumValue(length1,length2);
					if(edge1!=edge2 && edge1->ID==edge2->ID && edge1->ID!=i && edge1->typeOfEdge==edge2->typeOfEdge && length1<5000 && length2<5000 && llabs(length1-length2)<=10*max_length/100)
					{
						string1=stringInEdge(edge1,&gapCount,&nCount);
						string2=stringInEdge(edge2,&gapCount,&nCount);
						score=checkAlignment(1.0,-2.0,-1.0,string1,string2);
						free(string1); free(string2);
						if(score>=max_length*80/100)
						{
							edges_marked++;
							if(max_length==length2)
							{
								edge1->isReducible=2;
								edge1->twinEdge->isReducible=2;
								edge2->isReducible=-1;
								edge2->twinEdge->isReducible=-1;
							}
							else
							{
								edge1->isReducible=-1;
								edge1->twinEdge->isReducible=-1;
								edge2->isReducible=2;
								edge2->twinEdge->isReducible=2;
							}

						}
					}
				}
			}
		}
	}
	myfprintf("Edges marked: %10d\nFunction markSimilarEdges() in %4ld sec.\n",edges_marked,time(NULL)-second_s);
}

/* ============================================================================================
   This function checks similar string in edges (a,b) and (a,b). If they are very similar them
   remove one of them and add the flow to another.
   ============================================================================================ */
uint64_t deleteEdgesWithSimilarStrings()
{
	myfprintf("In function deleteEdgesWithSimilarStrings()\n");
	time_t seconds_s=time(NULL);
	uint64_t i,score,length1,length2,number_removed=0,max_length,gapCount,nCount;
	Edge *v,*edge1,*edge2;
	float in_flow,out_flow;
	char *string1, *string2;
	for(i=1;i<=numberOfUniqueReads;i++) // changed the marked similar edges
	{
		if(graph[i]!=NULL)
		{
			for(v=graph[i];v!=NULL;v=v->next)
			{
				if(v->isReducible<0)
				{
					v->isReducible=1;
					v->twinEdge->isReducible=1;
				}
			}
		}
	}
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graph[i]!=NULL)
		{
			in_flow=0;
			out_flow=0;
			for(v=graph[i];v!=NULL;v=v->next)
			{
				if(v->typeOfEdge==1 || v->typeOfEdge==0)
					in_flow+=v->flow;
				else
					out_flow+=v->flow;
			}
			if(in_flow!=out_flow)
				continue;
			for(edge1=graph[i];edge1!=NULL;edge1=edge1->next)
			{
				if(edge1->listOfReads[0]==0) continue;
				length1=edge1->lengthOfEdge;
				for(edge2=edge1->next;edge2!=NULL;edge2=edge2->next)
				{
					if(edge2->listOfReads[0]==0) continue;
					length2=edge2->lengthOfEdge;
					max_length=maximumValue(length1,length2);
					if(edge1!=edge2 && edge1->ID==edge2->ID && edge1->ID!=i && edge1->typeOfEdge==edge2->typeOfEdge && length1<5000 && length2<5000 && llabs(length1-length2)<=10*max_length/100)
					{
						string1=stringInEdge(edge1,&gapCount,&nCount);
						string2=stringInEdge(edge2,&gapCount,&nCount);
						score=checkAlignment(1.0,-2.0,-1.0,string1,string2);
						free(string1);free(string2);
						if(score>=max_length*80/100)
						{
							edge1->flow+=edge2->flow;
							edge1->isReducible=1;
							edge1->twinEdge->flow+=edge2->flow;
							edge1->twinEdge->isReducible=1;
							deleteEdge(edge2->twinEdge);
							deleteEdge(edge2);
							number_removed++;
							break;
						}
					}
				}
			}
		}
	}
	myfprintf("Edges deleted: %10"PRIu64"\nFunction deleteEdgesWithSimilarStrings() in %4ld sec.\n",number_removed,time(NULL)-seconds_s);
	return number_removed;
}


/* ============================================================================================
   Merge pair of edges according to supports.
   ============================================================================================ */

uint64_t resolvePairsOfEdges(Edge **supportedEdge1,Edge **supportedEdge2, uint64_t supported_pairs)
{
	myfprintf("\nIn function resolvePairsOfEdges().\n");
	time_t seconds_s=time(NULL);
	uint64_t i,j,k,*nodes_cannot_resolve,top_element=0,node_resolved=0;
	float flow1,flow2;
	if((nodes_cannot_resolve = (uint64_t *)malloc(defaultArraySize*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC, "nodes_cannot_resolve");
	Edge *v;
	pair_of_support=supported_pairs;
	myfprintf("\nTable: Partial list of merged pairs of edges in resolvePairsOfEdges() in iteration %d.\n-----------------------------------------------------------------------------------------------------------\n| # |      EDGE1      |  ADDRESS  | LENGTH | FLOW |      EDGE2      |  ADDRESS  | LENGTH | FLOW | SUPPORT |\n-----------------------------------------------------------------------------------------------------------\n",iteration);
	for(i=0;i<supported_pairs;i++)
	{
		if(supportedEdge1[i]==NULL || supportedEdge2[i]==NULL) // already deleted
			continue;
		for(j=0;j<top_element;j++) //not in independent set
			if(supportedEdge1[i]->ID==nodes_cannot_resolve[j])
				break;
		if(j<top_element)
			continue;
		if(support[i]>=5)
		{
			for(v=graph[supportedEdge1[i]->ID];v!=NULL;v=v->next) // Insert all neightbours in the list
			{
				for(j=0;j<top_element;j++)
				{
					if(nodes_cannot_resolve[j]==v->ID)
					{
						break;
					}
				}
				if(j==top_element && v->ID!=v->from_ID)
				{
					nodes_cannot_resolve[top_element++]=v->ID;
				}
			}
			if(node_resolved<=9) //Print only 10 pairs of deleted edges;
			{
				uint64_t a1=(uint64_t)supportedEdge1[i],a2=(uint64_t)supportedEdge1[i]->twinEdge,a3=(uint64_t)supportedEdge2[i],a4=(uint64_t)supportedEdge2[i]->twinEdge;
				a1=a1&0XFFFF;a2=a2&0XFFFF;a3=a3&0XFFFF;a3=a3&0XFFFF; a4=a4&0XFFFF; // We only print lower 16 bits of the address
				myfprintf("|%2"PRIu64" |%8"PRIu64",%8"PRIu64"|%5"PRIu64" %5"PRIu64"|%7"PRIu64" |%5.1f |%8"PRIu64",%8"PRIu64"|%5"PRIu64" %5"PRIu64"|%7"PRIu64" |%5.1f |%8d |\n",node_resolved+1,supportedEdge1[i]->from_ID, supportedEdge1[i]->ID,a1,a2,supportedEdge1[i]->lengthOfEdge,supportedEdge1[i]->flow,supportedEdge2[i]->from_ID,supportedEdge2[i]->ID,a3,a4,supportedEdge2[i]->lengthOfEdge,supportedEdge2[i]->flow,support[i]);
			}
			Edge *temp1=supportedEdge1[i],*temp1Twin=supportedEdge1[i]->twinEdge,*temp2=supportedEdge2[i],*temp2Twin=supportedEdge2[i]->twinEdge;
			flow1=supportedEdge1[i]->flow,flow2=supportedEdge2[i]->flow;
			if(mergeEdges(supportedEdge1[i],supportedEdge2[i],1))
			{
				node_resolved++;
				for(k=0;k<pair_of_support;k++) // Update the list of supported edges.
				{
					if(flow1<=1.0 && (supportedEdge1[k]==temp1 || supportedEdge1[k]==temp1Twin || supportedEdge2[k]==temp1 || supportedEdge2[k]==temp1Twin)) // first edge is deleted already
					{
						supportedEdge1[k]=NULL; supportedEdge2[k]=NULL;
					}
					if(flow2<=1.0 && (supportedEdge1[k]==temp2 || supportedEdge1[k]==temp2Twin || supportedEdge2[k]==temp2 || supportedEdge2[k]==temp2Twin)) // second edge is deleted
					{
						supportedEdge1[k]=NULL; supportedEdge2[k]=NULL;
					}
				}
			}
		}
		else break;
	}
	myfprintf("-----------------------------------------------------------------------------------------------------------\n\n");
	free(nodes_cannot_resolve);
	if(DEBUGGING)
	{
		char graph_file_name[10000],temp_number[10];
		itoa(iteration,temp_number);
		strcpy(graph_file_name,temp_number);
		printGraph(graph_file_name);
	}
	if(TEST)
	{
		printGraph("DUMMY");
		myfprintf("Function resolvePairsOfEdges() %4ld sec.\n",time(NULL)-seconds_s);
	}
	return node_resolved;
}


/* ============================================================================================
   This function will find all paths between mate_pair_1 and mate_pair_2 and will find
   if the all goes through (a,node), (node,b) for some a and b.
   ============================================================================================ */
int findPathsBetweenMatepairs(uint64_t mate_pair_1,uint64_t mate_pair_2,int type1,int type2, int library)
{
	Edge **edges,*edge1=NULL,*edge2=NULL;
	struct readToEdgeStructure *edge_ptr1,*edge_ptr2;
	uint64_t i,k,l,m,path_found=0,total_supported=0,x;
	int64_t *dist_on_edge1,*dist_on_edge2;
	int *support_flag;
	if(read_to_edge[mate_pair_1]==NULL || read_to_edge[mate_pair_2]==NULL)  // not in any edge
		return 0;
	if((edges = (Edge **)malloc(1000*sizeof(Edge *)))==NULL) printError(MEM_ALLOC, "edges");
	if((support_flag = (int *)malloc(1000*sizeof(int)))==NULL) printError(MEM_ALLOC, "support_flag");
	for(edge_ptr1=read_to_edge[mate_pair_1];edge_ptr1!=NULL;edge_ptr1=edge_ptr1->next)// find an edge contining mate_pair_1
	{
		for(edge_ptr2=read_to_edge[mate_pair_2];edge_ptr2!=NULL;edge_ptr2=edge_ptr2->next)// find an edge contining mate_pair_2
		{
			for(x=0;x<4;x++)
			{
				if(x==0)
				{
					edge1=edge_ptr1->edge;
					edge2=edge_ptr2->edge;
				}
				else if(x==1)
				{
					edge1=edge_ptr1->edge;
					edge2=edge_ptr2->edge->twinEdge;
				}
				else if(x==2)
				{
					edge1=edge_ptr1->edge->twinEdge;
					edge2=edge_ptr2->edge;
				}
				else if(x==3)
				{
					edge1=edge_ptr1->edge->twinEdge;
					edge2=edge_ptr2->edge->twinEdge;
				}
				dist_on_edge1=findDistanceOnEdge(edge1->twinEdge,mate_pair_1);
				dist_on_edge2=findDistanceOnEdge(edge2,mate_pair_2);
				if(dist_on_edge1[0]==0 || dist_on_edge2[0]==0)
				{
					free(dist_on_edge1);
					free(dist_on_edge2);
					continue;
				}
				int strt;
				strt=startingPosition(edge1->ID,edge2->from_ID);
				if(strt==-1) //not path possible
				{
					free(dist_on_edge1);
					free(dist_on_edge2);
					continue;
				}
				for(i=strt;i<totalPaths;i++)
				{
					if(allPairPaths[i][1]->from_ID>edge1->ID || (allPairPaths[i][1]->from_ID==edge1->ID && allPairPaths[i][0]->from_ID>edge2->from_ID))
						break;
					if(edge1->ID==allPairPaths[i][1]->from_ID && matechEdgeType(edge1,allPairPaths[i][1]) && allPairPaths[i][0]==edge2)
					{
						int flg=0;
						for(k=1;k<=(uint64_t)dist_on_edge1[0];k++)
						{
							for(l=1;l<=(uint64_t)dist_on_edge2[0];l++)
							{
								uint64_t distance=100000;
								//if(type1*type2*dist_on_edge1[k]*dist_on_edge2[l]>0)
								if(type1*(int)dist_on_edge1[k]<0 && type2*(int)dist_on_edge2[l]<0)
									distance=llabs(dist_on_edge1[k])+llabs(dist_on_edge2[l])+allPairPathLengths[i];
								if(distance<upperBoundOfInsert[library] && distance>lowerBoundOfInsert[library])
								{
									flg=1;
									break;
								}
							}
							if(flg==1)
							break;
						}
						if(flg==1)
						{
							edges[0]=edge1;
							for(m=1;allPairPaths[i][m]!=NULL;m++)
								edges[m]=allPairPaths[i][m];
							if(insertPath(edges,support_flag,m-1,path_found))
							{
								path_found++;
							}
							else
							{
								free(dist_on_edge1);
								free(dist_on_edge2);
								free(edges);
								free(support_flag);
								return 0;
							}
						}
					}
				}
				free(dist_on_edge1);
				free(dist_on_edge2);

			}
		}
	}
	if(path_found>0)
	{
		for(i=0;support_flag[i]>=0;i++) // move the supported pairs (flag=1) at the beginning of the stack
		{
			if(support_flag[i]==1)
			{
				supportedEdg1[total_supported]=supportedEdg1[i];
				supportedEdg2[total_supported]=supportedEdg2[i];
				total_supported++;
			}
		}
	}
	else
	{
		total_supported=0;
	}
	free(edges);
	free(support_flag);
	return total_supported;

}


/* ============================================================================================
   This function insert the path in the support pair list. if no pair of edge is supported it
   returns 0. Otherwise returns 1.
   ============================================================================================ */
uint64_t insertPath(Edge **edges,int *support_flag, uint64_t top,uint64_t path_found)
{
	uint64_t k,l;
	if(path_found==0) // if it is the first path found then store pair of edges (a,b) and (b,c) and mark them as supported
	{
		for(k=0;k<top+10;k++)
			support_flag[k]=-1;
		for(k=0;k<top;k++)
		{
			supportedEdg1[k]=edges[k];
			supportedEdg2[k]=edges[k+1];
			support_flag[k]=1;
		}
	}
	else // if not the first path
	{
		for(k=0;support_flag[k]>=0;k++) // unmark the pair of edges (a,b) and (b,c) not present in the new path
		{
			for(l=0;l<top;l++)
			{
				if((supportedEdg1[k]==edges[l] && supportedEdg2[k]==edges[l+1])||(supportedEdg1[k]==edges[l+1]->twinEdge && supportedEdg2[k]==edges[l]->twinEdge))
				{
					break;
				}
			}
			if(l==top) // not supported in new path
			{
				support_flag[k]=0;
			}
		}

	}
	for(k=0;support_flag[k]>=0;k++)
	{
		if(support_flag[k]==1)
		{
			break;
		}
	}
	if(support_flag[k]==-1) //no pair of edges supported. no need to continue more
	{
		return 0;
	}
	return 1;
}
/* ============================================================================================
   This function will simplify the graph.
   ============================================================================================ */
void updateGraph()
{
	removeTriangles();
	mapReadsToEdges();
	mapReadLocations();
}


/* ============================================================================================
   This function will delete simple edges in clique and diamonds. Only used at the end of
   this program.
   TODO: Complete this function.
   ============================================================================================ */
void cleanGraph(void)
{
	uint64_t i;
	Edge *v;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		for(v=graph[i];v!=NULL;v=v->next)
			if(v->listOfReads[0]!=0)
				break;
		if(v==NULL) // all simple edges.
		{
			while(graph[i]!=NULL) // delete such nodes with all simple edges.
			{
				deleteEdge(graph[i]->twinEdge);
				deleteEdge(graph[i]);
			}
		}
	}
}


/* ============================================================================================
   This function check if edge1 and edge2 have proper orientation.
   ============================================================================================ */
int matechEdgeType(Edge *edge1,Edge *edge2)
{
	if(((edge1->typeOfEdge==1 || edge1->typeOfEdge==3) && (edge2->typeOfEdge==2 || edge2->typeOfEdge==3))||((edge1->typeOfEdge==0 || edge1->typeOfEdge==2) && (edge2->typeOfEdge==0 || edge2->typeOfEdge==1)))
		return 1;
	return 0;
}

/* ============================================================================================
   This function checks if read is present in edge.
   ============================================================================================ */
int isUniqueInEdge(Edge *edge,uint64_t read)
{
	if(read_to_edge[read]==NULL) //not present;
		return 0;

	if((read_to_edge[read]->edge==edge || read_to_edge[read]->edge->twinEdge==edge) && read_to_edge[read]->next==NULL && graph[read]==NULL)
		return 1;
	return 0;
}

/* ============================================================================================
   This function finds distance on edge only for the first time. Otherwise use the previously
   calculated values.
   ============================================================================================ */

int64_t* findDistanceOnEdge(Edge *edge, uint64_t read)
{
	struct readToEdgeStructure *edge_ptr;
	uint64_t i;
	int64_t *distance;
	for(edge_ptr=read_to_edge[llabs(read)];edge_ptr!=NULL;edge_ptr=edge_ptr->next)
	{
		if(edge_ptr->edge==edge)
		{
			if(edge_ptr->locationForward!=NULL)
			{
				if((distance=(int64_t *) malloc((edge_ptr->locationForward[0]+1)*sizeof(int64_t)))==NULL)printError(MEM_ALLOC, "distance");
				for(i=0;i<=edge_ptr->locationForward[0];i++)
					distance[i]=edge_ptr->locationForward[i];
				return distance;
			}
		}
		else if(edge_ptr->edge==edge->twinEdge)
		{
			if(edge_ptr->locationReverse!=NULL)
			{
				if((distance=(int64_t *) malloc((edge_ptr->locationReverse[0]+1)*sizeof(int64_t)))==NULL)printError(MEM_ALLOC, "distance");
				for(i=0;i<=edge_ptr->locationReverse[0];i++)
					distance[i]=edge_ptr->locationReverse[i];
				return distance;
			}
		}
	}
	if((distance=(int64_t *)malloc(1*sizeof(int64_t)))==NULL)printError(MEM_ALLOC, "distance"); // Not found.
	distance[0]=0;
	return distance;
}

/* ============================================================================================
   This function will map each read on every edge and location of the read on the edges.
   ============================================================================================ */
void mapReadLocations(void)
{
	myfprintf("\nIn function mapReadLocations().\n");
	time_t seconds_s=time(NULL);
	uint64_t i,dist,index;
	struct readToEdgeStructure *edge_ptr;
	Edge *v;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		for(v=graph[i];v!=NULL;v=v->next)
		{
			if(v->listOfReads[0]==0) continue;
				dist=0;
			for(index=1;index<=(uint64_t)(v->listOfReads[0]);index++)
			{
				dist+=(v->listOfReads[index]&0X00000000003FF800)>>11;
				for(edge_ptr=read_to_edge[v->listOfReads[index]>>24];edge_ptr!=NULL;edge_ptr=edge_ptr->next)
				{
					if(edge_ptr->edge==v)
					{
						if(edge_ptr->locationForward==NULL) // No locationForward mapped yet. Should allocate the memory for the first time.
						{
							if((edge_ptr->locationForward=(uint64_t *) malloc((2)*sizeof(uint64_t)))==NULL)printError(MEM_ALLOC, "locationForward");
							edge_ptr->locationForward[0]=1;
							if(v->listOfReads[index]&0X0000000000800000)
								edge_ptr->locationForward[edge_ptr->locationForward[0]]=dist;
							else
								edge_ptr->locationForward[edge_ptr->locationForward[0]]=-dist;
						}
						else // Increase the allocate memory to make space for the new locationForward.
						{
							if((edge_ptr->locationForward=(uint64_t *)realloc(edge_ptr->locationForward,((edge_ptr->locationForward[0]+2))*sizeof(uint64_t)))==NULL)printError(MEM_ALLOC, "Reallocating locationForward");
							edge_ptr->locationForward[0]++;
							if(v->listOfReads[index]&0X0000000000800000)
								edge_ptr->locationForward[edge_ptr->locationForward[0]]=dist;
							else
								edge_ptr->locationForward[edge_ptr->locationForward[0]]=-dist;
						}
						break;
					}
					else if(edge_ptr->edge==v->twinEdge)
					{
						if(edge_ptr->locationReverse==NULL)// No locationReverse mapped yet. Should allocate the memory for the first time.
						{
							if((edge_ptr->locationReverse=(uint64_t *) malloc((2)*sizeof(uint64_t)))==NULL)printError(MEM_ALLOC, "locationReverse");
							edge_ptr->locationReverse[0]=1;
							if(v->listOfReads[index]&0X0000000000800000) // Check the orientation
								edge_ptr->locationReverse[edge_ptr->locationReverse[0]]=dist;
							else
								edge_ptr->locationReverse[edge_ptr->locationReverse[0]]=-dist;

						}
						else // Increase the allocate memory to make space for the new locationReverse.
						{
							if((edge_ptr->locationReverse=(uint64_t *)realloc(edge_ptr->locationReverse,((edge_ptr->locationReverse[0]+2))*sizeof(uint64_t)))==NULL)printError(MEM_ALLOC, "Reallocating locationReverse");
							edge_ptr->locationReverse[0]++;
							if(v->listOfReads[index]&0X0000000000800000) // Check the orientation
								edge_ptr->locationReverse[edge_ptr->locationReverse[0]]=dist;
							else
								edge_ptr->locationReverse[edge_ptr->locationReverse[0]]=-dist;
						}
						break;
					}
				}
			}
		}
	}
	myfprintf("Function mapReadLocations() in %4ld sec.\n",time(NULL)-seconds_s);
}


/* ============================================================================================
   This function will find all pair of feasible path in the graph starting from node.
   Backtracking used by recursion.
   ============================================================================================ */
uint64_t exploreGraph(uint64_t node,int level) // Recursive.
{
	static int allocationFlag=1;
	static uint64_t pathSaved,*pathLengths;
	static Edge **pathEdges;
	Edge *v;
	if(level>10) return 0; // do not go very deep.
	if(allocationFlag) // Allocated when called for the first time.
	{
		if((pathLengths=(uint64_t *)malloc(defaultArraySize*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC, "pathLengths");
		if((pathEdges=(Edge **)malloc(defaultArraySize*sizeof(Edge *)))==NULL) printError(MEM_ALLOC, "pathEdges");
		allocationFlag=0; pathSaved=0;
	}
	for(v=graph[node];v!=NULL;v=v->next)
	{
		if(v->isReducible<=0) // do not want to explore marked edges
			continue;
		if(level==1)
		{
			if(v->listOfReads[0]!=0 && v->flow>0) // only composite edges with flow
			{
				pathEdges[level]=v; pathLengths[level]=0;
				savePath(pathEdges,pathLengths[level],level);  pathSaved++;
				exploreGraph(v->ID,level+1);
			}
		}
		else if(matechEdgeType(pathEdges[level-1],v) && checkFlow(pathEdges,level-1,v) && (pathLengths[level-1]+pathEdges[level-1]->lengthOfEdge)<maximumUpperBoundOfInsert) // go deeper
		{
			pathEdges[level]=v; pathLengths[level]=pathLengths[level-1]+pathEdges[level-1]->lengthOfEdge;
			savePath(pathEdges,pathLengths[level],level); pathSaved++;
			exploreGraph(v->ID,level+1);
		}

	}
	if(level==1) // deallocate.
	{
		free(pathLengths);free(pathEdges);allocationFlag=1;
	}
	return pathSaved;
}

/* ============================================================================================
   This function will find paths from each node in the graph, one by one
   ============================================================================================ */
uint64_t findAllPairPaths(void)
{
	myfprintf("\nIn function findAllPairPaths().\n");
	time_t seconds_s=time (NULL);
	uint64_t mate_pair_1,mate_pair_2, tos=0,total_supported_pairs,i,j,l,index,returnValue=0;
	int64_t *dummy;
	struct matePairStructure *w;
	Edge *v,**supportedEdge1,**supportedEdge2;

	struct readToEdgeStructure *a;
	int *flag;//this is the flag whether we already found paths of the mate pairs or not.
	if((flag=(int *)malloc((numberOfUniqueReads+1)*sizeof(int)))==NULL) printError(MEM_ALLOC,"flag");
	for(i=1;i<=numberOfUniqueReads;i++)
		flag[i]=1;
	if((supportedEdge1=(Edge**)malloc(defaultArraySize*sizeof(Edge*)))==NULL) printError(MEM_ALLOC, "supportedEdge1");
	if((supportedEdge2=(Edge**)malloc(defaultArraySize*sizeof(Edge*)))==NULL) printError(MEM_ALLOC, "supportedEdge2");
	for(i=0;i<defaultArraySize;i++)
	{
		support[i]=0;
		supportedEdge1[i]=NULL;
		supportedEdge2[i]=NULL;
	}
	time_t second_start;
	uint64_t *list_of_reads,*list_of_explored,no_of_reads=0,last;
	if((list_of_reads=(uint64_t *)malloc(2000000*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC,"list_of_reads");
	if((list_of_explored=(uint64_t *)malloc(1000000*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC,"list_of_explored");
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graph[i]!=NULL) // for each node in the graph
		{
			arraySize=1000; totalPaths=0;last=0; no_of_reads=0;
			second_start=time(NULL);
			for(v=graph[i];v!=NULL;v=v->next)
			{
				uint64_t distance=0;
				for(index=1;index<=(uint64_t)(v->listOfReads[0]);index++)
				{
					distance+=(v->listOfReads[index]&0X00000000003FF800)>>11;
					if(distance>maximumUpperBoundOfInsert)
						break;
					if(flag[v->listOfReads[index]>>24]==1)
					{
						list_of_reads[no_of_reads++]=v->listOfReads[index]>>24;
						flag[v->listOfReads[index]>>24]=0;
					}
				}
				distance=0;
				for(index=1;index<=(uint64_t)(v->twinEdge->listOfReads[0]);index++)
				{
					distance+=(v->twinEdge->listOfReads[index]&0X00000000003FF800)>>11;
					if(distance>maximumUpperBoundOfInsert)
						break;
					if(flag[v->twinEdge->listOfReads[index]>>24]==1)
					{
						list_of_reads[no_of_reads++]=v->twinEdge->listOfReads[index]>>24;
						flag[v->twinEdge->listOfReads[index]>>24]=0;
					}
				}
			}
			for(l=0;l<no_of_reads;l++)
			{
				for(a=read_to_edge[list_of_reads[l]];a!=NULL;a=a->next)
				{
					for(j=0;j<last;j++)
					{
						if(list_of_explored[j]==a->edge->from_ID)
							break;
					}
					if(j==last)
					{
						list_of_explored[last++]=a->edge->from_ID;
						exploreGraph(a->edge->from_ID,1);
					}
					for(j=0;j<last;j++)
					{
						if(list_of_explored[j]==a->edge->ID)
							break;
					}
					if(j==last)
					{
						list_of_explored[last++]=a->edge->ID;
						exploreGraph(a->edge->ID,1);
					}
				}
			}
			if(totalPaths==0 || last==0 || no_of_reads==0)
				continue;
			quicksortPaths(0,totalPaths-1);
			indexPaths(last);
			for(l=0;l<no_of_reads;l++)
			{
				mate_pair_1=list_of_reads[l];
				for(w=mate_pair[mate_pair_1];w!=NULL;w=w->next)
				{
					mate_pair_2=w->ID;
					if(w->flag==1) //On different edge and do not find path between same pair of matepair twice
					{
						if(graph[mate_pair_1]==NULL && graph[mate_pair_2]==NULL && flag[mate_pair_2]==1) // both on edge
						{
							if(totalPaths==0)
							{
								total_supported_pairs=0;
							}
							else
							{
								total_supported_pairs=findPathsBetweenMatepairs(mate_pair_1,mate_pair_2,w->type1,w->type2,w->library);
							}
							if(total_supported_pairs>0)
							{
								for(j=0;j<total_supported_pairs;j++)
								{
									hashSupportInsert(supportedEdge1,supportedEdge2,supportedEdg1[j],supportedEdg2[j]);
									/*
									for(k=0;k<tos;k++)
									{

										if((supportedEdg1[j]==supportedEdge1[k] && supportedEdg2[j]==supportedEdge2[k]) || (supportedEdg2[j]->twinEdge==supportedEdge1[k] && supportedEdg1[j]->twinEdge==supportedEdge2[k]))
										{
											support[k]++;
											break;
										}
									}
									if(k==tos)
									{
										supportedEdge1[tos]=supportedEdg1[j];
										supportedEdge2[tos]=supportedEdg2[j];
										support[tos]=1;
										tos++;

									}
									*/
								}
							}
						}
					}
				}
			}
			freeIndexedPaths();
			if(totalPaths>0)
			{
				for(j=0;j<totalPaths;j++)
				{
					free(allPairPaths[j]);
				}
				free(allPairPaths);
				free(allPairPathLengths);
			}
			if(time(NULL)-second_start>100 || totalPaths>50000000) //print only if it takes long time
				myfprintf("node: %"PRIu64" totalPaths: %"PRIu64" in %4ld sec. reads: %"PRIu64" Last: %"PRIu64"\n",i, totalPaths,time(NULL)-second_start,no_of_reads,last);
		}

	}
	free(list_of_reads); free(list_of_explored); free(flag);
	for(i=0,j=0;i<defaultArraySize;i++)
	{
		if(supportedEdge1[i]!=NULL)
		{
			supportedEdge1[j]=supportedEdge1[i];
			supportedEdge2[j]=supportedEdge2[i];
			support[j]=support[i];
			j++;
			tos++;
		}
	}
	if(tos>0) //Edge pairs found
	{
		if((dummy=(int64_t *)malloc(tos*sizeof(int64_t)))==NULL) printError(MEM_ALLOC,"dummy");
		quicksortListOfLongEdges(supportedEdge1,supportedEdge2,support, dummy,0,tos-1);
		free(dummy);
		myfprintf("\nTable: Partial list of supported pairs of edges in findAllPairPaths() iteration %d.\n-----------------------------------------------------------------------------------------------------------\n| # |      EDGE1      |  ADDRESS  | LENGTH | FLOW |      EDGE2      |  ADDRESS  | LENGTH | FLOW | SUPPORT |\n-----------------------------------------------------------------------------------------------------------\n",iteration);
		for(i=0;i<tos;i++)
		{
			uint64_t a1=(uint64_t)supportedEdge1[i],a2=(uint64_t)supportedEdge1[i]->twinEdge,a3=(uint64_t)supportedEdge2[i],a4=(uint64_t)supportedEdge2[i]->twinEdge;
					a1=a1&0XFFFF;a2=a2&0XFFFF;a3=a3&0XFFFF;a3=a3&0XFFFF; a4=a4&0XFFFF; // We only print lower 16 bits of the address
					myfprintf("|%2"PRIu64" |%8"PRIu64",%8"PRIu64"|%5"PRIu64" %5"PRIu64"|%7"PRIu64" |%5.1f |%8"PRIu64",%8"PRIu64"|%5"PRIu64" %5"PRIu64"|%7"PRIu64" |%5.1f |%8d |\n",i+1,supportedEdge1[i]->from_ID, supportedEdge1[i]->ID,a1,a2,supportedEdge1[i]->lengthOfEdge,supportedEdge1[i]->flow,supportedEdge2[i]->from_ID,supportedEdge2[i]->ID,a3,a4,supportedEdge2[i]->lengthOfEdge,supportedEdge2[i]->flow,support[i]);
			if(i==9)
				break;
		}
		myfprintf("-----------------------------------------------------------------------------------------------------------\nPairs of supported edges: %8d\n\n",tos);
		myfprintf("\nFunction findAllPairPaths() in %4ld sec.\n",time(NULL)-seconds_s);
		returnValue=resolvePairsOfEdges(supportedEdge1,supportedEdge2,tos-1);
	}
	free(supportedEdge1);free(supportedEdge2);
	return returnValue;
}

int hashSupportInsert(Edge ** supportedEdge1, Edge ** supportedEdge2, Edge *edge1,Edge *edge2)
{
	uint64_t probe=(edge1->from_ID + edge1->ID + edge2->from_ID + edge2->ID ) %defaultArraySize;
	while(supportedEdge1[probe]!=NULL)
	{
		if(supportedEdge1[probe]==edge1 && supportedEdge2[probe]==edge2)
		{
			support[probe]++;
			return 1;
		}
		if(probe<defaultArraySize-1) probe++;
		else probe=0;
	}
	probe=(edge1->from_ID + edge1->ID + edge2->from_ID + edge2->ID ) %defaultArraySize;
	while(supportedEdge1[probe]!=NULL)
	{
		if(supportedEdge1[probe]==edge2->twinEdge && supportedEdge2[probe]==edge1->twinEdge)
		{
			support[probe]++;
			return 1;
		}
		if(probe<defaultArraySize-1) probe++;
		else probe=0;
	}

	probe=(edge1->from_ID + edge1->ID + edge2->from_ID + edge2->ID ) %defaultArraySize;
	while(supportedEdge1[probe]!=NULL)
	{
		if(probe<defaultArraySize-1) probe++;
		else probe=0;
	}
	supportedEdge1[probe]=edge1;
	supportedEdge2[probe]=edge2;
	support[probe]=1;
	return 1;

}

/* ============================================================================================
   This function will save the path found in the main memory.
   ============================================================================================ */
void savePath(Edge **edges,int distance,uint64_t top)
{
	uint64_t i;
	if(totalPaths==0)
	{
		if((allPairPaths=(Edge ***)malloc(arraySize*sizeof(Edge**)))==NULL) printError(MEM_ALLOC,"allPairPaths");
		if((allPairPathLengths=(uint64_t*)malloc(arraySize*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC,"allPairPathLengths");
	}
	if(totalPaths>1000000000)
	{
		myfprintf("Highest number of paths: %"PRIu64"\n    Longest path Length: %"PRIu64"\n",highestPaths,longestPathLength);
		printError(OTHERS,"Total path found >10000000 exiting program\n");
	}
	if(totalPaths>=arraySize) //increase the size of array if needed
	{
		arraySize=2*arraySize;
		if((allPairPaths = realloc(allPairPaths, arraySize*sizeof(Edge**)))==NULL)
			 printError(MEM_ALLOC,"reallocating allPairPaths");
		if((allPairPathLengths = (uint64_t *)realloc (allPairPathLengths, arraySize*sizeof(uint64_t)))==NULL)
			printError(MEM_ALLOC,"reallocating allPairPathLengths");
	}
	if((allPairPaths[totalPaths]=(Edge**)malloc((top+2)*sizeof(Edge*)))==NULL)
	{
		myfprintf("Highest number of paths: %"PRIu64"\n    Longest path Length: %"PRIu64"\n",highestPaths,longestPathLength);
		printError(MEM_ALLOC,"allPairPaths[totalPaths]");
	}
	if(longestPathLength<top)
		longestPathLength=top;
	for(i=1;i<=top;i++)
		allPairPaths[totalPaths][i]=edges[i];
	allPairPaths[totalPaths][0]=edges[top];
	allPairPaths[totalPaths][top+1]=NULL;
	allPairPathLengths[totalPaths]=distance;
	totalPaths++;
	if(totalPaths>highestPaths)
		highestPaths=totalPaths;
}

/* ============================================================================================
   This function checks if inserting next edge in the path violates the flow of next edge or not
   Here we added allow the edge to be used flow+2 times (as the flow in not always accurate.
   ============================================================================================ */
int checkFlow(Edge **edges, uint64_t top, Edge *next_edge)
{
	uint64_t i;
	float flow=0;
	for(i=1;i<=top;i++)
		if(edges[i]==next_edge || edges[i]==next_edge->twinEdge)
			flow+=1;
	if(next_edge->listOfReads[0]==0 && flow<=next_edge->flow) // 1 extra units of flow for simple edges
		return 1;
	if(flow<next_edge->flow)
		return 1;
	return 0;
}



/* ============================================================================================
   quicksortPaths() for sorting paths algorithm
   ============================================================================================ */
void quicksortPaths(int64_t left, int64_t right)
{
	uint64_t j;
	if( left < right )
	{
		j = partitionPaths(left,right);
		if(j>0)
		quicksortPaths(left,j-1);
		quicksortPaths(j+1,right);
	}
}


/* ============================================================================================
   This function is used to partition the array in quicksortPaths() algorithm
   ============================================================================================ */

uint64_t partitionPaths(int64_t left, int64_t right)
{
	int64_t i=left,j=right+1,temp_dist;
	Edge **temp_location;
	while(1)
	{
		do
			++i;
		while(i<=right && (allPairPaths[i][1]->from_ID<allPairPaths[left][1]->from_ID || (allPairPaths[i][1]->from_ID==allPairPaths[left][1]->from_ID && allPairPaths[i][0]->from_ID<allPairPaths[left][0]->from_ID)));
		do
			--j;
		while(allPairPaths[j][1]->from_ID>allPairPaths[left][1]->from_ID || (allPairPaths[j][1]->from_ID==allPairPaths[left][1]->from_ID && allPairPaths[j][0]->from_ID>allPairPaths[left][0]->from_ID));
		if( i >= j )
			break;
		temp_location=allPairPaths[i]; allPairPaths[i]=allPairPaths[j]; allPairPaths[j]=temp_location;
		temp_dist=allPairPathLengths[i]; allPairPathLengths[i]=allPairPathLengths[j];allPairPathLengths[j]=temp_dist;
	}
	temp_location=allPairPaths[left]; allPairPaths[left]=allPairPaths[j];allPairPaths[j]=temp_location;
	temp_dist=allPairPathLengths[left]; allPairPathLengths[left]=allPairPathLengths[j]; allPairPathLengths[j]=temp_dist;
  	return j;
}


/* ============================================================================================
   Index the string node of the paths for quick access.
   ============================================================================================ */
void indexPaths(uint64_t num)
{
	uint64_t i,firstReadOnEdge=-1,second_read=-1,count=-1;
	if((indx = (struct indexes **)malloc((num+2)*sizeof(struct indexes **)))==NULL) printError(MEM_ALLOC, "indx");
	for(i=0;i<num+2;i++)
		indx[i]=NULL;
	for(i=0;i<totalPaths;i++)
	{
		if(allPairPaths[i][1]->from_ID!=firstReadOnEdge)
		{
			struct indexes *v;
			if((v = (struct indexes *)malloc(sizeof(struct indexes)))==NULL) printError(MEM_ALLOC, "v");
			v->firstNode=allPairPaths[i][1]->from_ID;
			v->lastNode=allPairPaths[i][0]->from_ID;
			v->position=i;
			v->next=NULL;
			firstReadOnEdge=allPairPaths[i][1]->from_ID;
			second_read=allPairPaths[i][0]->from_ID;
			count++;
			if(indx[count]==NULL)
			{
				indx[count]=v;
			}
			else
			{
				v->next=indx[count]; indx[count]=v;
			}
		}
		else if(allPairPaths[i][0]->from_ID!=second_read)
		{
			struct indexes *vv;
			if((vv = (struct indexes *)malloc(sizeof(struct indexes)))==NULL) printError(MEM_ALLOC, "vv");
			vv->firstNode=allPairPaths[i][1]->from_ID;
			vv->lastNode=allPairPaths[i][0]->from_ID;
			vv->next=NULL;
			vv->position=i;
			second_read=allPairPaths[i][0]->from_ID;
			if(indx[count]==NULL)
			{
				indx[count]=vv;
			}
			else
			{
				vv->next=indx[count]; indx[count]=vv;
			}
		}
	}
}


/* ============================================================================================
   This function will free memory used by path index
   ============================================================================================ */
void freeIndexedPaths(void)
{
	uint64_t i;
	struct indexes *v,*v_next;
	for(i=0;indx[i]!=NULL;i++)
	{

		for(v=indx[i];v!=NULL;v=v_next)
		{
			v_next=v->next;
			free(v);

		}
	}
	free(indx);
}

/* ============================================================================================
   This function will return that appropriate location in the array.
   ============================================================================================ */

uint64_t startingPosition(uint64_t firstReadOnEdge,uint64_t second_read)
{
	uint64_t i;
	struct indexes *v;
	for(i=0;indx[i]!=NULL;i++)
	{
		if(indx[i]->firstNode==firstReadOnEdge)
		{
			for(v=indx[i];v!=NULL;v=v->next)
			{
				if(v->lastNode==second_read)
				{
					return v->position;
				}
			}
		}
	}
	return -1;
}


/* ============================================================================================
   This function generates flow in the graph.
   ============================================================================================ */
void computeMinCostFlow(void)
{
	time_t seconds_s=time(NULL);
	uint64_t edgesUsedMoreThanOnce=0, simple_edge=0,edgesUsedExactlyOnce=0,nodes_in_graph=0;
	uint64_t numberOfEdges=0, numberOfNodes=0, cnt, sinkToSourceCap;
	uint64_t *node_index, *inverse_node_index, *from, *to, *flow, i, j, k, delta, index;
	int xi, cost_scale_factor=10,flow_lb,lb_1=0,lb_2=0,lb_3=0,ub_1=0,ub_2=0,ub_3=0;
	int lb=0,ub=100,first_UB=3,second_UB=5,third_UB=1000;
	uint64_t n_1in,n_1out,n_2in,n_2out,node_1=0,node_2=0,node_3=0,node_4=0;
	uint64_t u_1in=0,u_1out=0,u_2in=0,u_2out=0,v_1in=0,v_1out=0,v_2in=0,v_2out=0;
	float cost,cost_1,cost_2,cost_3,a_statistics[25];
	double cost_of_cs2;
	char cs2InputFileName[10000],cs2OutputFileName[10000];
	strcpy(cs2InputFileName,directoryName); strcat(cs2InputFileName,"cs2_input.tmp"); // cs2 input file for flow
	strcpy(cs2OutputFileName,directoryName); strcat(cs2OutputFileName,"cs2_output.tmp"); // cs2 output file with flow
	FILE  *fp_flow,*fp_overlap1=fopen(cs2InputFileName,"w");
	myfprintf("\nIn function computeMinCostFlow().\n");
	Edge *v;
	if(fp_overlap1==NULL) printError(OPEN_FILE,cs2InputFileName);
	//for(i=1;i<=numberOfUniqueReads;i++)
	//	if(graph[i]!=NULL)
	//		numberOfNodes++;
	edgesUsedMoreThanOnce=0; count=0;
	simple_edge=0,edgesUsedExactlyOnce=0,nodes_in_graph=0;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graph[i]!=NULL)
			nodes_in_graph++;
		for(v=graph[i];v!=NULL;v=v->next)
		{
			if(i>=v->ID)
			{
				count++;
				if(v->listOfReads[0]==0)
					simple_edge++;
				else
				{
					delta=v->lengthOfEdge-1; k=0;
					for(index=1;index<=(uint64_t)(v->listOfReads[0]);index++)
						k+=frequency[v->listOfReads[index]>>24];
					for(j=1;j<=20;j++)
						a_statistics[j]=(float)((float)delta*(float)((float)numberOfReads/(float)genomeSize)-(float)((float)k*log((float)((float)(j+1)/(float)j))));
					if(a_statistics[1]>=aStatisticsThreshold && delta>=minDelta)
						edgesUsedExactlyOnce++;
					else
						edgesUsedMoreThanOnce++;
				}
			}
		}
	}
	myfprintf("          Number of edges: %10"PRIu64"\n             Simple edges: %10"PRIu64"\n  Edges used exactly once: %10"PRIu64"\nEdges used more than once: %10"PRIu64"\n          Number of nodes: %10"PRIu64"\n",count,simple_edge, edgesUsedExactlyOnce, edgesUsedMoreThanOnce,nodes_in_graph);
	numberOfEdges=1+6*nodes_in_graph+2*(count-(edgesUsedMoreThanOnce))+6*(edgesUsedMoreThanOnce);
	numberOfNodes=nodes_in_graph*4+2;
	if((from = (uint64_t *)malloc(numberOfEdges*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC, "from");
	if((to = (uint64_t *)malloc(numberOfEdges*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC, "to");
	if((flow = (uint64_t *)malloc(numberOfEdges*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC, "flow");

	if((node_index=(uint64_t *)malloc((numberOfUniqueReads+1)*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC, "node index");
	for(i=0;i<=numberOfUniqueReads;i++)
		node_index[i]=0;
	if((inverse_node_index=(uint64_t *)malloc((numberOfNodes+4)*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC, "inverse_node index");
	for(i=0;i<numberOfNodes+4;i++)
		inverse_node_index[i]=0;

	sinkToSourceCap = numberOfNodes * 100;

	fprintf(fp_overlap1,"p min %"PRIu64" %"PRIu64"\n",numberOfNodes,numberOfEdges);
	fprintf(fp_overlap1,"n 1 0\n");
	fprintf(fp_overlap1,"n 2 0\n");
	fprintf(fp_overlap1,"a 2 1 1 %"PRIu64" 0\n", sinkToSourceCap);
	cnt=3;
	for(i=1;i<=numberOfUniqueReads;i++) // edges for nodes
	{
		if(graph[i]!=NULL)
		{
			n_1in=cnt; n_1out=cnt+1; n_2in=cnt+2; n_2out=cnt+3;

			node_index[i]=cnt;
			inverse_node_index[cnt]=i;
			inverse_node_index[cnt+1]=i;
			inverse_node_index[cnt+2]=i;
			inverse_node_index[cnt+3]=i;

			lb=0; ub=100; cost=1000000;
			fprintf(fp_overlap1,"a 1 %"PRIu64" %d %d %f\n",n_1in,lb,ub,cost*cost_scale_factor);//source to node with lower bound 0 upper bound 100 and cost 1000000
			fprintf(fp_overlap1,"a %"PRIu64" 2 %d %d %f\n",n_1out,lb,ub,cost*cost_scale_factor );//node to sink with lower bound 0 upper bound 100 and cost 1000000
			fprintf(fp_overlap1,"a 1 %"PRIu64" %d %d %f\n",n_2in,lb,ub,cost*cost_scale_factor);//source to node with lower bound 0 upper bound 100 and cost 1000000
			fprintf(fp_overlap1,"a %"PRIu64" 2 %d %d %f\n",n_2out,lb,ub,cost*cost_scale_factor);//node to sink with lower bound 0 upper bound 100 and cost 1000000
			lb=0;ub=100;cost=1;// every read can be used any number of times for real reads. should be used at least once for corrected reads.
			fprintf(fp_overlap1,"a %"PRIu64" %"PRIu64" %d %d %f\n",n_1in,n_1out,lb,ub,cost*cost_scale_factor);//node_in to node_out with lower bound 0 upper bound 100 and cost 1
			fprintf(fp_overlap1,"a %"PRIu64" %"PRIu64" %d %d %f\n",n_2in,n_2out,lb,ub,cost*cost_scale_factor);//node_in to node_out with lower bound 0 upper bound 100 and cost 1
			cnt=cnt+4;
		}
	}
	for(i=1;i<=numberOfUniqueReads;i++) // edges for edges
	{
		for(v=graph[i];v!=NULL;v=v->next)
		{
			if(i>=v->ID)
			{
				j=node_index[v->from_ID];
				u_1in=j; u_1out=j+1; u_2in=j+2;	u_2out=j+3;
				j=node_index[v->ID];
				v_1in=j; v_1out=j+1;v_2in=j+2;v_2out=j+3;
				delta=v->lengthOfEdge-1; k=0;
				if(v->listOfReads[0]==0) // these are simple edges
				{
					k=0; a_statistics[1]=0;
				}
				else	//composite edges
				{
					for(index=1;index<=(uint64_t)(v->listOfReads[0]);index++)
					{
						k+=frequency[v->listOfReads[index]>>24];
					}
					for(j=1;j<=20;j++)
						a_statistics[j]=(float)((float)delta*(float)((float)numberOfReads/(float)genomeSize)-(float)((float)k*log((float)((float)(j+1)/(float)j))));

				}
				if(v->typeOfEdge==0)
				{
					node_1=v_1out;node_2=u_1in;node_3=u_2out;node_4=v_2in;
				}
				else if(v->typeOfEdge==1)
				{
					node_1=u_2out;node_2=v_1in;node_3=v_2out;node_4=u_1in;
				}
				else if(v->typeOfEdge==2)
				{
					node_1=u_1out;node_2=v_2in;node_3=v_1out;node_4=u_2in;
				}
				else if(v->typeOfEdge==3)
				{
					node_1=u_1out;node_2=v_1in;node_3=v_2out;node_4=u_2in;
				}
				if(v->listOfReads[0]==0) // simple edge: could be used any number of times
				{
					lb=0;ub=100;cost=10000;
					fprintf(fp_overlap1,"a %"PRIu64" %"PRIu64" %d %d %f\n",node_1,node_2,lb,ub,cost*cost_scale_factor);
					fprintf(fp_overlap1,"a %"PRIu64" %"PRIu64" %d %d %f\n",node_3,node_4,lb,ub,cost*cost_scale_factor);
				}
				else if(a_statistics[1]>=aStatisticsThreshold && delta>=minDelta) //these edges should traversed exactly once;
				{
					lb=1; ub=1; cost=1;
					fprintf(fp_overlap1,"a %"PRIu64" %"PRIu64" %d %d %f\n",node_1,node_2,lb,ub,cost*cost_scale_factor);
					fprintf(fp_overlap1,"a %"PRIu64" %"PRIu64" %d %d %f\n",node_3,node_4,lb,ub,cost*cost_scale_factor);
				}
				else // these are composite edges and sould be used 1 to any-number of times times
				{
					flow_lb=0; // default flow lower bound. To be in the safe side we can set it to 0 for real reads
					if(delta>minDelta)
					{
						for(j=1;j<=19;j++)
						{
							if (a_statistics[j]<-aStatisticsThreshold && a_statistics[j+1]>aStatisticsThreshold)
							{
								flow_lb=j+1;
								break;
							}
						}
					}
					else
					{
						if(v->listOfReads[0]>35) // composite edges with many reads should have flow at least 1
							flow_lb=1;
					}
					cost_1=0;cost_2=0;cost_3=0;
					uint64_t n=numberOfReads, N=genomeSize;
					for(index=1;index<=(uint64_t)(v->listOfReads[0]);index++)
					{
						xi=frequency[v->listOfReads[index]>>24];
						cost_1=cost_1+(-(xi*log((float)(first_UB)))-(n-xi)*log(N-first_UB)-(-(xi*log((float)(1)))-(n-xi)*log((float)(N-1))));
						cost_2=cost_2+(-(xi*log((float)(second_UB)))-(n-xi)*log((float)(N-second_UB))-(-(xi*log((float)(first_UB)))-(n-xi)*log((float)(N-first_UB))));
						cost_3=cost_3+(-(xi*log((float)(third_UB)))-(n-xi)*log(N-third_UB)-(-(xi*log((float)(second_UB)))-(n-xi)*log((float)(N-second_UB))));
					}
					if(flow_lb>1)
					{
						myfprintf("Flow bound on edge (%"PRIu64",%"PRIu64") set to %d. Length: %"PRIu64"\n",v->from_ID,v->ID,flow_lb,v->lengthOfEdge);
						for(index=1;index<=(uint64_t)(v->listOfReads[0]);index++)
						{
							myfprintf("%d, ",frequency[v->listOfReads[index]>>24]);
						}
						myfprintf("\n");
					}

					lb_1=minimumValue(flow_lb,first_UB);
					lb_2=maximumValue(0,minimumValue(flow_lb,second_UB)-first_UB);
					lb_3=maximumValue(0,minimumValue(flow_lb,third_UB)-second_UB);
					ub_1=first_UB;
					ub_2=second_UB-first_UB;
					ub_3=third_UB-second_UB;

					cost_1=cost_1/(ub_1-1);
					cost_2=cost_2/ub_2;
					cost_3=cost_3/ub_3;

					fprintf(fp_overlap1,"a %"PRIu64" %"PRIu64" %d %d %f\n",node_1,node_2,lb_1,ub_1,cost_1*cost_scale_factor);
					fprintf(fp_overlap1,"a %"PRIu64" %"PRIu64" %d %d %f\n",node_3,node_4,lb_1,ub_1,cost_1*cost_scale_factor);

					fprintf(fp_overlap1,"a %"PRIu64" %"PRIu64" %d %d %f\n",node_1,node_2,lb_2,ub_2,cost_2*cost_scale_factor);
					fprintf(fp_overlap1,"a %"PRIu64" %"PRIu64" %d %d %f\n",node_3,node_4,lb_2,ub_2,cost_2*cost_scale_factor);

					fprintf(fp_overlap1,"a %"PRIu64" %"PRIu64" %d %d %f\n",node_1,node_2,lb_3,ub_3,cost_3*cost_scale_factor);
					fprintf(fp_overlap1,"a %"PRIu64" %"PRIu64" %d %d %f\n",node_3,node_4,lb_3,ub_3,cost_3*cost_scale_factor);
				}
			}
		}
	}
	fclose(fp_overlap1);
	time_t seconds_start_CS2=time(NULL);
	myfprintf("Calling CS2...\n");
	cost_of_cs2=main_cs2(cs2InputFileName,cs2OutputFileName); // cs2 is called and the output is written to file
	myfprintf("             Cost of CS2 : %.0f\nFunction CS2() in %4ld sec.\n",cost_of_cs2,time(NULL)-seconds_start_CS2);
	if((fp_flow=fopen(cs2OutputFileName,"r"))==NULL) printError(OPEN_FILE,cs2OutputFileName);
	for(i=0;i<(uint64_t)numberOfEdges;i++)
		if(fscanf(fp_flow,"%"PRIu64" %"PRIu64" %"PRIu64"\n",&from[i],&to[i],&flow[i])<0) printError(OTHERS,"Unable to scan from");
	fclose(fp_flow);

	for(i=0;i<(uint64_t)numberOfEdges;i++) // for flow in the edges.
	{
		uint64_t from_node,to_node;
		if((from[i]==1 && to[i]==2) || (from[i]==2 && to[i]==1)) // between sourse and sink
			myfprintf("         Flow from T to S: %10"PRIu64"\n", flow[i]);
		else if( from[i]>2 && to[i]>2) // other nodes
		{
			from_node=inverse_node_index[from[i]];
			to_node=inverse_node_index[to[i]];
			if(from_node!=to_node && to_node>=1 && from_node>=1)
			{
				for(v=graph[from_node];v!=NULL;v=v->next)
				{
					if(v->ID==to_node)
					{
						v->flow+=flow[i];
						break;
					}
				}
			}
		}
	}
	uint64_t flow_different=0;
	for(i=1;i<=numberOfUniqueReads;i++) //take average of flow in the twin edges
	{
		for(v=graph[i];v!=NULL;v=v->next)
		{
			float averageFlow=(v->flow+v->twinEdge->flow)/2;
			if(v->flow!=v->twinEdge->flow)
			{
				flow_different++;
			}
			v->flow=averageFlow;
			v->twinEdge->flow=averageFlow;
		}
	}
	free(from); free(to);free(flow);free(node_index);free(inverse_node_index);
	myfprintf("  Flow different in edges: %10"PRIu64"\nFunction computeMinCostFlow() in %4ld sec.\n",flow_different,time(NULL)-seconds_s);
	if(!DEBUGGING)
	{
		deleteTemporaryFiles();
	}
}



/* ============================================================================================
   This function will free all memory at the very end.
   ============================================================================================ */
void freeAllocatedMemory(void)
{
	uint64_t i;
	Edge *edge_v,*edge_v_next;
	struct readToEdgeStructure *read_to_edge_v,*read_to_edge_v_next;
	struct matePairStructure *mate_pair_v,*mate_pair_v_next;
	for(i=0;i<=numberOfUniqueReads;i++)
	{
		free((uint8_t *)readsInt[i]);
		free((uint8_t *)readsReverseInt[i]);
		if(graph[i]!=NULL)
		{
			for(edge_v=graph[i];edge_v!=NULL;edge_v=edge_v_next)
			{
				edge_v_next=edge_v->next;
				deleteEdge(edge_v); // delete edges
			}
		}
	}
	free((uint8_t **)readsInt);
	free((uint8_t *)readsReverseInt);
	free(graph); //free the whole graph
	// free all the global variable used so far
	free(supportedEdg1); free(supportedEdg2);
	free(support); free(frequency);
	for(i=0;i<=numberOfUniqueReads;i++)
	{
		if(mate_pair[i]!=NULL)
		{

			for(mate_pair_v=mate_pair[i];mate_pair_v!=NULL;mate_pair_v=mate_pair_v_next)
			{
				mate_pair_v_next=mate_pair_v->next;
				free(mate_pair_v);
			}
		}
	}
	free(mate_pair); // free mate_pair structure
	for(i=0;i<=numberOfUniqueReads;i++) // delete previous mapping of read to edge
	{
		if(read_to_edge[i]!=NULL)
		{
			for(read_to_edge_v=read_to_edge[i];read_to_edge_v!=NULL;read_to_edge_v=read_to_edge_v_next)
			{
				read_to_edge_v_next=read_to_edge_v->next;
				free(read_to_edge_v->locationForward);
				free(read_to_edge_v->locationReverse);
				free(read_to_edge_v);
			}
		}
	}
	free(read_to_edge); //free mapping of reads
	free(Mean); free(standardDeviation); free(upperBoundOfInsert); free(lowerBoundOfInsert);
}

/* ============================================================================================
   This function will return type of edge if you merge edge1 and edge2
   ============================================================================================ */
int combinedEdgeType(Edge *edge1, Edge* edge2)
{
	     if((edge1->typeOfEdge==0 && edge2->typeOfEdge==0) || (edge1->typeOfEdge==1 && edge2->typeOfEdge==2)) return 0;
	else if((edge1->typeOfEdge==0 && edge2->typeOfEdge==1) || (edge1->typeOfEdge==1 && edge2->typeOfEdge==3)) return 1;
	else if((edge1->typeOfEdge==2 && edge2->typeOfEdge==0) || (edge1->typeOfEdge==3 && edge2->typeOfEdge==2)) return 2;
	else if((edge1->typeOfEdge==2 && edge2->typeOfEdge==1) || (edge1->typeOfEdge==3 && edge2->typeOfEdge==3)) return 3;
	else return -1;
}

/* ============================================================================================
   This function checks if any reads are missing in the graph or not. For debugging only.
   ============================================================================================ */

int checkIfAllReadsPresent(void)
{
	myfprintf("\nIn function checkIfAllReadsPresent().\n");
	time_t seconds_s=time(NULL);
	int *flag;
	uint64_t i,counter=0,index, readsPresent=0;
	Edge *v;
	if((flag=(int*)malloc((numberOfUniqueReads+1)*sizeof(int)))==NULL) printError(MEM_ALLOC, "Flag");
	for(i=1;i<=numberOfUniqueReads;i++)
		flag[i]=0;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graph[i]!=NULL)
		{
			flag[i]=1;
			for(v=graph[i];v!=NULL;v=v->next)
			{
				for(index=1;index<=(uint64_t)(v->listOfReads[0]);index++)
				{
					flag[v->listOfReads[index]>>24]=1;
				}
				if(DEBUGGING)
				{
					if(v->twinEdge->twinEdge!=v)
						myfprintf("Check twin edge (%"PRIu64",%"PRIu64")\n",v->from_ID,v->ID);
				}
			}
		}
	}
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(flag[i])
		{
			readsPresent+=frequency[i];//fprintf(fpo_readsRemaining,">read_%d\n%s\n",i,reads[i]);

		}
		else
		{
			counter+=frequency[i];
		}

	}
	myfprintf("Number of missing reads: %10"PRIu64"\n",counter);
	if(DEBUGGING)
	{
		for(i=1;i<=numberOfUniqueReads;i++)
		{
			if(graph[i]!=NULL)
			{
				for(v=graph[i];v!=NULL;v=v->next)
				{
					if(v->flow!=v->twinEdge->flow && v->isReducible)
						myfprintf("flow different in edge (%"PRIu64",%"PRIu64") edge flow: %.2f and twinEdge_flow: %.2f\n",v->from_ID, v->ID, v->flow, v->twinEdge->flow);
				}
			}
		}
	}
	free(flag);
	myfprintf("Function checkIfAllReadsPresent() in %4ld sec.\n",time(NULL)-seconds_s);
	return counter;
}


/* ============================================================================================
   Merge edges (not necessarily connected) in the graph that are supported.
   ============================================================================================ */
void mergeContigs(void)
{
	int i;
	time_t seconds_s=time(NULL);
	myfprintf("\nIn function mergeContigs().\n");
	for(i=0;i<20;i++)
	{
		int merged=mergeFinal();
		if(merged==0)
		{
			myfprintf("Function mergeContigs() in %4ld sec.\n",time(NULL)-seconds_s);
			return;
		}
	}
	myfprintf("Function mergeContigs() in %4ld sec.\n",time(NULL)-seconds_s);
	return;
}


/* ============================================================================================
   This function merges long contigs that are not connected by overlaps.
   TODO: Clean this function.
   ============================================================================================ */

uint64_t mergeFinal(void)
{
	static int counter=1;
	myfprintf("\nIn function mergeFinal().\nCounter: %10d\n",counter++);
	time_t seconds_s=time(NULL);
	uint64_t i,j,k,m,n,tos=0,number_merged=0,mate_pair_1,mate_pair_2,support_threshold=2,dist,total_edge_pair_to_merge=0,numberOfFeasibleEdges,arraySize=10000, flag_distance,index;
	int64_t *gapDistance, gapSum,support,distance_all,*distance_1,*distance_2;
	int *pair_support;
	Edge **listOfLongEdges,*edge1=NULL,*edge2=NULL,**edge_1,**edge_2,**feasibleListOfEdges;
	struct matePairStructure *w;
	updateGraph();
	listOfLongEdges=getListOfCompositeEdges(&tos);
	myfprintf("Total edges to consider: %10"PRIu64"\n",tos);
	if((edge_1=(Edge **)malloc(arraySize*sizeof(Edge *)))==NULL) printError(MEM_ALLOC, "edge_1");
	if((edge_2=(Edge **)malloc(arraySize*sizeof(Edge *)))==NULL) printError(MEM_ALLOC, "edge_2");
	if((pair_support=(int *)malloc(arraySize*sizeof(int)))==NULL) printError(MEM_ALLOC, "pair_support");
	if((gapDistance=(int64_t *)malloc(arraySize*sizeof(int64_t)))==NULL) printError(MEM_ALLOC, "gapDistance");
	for(i=0;i<tos;i++)
	{
		feasibleListOfEdges=getListOfFeasibleEdges(listOfLongEdges[i],&numberOfFeasibleEdges);
		for(j=0;j<numberOfFeasibleEdges;j++)
		{
			int64_t len1=listOfLongEdges[i]->lengthOfEdge,len2=feasibleListOfEdges[j]->lengthOfEdge;
			//if(len1>minimumUpperBoundOfInsert || len2>minimumUpperBoundOfInsert) // at least one edge should be longer than the bound
			//if((listOfLongEdges[i]->flow>0 || len1>minimumUpperBoundOfInsert) && (feasibleListOfEdges[j]->flow>0 || len2>minimumUpperBoundOfInsert))
			if((listOfLongEdges[i]->flow>0 || len1>contigLengthThreshold) && (feasibleListOfEdges[j]->flow>0 || len2>contigLengthThreshold))
			{
				for(k=0;k<4;k++)
				{
					if(k==0){ edge1=listOfLongEdges[i]; edge2=feasibleListOfEdges[j];}
					else if(k==1){edge1=listOfLongEdges[i]; edge2=feasibleListOfEdges[j]->twinEdge;}
					else if(k==2){edge1=listOfLongEdges[i]->twinEdge; edge2=feasibleListOfEdges[j];}
					else if(k==3){edge1=listOfLongEdges[i]->twinEdge; edge2=feasibleListOfEdges[j]->twinEdge;}
					dist=0;support=0;gapSum=0;
					if(edge1>edge2) continue; // do not double count the pair of edges
					for(index=1;index<=(uint64_t)(edge1->listOfReads[0]);index++)
					{
						dist+=(edge1->listOfReads[index]&0X00000000003FF800)>>11;
						if(dist>maximumUpperBoundOfInsert) break;
						mate_pair_1=edge1->listOfReads[index]>>24;
						if(isUniqueInEdge(edge1,mate_pair_1)) // mate pair1 is only present in this edge
						{
							for(w=mate_pair[mate_pair_1];w!=NULL;w=w->next)
							{
								mate_pair_2=w->ID;
								int lib=w->library;
								//if(len1>upperBoundOfInsert[lib] || len2>upperBoundOfInsert[lib])
								{
									if(isUniqueInEdge(edge2,mate_pair_2)) // matepair 2 is only present in this edge
									{
										distance_1=findDistanceOnEdge(edge1,mate_pair_1);
										distance_2=findDistanceOnEdge(edge2,mate_pair_2);
										if(distance_1[0]>0 && distance_2[0]>0)
										{
											flag_distance=0;
											for(m=1;m<=(uint64_t)distance_1[0];m++)
											{
												for(n=1;n<=(uint64_t)distance_2[0];n++)
												{
													distance_all=10000000;
													if(w->type1*(int)distance_1[m]<0 && w->type2*(int)distance_2[n]<0)
														distance_all=llabs(distance_1[m])+llabs(distance_2[n]);
													if(distance_all+(uint64_t)readLength<(uint64_t)(upperBoundOfInsert[lib]+minOverlap))
													{
														flag_distance=1; break;
													}
												}
												if(flag_distance==1)
													break;
											}
											if(flag_distance==1)
											{
												support++;
												gapSum+=Mean[lib]-(distance_all+readLength);
											}
										}
										free(distance_1);
										free(distance_2);
									}
								}
							}
						}
					}
					if(support>0)
					{
						edge_1[total_edge_pair_to_merge]=edge1->twinEdge;
						edge_2[total_edge_pair_to_merge]=edge2;
						pair_support[total_edge_pair_to_merge]=support;
						gapDistance[total_edge_pair_to_merge++]=gapSum/support;
						if(total_edge_pair_to_merge>=arraySize-5)
						{
							arraySize=arraySize<<1;
							if((edge_1=(Edge **)realloc(edge_1,arraySize*sizeof(Edge *)))==NULL) printError(MEM_ALLOC, "realloc edge_1");
							if((edge_2=(Edge **)realloc(edge_2,arraySize*sizeof(Edge *)))==NULL) printError(MEM_ALLOC, "realloc edge_2");
							if((pair_support=(int *)realloc(pair_support,arraySize*sizeof(int)))==NULL) printError(MEM_ALLOC, "realloc pair_support");
							if((gapDistance=(int64_t *)realloc(gapDistance,arraySize*sizeof(int64_t)))==NULL) printError(MEM_ALLOC, "realloc gapDistance");
						}
					}
				}
			}
		}
		free(feasibleListOfEdges);
	}
	if(total_edge_pair_to_merge>0)
	{
		quicksortListOfLongEdges(edge_1,edge_2,pair_support,gapDistance, 0,total_edge_pair_to_merge-1); // Sort edges according to the support.
		myfprintf("\nTable: Partial list of supported pairs of edges in mergeFinal().\n------------------------------------------------------------------------------------------\n| # |      EDGE1      | LENGTH | FLOW |      EDGE2      | LENGTH | FLOW | SUPPORT | DIST |\n------------------------------------------------------------------------------------------\n");
		for(i=0;i<total_edge_pair_to_merge;i++) //print according to support
		{
			if(i==10) break;
			myfprintf("|%2"PRIu64" |%8"PRIu64",%8"PRIu64"|%7"PRIu64" |%5.1f |%8"PRIu64",%8"PRIu64"|%7"PRIu64" |%5.1f |%8d | %4"PRId64" |\n",i+1, edge_1[i]->from_ID,edge_1[i]->ID,edge_1[i]->lengthOfEdge,edge_1[i]->flow,edge_2[i]->from_ID,edge_2[i]->ID,edge_2[i]->lengthOfEdge,edge_2[i]->flow,pair_support[i],gapDistance[i]);
		}
		myfprintf("------------------------------------------------------------------------------------------\nTotal pairs of edges supported: %10d\n\n",total_edge_pair_to_merge);
		number_merged=mergeAccordingToSupport(edge_1,edge_2,pair_support,gapDistance,total_edge_pair_to_merge,support_threshold);
	}
	free(edge_1); free(edge_2); free(pair_support); free(listOfLongEdges); free(gapDistance);
	if(DEBUGGING)
	{
		char graph_file_name[10000],temp_number[10];
		itoa(counter-1,temp_number);
		strcpy(graph_file_name,"merge-final-");
		strcat(graph_file_name,temp_number);
		printGraph(graph_file_name);
	}
	if(TEST)
	{
		printGraph("DUMMY");
	}
	myfprintf("Function mergeFinal() in %4ld sec.\n",time(NULL)-seconds_s);
	return number_merged;
}


/* ============================================================================================
   Merge edge pairs according to support.
   ============================================================================================ */
uint64_t mergeAccordingToSupport(Edge **edge_1,Edge **edge_2,int *pair_support,int64_t *gapDistance, int total_edge_pair_to_merge,int support_threshold)
{
	uint64_t number_merged=0;
	int i,j;
	float flow1,flow2;
	Edge *edge1_p,*edge2_p,*edge1_t_p,*edge2_t_p;
	myfprintf("\nTable: Partial list of merged pairs of edges in mergeFinal().\n------------------------------------------------------------------------------------------\n| # |      EDGE1      | LENGTH | FLOW |      EDGE2      | LENGTH | FLOW | SUPPORT | DIST |\n------------------------------------------------------------------------------------------\n");
	for(i=0;i<total_edge_pair_to_merge;i++) // merge according to support
	{
		if(pair_support[i]<support_threshold)
			break;
		if(edge_1[i]!=NULL && edge_2[i]!=NULL)
		{
			number_merged++;
			if(number_merged<=10)
					myfprintf("|%2"PRIu64" |%8"PRIu64",%8"PRIu64"|%7"PRIu64" |%5.1f |%8"PRIu64",%8"PRIu64"|%7"PRIu64" |%5.1f |%8d | %4"PRId64" |\n",number_merged, edge_1[i]->from_ID,edge_1[i]->ID,edge_1[i]->lengthOfEdge,edge_1[i]->flow, edge_2[i]->from_ID,edge_2[i]->ID,edge_2[i]->lengthOfEdge,edge_2[i]->flow,pair_support[i],gapDistance[i]);
			edge1_p=edge_1[i]; edge2_p=edge_2[i]; edge1_t_p=edge_1[i]->twinEdge; edge2_t_p=edge_2[i]->twinEdge;
			flow1=edge_1[i]->flow; flow2=edge_2[i]->flow;
			if(mergeDisconnectedEdges(edge_1[i],edge_2[i],gapDistance[i]))
			{

				for(j=0;j<total_edge_pair_to_merge;j++)
				{
					if(i!=j && edge_1[j]!=NULL && edge_2[j]!=NULL)
					{
						if(flow1<=1 && (edge_1[j]==edge1_p || edge_1[j]==edge1_t_p || edge_2[j]==edge1_p || edge_2[j]==edge1_t_p)) // edge 1 is deleted
						{
							edge_1[j]=NULL; edge_2[j]=NULL;
						}
						if(flow2<=1 && (edge_1[j]==edge2_p || edge_1[j]==edge2_t_p || edge_2[j]==edge2_p || edge_2[j]==edge2_t_p)) // edge 2 is deleted
						{
							edge_1[j]=NULL; edge_2[j]=NULL;
						}
					}
				}
				edge_1[i]=NULL;
				edge_2[i]=NULL;
			}
		}
	}
	myfprintf("------------------------------------------------------------------------------------------\nTotal pairs of edges merged: %10d\n\n",number_merged);
	return number_merged;
}

/* ============================================================================================
   This function returns the list of all composite edges in the graph.
   ============================================================================================ */
Edge** getListOfCompositeEdges(uint64_t *number)
{
	Edge **listOfLongEdges,*v,*u;
	uint64_t arraySize=100,tos=0,i;
	if((listOfLongEdges=(Edge **)malloc(arraySize*sizeof(Edge *)))==NULL) printError(MEM_ALLOC, "listOfLongEdges");
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graph[i]!=NULL)
		{
			for(v=graph[i];v!=NULL;v=v->next)
			{
				if( i<=v->ID && v->listOfReads[0]!=0) //Take all composite edges.
				{
					if(i<v->ID)
						listOfLongEdges[tos++]=v;
					else
					{
						for(u=graph[i];u!=v;u=u->next)
						{
							if(u==v->twinEdge)// Twin edge already considered
								break;
						}
						if(u==v) // Not considered already
							listOfLongEdges[tos++]=v;
					}
					if(tos>=arraySize-5)
					{
						arraySize=arraySize<<1;
						if((listOfLongEdges=(Edge **)realloc(listOfLongEdges,arraySize*sizeof(Edge *)))==NULL) printError(MEM_ALLOC, "reallocating listOfLongEdges");
					}
				}
			}
		}
	}
	*number=tos;
	return listOfLongEdges;
}


/* ============================================================================================
   This function returns the list of edges that can possibly be connected to 'edge' by matepairs.
   ============================================================================================ */

Edge** getListOfFeasibleEdges(Edge *edge,uint64_t *number)
{
	Edge** feasibleListOfEdges;
	struct matePairStructure *w;
	uint64_t k,mate_pair_1,mate_pair_2,arraySize=10,numberOfFeasibleEdges=0,dist=0,index;
	if((feasibleListOfEdges=(Edge **)malloc(arraySize*sizeof(Edge *)))==NULL) printError(MEM_ALLOC, "feasibleListOfEdges");
	for(index=1;index<=(uint64_t)(edge->listOfReads[0]);index++)
	{
		dist+=(edge->listOfReads[index]&0X00000000003FF800)>>11;
		if(dist>maximumUpperBoundOfInsert)
			 break;
		mate_pair_1=edge->listOfReads[index]>>24;
		if(isUniqueInEdge(edge,mate_pair_1))
		{
			for(w=mate_pair[mate_pair_1];w!=NULL;w=w->next)
			{
				mate_pair_2=w->ID;
				if(read_to_edge[mate_pair_2]==NULL) continue;
				if(read_to_edge[mate_pair_2]->edge==edge) continue;
				if(read_to_edge[mate_pair_2]->edge==edge->twinEdge) continue;
				if(isUniqueInEdge(read_to_edge[mate_pair_2]->edge,mate_pair_2))
				{
					for(k=0;k<numberOfFeasibleEdges;k++)
						if(feasibleListOfEdges[k]==read_to_edge[mate_pair_2]->edge || feasibleListOfEdges[k]==read_to_edge[mate_pair_2]->edge->twinEdge)
							break;
					if(k==numberOfFeasibleEdges)
					{
						feasibleListOfEdges[numberOfFeasibleEdges++]=read_to_edge[mate_pair_2]->edge;
						if(numberOfFeasibleEdges>=arraySize-5)
						{
							arraySize=arraySize<<1;
							if((feasibleListOfEdges=(Edge **)realloc(feasibleListOfEdges,arraySize*sizeof(Edge *)))==NULL) printError(MEM_ALLOC, "realloc feasibleListOfEdges");
						}
					}
				}
			}
		}
	}
	dist=0;
	for(index=1;index<=(uint64_t)(edge->twinEdge->listOfReads[0]);index++)
	{
		dist+=(edge->twinEdge->listOfReads[index]&0X00000000003FF800)>>11;
		if(dist>maximumUpperBoundOfInsert)
			 break;
		mate_pair_1=edge->twinEdge->listOfReads[index]>>24;
		if(isUniqueInEdge(edge->twinEdge,mate_pair_1))
		{
			for(w=mate_pair[mate_pair_1];w!=NULL;w=w->next)
			{
				mate_pair_2=w->ID;
				if(read_to_edge[mate_pair_2]==NULL) continue;
				if(read_to_edge[mate_pair_2]->edge==edge) continue;
				if(read_to_edge[mate_pair_2]->edge==edge->twinEdge) continue;
				if(isUniqueInEdge(read_to_edge[mate_pair_2]->edge,mate_pair_2))
				{
					for(k=0;k<numberOfFeasibleEdges;k++)
					{
						if(feasibleListOfEdges[k]==read_to_edge[mate_pair_2]->edge || feasibleListOfEdges[k]==read_to_edge[mate_pair_2]->edge->twinEdge)
						{
							break;
						}
					}
					if(k==numberOfFeasibleEdges)
					{
						feasibleListOfEdges[numberOfFeasibleEdges++]=read_to_edge[mate_pair_2]->edge;
						if(numberOfFeasibleEdges>=arraySize-5)
						{
							arraySize=arraySize<<1;
							if((feasibleListOfEdges=(Edge **)realloc(feasibleListOfEdges,arraySize*sizeof(Edge *)))==NULL) printError(MEM_ALLOC, "realloc feasibleListOfEdges");
						}
					}
				}
			}
		}
	}
	*number=numberOfFeasibleEdges;
	return feasibleListOfEdges;
}

/* ============================================================================================
   Sort list of long edges according to support. If the supports are same then sort by length.
   ============================================================================================ */
void quicksortListOfLongEdges(Edge **edge_1,Edge **edge_2,int * pair_support, int64_t *gapDistance, int64_t left, int64_t right)
{
	int64_t i=left,j=right,pivot=*(pair_support+((left+right)>>1)),temporary,lengthSum=edge_1[(left+right)>>1]->lengthOfEdge+edge_2[(left+right)>>1]->lengthOfEdge;
	Edge *temp;
	while(i<j)
	{
		while((int64_t)*(pair_support+i)>pivot || ((int64_t)*(pair_support+i)==pivot && ((int64_t)(edge_1[i]->lengthOfEdge+edge_2[i]->lengthOfEdge))>lengthSum))
			i++;
		while((int64_t)*(pair_support+j)<pivot || ((int64_t)*(pair_support+j)==pivot && ((int64_t)(edge_1[j]->lengthOfEdge+edge_2[j]->lengthOfEdge))<lengthSum))
			j--;
		if (i<=j)
		{
			temp=*(edge_1+i); *(edge_1+i)=*(edge_1+j); *(edge_1+j)=temp;
			temp=*(edge_2+i); *(edge_2+i)=*(edge_2+j); *(edge_2+j)=temp;
			temporary=*(pair_support+i); *(pair_support+i)=*(pair_support+j); *(pair_support+j)=temporary;
			temporary=*(gapDistance+i); *(gapDistance+i)=*(gapDistance+j); *(gapDistance+j)=temporary;
			i++; j--;
		}
	}
	if (left < j )
		quicksortListOfLongEdges(edge_1,edge_2,pair_support,gapDistance,left,j);
    	if (i < right)
		quicksortListOfLongEdges(edge_1,edge_2,pair_support,gapDistance,i,right);
}


/* ============================================================================================
   This function will remove all temporary files.
   ============================================================================================ */
int deleteTemporaryFiles(void)
{
	myfprintf("\nRemoving all temp files.\n");
	int systemRet;
	char clean_directoryName[10000];
	strcpy(clean_directoryName,"rm -f "); strcat(clean_directoryName,directoryName); strcat(clean_directoryName,"*.tmp");
	systemRet=system(clean_directoryName);
	myfprintf("All temp files removed.\n");
	return systemRet;
}


/* ============================================================================================
   This function will mark edges in the graph that are possibly for errors related to all single
   A's or C's or G's or T's in the read.
   ============================================================================================ */
void getDegreeStatistics(void)
{
	uint64_t i,highestDegreeNode=-1, totalDegree=0,numberOfNodes=0,maxDegree=0,averageDegree,nodesWithHighDegree=0;
	time_t seconds_s=time(NULL);
	myfprintf("\nIn function getDegreeStatistics().\n");
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graphEconomy[i]!=NULL)
		{
			numberOfNodes++;
			totalDegree+=graphEconomy[i][0];
			if((uint64_t)(graphEconomy[i][0])>maxDegree)
			{
				maxDegree=graphEconomy[i][0];
				highestDegreeNode=i;
			}
		}
	}
	averageDegree=totalDegree/numberOfNodes;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graphEconomy[i]!=NULL)
		{
			if((uint64_t)(graphEconomy[i][0])>3000*averageDegree) // If a node has more than 20 times the average number of edges
			{
				selfOverlappingReads[i]=2;
				nodesWithHighDegree++;
			}
		}
	}
	myfprintf("    Total degrees: %10"PRIu64" (Double counted each edge)\n      Total nodes: %10"PRIu64"\n   Average degree: %10"PRIu64"\nHigh Degree nodes: %10"PRIu64" (Degree more than 5*averageDegree)\n",totalDegree,numberOfNodes,averageDegree,nodesWithHighDegree);
	myfprintf("   Highest degree: %10"PRIu64" (%10"PRIu64")\nFunction getDegreeStatistics() in %4ld sec.\n",highestDegreeNode,maxDegree,time(NULL)-seconds_s);
}


/* ============================================================================================
   Sort the edges. First by length, then id then type.
   ============================================================================================ */

void sortEdgesEconomy(int type)
{
	uint64_t i;
	myfprintf("\nIn function sortEdgesEconomy()\n");
	time_t second_s=time(NULL);
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graphEconomy[i]!=NULL)
		{
			if(graphEconomy[i][0]>1) //no need to sort node with one neighbor.
			{
				quicksortEdgesEconomy(graphEconomy[i],1,graphEconomy[i][0],type);
			}
		}
	}
	myfprintf("Function sortEdgesEconomy() in %4ld sec.\n",time(NULL)-second_s);
}


/* ============================================================================================
   Sort the edges. First by length, then id then type.
   ============================================================================================ */
void quicksortEdgesEconomy(uint64_t *listForSorting,int64_t left, int64_t right,int type) // quiksort to sort the list of neighbouts of a read
{
	int64_t i=left,j=right;
	uint64_t pivot=listForSorting[(left+right)>>1], temp;
	while(i<j)
	{
		while(compareEdges(pivot,listForSorting[i],type)<0)
			i++;
		while(compareEdges(pivot,listForSorting[j],type)>0)
			j--;
		if (i<=j)
		{
			temp=*(listForSorting+i);
			*(listForSorting+i)=*(listForSorting+j);
			*(listForSorting+j)=temp;
			i++; j--;
		}
	}
	if (left < j)
		quicksortEdgesEconomy(listForSorting,left,j,type);
    	if (i < right)
		quicksortEdgesEconomy(listForSorting,i,right,type);
}

/* ============================================================================================
   Sort the edges. First by length, then id then type.
   ============================================================================================ */
int compareEdges(uint64_t a,uint64_t b, int type)
{
	uint64_t tempa=0,tempb=0;
	if(type==1) // first length,then id, then type
	{
		tempa=((a<<40)&0XC000000000000000) | ((a>>2)&0X3FFFFFFFFFC00000) | (a&0X00000000000FFFFF);
		tempb=((b<<40)&0XC000000000000000) | ((b>>2)&0X3FFFFFFFFFC00000) | (b&0X00000000000FFFFF);
	}
	if(type==2) //first id, then type then length
	{
		tempa=a;
		tempb=b;
	}
	if(tempa<tempb)
		return -1;
	if(tempa>tempb)
		return 1;
	return 0;
}


void sortEdgeListOfOverlapGraph(void)
{
	uint64_t i;
	myfprintf("\nIn function sortEdgesOverlapGraph()\n");
	time_t second_s=time(NULL);
	Edge *v;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(graph[i]!=NULL)
		{
			int64_t counter=0;
			for(v=graph[i];v!=NULL;v=v->next )
			{
				counter++;
			}
			if(counter>1) //no need to sort node with one neighbor.
			{
				sortListOfNeighborsOverlapGraph(graph[i],counter,i);
			}
		}
	}
	myfprintf("Function sortEdgesOverlapGraph() in %4ld sec.\n",time(NULL)-second_s);
}


/* ============================================================================================
   Quicksort neighbours of the node from_ID.
   ============================================================================================ */
void sortListOfNeighborsOverlapGraph(Edge *edge,uint64_t number,uint64_t from_ID) // quicksort.
{
	Edge *v,**listForSorting;
	uint64_t i=0;
	if((listForSorting =(Edge **)malloc(number*sizeof(Edge *)))==NULL) printError(MEM_ALLOC, "listForSorting");
	for(v=edge;v!=NULL;v=v->next)
		listForSorting[i++]=v;
	quicksortEdges(listForSorting,0,number-1);
	for(i=0;i<number;i++) // update the sorted list
	{
		if(i==0)
		{
			listForSorting[i]->next=NULL;
			listForSorting[i]->previous=NULL;
			graph[from_ID]=listForSorting[i];
			v=graph[from_ID];
		}
		else
		{
			v->next=listForSorting[i];
			v->next->previous=v;
			v=v->next;
			v->next=NULL;
		}
	}
	free(listForSorting);
}

/* ============================================================================================
   Quicksort neighbours of the node from_ID.
   ============================================================================================ */
void quicksortEdges(Edge ** listForSorting,int64_t left, int64_t right) // quiksort to sort the list of neighbouts of a read
{
	int64_t i=left,j=right;
	int64_t pivot=listForSorting[(left+right)>>1]->lengthOfEdge;
	Edge *temp;
	while(i<j)
	{
		while((int64_t)(listForSorting[i]->lengthOfEdge)<pivot)
			i++;
		while((int64_t)(listForSorting[j]->lengthOfEdge)>pivot)
			j--;
		if (i<=j)
		{
			temp=*(listForSorting+i);
			*(listForSorting+i)=*(listForSorting+j);
			*(listForSorting+j)=temp;
			i++; j--;
		}
	}
	if (left < j)
		quicksortEdges(listForSorting,left,j);
    	if (i < right)
		quicksortEdges(listForSorting,i,right);
}





/* ============================================================================================
   This function inserts 64 bit integers in the hash table.
   MSB 40 bits to represent read number.
   LSB 2 bits for orientation.
   other bits are unused.
   ============================================================================================ */
uint64_t hashTableInsert(uint64_t *value,uint64_t readNumber,uint8_t type)
{
	static uint64_t hashMiss=0;
	uint64_t probe = getHashValue(value),*returnValue=NULL;
	uint64_t read, type64=type, type2;
	while(listOfHashReadsInt[probe]!=NULL) // Wile value not found or an empty space not found.
	{
		read=listOfHashReadsInt[probe][1]>>24;
		type2=listOfHashReadsInt[probe][1]&0X0000000000000003;
		if(type2==0)
			returnValue=get64Bit2Int(readsInt[read],0,hashStringLength);
		else if(type2==1)
			returnValue=get64Bit2Int(readsInt[read],readLength-hashStringLength,hashStringLength);
		else if(type2==2)
			returnValue=get64Bit2Int(readsReverseInt[read],0,hashStringLength);
		else if(type2==3)
			returnValue=get64Bit2Int(readsReverseInt[read],readLength-hashStringLength,hashStringLength);
		if(value[1]==returnValue[1] && value[0]==returnValue[0]) // check less significant variable first
		{
			free((uint64_t *) returnValue);
			break;
		}
		free((uint64_t *) returnValue);
		if(probe<(sizeOfHashTable-1)) // Increase linear probe.
			probe++;
		else	// Reached at the end of hash table.
			probe=0;
		hashMiss++; // Number of hash miss while inserting elements in the hash table
	}
	if(listOfHashReadsInt[probe]==NULL) // Inserting read for the first time
	{
		if((listOfHashReadsInt[probe]=(uint64_t *)malloc((2)*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC, "listOfHashReadsInt[probe]");
		listOfHashReadsInt[probe][0]=1;
		listOfHashReadsInt[probe][1]=(readNumber<<24)|type64;
	}
	else // Increasing the array of reads.
	{
		if((listOfHashReadsInt[probe]=(uint64_t *)realloc(listOfHashReadsInt[probe],(listOfHashReadsInt[probe][0]+2)*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC,"reallocating listOfHashReadsInt[probe] failed");
		listOfHashReadsInt[probe][0]++;
		listOfHashReadsInt[probe][listOfHashReadsInt[probe][0]]=(readNumber<<24)|type64;
	}
	return hashMiss; // Return the number of hash miss while searching in total.
}

uint64_t getHashValue(uint64_t *value)
{
	uint64_t returnValue=((value[1]%sizeOfHashTable)+(value[0]%sizeOfHashTable)*precomputeHash)%sizeOfHashTable;
	return returnValue;
}

/* ============================================================================================
   This function searches 64 bit integers in the hash table.
   MSB 40 bits to represent read number.
   LSB 2 bits for orientation.
   other bits are unused.
   ============================================================================================ */
uint64_t hashTableSearch(uint64_t *value)
{
	uint64_t probe = getHashValue(value),*returnValue=NULL;
	uint64_t i, read, type2;
	for(i=0;i<(uint64_t)sizeOfHashTable;i++)
	{
		if(listOfHashReadsInt[probe]==NULL) // Empty space found. Means not in the hash table.
			return -1; // Not found because the space is empty.
		else
		{
			read=listOfHashReadsInt[probe][1]>>24;
			type2=listOfHashReadsInt[probe][1]&0X0000000000000003;
			if(type2==0)
				returnValue=get64Bit2Int(readsInt[read],0,hashStringLength);
			else if(type2==1)
				returnValue=get64Bit2Int(readsInt[read],readLength-hashStringLength,hashStringLength);
			else if(type2==2)
				returnValue=get64Bit2Int(readsReverseInt[read],0,hashStringLength);
			else if(type2==3)
				returnValue=get64Bit2Int(readsReverseInt[read],readLength-hashStringLength,hashStringLength);

			if(value[1]==returnValue[1] && value[0]==returnValue[0]) // check less significant variable first
			{
				free((uint64_t *) returnValue);
				return probe; // Return the index in the hash table.
			}
			free((uint64_t *) returnValue);

		}
		if(probe<(sizeOfHashTable-1)) // Increase linear probe.
			probe++;
		else //reached at the end of the hash table.
			probe=0;
		numberOfHashMiss++; // Count the number of hash miss while searching.
	}
	return -1; // Not found searched the entire table
}

/* ============================================================================================
   This function creates the hash table for each 28 bp prefix and suffix of the unique reads and
   their reverse complements.
   ============================================================================================ */
void hashPrefixesAndSuffix(void)
{
	time_t second_s=time(NULL);
	myfprintf("\nIn function hashPrefixesAndSuffix().\n");
	uint64_t i;
	uint64_t hashMiss=0;
	if(minOverlap>64)
		hashStringLength=64; // string to hash in hash table
	else
		hashStringLength=minOverlap;
	myfprintf("    Hash string length: %10d\n",hashStringLength);
	sizeOfHashTable=findNextPrimeFromFile(8*numberOfUniqueReads); // twice the number of elements to put in the table
	precomputeHash=((0XFFFFFFFFFFFFFFFF)%sizeOfHashTable+1)%sizeOfHashTable;
	myfprintf("Hash table size set to: %10"PRIu64"\n",sizeOfHashTable);
	if((listOfHashReadsInt=(uint64_t  **)malloc((sizeOfHashTable)*sizeof(uint64_t *)))==NULL) printError(MEM_ALLOC, "listOfHashReadsInt");
	for(i=0;i<(uint64_t)sizeOfHashTable;i++)
		listOfHashReadsInt[i]=NULL;
	uint64_t *prefixRead,*suffixRead,*prefixReadReverseComplement,*suffixReadReverseComplement;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		prefixRead=get64Bit2Int(readsInt[i],0,hashStringLength);
		suffixRead=get64Bit2Int(readsInt[i],readLength-hashStringLength,hashStringLength);
		prefixReadReverseComplement=get64Bit2Int(readsReverseInt[i],0,hashStringLength);
		suffixReadReverseComplement=get64Bit2Int(readsReverseInt[i],readLength-hashStringLength,hashStringLength);
		hashTableInsert(prefixRead,i,0); hashTableInsert(suffixRead,i,1);
		hashTableInsert(prefixReadReverseComplement,i,2); hashMiss=hashTableInsert(suffixReadReverseComplement,i,3);
		free((uint64_t *) prefixRead);free((uint64_t *) suffixRead);
		free((uint64_t *) prefixReadReverseComplement);free((uint64_t *) suffixReadReverseComplement);
	}
	uint64_t max=0;
	for(i=0;i<(uint64_t)sizeOfHashTable;i++)
	{
		if(listOfHashReadsInt[i]!=NULL)
		{
			if(listOfHashReadsInt[i][0]>max)
			{
				max=listOfHashReadsInt[i][0];
			}
		}
	}
	//for(i=0;i<(uint64_t)sizeOfHashTable;i++)
	//{
	//	if(listOfHashReadsInt[i]!=NULL)
	//	{
	//		if(listOfHashReadsInt[i][0]>1000)
	//		{
	//			listOfHashReadsInt[i][0]=1000;
	//		}
	//	}
	//}
	myfprintf("Highest entry in hash table: %"PRIu64"\n",max);
	myfprintf("Total hash miss while building the hash table: %"PRIu64"\nFunction hashPrefixesAndSuffix() in %4ld sec.\n",hashMiss, time(NULL)-second_s);
}

/* ============================================================================================
   This function returns the smallest prime larger than value.
   The size of the hash table should prime number to minimize hash collision.
   ============================================================================================ */
int64_t findNextPrimeFromFile(int64_t value)
{
	/* Create and initialize array to store hash table sizes. All values are prime numbers. */
	int64_t hashTableSizes[450]={100003, 200003, 300007, 400009, 500009, 600011, 700001, 800011, 900001, 1000003, 1769627, 1835027, 1900667, 1966127, 2031839, 2228483, 2359559, 2490707, 2621447, 2752679, 2883767, 3015527, 3145739, 3277283, 3408323, 3539267, 3670259, 3801143, 3932483, 4063559, 4456643, 4718699, 4980827, 5243003, 5505239, 5767187, 6029603, 6291563, 6553979, 6816527, 7079159, 7340639, 7602359, 7864799, 8126747, 8913119, 9437399, 9962207, 10485767, 11010383, 11534819, 12059123, 12583007, 13107923, 13631819, 14156543, 14680067, 15204467, 15729647, 16253423, 17825999, 18874379, 19923227, 20971799, 22020227, 23069447, 24117683, 25166423, 26214743, 27264047, 28312007, 29360147, 30410483, 31457627, 32505983, 35651783, 37749983, 39845987, 41943347, 44040383, 46137887, 48234623, 50331707, 52429067, 54526019, 56623367, 58720307, 60817763, 62915459, 65012279, 71303567, 75497999, 79691867, 83886983, 88080527, 92275307, 96470447, 100663439, 104858387, 109052183, 113246699, 117440699, 121635467, 125829239, 130023683, 142606379, 150994979, 159383759, 167772239, 176160779, 184549559, 192938003, 201327359, 209715719, 218104427, 226493747, 234882239, 243269639, 251659139, 260047367, 285215507, 301989959, 318767927, 335544323, 352321643, 369100463, 385876703, 402654059, 419432243, 436208447, 452986103, 469762067, 486539519, 503316623, 520094747, 570425399, 603979919, 637534763, 671089283, 704643287, 738198347, 771752363, 805307963, 838861103, 872415239, 905971007, 939525143, 973079279, 1006633283, 1040187419, 1140852767, 1207960679, 1275069143, 1342177379, 1409288183, 1476395699, 1543504343, 1610613119, 1677721667, 1744830587, 1811940419, 1879049087, 1946157419, 2013265967, 2080375127, 2281701827, 2415920939, 2550137039, 2684355383, 2818572539, 2952791147, 3087008663, 3221226167, 3355444187, 3489661079, 3623878823, 3758096939, 3892314659, 4026532187, 4160749883, 4563403379, 4831838783, 5100273923, 5368709219, 5637144743, 5905580687, 6174015503, 6442452119, 6710886467, 6979322123, 7247758307, 7516193123, 7784629079, 8053065599, 8321499203, 9126806147, 9663676523, 10200548819, 10737418883, 11274289319, 11811160139, 12348031523, 12884902223, 13421772839, 13958645543, 14495515943, 15032386163, 15569257247, 16106127887, 16642998803, 18253612127, 19327353083, 20401094843, 21474837719, 22548578579, 23622320927, 24696062387, 25769803799, 26843546243, 27917287907, 28991030759, 30064772327, 31138513067, 32212254947, 33285996803, 36507222923, 38654706323, 40802189423, 42949673423, 45097157927, 47244640319, 49392124247, 51539607599, 53687092307, 55834576979, 57982058579, 60129542339, 62277026327, 64424509847, 66571993199, 73014444299, 77309412407, 81604379243, 85899346727, 90194314103, 94489281203, 98784255863, 103079215439, 107374183703, 111669150239, 115964117999, 120259085183, 124554051983, 128849019059, 133143986399, 146028888179, 154618823603, 163208757527, 171798693719, 180388628579, 188978561207, 197568495647, 206158430447, 214748365067, 223338303719, 231928234787, 240518168603, 249108103547, 257698038539, 266287975727, 292057776239, 309237645803, 326417515547, 343597385507, 360777253763, 377957124803, 395136991499, 412316861267, 429496730879, 446676599987, 463856468987, 481036337207, 498216206387, 515396078039, 532575944723, 584115552323, 618475290887, 652835029643, 687194768879, 721554506879, 755914244627, 790273985219, 824633721383, 858993459587, 893353198763, 927712936643, 962072674643, 996432414899, 1030792152539, 1065151889507, 1168231105859, 1236950582039, 1305670059983, 1374389535587, 1443109012607, 1511828491883, 1580547965639, 1649267441747, 1717986918839, 1786706397767, 1855425872459, 1924145348627, 1992864827099, 2061584304323, 2130303780503, 2336462210183, 2473901164367, 2611340118887, 2748779070239, 2886218024939, 3023656976507, 3161095931639, 3298534883999, 3435973836983, 3573412791647, 3710851743923, 3848290698467, 3985729653707, 4123168604483, 4260607557707, 4672924419707, 4947802331663, 5222680234139, 5497558138979, 5772436047947, 6047313952943, 6322191860339, 6597069767699, 6871947674003, 7146825580703, 7421703488567, 7696581395627, 7971459304163, 8246337210659, 8521215117407, 9345848837267, 9895604651243, 10445360463947, 10995116279639, 11544872100683, 12094627906847, 12644383722779, 13194139536659, 13743895350023, 14293651161443, 14843406975659, 15393162789503, 15942918604343, 16492674420863, 17042430234443, 18691697672867, 19791209300867, 20890720927823, 21990232555703, 23089744183799, 24189255814847, 25288767440099, 26388279068903, 27487790694887, 28587302323787, 29686813951463, 30786325577867, 31885837205567, 32985348833687, 34084860462083, 37383395344739, 39582418600883, 41781441856823, 43980465111383, 46179488367203, 48378511622303, 50577534878987, 52776558134423, 54975581392583, 57174604644503, 59373627900407, 61572651156383, 63771674412287, 65970697666967, 68169720924167, 74766790688867, 79164837200927, 83562883712027, 87960930223163, 92358976733483, 96757023247427, 101155069756823, 105553116266999, 109951162779203, 114349209290003, 118747255800179, 123145302311783, 127543348823027, 131941395333479, 136339441846019, 149533581378263, 158329674402959, 167125767424739, 175921860444599, 184717953466703, 193514046490343, 202310139514283, 211106232536699, 219902325558107, 228698418578879, 237494511600287, 246290604623279, 255086697645023, 263882790666959, 272678883689987, 299067162755363, 316659348799919, 334251534845303, 351843720890723, 369435906934019, 387028092977819, 404620279022447, 422212465067447, 439804651111103, 457396837157483, 474989023199423, 492581209246163, 510173395291199, 527765581341227, 545357767379483, 598134325510343, 633318697599023, 668503069688723, 703687441776707, 738871813866287, 774056185954967, 809240558043419, 844424930134187, 879609302222207, 914793674313899, 949978046398607, 985162418489267, 1020346790579903, 1055531162666507, 1090715534754863};

	/* Set the initial size of the hash table to 3 times the size of the genome. */
	int n = 0;
	for (n=0; n<449; ++n)
		if (hashTableSizes[n] > value)
			return hashTableSizes[n];

	return hashTableSizes[n];
}


/* ============================================================================================
   Converts the economy graph to overlap graph.
   ============================================================================================ */
void convertGraph(void)
{
	myfprintf("\nIn function convertGraph().\n");
	time_t second_s=time(NULL);
	uint64_t i,j,nodea,nodeb;
	int type,length;
	sortEdgesEconomy(2); // sort first by id, then type then length. we will not add multiple edges
	for(i=1;i<=numberOfUniqueReads;i++) // Total execution time O(nd)
	{
		nodea=i;
		if(graphEconomy[i]!=NULL) // Delete the edges of a node in O(d) time.
		{
			for(j=1;j<=(uint64_t)(graphEconomy[i][0]);j++)
			{
				if(graphEconomy[i][j]!=0)
				{
					nodeb=(graphEconomy[i][j]>>24);
					if(j>1 && (graphEconomy[i][j-1]&0XFFFFFFFFFFC00000)==(graphEconomy[i][j]&0XFFFFFFFFFFC00000)) //same endpoints, same type of edge
						continue;
					if(nodea<nodeb && selfOverlappingReads[nodea]!=2 && selfOverlappingReads[nodeb]!=2)
					{
						  type=(graphEconomy[i][j]&0X0000000000C00000)>>22;
						length=(graphEconomy[i][j]&0X00000000000FFFFF);
						insertEdgeInGraph(nodea,nodeb,length,type);
					}
				}
			}
			free(graphEconomy[i]);
		}
	}
	int64_t numberOfselfOverlappingReads=0;
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		if(selfOverlappingReads[i]==1)
			numberOfselfOverlappingReads++;
	}
	free(selfOverlappingReads);
	free(graphEconomy);
	myfprintf("Self overlapping reads: %10"PRIu64"\nFunction convertGraph() in %4ld sec.\n",numberOfselfOverlappingReads,time(NULL)-second_s);
}


/* ============================================================================================
   Insert a good read from the file to the list of reads.
   ============================================================================================ */
uint64_t insertReadIntoList(char *read1, char *read2,uint64_t *tableSize)
{
	static uint64_t sizeOfList,currentPosition=1;
	char *read1Reverse,*read2Reverse;
	if(readsInt==NULL)
	{
		currentPosition=1;
		sizeOfList=100000;
		if((readsInt=(uint8_t **)malloc(sizeOfList*sizeof(uint8_t *)))==NULL) printError(MEM_ALLOC, "readsInt");
	}
	if(numberOfReads>=sizeOfList-10) // list is full
	{
		sizeOfList=sizeOfList+100000;
		if((readsInt=(uint8_t **)realloc(readsInt,sizeOfList*sizeof(uint8_t *)))==NULL) printError(MEM_ALLOC,"reallocating readsInt failed");
	}
	*tableSize=sizeOfList;
	read1Reverse=reverseComplement(read1);
	read2Reverse=reverseComplement(read2);
	if(strcmp(read1,read1Reverse)<0) // Take the lexicographically smaller one.
		readsInt[currentPosition++]=charsToBytes(read1);
	else
		readsInt[currentPosition++]=charsToBytes(read1Reverse);
	if(strcmp(read2,read2Reverse)<0) // Take the lexicographically smaller one.
		readsInt[currentPosition++]=charsToBytes(read2);
	else
		readsInt[currentPosition++]=charsToBytes(read2Reverse);
	free(read1Reverse);free(read2Reverse);
	return currentPosition-1;
}

/* ============================================================================================
   This function computes the bit representation of the reads.
   ============================================================================================ */
uint8_t *charsToBytes(char *read)
{
	int i,j=0,shift=2*(4*readArrayLength-readLength);
	uint8_t x,*returnArray;
	if((returnArray=(uint8_t *)malloc((readArrayLength)*sizeof(uint8_t)))==NULL) printError(MEM_ALLOC, "returnArray");
	for(i=0;i<(int)readArrayLength;i++)
		returnArray[i]=0;
	for(i=0;i<(int)(readLength);i++)
	{
		j=i/4;
		x=0;
		if(read[i]=='A')	x=0;
		else if(read[i]=='C')	x=1;
		else if(read[i]=='G')	x=2;
		else if(read[i]=='T')	x=3;

		returnArray[j]=(returnArray[j]<<2)|x;
	}
	returnArray[j]=returnArray[j]<<shift;
	return returnArray;
}


/* ============================================================================================
   This function returns the read from the bit representation.
   ============================================================================================ */
char *bytesToChars(uint8_t *readInt)
{
	char *array;
	int i,x;
	if((array=(char *)malloc((readLength+1)*sizeof(char)))==NULL) printError(MEM_ALLOC, "array");
	for(i=0;i<(int)readLength;i++)
	{
		x=(readInt[i>>2]>>(8-2*(i%4+1)))&0X03;
		if(x==0) array[i]='A';
		else if(x==1) array[i]='C';
		else if(x==2) array[i]='G';
		else if(x==3) array[i]='T';
	}
	array[readLength]='\0'; // String Terminating Character;
	return array;
}


/* ============================================================================================
   Fprintf function to write output in the File.
   ============================================================================================ */
void myfprintf(char * format, ...)
{
	va_list args;
	va_start (args, format);
	vfprintf (fpo_output_file, format, args);
	va_end (args);
	fflush(fpo_output_file);

}

/* ============================================================================================
   Checks if a read is good or not.
   Good read = No base pairs other than A,C,G and T and does not contain more than 80% of the
   same basepair.
   ============================================================================================ */
int isGoodRead(char *read)
{
	int i,A=0,C=0,G=0,T=0;
	for(i=0;i<(int)readLength;i++)
	{
		if(read[i]!='A' && read[i]!='C' && read[i]!='G' && read[i]!='T')
			break;
		if(read[i]=='A') A++; if(read[i]=='C') C++;if(read[i]=='G') G++; if(read[i]=='T') T++;
	}
	if(i==(int)readLength)
	{
		if(A>=numberOfSameBasepairThreshold || C>=numberOfSameBasepairThreshold || G>=numberOfSameBasepairThreshold || T>=numberOfSameBasepairThreshold)
			return 0;
		else
			return 1;
	}
		return 0;
}

/* ============================================================================================
   Saves the edges of the overlap graph in a file. Each edge is represented by three uint64_t.
   uint64_t [0] = source node ID.
   uint64_t [1] = destination node ID.
   uint64_t [2] = type and length. Type (MSB 32 bits). Length (LSB 32 bits)
   ============================================================================================ */
void saveOverlapGraph(void)
{
	if(!SAVEGRAPH) return;
	time_t second_s=time(NULL);
	char fileName[10000];
	strcpy(fileName,directoryName);
	strcat(fileName,"graph.dat");
	myfprintf("\nIn function saveOverlapGraph().\n");
	uint64_t *edgesToSave,edgesSaved=0,i;
	Edge *v;
	if((edgesToSave=(uint64_t *)malloc(3*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC, "edgesToSave");
	FILE *fpo;
	if((fpo=fopen(fileName,"w"))==NULL) printError(OPEN_FILE,fileName);
	for(i=1;i<=numberOfUniqueReads;i++)
	{
		for(v=graph[i];v!=NULL;v=v->next)
		{
			if(v->from_ID<v->ID)
			{
				edgesToSave[0]=v->from_ID;
				edgesToSave[1]=v->ID;
				uint64_t type64 = v->typeOfEdge;
				edgesToSave[2]=(type64<<32)|v->lengthOfEdge;
				fwrite (edgesToSave,sizeof(uint64_t),3,fpo);
				edgesSaved++;
			}
		}
	}
	fclose(fpo);
	myfprintf("Edges saved: %10"PRIu64"\nFunction saveOverlapGraph() in %4ld sec.\n",edgesSaved,time(NULL)-second_s);
}

/* ============================================================================================
   Reads the overlap graph from file.
   ============================================================================================ */
void readOverlapGraph(void)
{
	char fileName[10000];
	strcpy(fileName,directoryName);
	strcat(fileName,"graph.dat"); //  Previously saved graph.
	time_t second_s=time(NULL);
	myfprintf("\nIn function readOverlapGraph().\n");
	uint64_t *edgesToLoad,edgesLoaded=0;
	if((edgesToLoad=(uint64_t *)malloc(3*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC, "readsload");
	FILE *fpi;
	if((fpi=fopen(fileName,"r"))==NULL) printError(OPEN_FILE,fileName);
	while(fread(edgesToLoad,sizeof(uint64_t),3,fpi)!=0)
	{
		insertEdgeInGraph(edgesToLoad[0],edgesToLoad[1],edgesToLoad[2]&0X00000000FFFFFFFF,edgesToLoad[2]>>32);
		edgesLoaded++;
	}
	fclose(fpi);
	myfprintf("Edges loaded: %10"PRIu64"\nFunction readOverlapGraph() in %4ld sec.\n",edgesLoaded,time(NULL)-second_s);
}

/* ============================================================================================
   This function compare two strings by comparing their bit representation.
   ============================================================================================ */

int compareStringInReads(char *string1,char *string2,int start)
{
	int i;
	for(i=0;i<(int)(readLength-start);i++)
		if (string1[i+start]!=string2[i])
			return 0;
	return 1;
}

/* ============================================================================================
   Returns 64bit integers that represent a substring of 32bps.
   ============================================================================================ */
uint64_t * get64Bit2Int(uint8_t *read,int start,int length)
{
	uint64_t *number;
	if((number=(uint64_t *)malloc(2*sizeof(uint64_t)))==NULL) printError(MEM_ALLOC, "number");
	number[0]=0;number[1]=0;
	if(length<=32)
	{
		number[1]=get64BitInt(read,start,length);
	}
	else
	{
		number[0]=get64BitInt(read,start,length-32);
		number[1]=get64BitInt(read,start+length-32,32);
	}
	return number;
}


uint64_t get64BitInt(uint8_t *read,int start,int length)
{
	uint64_t number=0;
	int byte,fraction1=(start&0X3)<<1,fraction2=((start+length)&0X3)<<1;
	if(start>>2==(start+length)>>2) // Start and End in the same byte.
	{
		number=(read[start>>2]&(0XFF>>fraction1))>>(8-fraction2);
		return number;
	}
	for(byte=start>>2;byte<(start+length)>>2;byte++)
		if(byte==start>>2)
			number=read[byte]&(0XFF>>fraction1);
		else
			number=(number<<8)|read[byte];
	number=(number<<fraction2)|(read[byte]>>(8-fraction2));
	return number;
}
/* ============================================================================================
   Performs string comparison using 64 bit ints.
   ============================================================================================ */
int compareStringInBytes(uint8_t *read2bits1,uint8_t *read2bits2,int start)
{
	uint64_t *inta,*intb;
	int start1=start+hashStringLength,start2=hashStringLength,length;
	while(start1<(int)readLength)
	{
		length=minimumValue(64,readLength-start1);
		inta=get64Bit2Int(read2bits1,start1,length);
		intb=get64Bit2Int(read2bits2,start2,length);
		if(inta[0]!=intb[0] || inta[1]!=intb[1] )
		{
			free((uint64_t *) inta); free((uint64_t *) intb);
			return 0;
		}
		free((uint64_t *) inta); free((uint64_t *) intb);
		start1+=length; start2+=length;
	}
	return 1;
}


/* ============================================================================================
   Make a directory called directoryName. And delete any previous files in the directory excetp
   graph.dat file (previously saved overlap graph).
   ============================================================================================ */

int makeDirectory(char *directoryName)
{
	char deleteFolder[10000];
	int systemRet;
	mkdir(directoryName,0777); // creat the directoryName
	strcpy(deleteFolder,"rm -f "); strcat(deleteFolder,directoryName); strcat(deleteFolder,"*.txt"); // Remove any *.txt file from the directory.
	systemRet=system(deleteFolder); //clean the directoryName
	strcpy(deleteFolder,"rm -f "); strcat(deleteFolder,directoryName); strcat(deleteFolder,"*.gdl"); // Remove any *.gdl file from the directory.
	systemRet=system(deleteFolder); //clean the directoryName
	strcpy(deleteFolder,"rm -f "); strcat(deleteFolder,directoryName); strcat(deleteFolder,"*.fasta"); // Remove any *.fasta file from the directory.
	systemRet=system(deleteFolder); //clean the directoryName
	return systemRet;
}


/* ============================================================================================
   This function returns the string in an edge
   ============================================================================================ */
char *stringInEdge(Edge *edge,uint64_t *gapCount,uint64_t *nCount)
{
	char *returnString,*string; *nCount=0; *gapCount=0;
	uint64_t length=edge->lengthOfEdge+readLength+1,nextDistance=0,index,counter=0,i,nextReadDistance;
	string=bytesToChars(readsInt[edge->ID]);
	if((returnString=(char *)malloc(length*sizeof(char)))==NULL) printError(MEM_ALLOC, "returnString");
	if(edge->typeOfEdge==0 || edge->typeOfEdge==1)
		string=bytesToChars(readsReverseInt[edge->from_ID]);
	else
		string=bytesToChars(readsInt[edge->from_ID]);
	for(i=0;i<(uint64_t)(readLength);i++)
		returnString[counter++]=string[i]; // Copy the first read
	free(string);

	if(edge->listOfReads[0]==0) nextDistance=edge->lengthOfEdge;
	else
	{
		for(index=1;index<=(uint64_t)(edge->listOfReads[0]);index++) // Copy intermediate reads.
		{
			if(edge->listOfReads[index]&0X0000000000800000) // Forward orientation
				string=bytesToChars(readsInt[edge->listOfReads[index]>>24]);
			else // Reverse orientation.
				string=bytesToChars(readsReverseInt[edge->listOfReads[index]>>24]);
			nextReadDistance=((edge->listOfReads[index]&0X00000000003FF800)>>11);
			if(nextReadDistance>(uint64_t)readLength) // there is a gap
			{
				*gapCount+=1; *nCount+=nextReadDistance-readLength;
				for(i=0;i<(uint64_t)(nextReadDistance-readLength);i++) // fill in the gaps with N's
					returnString[counter++]='N';
				for(i=0;i<(uint64_t)readLength;i++) // add the read
					returnString[counter++]=string[i];
			}
			else//Not a gap, append the strting from read
			{
				for(i=(readLength-((edge->listOfReads[index]&0X00000000003FF800)>>11));i<(uint64_t)(readLength);i++)
					returnString[counter++]=string[i];
			}
			free(string);
			nextDistance=(edge->listOfReads[index]&0X00000000000007FF);
		}
	}
	if(edge->typeOfEdge==1 || edge->typeOfEdge==3)
		string=bytesToChars(readsInt[edge->ID]);
	else
		string=bytesToChars(readsReverseInt[edge->ID]);
	for(i=readLength-nextDistance;i<(uint64_t)(readLength);i++)
		returnString[counter++]=string[i];// Copy the last read.
	free(string);
	returnString[counter]='\0';
	return returnString;
}


/* ============================================================================================
   This function returns the reverse complement of a given string
   If the string is GAAAAAACCCCCCG
   It will retrun   CGGGGGGTTTTTTC
   ============================================================================================ */

char* reverseComplement(char *string)
{
	char *str;
	int stringLength1=strlen(string);
	if((str=(char *)malloc((stringLength1+1)*sizeof(char)))==NULL) printError(102, "str");
	if(stringLength1<(int)readLength)
	{
		strcpy(str,string);
		return str;
	}
	int i;
	int j=0;
	for(i=stringLength1-1;i>=0;i--)
	{
		if(string[i]=='A' || string[i]=='a')
			str[j]='T';
		else if(string[i]=='C' || string[i]=='c')
			str[j]='G';
		else if(string[i]=='G'|| string[i]=='g')
			str[j]='C';
		else if(string[i]=='T' || string[i]=='t')
			str[j]='A';
		j++;
	}
	str[j]='\0';
	return str;
}

void whereIs(uint64_t read)
{
	myfprintf("In function whereIs()\n");
	struct readToEdgeStructure *edge;
	if(read_to_edge[read]!=NULL) // Nodes
	{
		for(edge=read_to_edge[read];edge!=NULL;edge=edge->next)
		{
			myfprintf("Read %"PRIu64" in Edge (%"PRIu64",%"PRIu64") length: %"PRIu64" flow: %f Number of Reads: %"PRIu64"\n",read,edge->edge->from_ID,edge->edge->ID,edge->edge->lengthOfEdge,edge->edge->flow,edge->edge->listOfReads[0]);
		}
	}
	myfprintf("\n");
}

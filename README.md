# SAGE
String-overlap Assembly of GEnomes
***************************************************
# NOTE - Please use SAGE2 found [here](https://github.com/lucian-ilie/SAGE2); SAGE is no longer maintained. 
***************************************************

SAGE is a new string-overlap graph-based de novo genome assembler. To run SAGE, first correct the input dataset using RACER, then use the command: 

     SAGE [inputFile(s)] [outputDir] [minOverlapLength] 

where 
 - SAGE is the appropriate binary used, e.g., SAGE_Linux
 - [inputFile(s)] is the list of input files with (corrected) reads
 - [outputDir] contains the assembly produced
 - [overlapLength] is the length of the minimum overlap between reads

SAGE assumes all reads have the same length and the paired reads are interleaved. 

The assemblies for the datasets in the paper can be found below:
 - [B.subtilis](http://www.csd.uwo.ca/~ilie/SAGE/B.subtilis.fasta)
 - [C.trachomatis](http://www.csd.uwo.ca/~ilie/SAGE/C.trachomatis.fasta)
 - [S.pseudopneumonie](http://www.csd.uwo.ca/~ilie/SAGE/S.pseudopneumonie.fasta)
 - [F.tularensis](http://www.csd.uwo.ca/~ilie/SAGE/F.tularensis.fasta)
 - [L.interrogans](http://www.csd.uwo.ca/~ilie/SAGE/L.interrogans.fasta)
 - [P.gingivalis](http://www.csd.uwo.ca/~ilie/SAGE/P.gingivalis.fasta)
 - [E.coli](http://www.csd.uwo.ca/~ilie/SAGE/E.coli.fasta)
 - [C.thermocellum](http://www.csd.uwo.ca/~ilie/SAGE/C.thermocellum.fasta)
 - [C.elegans](http://www.csd.uwo.ca/~ilie/SAGE/C.elegans.fasta)

# CITE
If you use SAGE, please cite:

L. Ilie, B. Haider, M. Molnar, R. Solis-Oba, [SAGE: String-graph Assembly of GEnomes, BMC Bioinformatics](http://www.biomedcentral.com/1471-2105/15/302/abstract) 15 (2014) 302.

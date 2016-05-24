Introduction
------------

This program optimizes access files and parallelizes the jobs alignment with BWA-MEM alignment v0.7.12.
The input are fasta files of pair reads sequenced with Illumina technology (see sample files).  
Batch of 100M pair bases are loaded and aligned assuring the result is identical to classic BWA-MEM. 
This program optimize access file and parallelized the jobs for BWA-MEM alignment.
The input are fastq files of pair reads sequenced with Illumina technology. 
Batch of 100M bases are loaded and aligned assuring the result is identical to classic BWA. 

Requirements
------------

You need a C compiler as required for classic BWA program.
You need a mpi compiler too. to check your mpi installation tell in a command window whereis mpirun normally it is installed in /usr/bin/mpirun.
This program runs on supercomputer architecture and supports also NFS file system. 
A classic 1Gb or 10Gb network is sufficient.

Your reads should be paired or single but not trimmed.

Known issues:
--------------------------

Primary hits are reproduced between the serial version and the parallel but you can see differences in mapping position for alternate contigs.
this issues is under investigation.  

Compilation 
-----------

To compile the source open the makefile tell where mpicc is installed (in CC directive of makefile) and make it.
Results are 2 executables pidx and pbwa7.

Build a reference
-----------------
After the creation of the reference file with BWA, you need to build a mapped reference genome. 
To do that: pidx my_ref.fa (where my_ref.fa has been build with BWA).
Pidx build a reference with the extension .map (my_ref.fa.map). 
This reference will be mapped in share memory.

If you want to speed-up the mapping of the reference in shared memory you can copy the reference in /tmp of each servers before launching the alignment. 

Start alignment
---------------

Example of bash to run the pbwa7 on Torque scheduler

SPLITTED_READS_DIR=/data/tmp/WholeGenomeTestSplittedReads

MAIL="mymail@toto.fr"

PBS_OUTPUT=/data/Test/OUTPUT

PBS_ERROR=/data/Test/ERROR

pBWA_BIN_DIR=/data/Test/mpiBWA

BWA_REF_TMP=/data/tmp/fjarlier/BWA_reference_idx/hg19.fasta

SAMPLENAME="MPIBWA"

FILE_TO_ALIGN_R1=/data/HCC1187C/HCC1187C_R1_10X.fastq

FILE_TO_ALIGN_R2=/data/HCC1187C/HCC1187C_R2_10X.fastq

TOTAL_PROC=600

OUTPUT_DIR=/data/Test/RESULTS/

FILE_TO_WRITE=/data/Test/RESULTS/test.sam

Only for Lustre usage tell the striping of the results

lfs setstripe -c -1 -s 3g $OUTPUT_DIR 

launch the jobs with torque

echo " mpirun -n $TOTAL_PROC $pBWA_BIN_DIR/pbwa7 mem -t 1 -o $FILE_TO_WRITE $BWA_REF_TMP $FILE_TO_ALIGN_R1 $FILE_TO_ALIGN_R2" | qsub -o $PBS_OUTPUT -e $PBS_ERROR -N ${SAMPLENAME} -q batch  -l nodes=40:ppn=15

or launch the jobs with mpirun

mpirun -n $TOTAL_PROC $pBWA_BIN_DIR/pbwa7 mem -t 1 -o $FILE_TO_WRITE $BWA_REF_TMP $FILE_TO_ALIGN_R1 $FILE_TO_ALIGN_R2

Remarks
-------

1) Do not type the .map extension when you give the reference to pbwa7

2) For Lustre or parallel file system users. If you intend to run the mpiSort after the alignment you have to tell the striping of the results. 
This is done according to the striping you set in the mpiSort program with lfs setstripe command (lfs setstripe -c -1 -s 2gb .). 


Authors
-------

This program has been developed by 

Frederic Jarlier from Institut Curie and 
Nicolas Joly from Institut Pasteur

and supervised by

Philippe Hupe from Institut Curie

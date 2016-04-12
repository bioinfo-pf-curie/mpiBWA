Introduction
------------


This program optimizes access files and parallelized the jobs for BWA-MEM alignment v7.12.
The input are fasta files of pair reads sequenced with Illumina technology. 
Batch of 100M pair bases are loaded and aligned assuring the result is identical to classic BWA-MEM. 
=======
This program optimize access file and parallelized the jobs for BWA-MEM alignment.
The input are fastq files of pair reads sequenced with Illumina technology. 
Batch of 100M bases are loaded and aligned assuring the result is identical to classic BWA. 

Requirements
------------

You need a C compiler as required for classic BWA program.
You need a mpi compiler too. This program runs on supercomperter architecture and supports also NFS file system. 
=======
This program runs on supercomputer architecture but supports also NFS file system. 

A classic 1Gb or 10Gb network is sufficient.

Compilation 
-----------

To compile the source open the makefile tell where mpicc is installed and make it.
Results are 2 executables pidx and pbwa7.

Build a reference
-----------------

You need to build a mapped reference genome. 
To do that: pidx my_reg.fa 

pidx build a reference with the extension .map: my_ref.fa.map 

This reference will be mapped in share memory.

To speed-up the mapping of the reference copy it in /tmp of each servers before launching the alignment. 

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

tell the striping of the results

lfs setstripe -c 128 -s 3g $OUTPUT_DIR 

launch the job with torque

echo " mpirun -n $TOTAL_PROC $pBWA_BIN_DIR/pbwa7 mem -t 1 -o $FILE_TO_WRITE $BWA_REF_TMP $FILE_TO_ALIGN_R1 $FILE_TO_ALIGN_R2" | qsub -o $PBS_OUTPUT -e $PBS_ERROR -N ${SAMPLENAME} -q batch  -l nodes=40:ppn=15

launch the job with mpirun

mpirun -n $TOTAL_PROC $pBWA_BIN_DIR/pbwa7 mem -t 1 -o $FILE_TO_WRITE $BWA_REF_TMP $FILE_TO_ALIGN_R1 $FILE_TO_ALIGN_R2

Remarks
-------

If you intend to run the mpiSort after the alignment you have to tell the striping of the results. 
This is done according to the striping you set in the mpiSort program with lfs setstripe command. 

Authors
-------

This program has been developed by 

Frederic Jarlier from Institut Curie
Nicolas Joly from Institut Pasteur

and supervised by

Philippe Hupe from Institut Curie

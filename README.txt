Release notes
------------

Release 1.0 from 07/11/2016

1) Support for trimmed reads. Works with an even number of jobs. Tested with toy data up to 10 jobs.
Need more test for scalability, and load balancing of chunks sizes.

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

You need to install a version of MPI.

You need a mpi compiler too. to check your mpi installation tell in a command window whereis mpirun, normally it is installed in /usr/bin/mpirun.
This program runs on supercomputer architecture and supports also NFS file system. 
A classic 1Gb or 10Gb network is sufficient.

Your reads should be paired or single.

Options
------

All the BWA-MEM options are available in this version, of course according to the bwa release wrapped. 

Known issues:
-------------

1) Primary hits are reproduced between the serial version and the parallel but you can see differences in mapping position for alternate contigs. 
This problem stems from the randomization of multi-hits reads. When running with the same number of MPI jobs alternative positions are reproduced but when the number of jobs varies the positions can switch for secondary alignments.

2) The reference genome is an image in binary format of the .fa reference file and its extension is .map (from the hgXX.fa we create a hgXX.fa.map). On some architecture the mmaping of the references .map can be slow. Ordinary mmaping of the should 1 or 2 minutes. To solve this issue before the mapping copy with cp or rsync the reference file .map in the /tmp of each computers nodes before aligning. For instance with the command "mprun -n $NUM_CPU cp $BWA_REF /tmp".

How to integrate further version
--------------------------------

This version of mpiBWA has been build with 0.7.12 BWA version.
To integrate the 0.7.13 or 0.7.14:

1) git clone the 0.7.14 of BWA.

2) in the folder of bwa copy-pass the following function from mpiBWA:

makefile
main_parallel_version.c 
pidx.c

then make.

Compilation 
-----------

You need automake 1.15 for the installation.
You can install automake and autoconf in differents directories and export the path like this:
export PATH=../automake-1.15/bin:../autoconf-2.69/bin:$PATH

Download from git. In the folder mpiBWA type:
./configure && make install && make

or for distribution:
make dist
tar xzf .tar.gz
cd pbwa7-1.0
./configure && make install && make

for passing mpi path:
./configure CC=mpi_bin_path
add --prefix in configure if you need 

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

1) Do not type the .map extension when you give the reference to pbwa7.

2) Hybrid mode is possible with -t options. MPI rank fix the number of servers and -t options the number of threads per job. 

3) For Lustre or parallel file system users: If you intend to run the mpiSort after the alignment you can tell the striping of the results. 
This is done according to the striping you set in the mpiSort program with the Lustre "lfs setstripe" command (lfs setstripe -c -1 -s 2gb .). 
The "-c -1" option tells Lustre to use all the file system servers and and "-s 2gb" is the size of contigues data blocks. 
For speed purpose reading commands (particularly MPI commands) are aligned on those blocks.

Notes: The striping (like with Hadoop) is the way your data is distributed among servers. This technic accelerates access files and meta file information.

Authors
-------

This program has been developed by 

Frederic Jarlier from Institut Curie and 
Nicolas Joly from Institut Pasteur

and supervised by

Philippe Hupe from Institut Curie

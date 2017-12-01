Your are in the master branch of the the mpiBWA project.

This version is dedicated to fit a large variety of cluster.

If you are equiped with low latency network and network file system (NFS) use this branch.

If your are equiped with Lustre and Low latency network use the FULLMPI branch of the project.
  
If your are equiped with NFS and high latency network and want to go very fast but less precise use the LAZYCHUNK branch.

Add an experimental branch for higher scalability, full reproducibility and better accuracy.

Release notes
------------

Release 1.0 from 01/12/2017

Changes in the branch Experimental. <br />
Fix the mpi_read_at buffer limit size. Now we load the buffer in chunk of 1gb. <br />

Release 1.0 from 29/11/2017

Add a new branch called Experimental. <br />
Warning: This is experimental work do not use in production. But test it and send us reports.<br />

Rationnal:<br />

The master and FULLMPI branches are made for full reproducibility (independant to the number of jobs) and accuracy but they reach the Amdah'ls law point. 
Indeed the locking file RMA implementation (the serialization when computing offsets) introduces a bottle neck we are not able to overpass with RMA technics. <br />
 
This is why we have implemented a new algorithm. Now a lot more master jobs are responsible for chuncking the data the way bwa does and present it to bwa-mem aligner. 
And instead of doing it linearly on the fastq now they do it independently and in parallel. With a little inter communication they adjust the chunks offets and sizes. 
This method removes the serialization bottle neck. <br />

As in the master and FULLMPI branches we obtain a full reproducibility and with a better efficiency and scalability. <br />

When testing this branch make sure the total number of jobs you take is (master jobs) * 8. <br />
8 is the number of aligner threads used by bwa-mem. <br />
According to bwa-mem policy all chunks are 10e6 by the number of threads nucleotide bases big. <br />

Remark: The initial buffer of each master jobs is limited to 2gb (due to mpi_read_at buffer size). <br />
This version does not work on trimmed reads. <br />

First results test on broadwell. <br />

Sample: <br />

NA12878 Illumina 300X 2x150 WGS from GIAB chinese trio.

The alignement is done with 352*8 = 2816 cpu<br />
352 = (Forward Fastq size in gb) / 2g,  8 is thenumber of bwa-mem aligner jobs (8 per master jobs)<br />

MPI parameters :<br />
mpi_run -n 352 -c 8<br />
or : <br />
MSUB -n 352<br />
MSUB -c 8<br />

alignment time: 26 mn<br />
Time to compute chunks: 8s<br />

Reproducibility: the pipeline has been tested tested with 5632 and 2816 cpu results are the same. <br />

Next step:<br />

1) remove initial buffer size limitation<br />
2) support trimmed reads<br />
3) need tests on low throughput infrastructures<br />

Release 1.0 from 21/11/2017

Changes in  LAZYCHUNCK branch:

1) forget a end condition

Release 1.0 from 20/11/2017

Changes in  LAZYCHUNCK branch:
1) remove memory leaks
2) update results section
3) test with 600Gb (x2) fastq files and 100 jobs: ok. (20 mn to approximate chuncks offset) 
4) reproducibility with constant number of jobs: ok.
 
Notes:
If want 100% reproducibility whatever the number of jobs use the master branch or the FULLMPI branch.
But if you want go faster use LAZYCHUNK branch. 
 
Release 1.0 from 17/11/2017

Changes in  LAZYCHUNCK branch

1) Improvement of the algorithm for window sizing.
2) Change the number of bases for the estimation line 541. 

Release 1.0 from 15/11/2017

Changes in LAZYCHUNK branch

1) Due to a mmap in the previous version the virtual memory may be big when the fastq file is large.
Some scheduler like Torque doesn't like this. We review the algorithm of computing chuncks.
first we divide the file in window of size (number of jobs) * 1Gb. On each window jobs approximate chuncks of 10 mega bases.
this way the virtual memory stay low. 
 
2) We also saw problem on some architecture with MPI file read shared we replace it with MPI file read at. 

3) We have also other projects of parallelization (sorting, marking duplicates, clustering... ) and we are looking for people willing to help us for developments and tests. Don't hesitate to contact me (frederic.jarlier@curie.fr). 


Release 1.0 from 30/06/2017

Major changes:

1) new algorithm. 

We mimic original BWA algorithm when spliting chunk of reads.

In this version 2 workers are responsible for computing offsets chuncks. 
Then offset workers sends the offsets to aligning workers.

To change the number of offset workers modify the line 65 (we will pass it in parameters in the next release).
Adapt the number of offset workers according to the number of aligners.
We have tested 2 offset workers for 80 aligners with a NFS and low ethernet with a good effciency.
 
This version has a low adherence with MPI in order to be compatible with NFS and low ethernet connectors.
Tested with openMPI 2.1.1.


Release 1.0 from 13/06/2017

Update version of bwa 7.15

Release 1.0 from 30/05/2017

To overcome lock contention problem we use a RMA MCS lock inspired from Latham R. et al. (2007). 


Release 1.0 from 17/05/2017

To enhance reproducibility and accuracy we have implemented a new method for read's chunk evaluation. 
The reads chunk sizes are now identical compare to serial BWA-MEM for both trimmed reads or not.

But the algorithm suffers from scalability. The algorithm scales up to 300 cpu (test on Cobalt TGCC) with same speed up than previous version. 
Above that number we encounter lock contention. We are working on a solution now. 
Go back to february 2017 version for more scalability. 
 
Release 1.0 from 20/01/2017

No need to copy the reference genome in /tmp before the mapping. This is done automatically by MPI. 
To change the target directory set $TMPDIR to the new location.

Release 1.0 from 07/11/2016

1) Support for trimmed reads. Works with an even number of jobs. Tested with sample data up to 10 jobs.
Need more test for scalability, load balancing, and evaluation of the chunks sizes.

Introduction
------------

This program optimizes access files and parallelizes the jobs alignment with BWA-MEM alignment v0.7.12.
The input are fasta files of pair reads sequenced with Illumina technology (see sample files). 
Batch of 100M pair bases are loaded and aligned assuring the result is identical to classic BWA-MEM. 

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

Primary hits are reproduced between the serial version and the parallel but you can see differences in mapping position for alternate contigs. 
This problem stems from the randomization of multi-hits reads. When running with the same number of MPI jobs alternative positions are reproduced but when the number of jobs varies the positions can switch for secondary alignments.

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

The pidx needs the following files my_ref.fa.sa, my_ref.fa.bwt, my_ref.fa.ann, my_ref.fa.pac, my_ref.fa.amb to construct the fa.map.
 Those files are generated with bwa index. It works also with the fasta.

How to manage the multithreading
----------------

You have the option to use the multithreading options. 
To do that you fix the number of nodes with mpirun and then the number of threads per nodes with the -t.

Examples:

1 server with 8 threads per server:
mpirun -n 1  pbwa7 mem -t 8... 

2 servers with 10 threads per servers
mpirun -n 2  pbwa7 mem -t 10... 

And so on.

Alignment with a scheduler
------------

Example of bash to run the pbwa7 on Torque scheduler

SPLITTED_READS_DIR=../WholeGenomeSampleReads
MAIL="mymail@toto.fr"
PBS_OUTPUT=../OUTPUT
PBS_ERROR=../ERROR
pBWA_BIN_DIR=../mpiBWA
BWA_REF_TMP=../hg19.fasta
SAMPLENAME="MPIBWA"
FILE_TO_ALIGN_R1=../myread_forward.fastq
FILE_TO_ALIGN_R2=../myread_backward.fastq
TOTAL_PROC=600
OUTPUT_DIR=../RESULTS/
FILE_TO_WRITE=../RESULTS/test.sam

//launch the jobs with torque and one job per core (-t 1)
echo " mpirun -n $TOTAL_PROC $pBWA_BIN_DIR/pbwa7 mem -t 1 -o $FILE_TO_WRITE $BWA_REF_TMP $FILE_TO_ALIGN_R1 $FILE_TO_ALIGN_R2" | qsub -o $PBS_OUTPUT -e $PBS_ERROR -N ${SAMPLENAME} -q batch  -l nodes=40:ppn=15

or launch the jobs with mpirun

mpirun -n $TOTAL_PROC $pBWA_BIN_DIR/pbwa7 mem -t 1 -o $FILE_TO_WRITE $BWA_REF_TMP $FILE_TO_ALIGN_R1 $FILE_TO_ALIGN_R2

Results
-------

See the results here:

http://devlog.cnrs.fr/jdev2017/posters

in the section HPC@NGS


Remarks
-------

1) Do not type the .map extension when you give the reference to pbwa7.

2) Hybrid mode is possible with -t options. MPI rank fix the number of servers and -t options the number of threads per job. 

3) For Lustre or parallel file system users: If you intend to run the mpiSort after the alignment you can tell the striping of the results. 
This is done according to the striping you set in the mpiSort program with the Lustre "lfs setstripe" command (lfs setstripe -c -1 -s 2gb .). 
The "-c -1" option tells Lustre to use all the file system servers and and "-s 2gb" is the size of contigues data blocks. 
For speed purpose reading commands (particularly MPI commands) are aligned on those blocks.

Notes: The striping (like with Hadoop) is the way your data is distributed among servers. This technic accelerates access files and meta file information.

4) Make sure you have the same version of MPI for compiling and running.

5) From our experiences some scheduler do not interpret shared memory and jobs ends up stuck for overloading memory.
In this case take 6gb per job. This is approximately the total reference plus the chunk size.


Future work:
----------
1) Test the reproducibility after randomization of the fastq.

2) Manage the randomization of alternate contigs. To mimic original algorithm.

3) Manage the insert size statistics between jobs.

References:
---------

This work is based on the original bwa aligner written by Li et al.

Li H. and Durbin R. (2010) Fast and accurate long-read alignment with Burrows-Wheeler transform. Bioinformatics, 26, 589-595. [PMID: 20080505] 

Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v1 [q-bio.GN]

Latham R. et al. (2007)  Implementing MPI-IO Atomic Mode and Shared File Pointers Using MPI One-Sided Communication Authors 

-------

The program has been developed by 

Frederic Jarlier from Institut Curie and 
Nicolas Joly from Institut Pasteur

and supervised by

Philippe Hupe from Institut Curie

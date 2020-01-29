
You are in the Master branch of the mpiBWA project
In this branch we implement new algorithm for chunking the data before sending them to bwa-mem.
We want to be fast, scalable, accurate and reproducible whatever the number of jobs you chose.

Rationnal:<br />
 
We have implemented a new algorithm. A lot more master jobs are responsible for chuncking the data the way bwa does and present it to bwa-mem aligner. 
And instead of doing it linearly on the fastq now they do it independently and in parallel. With a little inter communication they adjust the chunks offets and sizes. 
This method removes the serialization bottle neck. <br />

When testing this branch make sure the total number of jobs you take is (master jobs) * 8. <br />
8 is the number of aligner threads used by bwa-mem. <br />
According to bwa-mem policy all chunks are 10e6 by the number of threads nucleotide bases big. <br />

First results test on broadwell. <br />

Sample: <br />

NA12878 Illumina 300X 2x150 WGS from GIAB chinese trio.

The alignement is done with 352*8 = 2816 cpu<br />
352 = (Forward Fastq size in gb) / 2g,  8 is thenumber of bwa-mem aligner jobs (8 per master jobs)<br />
MPI parameters :<br />
mpi_run -n 352 -c 8<br />
or : <br />
MSUB -n 352 <br />
MSUB -c 8 <br />

alignment time: 26 mn<br />
Time to compute chunks: 8s<br />

Reproducibility: the pipeline has been tested tested with 5632 and 2816 cpu results are the same. <br />


Installation
---------
 
see INSTALL.md


Requirements
------------

You need a C compiler as required for classic BWA program. <br />
You need to install a version of openMPI. (see: https://www.open-mpi.org/). <br /> 
You need a mpi compiler too. to check your mpi installation tell in a command window whereis mpirun, normally it is installed in /usr/bin/mpirun. <br />
You need a low latency network and a parallel file system to use this branch. <br />
Your reads can be single, paired, trimmed or not. <br />

Options
------

All the BWA-MEM options are available in this version, of course according to the bwa release wrapped. 

How to integrate further version
--------------------------------

This version of mpiBWA has been build with 0.7.15 BWA version.
To integrate the 0.7.+:

1) git clone the 0.7.+ of BWA. <br />

2) in the folder of bwa copy-pass the following function from mpiBWA: <br />

makefile.am <br />
configure.ac <br />
install-sh <br />
main_parallel_version.c <br /> 
pidx.c <br />

then <br />

see INSTALL.md

Compilation 
-----------

This version ask for automake 1.15 during installation. <br />
If you don't have 1.15 change in the configure.ac the line  <br />
AM_INIT_AUTOMAKE([1.15 foreign -Wall]) with AM_INIT_AUTOMAKE([1.13 foreign -Wall])  <br />

You can install automake and autoconf in differents directories and export the path like this: <br />
export PATH=../automake-1.15/bin:../autoconf-2.69/bin:$PATH <br />

Download from git. In the folder mpiBWA type: <br />
./configure && make install && make <br />

or for distribution: <br />
make dist <br />
tar xzf .tar.gz <br />
cd pbwa7-1.0 <br />
./configure && make install && make <br />

for passing mpi path: <br />
./configure CC=mpi_bin_path <br />
add --prefix in configure if you need <br /> 

Results are 2 executables pidx and pbwa7. <br />

Build a reference
-----------------
After the creation of the reference file with BWA, you need to build a mapped reference genome. 
To do that: <br />

pidx my_ref.fa (where my_ref.fa has been build with BWA).

<br />

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

Example of command lines
---------------------

Example of bash to run the pbwa7 with openmpi 1.10.7

We launch 30 master jobs with 8 threads each it make a total of 240 jobs 

ERROR=$PATH/error.txt <br />
SAM=$PATH/mysample.sam <br />
TASKS=30 <br />
PPN=8 <br />
mpirun -np $TASKS --map-by ppr:$PPN:socket:pe=$PE $PBWA mem -t 8 -o $SAM $REF $FASTQ1 $FASTQ2 &> $ERROR <br />

Results
-------

Here results wo obtain from test realized with TGCC (Très Grand Centre de Calcul - Bruyères le Chatel - France).

![img](Results_TGCC_Broadwell.jpg)

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
In this case take less than 1gb per job per threads (if you take 8 threads). This is approximately the total reference plus the chunk size.

Questions:
--------

We have created a google forum for questions.

https://groups.google.com/forum/#!forum/hpc-bioinformatics-pipeline

Feel free to ask any question

Future work:
----------

1) Overlap writing, reading and aligning.

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

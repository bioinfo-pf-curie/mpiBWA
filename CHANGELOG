
Release notes
------------

Release 1.0 from the 15/01/2020<br />

1) Add a version of the main that align and split the result by chromosom. <br />
We create the output file with the header and SAM files by chromosom, a SAM for discordant reads and unmapped. <br /> 

Rational:

The current version of mpiBWA creates one big SAM file. In some cases this SAM is big of several tera bytes. <br />
This is difficult for the sorting to deal with files that big. The idea with this version is to create a SAM file by chromosom. <br />
Each chromosom's name come from the header line of the genome reference. This way the sam file to sort and markdup is much smaller so we need less RAM and CPU by chromosom. <br />
The extra overhead of the splitting is negligible compare with previous version.<br />
Now we can sort individual chromosom in parallel, for instance the chr1 of 300X WGS is equivalent to a 30X. <br />

Warning: 

This version is under construction. It has been tested in the case the fastq files are trimmed.
Another tests are conducted. The sorting program is not updated to deal with chromosom independantly.<br />
This is under construction too.<br />

How to use: 

In the Makefile.am replace the line <br />

pbwa7_SOURCES = main_parallel_version.c <br />
with <br />
pbwa7_SOURCES = main_parallel_version_split_by_chr.c <br />

What's next:

1) Shall we include secondary alignment
2) Shall we include discordant reads in the chromosom it belongs 


Release 1.0 from the 12/12/2019<br />

1) remove MPI call after finalyze

Release 1.0 from the 3/04/2019<br />

Changes in Master branch 

1) Add support of single read trimmed or not

Release 1.0 from the 10/07/2018<br />

Changes in Master branch 

1) Fix a bug during the mapping in shared memory of the reference genome
This bug didn't appear with openMPI version but Intel compiler complains.

2) Creation of a google group <br />
	https://groups.google.com/forum/#!forum/hpc-bioinformatics-pipeline <br />

3) To improve performances on Lustre file system removed the “suid” mount option (rw,nosuid,flock,lazystatfs) <br />
With Beegfs set the "flock" to "on" for reproducibility.<br />

Release 1.0 from the 30/04/2018 <br />

Changes in Master branch <br />

1) Add support for trimmed reads.  <br />
2) Be aware of the flock mode on parallel file system (Lustre, beegfs): flock must be on.
Otherwise reproducibility is not guaranteed. <br />

Release 1.0 from the 21/03/2018 <br />

1) Fix an inversion in file handle (line 1039 and 1056) <br />

Release 1.0 from 23/12/2017

changes in Master branch <br />

1) 100% reproducibility with the pipeline control (bwa mem -t 8) <br />
2) remove memory limits now offset are computed on 1gb buffer <br />
3) remove some memory leaks <br />
4) optimization of code <br />

Release 1.0 from 04/12/2017

Changes in the branch Master.

1) Fix a memory leak.

Release 1.0 from 01/12/2017

Changes in the branch Master. <br />
Fix the mpi_read_at buffer limit size.<br />


Release 1.0 from 29/11/2017

Add a new branch called Experimental. <br />
Warning: This is experimental work do not use in production. But test it and send us reports.<br />

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
=======

First release

We have implemented a new algorithm.
In this algorithm you have master jobs and aligner jobs.
Master jobs are responsible for computing offset and chunk sizes.
This way all the chuncks have the same number of bases. 
Then master jobs send chunks to bwa-mem and are aligned.
Finally master jobs write ni the result sam file.

the total number of jobs muste be: (number of master) * 8 
8 is the number of threads used by bwa-mem.

First result test on broadwell of TGCC

Condition: Due to mpi_read_at buffer limit of 2g. 
You should limit the initial buffer read to 2gb per master jobs.
Futur release will remove this condition. 

Tested on the SRR2052 WGS with 352*8 cpu 
alignment time: 26 mn
Time to compute chunks: 8s
Scalability : ok
Reproducibility : ok 
Command line :mpi_run -n 352 -c 8
or : MSUB -n 352
MSUB -c 8



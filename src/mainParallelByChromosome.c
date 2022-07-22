/*
This file is part of mpiBWA

The project was developped by Frederic Jarlier from Institut Curie and Nicolas Joly from Institut Pasteur

NGS aligner inspired by BWA-MEM 

Copyright (C) 2016-2021  Institut Curie / Institut Pasteur

You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
   
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
   
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

*/

#define _GNU_SOURCE

#include <sys/mman.h>
#include <sys/stat.h>

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <inttypes.h> /* For PRIu64 */
#include <libgen.h>
#include <limits.h>   /* For PATH_MAX */
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <pthread.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bwa.h"
#include "bwamem.h"
#include "utils.h"
#include "tokenizer.h"
#include "fixmate.h"
#include "parallel_aux.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

//due to mpi_read_at limit buffer 1gb
#define DEFAULT_INBUF_SIZE (1024*1024*1024)
//#define DEFAULT_INBUF_SIZE 1024

/* We require 64bit offsets */
#ifndef MPI_OFFSET
#define MPI_OFFSET MPI_LONG_LONG
#endif

#ifndef OPEN_MPI
#define OPEN_MPI 0
#endif

#ifdef TIMIMG
#define xfprintf fprintf
#else
#define xfprintf(...) /**/
#endif

#define STRIPING_FACTOR "6"
//#define STRIPING_UNIT "4194304"   // 4 MB
//#define STRIPING_UNIT "268435456"   // 256 MB
//#define STRIPING_UNIT "536870912"   // 500 MB
#define STRIPING_UNIT "1073741824"  // 1GB
//#define STRIPING_UNIT "1610612736"  // 1.5GB
//#define STRIPING_UNIT "2147483648"  // 2GB
//#define STRIPING_UNIT "2684354560"  // 2.5GB
//#define STRIPING_UNIT "3221225472"  // 3GB
//#define STRIPING_UNIT "3758096384"  // 3.5GB


#define NB_PROC  "16" //numer of threads for writing
#define CB_NODES "2" //numer of server for writing
#define CB_BLOCK_SIZE  "268435456" /* 256 MBytes - should match FS block size */
#define CB_BUFFER_SIZE  "3758096384" /* multiple of the block size by the number of proc*/
#define DATA_SIEVING_READ "enable"
#define min(a,b) (a>=b?b:a)
#define MAX_CHAR_SIZE 2048
#define MAX_CHR_NAME_SIZE 200
#define SMALL_STACK (1024*1024)
#define BIG_STACK (1024*1024*512)

int main(int argc, char *argv[]) {
	

	const char *mode = NULL;
	char *progname = basename(argv[0]);
	char *file_r1 = NULL, *file_r2 = NULL;
	char *buffer_r1, *buffer_r2;
	char *file_out = NULL;
	char *shared_mem = NULL;
	char *file_ref = NULL;
	char *rg_line = NULL, *hdr_line = NULL, *pg_line = NULL;
	char file_map[PATH_MAX];
	char *p, *q, *s, *e;

	int fd_in1;
	int proc_num, rank_num;
	int fixed_chunk_size = 0;
	int res, count;
	int files, nargs;
	int dofixmate = 0;
	int c, copy_comment = 0;
	int ignore_alt = 0;
	int write_format = 2; //default is SAM
    	int compression_level = 3;
	
	double bef, aft;
	size_t n = 0;
	off_t maxsiz;
	struct stat stat_r1, stat_r2, stat_map;

	MPI_File fh_r1, fh_r2;
	MPI_Status status;
	MPI_Win win_shr;

	mem_opt_t *opt, opt0;
	mem_pestat_t pes[4], *pes0 = NULL;
	bwaidx_t indix;
	bseq1_t *seqs;

	if (argc < 2) {

		fprintf(stderr, "program: %s is a MPI version of BWA MEM\n"
			"version: %s\n"
			"\nusage : mpirun -n TOTAL_PROC %s mem -t 8 [-f] [-b] -o SAM_FILE REFERENCE_GENOME FASTQ_R1 [FASTQ_R2]\n"
            "\n\tTOTAL_PROC tells how many cores will be used by MPI to parallelize the computation.\n"
			"\nrequirements : from the reference genome index file generated with the command 'bwa index'\n"
            "\tyou need to create a reference genome map file with 'mpiBAWIdx' that comes along\n"
            "\twith this program as follows:\n"
            "\n\t\tmpiBWAIdx myReferenceGenome.fa\n\n"
			"\tIt creates a .map file that will be used in shared memory as reference genome.\n"
            "\ninput:\n"
	    	"\t-f to fix the mate on the fly for compatibility with samtools markdup (optional)\n"
	    	"\t-b to write output in bam format (optional)\n"
		"\t-g to write output in BGZF format (optional)\n"
		"\t-z to tell where to place the reference genome (available with openMPI: socket, numa, l1, l2, l3, shared) (optional)\n"		
            "\tREFERENCE_GENOME: reference genome name (e.g. myReferenceGenome.fa).\n"
            "\t\tDo not provide the '.map' extension of the file genareted with 'mpiBWAIdx'\n"
            "\n\tFASTQ_R1: fastq file for R1\n"
            "\n\tFASTQ_R2: fastq file for R2 is the data come from paired-end sequencing (optional)\n"
            "\noutput:\n"
            "\tSAM_FILE\n"
            "\n\tIndividual chrN.sam files with aligned reads on each chrN (ChrN are the chromosome name\n"
            "\tfrom the header). The chrN.sam contains the header for the chrN but the positions\n"
            "\tare not sorted. They can be sorted independently with mpiSORT.\n"
            "\n\tThe file discordant.sam contains chimeric alignments.\n"
            "\n\tThe unmapped.sam contains unmapped reads.\n"
            "\noptions: 'bwa mem' options can be passed from command line (e.g. mpiBWA mem -t 8 -k 18)\n"
            "\nFor more detailed documentation visit:\n"
            "\thttps://github.com/bioinfo-pf-curie/mpiBWA\n"
            "\nCopyright (C) 2022  Institut Curie <http://www.curie.fr> \n"
            "\nThis program comes with ABSOLUTELY NO WARRANTY. \n"
            "This is free software, and you are welcome to redistribute it \n"
            "under the terms of the CeCILL License. \n"
	    	"\ncontact: Frederic Jarlier (frederic.jarlier@curie.fr) \n",
			progname, VERSION, progname);
			return 1; }
	/* Validate provided command (first argument) */
	if (strcmp(argv[1], "mem") != 0) {
		fprintf(stderr, "%s: unsupported %s command\n", progname, argv[1]);
		return 1; }
	else {
		 //create pg_line for create_sam_header
		 asprintf(&pg_line, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", VERSION, argv[0]);
		 int i = 0;
		 for (i = 1; i < argc; ++i) asprintf(&pg_line, "%s %s", pg_line, argv[i]);
	}

	/* initialize the BWA-MEM parameters to the default values */
	opt = mem_opt_init();
	memset(&opt0, 0, sizeof(opt0));
	while ((c = getopt(argc-1, argv+1, "bg51qpaMCSPVYjk:K:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:W:x:G:h:y:K:X:H:o:z:f")) >= 0) {
		if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
		else if (c == '1') ; /* FIXME: unsupported */
		else if (c == 'x') mode = optarg;
		else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;
		else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1;
		else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1;
		else if (c == 'b') write_format = 1;
        	else if (c == 'g') write_format = 0;
		else if (c == 'z') shared_mem = optarg;
		else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1;
		else if (c == 'U') opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1;
		else if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		else if (c == 'p') opt->flag |= MEM_F_PE | MEM_F_SMARTPE;
		else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
		else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
		else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
		else if (c == 'V') opt->flag |= MEM_F_REF_HDR;
		else if (c == '5') opt->flag |= MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ;
		else if (c == 'q') opt->flag |= MEM_F_KEEP_SUPP_MAPQ;
		else if (c == 'c') opt->max_occ = atoi(optarg), opt0.max_occ = 1;
		else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;
		else if (c == 'v') bwa_verbose = atoi(optarg);
		else if (c == 'j') ignore_alt = 1;
		else if (c == 'r') opt->split_factor = atof(optarg), opt0.split_factor = 1.;
		else if (c == 'D') opt->drop_ratio = atof(optarg), opt0.drop_ratio = 1.;
		else if (c == 'm') opt->max_matesw = atoi(optarg), opt0.max_matesw = 1;
		else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;
		else if (c == 'G') opt->max_chain_gap = atoi(optarg), opt0.max_chain_gap = 1;
		else if (c == 'N') opt->max_chain_extend = atoi(optarg), opt0.max_chain_extend = 1;
		else if (c == 'W') opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1;
		else if (c == 'y') opt->max_mem_intv = atol(optarg), opt0.max_mem_intv = 1;
		else if (c == 'C') copy_comment = 1;
		else if (c == 'K') fixed_chunk_size = atoi(optarg);
		else if (c == 'X') opt->mask_level = atof(optarg);
		else if (c == 'h') {
			opt0.max_XA_hits = opt0.max_XA_hits_alt = 1;
			opt->max_XA_hits = opt->max_XA_hits_alt = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->max_XA_hits_alt = strtol(p+1, &p, 10);
		}
		else if (c == 'Q') {
			opt0.mapQ_coef_len = 1;
			opt->mapQ_coef_len = atoi(optarg);
			opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;
		}
		else if (c == 'O') {
                        opt0.o_del = opt0.o_ins = 1;
                        opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
                        if (*p != 0 && ispunct(*p) && isdigit(p[1]))
                                opt->o_ins = strtol(p+1, &p, 10);
                } else if (c == 'E') {
                        opt0.e_del = opt0.e_ins = 1;
                        opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
                        if (*p != 0 && ispunct(*p) && isdigit(p[1]))
                                opt->e_ins = strtol(p+1, &p, 10);
                } else if (c == 'L') {
                        opt0.pen_clip5 = opt0.pen_clip3 = 1;
                        opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
                        if (*p != 0 && ispunct(*p) && isdigit(p[1]))
                                opt->pen_clip3 = strtol(p+1, &p, 10);
                }
		else if (c == 'R') {
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1;
		}
		else if (c == 'H') {
			if (optarg[0] != '@') {
				FILE *fp;
				if ((fp = fopen(optarg, "r")) != 0) {
					char *buf;
					buf = calloc(1, 0x10000);
					while (fgets(buf, 0xffff, fp)) {
						size_t i = strlen(buf);
						assert(buf[i-1] == '\n'); // a long line
						buf[i-1] = 0;
						hdr_line = bwa_insert_header(buf, hdr_line);
					}
					free(buf);
					fclose(fp);
				}
			} else hdr_line = bwa_insert_header(optarg, hdr_line);
		}
		else if (c == 'I') { // specify the insert size distribution
			pes0 = pes;
			pes[1].failed = 0;
			pes[1].avg = strtod(optarg, &p);
			pes[1].std = pes[1].avg * .1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].std = strtod(p+1, &p);
			pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);
			pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);
			if (pes[1].low < 1) pes[1].low = 1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].high = (int)(strtod(p+1, &p) + .499);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].low  = (int)(strtod(p+1, &p) + .499);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] mean insert size: %.3f, stddev: %.3f, max: %d, min: %d\n",
					__func__, pes[1].avg, pes[1].std, pes[1].high, pes[1].low);
		}
		/* Tool specific options */
		else if (c == 'o') file_out = optarg;
		else if (c == 'f') dofixmate = 1;
		else return 1; }

	if (mode) {
		if (strcmp(mode, "intractg") == 0) {
                        if (!opt0.o_del) opt->o_del = 16;
			if (!opt0.o_ins) opt->o_ins = 16;
			if (!opt0.b) opt->b = 9;
			if (!opt0.pen_clip5) opt->pen_clip5 = 5;
			if (!opt0.pen_clip3) opt->pen_clip3 = 5;
		} else if (strcmp(mode, "pacbio") == 0 || strcmp(mode, "pbref") == 0 || strcmp(mode, "ont2d") == 0) {
                        if (!opt0.o_del) opt->o_del = 1;
			if (!opt0.e_del) opt->e_del = 1;
			if (!opt0.o_ins) opt->o_ins = 1;
			if (!opt0.e_ins) opt->e_ins = 1;
			if (!opt0.b) opt->b = 1;
			if (opt0.split_factor == 0.) opt->split_factor = 10.;
			if (strcmp(mode, "ont2d") == 0) {
				if (!opt0.min_chain_weight) opt->min_chain_weight = 20;
				if (!opt0.min_seed_len) opt->min_seed_len = 14;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			} else {
				if (!opt0.min_chain_weight) opt->min_chain_weight = 40;
				if (!opt0.min_seed_len) opt->min_seed_len = 17;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			}
		} else {
			fprintf(stderr, "[E::%s] unknown read type '%s'\n", __func__, mode);
			return 1; // FIXME memory leak
		}
	} else {
		//update_a(opt, &opt0);
		if (opt0.a) { // matching score is changed
			if (!opt0.b) opt->b *= opt->a;
			if (!opt0.T) opt->T *= opt->a;
			if (!opt0.o_del) opt->o_del *= opt->a;
			if (!opt0.e_del) opt->e_del *= opt->a;
			if (!opt0.o_ins) opt->o_ins *= opt->a;
			if (!opt0.e_ins) opt->e_ins *= opt->a;
			if (!opt0.zdrop) opt->zdrop *= opt->a;
			if (!opt0.pen_clip5) opt->pen_clip5 *= opt->a;
			if (!opt0.pen_clip3) opt->pen_clip3 *= opt->a;
			if (!opt0.pen_unpaired) opt->pen_unpaired *= opt->a;
		}
	}

	nargs = argc - 1 - optind;
	if (opt->n_threads < 1) opt->n_threads = 1;
	if (nargs != 2 && nargs != 3) {
                fprintf(stderr, "\n");
		fprintf(stderr, "Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT        number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "       -k INT        minimum seed length [%d]\n", opt->min_seed_len);
		fprintf(stderr, "       -w INT        band width for banded alignment [%d]\n", opt->w);
		fprintf(stderr, "       -d INT        off-diagonal X-dropoff [%d]\n", opt->zdrop);
		fprintf(stderr, "       -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
		fprintf(stderr, "       -y INT        seed occurrence for the 3rd round seeding [%ld]\n", (long)opt->max_mem_intv);
		fprintf(stderr, "       -c INT        skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
		fprintf(stderr, "       -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [%.2f]\n", opt->drop_ratio);
		fprintf(stderr, "       -W INT        discard a chain if seeded bases shorter than INT [0]\n");
		fprintf(stderr, "       -m INT        perform at most INT rounds of mate rescues for each read [%d]\n", opt->max_matesw);
		fprintf(stderr, "       -S            skip mate rescue\n");
		fprintf(stderr, "       -P            skip pairing; mate rescue performed unless -S also in use\n");
		fprintf(stderr, "\nScoring options:\n\n");
		fprintf(stderr, "       -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [%d]\n", opt->a);
		fprintf(stderr, "       -B INT        penalty for a mismatch [%d]\n", opt->b);
		fprintf(stderr, "       -O INT[,INT]  gap open penalties for deletions and insertions [%d,%d]\n", opt->o_del, opt->o_ins);
		fprintf(stderr, "       -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [%d,%d]\n", opt->e_del, opt->e_ins);
		fprintf(stderr, "       -L INT[,INT]  penalty for 5'- and 3'-end clipping [%d,%d]\n", opt->pen_clip5, opt->pen_clip3);
		fprintf(stderr, "       -U INT        penalty for an unpaired read pair [%d]\n\n", opt->pen_unpaired);
		fprintf(stderr, "       -x STR        read type. Setting -x changes multiple parameters unless overriden [null]\n");
		fprintf(stderr, "                     pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)\n");
		fprintf(stderr, "                     ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)\n");
		fprintf(stderr, "                     intractg: -B9 -O16 -L5  (intra-species contigs to ref)\n");
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -p            smart pairing (ignoring in2.fq)\n");
		fprintf(stderr, "       -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n");
		fprintf(stderr, "       -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]\n");
		fprintf(stderr, "       -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "       -T INT        minimum score to output [%d]\n", opt->T);
		fprintf(stderr, "       -h INT[,INT]  if there are <INT hits with score >80%% of the max score, output all in XA [%d,%d]\n", opt->max_XA_hits, opt->max_XA_hits_alt);
		fprintf(stderr, "       -a            output all alignments for SE or unpaired PE\n");
		fprintf(stderr, "       -C            append FASTA/FASTQ comment to SAM output\n");
		fprintf(stderr, "       -V            output the reference FASTA header in the XR tag\n");
		fprintf(stderr, "       -Y            use soft clipping for supplementary alignments\n");
		fprintf(stderr, "       -M            mark shorter split hits as secondary\n\n");
		fprintf(stderr, "       -I FLOAT[,FLOAT[,INT[,INT]]]\n");
		fprintf(stderr, "                     specify the mean, standard deviation (10%% of the mean if absent), max\n");
		fprintf(stderr, "                     (4 sigma from the mean if absent) and min of the insert size distribution.\n");
		fprintf(stderr, "                     FR orientation only. [inferred]\n");
		fprintf(stderr, "\nExtra options:\n\n");
		fprintf(stderr, "       -K INT        process INT input bases in each batch regardless of nThreads (for reproducibility) []\n");
		fprintf(stderr, "       -o STR     the output file name\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "Note: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		free(opt);
		return 1;
	}
	bwa_fill_scmat(opt->a, opt->b, opt->mat);

	files = 0;

	file_ref = argv[optind+1+0];
	if (nargs > 1) {
		file_r1 = argv[optind+1+1]; files += 1;
	}
	if (nargs > 2) {
		file_r2 = argv[optind+1+2]; files += 1;
		opt->flag |= MEM_F_PE;
	}

	/* Derived file names */
	sprintf(file_map, "%s.map", file_ref);

	/* start up MPI */
	int threads_ok;
        int provided;
        res = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
        assert(res == MPI_SUCCESS);
        threads_ok = provided >= MPI_THREAD_FUNNELED;
        assert(res == MPI_SUCCESS);
	res = MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	assert(res == MPI_SUCCESS);
	res = MPI_Comm_rank(MPI_COMM_WORLD, &rank_num);
	assert(res == MPI_SUCCESS);
        if (rank_num == 0) {
		fprintf(stderr, "rank %d: is multithread ok ? %d\n", rank_num, threads_ok );
		fprintf(stderr, "all ranks will use %d threads per mpi job \n", opt->n_threads);
	}
	// some internal structures
	char *p1, *q1, *e1, *p2, *q2, *e2;
	int line_number;
	int NUM_THREADS = opt->n_threads;
	//threads for reading
	int NUM_THREADS_1 = 1;
	//threads for writing
	//int NUM_THREADS_2 = 1;

	int64_t bases;
	double before_local_mapping = 0;
	double after_local_mapping	= 0;
	double total_time_local_mapping  = 0;
	double total_time_mapping = 0;
	double total_time_reading_seq = 0;
	double grand_total_time_reading_seq = 0;
	double total_time_parsing = 0;
	double grand_total_time_parsing = 0;
	double total_time_writing = 0;
	double grand_total_time_writing = 0;
	size_t reads_r1, reads_r2, reads;	
	size_t offset_chunk;
	size_t offset_chunk_2;
	size_t size_chunk;
	size_t size_chunk_2;
	size_t total_local_reads_aligned= 0;
	size_t total_reads_check 		= 0;
	size_t grand_total_num_reads = 0;
	size_t *begin_offset_chunk 	= NULL;
	size_t *chunk_size 		        = NULL;
	size_t *reads_in_chunk 		= NULL;

	size_t *all_begin_offset_chunk  = NULL;
    	size_t *all_chunk_size          = NULL;
    	size_t *all_reads_in_chunk      = NULL;

    	size_t *all_begin_offset_chunk_2  = NULL;
    	size_t *all_chunk_size_2          = NULL;
    	size_t *all_reads_in_chunk_2      = NULL;

	
	pthread_t threads[NUM_THREADS];
        //pthread_t writter_threads[NUM_THREADS_2]; 
	
	//struct thread_data *td;//[NUM_THREADS];
	//td = malloc (NUM_THREADS * sizeof(struct thread_data));

	//MPI_Info finfo;
	//MPI_Info_create(&finfo);
	/*
	 * In this part you shall adjust the striping factor and unit according
	 * to the underlying filesystem.
	 * Harmless for other file system.
	 *
	 */
	//MPI_Info_set(finfo,"striping_factor", STRIPING_FACTOR);
	//MPI_Info_set(finfo,"striping_unit", STRIPING_UNIT); //2G striping
	//MPI_Info_set(finfo,"ind_rd_buffer_size", STRIPING_UNIT); //2gb buffer
	//MPI_Info_set(finfo,"romio_ds_read",DATA_SIEVING_READ);

	/*
	 * for collective reading and writing
	 * should be adapted too and tested according to the file system
	 * Harmless for other file system.
	 */
	//MPI_Info_set(finfo,"nb_proc", NB_PROC);
	//MPI_Info_set(finfo,"cb_nodes", CB_NODES);
	//MPI_Info_set(finfo,"cb_block_size", CB_BLOCK_SIZE);
	//MPI_Info_set(finfo,"cb_buffer_size", CB_BUFFER_SIZE);
	//MPI_Info_set(finfo, "collective_buffering", "true");
	//MPI_Info_set(finfo,"romio_cb_write","enable");
	//MPI_Info_set(finfo,"romio_cb_read","enable");


	/* Check that output file (-o) is not null ... */
	
	if (file_out == NULL) {
		fprintf(stderr, "missing mandatory output file (-o)\n");
		res = MPI_Finalize();
		assert(res == MPI_SUCCESS);
		exit(2);
	}
	
	/* here we derive the output path to construct a file by chr*/
	size_t f_out_sz = strlen(file_out);
	size_t mi = f_out_sz;
	char *k = file_out + f_out_sz;
	while ( (mi > 0) && (*k-- != '/')) mi--;
        //char output_path[FILENAME_MAX];

	if (mi == 0){
		fprintf(stderr, "We have a problem. the output file must be a path and a file name.Make sure it is a valid path. \n");
                res = MPI_Finalize();
                assert(res == MPI_SUCCESS);
                exit(2);

	}

	fixed_chunk_size = fixed_chunk_size > 0? fixed_chunk_size : opt->chunk_size * opt->n_threads;


	char *output_path = malloc (mi * sizeof(char) + 1);
	output_path[mi] = 0;
    	char *h = output_path;
    	memmove(h, file_out, mi);
    	assert(output_path);


	//if ( mi == 0) getcwd( output_path, FILENAME_MAX );
	/*else{ 
		output_path[mi] = 0;
		char *h = output_path;
		memmove(h, file_out, mi);
        	assert(output_path);
	}
	*/
	/* Check the map file is present otherwise send a message ... */
        if (stat(file_map, &stat_map) == -1) {
                fprintf(stderr, "There is a problem with the map file %s: %s. It is not present or you have not generate it with mpiBWAIdx \n", file_map, strerror(errno));
	        res = MPI_Finalize();
                assert(res == MPI_SUCCESS);
                exit(2);
        }


	/* Check R1 & R2 file sizes */
	if (file_r1 != NULL && stat(file_r1, &stat_r1) == -1) {
		fprintf(stderr, "%s: %s\n", file_r1, strerror(errno));
		res = MPI_Finalize();
		assert(res == MPI_SUCCESS);
		exit(2);
	}

	if (file_r2 != NULL && stat(file_r2, &stat_r2) == -1) {
		fprintf(stderr, "%s: %s\n", file_r2, strerror(errno));
		res = MPI_Finalize();
		assert(res == MPI_SUCCESS);
		exit(2);
	}

	 if (shared_mem != NULL){
                 if (OPEN_MPI) {
                        if ((strcmp(shared_mem, "numa") == 0 ||
                                strcmp(shared_mem, "l1") == 0 || strcmp(shared_mem, "l2") == 0 ||
                                        strcmp(shared_mem, "l3") == 0  || strcmp(shared_mem, "socket") == 0 || strcmp(shared_mem, "shared") == 0)){

                                        if (rank_num == 0)
                                                fprintf(stderr, "%s: shared memory is set to : %s \n", __func__, shared_mem);
                        }
                        else{
                                if (rank_num == 0){
                                        fprintf(stderr, "%s: shared memory options %s not recognized. \n", __func__, shared_mem);
                                        fprintf(stderr, "%s: available options are numa, l1, l2, l3, socket, shared(default) \n", __func__);
                                        res = MPI_Finalize();
                			assert(res == MPI_SUCCESS);
                			exit(2);

                                }
                        }
                }
                else{
			 if ((strcmp(shared_mem, "numa") == 0 ||
                                strcmp(shared_mem, "l1") == 0 || strcmp(shared_mem, "l2") == 0 ||
                                        strcmp(shared_mem, "l3") == 0  || strcmp(shared_mem, "socket") == 0 )){
                                
                                if (rank_num == 0){
                                        fprintf(stderr, "%s: shared memory options %s not recognized. \n", __func__, shared_mem);
                                        fprintf(stderr, "%s: you are not using openMPI only available options are shrared or core \n", __func__);
                                        res = MPI_Finalize();
                                        assert(res == MPI_SUCCESS);
                                        exit(2);
                                
                                }
                        }
                        if (strcmp(shared_mem, "core") == 0){
                                if (rank_num == 0){
                                        fprintf(stderr, "%s: shared memory is set to core \n", __func__);
                                        shared_mem = "core";
                                }
                        }

                        if (strcmp(shared_mem, "shared") == 0){
                                if (rank_num == 0){
                                        fprintf(stderr, "%s: shared memory is set to shared \n", __func__);
                                        shared_mem = "shared";
                                }
                        }
                }
        }
	else{
		fprintf(stderr, "%s: shared memory is set to core \n", __func__);
		shared_mem = "core";	
	} 

	if (rank_num == 0)
		fprintf(stderr, "%s: controls are done. Start analyzing fastqs it could take few minutes...\n", __func__);	


	 //we create a RMA window to hold index chunk
	MPI_Win win;
        uint64_t *index_chunk;
        uint64_t incr = 1;
        uint64_t u1 = 0;
        int rank_target=0;
	if ( proc_num > 1) 
        	MPI_Win_allocate(sizeof(uint64_t),  1, MPI_INFO_NULL, MPI_COMM_WORLD, &index_chunk, &win);
        
        else u1 = 0;

	if ( (file_r1 != NULL && file_r2 != NULL  && (stat_r1.st_size == stat_r2.st_size)))  {
	
		/*
	 	 * 
	     * We are in case we are paired and not trimmed
	     *
	 	 */

		/* Work around build warning in non timing case */
		aft = 0; aft++;
		bef = 0; bef++;
		/*
		 * Rank 0 estimate the size of a read
		 *
		 */
		size_t blen;
		off_t tmp_sz = 1024;
		if (rank_num == 0){
			// 512 Mo are
			int fd_tmp = open(file_r1, O_RDONLY, 0666);
			char *buffer = malloc(tmp_sz  + 1);
			buffer[tmp_sz] = '\0';
			size_t read_out = read(fd_tmp, buffer, tmp_sz);
			assert(read_out);
			assert(strlen(buffer) == tmp_sz);
			assert( *buffer == '@');
			/* Estimate sequence size (bases and bytes) */
			/* With this estimation we compute approximately the chuncks size and offset*/

			s = buffer;
			e = buffer + tmp_sz;
			p = q = s;
			while (q < e && *q != '\n') q++; 
			p = ++q;
			while (q < e && *q != '\n') q++; 
			blen = q - p; p = ++q;

			/* Split local buffer in chunks of 100000000 bases */
			free(buffer);
		}

		//Rank O broadcast the size of a read
	 	res = MPI_Bcast(&blen, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
		assert(res == MPI_SUCCESS);

		/*
		 * Split sequence files in chunks
		 */
		MPI_File mpi_fd_in1, mpi_fd_in2;
		bef = MPI_Wtime();

		res = MPI_File_open(MPI_COMM_WORLD, file_r1,  MPI_MODE_RDONLY , MPI_INFO_NULL, &mpi_fd_in1);
		assert(res == MPI_SUCCESS);
		
		res = MPI_File_open(MPI_COMM_WORLD, file_r1,  MPI_MODE_RDONLY , MPI_INFO_NULL, &mpi_fd_in2);
		assert(res == MPI_SUCCESS);
		
		/*
		 * first we parse the buffer and see
		 * how many reads we have
		 */

		assert(fd_in1 != -1);
		size_t *goff = NULL; //global offset contain the start offset in the fastq
		goff = malloc((proc_num * NUM_THREADS_1 + 1) * sizeof(size_t));
		size_t *goff_inter = calloc( (proc_num * NUM_THREADS_1 + 1) , sizeof(size_t));
	
		//we shall call 
		bef = MPI_Wtime();
                find_process_starting_offset_mt(goff, stat_r1.st_size, file_r1, proc_num, rank_num, NUM_THREADS_1);
                aft = MPI_Wtime();
		if (rank_num == 0)
                	fprintf(stderr, "%s: rank %d time spend in finding process start offset = (%.02f) \n", __func__, rank_num, aft - bef);

                int i12=0;
                for ( i12 = 0; i12 < proc_num * NUM_THREADS_1 + 1; i12++ )  goff_inter[i12] = goff[i12];
	
		//now we exchange the goff buffer between all proc
		//rank 0 gather the vector
		res = MPI_Allgather(&goff_inter[rank_num*NUM_THREADS_1], NUM_THREADS_1, MPI_LONG_LONG_INT, goff , NUM_THREADS_1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
		assert(res == MPI_SUCCESS);

		free(goff_inter);

		//we compute the new size according to the shift
		//We calculate the size to read for each process
		MPI_Barrier(MPI_COMM_WORLD);
		
		size_t local_num_reads          = 0;		
		size_t total_num_reads          = 0;
		size_t *local_read_offsets 	= NULL;
		size_t *local_read_bytes    	= NULL;
		int *local_read_size    	= NULL;
	
		bef = MPI_Wtime();

		pthread_attr_t attr;
                pthread_attr_init(&attr);
                //pthread_attr_setstacksize(&attr, SMALL_STACK);
                //pthread_attr_setdetachstate(&attr, 0);

                pthread_t threads_1[NUM_THREADS_1];

                struct struct_data_thread_1 *td_1 =  malloc(NUM_THREADS_1 * sizeof(struct struct_data_thread_1));

                size_t *local_num_reads_t            = calloc(NUM_THREADS_1, sizeof(size_t));
                size_t *total_num_reads_t            = calloc(NUM_THREADS_1, sizeof(size_t));
                size_t **local_read_offsets_t        = calloc(NUM_THREADS_1, sizeof(size_t*));
                size_t **local_read_bytes_t          = calloc(NUM_THREADS_1, sizeof(size_t*));
                int **local_read_size_t              = calloc(NUM_THREADS_1, sizeof(int*));

                
                int goff_idx = 0;
                for ( n = 0; n < NUM_THREADS_1; n++){

                        goff_idx = (rank_num * NUM_THREADS_1) + n;
                        td_1[n].offset_in_file_mt         = goff[goff_idx];
                        td_1[n].size2read_mt              = goff[goff_idx + 1] - goff[goff_idx];
                        td_1[n].file_r1_mt                = file_r1;
                        td_1[n].local_num_reads_mt        = &local_num_reads_t[n];
                        td_1[n].total_num_reads_mt        = &total_num_reads_t[n];
                        td_1[n].local_read_offsets_mt     = &local_read_offsets_t[n];
                        td_1[n].local_read_size_mt        = &local_read_size_t[n];
                        td_1[n].local_read_bytes_mt       = &local_read_bytes_t[n];
                        td_1[n].proc_num_mt               = proc_num;
                        td_1[n].rank_num_mt               = rank_num;
                        td_1[n].thread_num_mt             = n;
                        td_1[n].previous_read_num         = 0;
                        pthread_create(&threads_1[n], &attr, find_reads_size_and_offsets_mt, (void *)(&td_1[n]));

                }

                total_num_reads = 0;
                for (n = 0; n < NUM_THREADS_1; n++){
                        pthread_join(threads_1[n], (void *)(&td_1[n]));
                        total_num_reads += *(td_1[n].total_num_reads_mt);
                }

		local_read_offsets  = calloc(total_num_reads, sizeof(size_t));
                local_read_size     = calloc(total_num_reads, sizeof(int));
                local_read_bytes    = calloc(total_num_reads, sizeof(size_t));

                assert(local_read_offsets);
                assert(local_read_size);
                assert(local_read_bytes);

                size_t tmp_var = 0;
                for (n = 0; n < NUM_THREADS_1; n++){
                        td_1[n].local_read_offsets     = local_read_offsets;
                        td_1[n].local_read_size        = local_read_size;
                        td_1[n].local_read_bytes       = local_read_bytes;
                        td_1[n].previous_read_num      = tmp_var;
                        tmp_var                        += *(td_1[n].total_num_reads_mt);
		
			pthread_create(&threads_1[n], &attr, copy_local_read_info_mt, (void *)(&td_1[n]));

                }

                for (n = 0; n < NUM_THREADS_1; n++){
                        pthread_join(threads_1[n], (void *)(&td_1[n]));

                        free(local_read_offsets_t[n]);
                        free(local_read_bytes_t[n]);
                        free(local_read_size_t[n]);
                }

                free(local_num_reads_t);
                free(total_num_reads_t);

                free(local_read_offsets_t);
                free(local_read_bytes_t);
                free(local_read_size_t);

                pthread_attr_destroy(&attr);
                free(td_1);
		/*
		find_reads_size_and_offsets(goff[ind],
                                                siz2read,
                                                file_r1,
                                                &local_num_reads,
                                                &total_num_reads,
                                                &local_read_offsets,
                                                &local_read_size,
                                                &local_read_bytes,
                                                proc_num,
                                                rank_num);

		*/

		MPI_Barrier(MPI_COMM_WORLD);
		aft = MPI_Wtime();
		if ( rank_num == 0 )
		fprintf(stderr, "%s: rank %d num reads parsed: %zu ::: time spend reading and parsing entire buffer = (%.02f) \n", __func__, rank_num, total_num_reads, aft - bef);
			
		if (goff) free(goff);
	
		MPI_File_close(&mpi_fd_in2);			
		res = MPI_Reduce(&total_num_reads, &grand_total_num_reads, 1,MPI_LONG_LONG_INT, MPI_SUM, 0,MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);	
		res = MPI_Bcast(&grand_total_num_reads, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);	

		size_t num_reads_by_proc[proc_num];
		res = MPI_Allgather(&total_num_reads, 1, MPI_LONG_LONG_INT, num_reads_by_proc, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
		//we dispatch the previous result
	
		if ( rank_num == 0 || rank_num ==1 )
			fprintf(stderr, "rank %d ::: total_num_reads = %zu \n", rank_num, grand_total_num_reads);

		/*
		 * Now each rank compute buffers offset start and end.
		 */

		local_num_reads = total_num_reads;
		//now we estimate the number of chunk per rank
		size_t chunck_num = (local_num_reads * blen) / ( fixed_chunk_size / 2);
		chunck_num += 2; //the last chunk hold the remain bases

		size_t h=0;
		for ( h = 0; h < total_num_reads; h++){
			assert( local_read_size[h] == blen );
			assert( local_read_offsets[h] >= 0);
		}

		// we allocate vector for chunks offset
		begin_offset_chunk 	= calloc(chunck_num, sizeof(size_t));
		chunk_size      	= calloc(chunck_num, sizeof(size_t));
		reads_in_chunk 		= calloc(chunck_num, sizeof(size_t));

		assert( begin_offset_chunk != NULL );
		assert( chunk_size != NULL );
		assert( reads_in_chunk != NULL );
	
		size_t chunk_count = 0;
		//fprintf(stderr,"rank %d ::: chunck size = %zu \n", rank_num, chunk_size);
		maxsiz = ( fixed_chunk_size ) / 2; 
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank_num == 0)
			fprintf(stderr,"rank %d ::: Call find_chunks_info \n", rank_num);
		// the detail of he paramters is at the function definition
		// fprintf(stderr, "rank %d ::: begin_offset_chunk = %zu \n", rank_num, begin_offset_chunk[0]);
		// fprintf(stderr, "rank %d ::: begin_offset_chunk_2 = %zu \n", rank_num, begin_offset_chunk_2[0]);
		// fprintf(stderr, "rank %d ::: chunk_size = %zu \n", rank_num, chunk_size);
		// fprintf(stderr, "rank %d ::: chunk_size_2 = %zu \n", rank_num, chunk_size_2);

		bef = MPI_Wtime();
		find_chunks_info(begin_offset_chunk,
				 chunk_size,
				 reads_in_chunk,
				 local_read_size,
				 local_read_bytes,
				 local_read_offsets,
				 rank_num,
				 proc_num,
				 local_num_reads,
				 grand_total_num_reads,
				 maxsiz,
				 &chunk_count
				 );

		aft = MPI_Wtime();
		fprintf(stderr, "%s ::: rank %d ::: evaluating offsets chuncks and sizes (%.02f) found %zu chuncks \n", __func__, rank_num, aft - bef, chunk_count);
	
		free(local_read_offsets);
		free(local_read_size);
		free(local_read_bytes);	

		bef = MPI_Wtime();
		size_t total_chunks = 0;
                res = MPI_Reduce(&chunk_count, &total_chunks, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
                assert ( res == MPI_SUCCESS);
                res = MPI_Bcast(&total_chunks, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
                assert ( res == MPI_SUCCESS);

                size_t chunks_per_rank[proc_num];
                res = MPI_Allgather(&chunk_count, 1, MPI_LONG_LONG_INT, chunks_per_rank, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);

		all_begin_offset_chunk = malloc(total_chunks * sizeof(size_t));
                all_begin_offset_chunk[0]=0;
                all_chunk_size = malloc(total_chunks * sizeof(size_t));
                all_chunk_size[0] = 0;
                all_reads_in_chunk =  malloc(total_chunks * sizeof(size_t));
                all_reads_in_chunk[0] = 0;

                int displ_chunk[proc_num];
                displ_chunk[0] = 0;
		int i = 0;
                for (i = 1; i < proc_num; i++) displ_chunk[i] = (displ_chunk[i-1] + chunks_per_rank[i-1]);

                int indx=displ_chunk[rank_num];
                for (i = 0; i <  chunks_per_rank[rank_num]; i++) {
                        all_chunk_size[indx+i]=chunk_size[i];
                        all_begin_offset_chunk[indx+i]=begin_offset_chunk[i];
                        all_reads_in_chunk[indx+i]=reads_in_chunk[i];
                }

                if (rank_num > 0){
                res=MPI_Send(chunk_size, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
                assert(res == MPI_SUCCESS);
                }
                else{
                for (i = 1; i < proc_num; i++){
                        res=MPI_Recv(&(all_chunk_size[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        assert(res == MPI_SUCCESS);
                }
                }

                if (rank_num > 0){
                        res=MPI_Send(begin_offset_chunk, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
                        assert(res == MPI_SUCCESS);
                }
                else{
                for (i = 1; i < proc_num; i++){
                        res=MPI_Recv(&(all_begin_offset_chunk[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        assert(res == MPI_SUCCESS);
                        }
                }

		if (rank_num > 0){
                        res=MPI_Send(reads_in_chunk, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
                        assert(res == MPI_SUCCESS);
                }
                else{
                for (i = 1; i < proc_num; i++){
                        res=MPI_Recv(&(all_reads_in_chunk[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        assert(res == MPI_SUCCESS);
                        }
                }

                res=MPI_Bcast(all_chunk_size, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
                res=MPI_Bcast(all_reads_in_chunk, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
                res=MPI_Bcast(all_begin_offset_chunk, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

                assert(res == MPI_SUCCESS);

                for ( i = 0; i < total_chunks; i++)
                assert(all_chunk_size[i] != 0);

                free(chunk_size);
                free(begin_offset_chunk);
                free(reads_in_chunk);

		aft = MPI_Wtime();
		if ( rank_num == 0 )
                	fprintf(stderr, "%s: rank %d time spend gathering chunks info = (%.02f) \n", __func__, rank_num, aft - bef);

		/*
		 * Map reference genome indexes in shared memory (by host)
		 */
		bef = MPI_Wtime();
                map_indexes(file_map, &count, &indix, &ignore_alt, &win_shr, shared_mem);
                aft = MPI_Wtime();
		if ( rank_num == 0 )
			fprintf(stderr, "%s: mapped indexes (%.02f)\n", __func__, aft - bef);

		//we create a vector with chromosom names 
		int s;
		int incrmnt = 2;
                if (dofixmate) incrmnt = 1;

		MPI_File *fh_out          = malloc( (indix.bns->n_seqs + incrmnt) * sizeof(MPI_File));
		char *files_out_sam_name[indix.bns->n_seqs + incrmnt];
		char *file_map_by_chr[(indix.bns->n_seqs + incrmnt)];

		for (s = 0; s < indix.bns->n_seqs; ++s){
			files_out_sam_name[s] = malloc( (strlen( indix.bns->anns[s].name) + 1 ) * sizeof(char));
		        files_out_sam_name[s][strlen( indix.bns->anns[s].name)] = 0;
		        char *p0 = files_out_sam_name[s];
		        memmove(p0, indix.bns->anns[s].name, strlen( indix.bns->anns[s].name));
		}
		
		if (!dofixmate){
                        char DISCORDANT[] = "discordant";
                        files_out_sam_name[s++] = strdup(DISCORDANT);
                }

		char UNMAPPED[]   = "unmapped";
		files_out_sam_name[s++] = strdup(UNMAPPED);
		int file_name_len = 0;
		
		for (s = 0; s < (indix.bns->n_seqs + incrmnt); ++s){
			/* Derived file names */
		        file_name_len = strlen(output_path) + strlen(files_out_sam_name[s]) + 6;
		        file_map_by_chr[s] = calloc( file_name_len, sizeof(char));
			if ( write_format == 2	)	        
		        	sprintf(file_map_by_chr[s], "%s/%s.sam", output_path, files_out_sam_name[s]);
			if ( write_format == 1  )
                                sprintf(file_map_by_chr[s], "%s/%s.bam", output_path, files_out_sam_name[s]);
			if ( write_format == 0  )
                                sprintf(file_map_by_chr[s], "%s/%s.gz", output_path, files_out_sam_name[s]);


		}
		
		if ( write_format == 2  )
			create_sam_header_by_chr_file(file_map_by_chr, &indix, &count, hdr_line, rg_line, pg_line, rank_num);
		
		if ( write_format == 1 || write_format == 0)
			create_bam_header_by_chr_file(file_map_by_chr, &indix, &count, hdr_line, rg_line, pg_line, rank_num, compression_level, dofixmate);

		//we create a file which contains only the header
		//this file will be used by the sorting to get individual chromosom file
		//create_sam_header(file_out, &indix, &count, hdr_line, rg_line, rank_num);
		
		//we add the header in the discordant SAM
		//if ( write_format == 2  )
		//	create_sam_header(file_map_by_chr[indix.bns->n_seqs], &indix, &count, hdr_line, rg_line, rank_num);
                //if ( write_format == 1 || write_format == 0)
		//	 create_bam_header(file_map_by_chr[indix.bns->n_seqs], &indix, &count, hdr_line, rg_line, rank_num, compression_level);
		

		for (s = 0; s < (indix.bns->n_seqs + incrmnt); ++s){
                        res = MPI_File_open(MPI_COMM_WORLD, file_map_by_chr[s], MPI_MODE_CREATE|MPI_MODE_WRONLY|MPI_MODE_APPEND, MPI_INFO_NULL, &fh_out[s]);
                        assert(res == MPI_SUCCESS);
                }

		//bef = MPI_Wtime();
		res = MPI_Barrier(MPI_COMM_WORLD);
		assert(res == MPI_SUCCESS);
		//aft = MPI_Wtime();
		//xfprintf(stderr, "%s: synched processes (%.02f)\n", __func__, aft - bef);

		if (file_r1 != NULL) {
			//fprintf(stderr, "rank %d ::: open file %s \n",rank_num, file_r1);
			res = MPI_File_open(MPI_COMM_WORLD, file_r1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_r1);
			assert(res == MPI_SUCCESS);
		}
		if (file_r2 != NULL) {
			//fprintf(stderr, "rank %d ::: open file %s \n",rank_num, file_r2);
			res = MPI_File_open(MPI_COMM_WORLD, file_r2, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_r2);
			assert(res == MPI_SUCCESS);
		}

		buffer_r1 = buffer_r2 = NULL; seqs = NULL;
		//localsize_vec contain the size of the buffer to write in the sam file
		int *chr_buff_size  = calloc ( (indix.bns->n_seqs + incrmnt), sizeof(int) );
		char *buffer_out_vec[indix.bns->n_seqs + incrmnt];
		
		bef = MPI_Wtime();

		if ( proc_num > 1) {

                	MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank_target, 0, win);
                	MPI_Request req;
                	MPI_Rget(&index_chunk, 1, MPI_UINT64_T, rank_target, 0, 1, MPI_UINT64_T, win, &req);
                	MPI_Wait(&req, MPI_STATUS_IGNORE);
                	u1 = (uint64_t)index_chunk;
                	MPI_Fetch_and_op(&incr, &index_chunk, MPI_UINT64_T, rank_target, 0, MPI_SUM, win);
                	MPI_Win_flush(rank_target, win);
                	MPI_Win_unlock(rank_target, win);

                	aft = MPI_Wtime();
                	fprintf(stderr, "rank %d ::: initial u1 = %zu :: time %.02f \n",rank_num, u1, aft-bef);
		}
		// here we loop until there's nothing to read
		// in the offset and size file
		
		before_local_mapping = MPI_Wtime();

		//we loop the chunck_count
		
		while ( u1 < total_chunks ){

			offset_chunk = all_begin_offset_chunk[u1];
			size_chunk   = all_chunk_size[u1];
			assert(size_chunk != 0);
			/*
			 * Read sequence datas ...
			 *
			 */
			bef = MPI_Wtime();

			buffer_r1 = malloc(size_chunk+1);
			assert(buffer_r1 != NULL);
			buffer_r1[size_chunk] = 0;

			buffer_r2 = malloc(size_chunk+1);
                        assert(buffer_r2 != NULL);
                        buffer_r2[size_chunk]=0;

			struct struct_pread_fastq *td_pread1;
                        td_pread1 = malloc (NUM_THREADS_1 * sizeof(struct struct_pread_fastq));
                        bef = MPI_Wtime();
                        pthread_attr_t attr4;
                        pthread_attr_init(&attr4);
                        pthread_attr_setstacksize(&attr4, BIG_STACK);
                        pthread_attr_setdetachstate(&attr4, 0);

                        for( n = 0; n < NUM_THREADS_1; n++ ){
                                td_pread1[n].total_thread = NUM_THREADS_1;
                                td_pread1[n].thread_id = n;
                                td_pread1[n].job_rank = rank_num;
                                td_pread1[n].offset= offset_chunk;
                                td_pread1[n].size = size_chunk;
                                td_pread1[n].buffer = buffer_r1;
                                td_pread1[n].fd  = fh_r1;
                                int ret_code = pthread_create(&threads[n], &attr4, pread_fastq_chunck, (void *)(&td_pread1[n]));
                                if (ret_code) {
                                        fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                                }
                        }
                        for(n=0; n<NUM_THREADS_1; n++)
                        pthread_join(threads[n], (void *)(&td_pread1[n]));

                        pthread_attr_destroy(&attr4);
                        free(td_pread1);

                        assert(strlen(buffer_r1) == size_chunk);
                        assert(*buffer_r1 == '@');

			struct struct_pread_fastq *td_pread2;
                        td_pread2 = malloc (NUM_THREADS_1 * sizeof(struct struct_pread_fastq));
                        bef = MPI_Wtime();
                        pthread_attr_t attr5;
                        pthread_attr_init(&attr5);
                        pthread_attr_setstacksize(&attr5, BIG_STACK);
                        pthread_attr_setdetachstate(&attr5, 0);

                        for( n = 0; n < NUM_THREADS_1; n++ ){
                                td_pread2[n].total_thread = NUM_THREADS_1;
                                td_pread2[n].thread_id = n;
                                td_pread2[n].job_rank = rank_num;
                                td_pread2[n].offset= offset_chunk;
                                td_pread2[n].size = size_chunk;
                                td_pread2[n].buffer = buffer_r2;
                                td_pread2[n].fd  = fh_r2;

                                int ret_code = pthread_create(&threads[n], &attr5, pread_fastq_chunck, (void *)(&td_pread2[n]));
                                if (ret_code) {
                                        fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                                }
                        }
                        for(n=0; n<NUM_THREADS_1; n++)
                                pthread_join(threads[n], (void *)(&td_pread2[n]));

                        pthread_attr_destroy(&attr5);
                        free(td_pread2);

                        assert(strlen(buffer_r2) == size_chunk);
                        assert(*buffer_r2 == '@');

			if (u1 < total_chunks){
                                MPI_File_seek(fh_r1, (MPI_Offset)all_begin_offset_chunk[u1 + proc_num], MPI_SEEK_CUR ) ;
                                MPI_File_seek(fh_r2, (MPI_Offset)all_begin_offset_chunk[u1 + proc_num], MPI_SEEK_CUR ) ;
                        }
		

			aft = MPI_Wtime();
			fprintf(stderr, "%s: read sequences (%.02f)\n", __func__, aft - bef);
			total_time_reading_seq += (aft - bef);

			bef = MPI_Wtime();
			reads_r1 = all_reads_in_chunk[u1];
			
			if (file_r2 != NULL) reads_r2 = all_reads_in_chunk[u1];

			reads = reads_r1 + reads_r2; bases = 0;
			total_local_reads_aligned += reads/2;
			assert(reads <= INT_MAX);
			fprintf(stderr, "%s: num_rank = %d :: number of reads = %zu \n", __func__, rank_num, reads);

			/* Parse sequences ... */
			seqs = malloc(reads * sizeof(*seqs));
			assert(seqs != NULL);

			//case one we are paired
			if (file_r1 != NULL && file_r2 !=NULL) {
				p1 = q1 = buffer_r1; e1 = buffer_r1 + size_chunk; line_number = 0;
				p2 = q2 = buffer_r2; e2 = buffer_r2 + size_chunk; 
				
				while (q1 < e1) {
					if (*q1 != '\n') { q1++; q2++; continue; }
					/* We have a full line ... process it */
					*q1 = '\0'; 
					*q2 = '\0';
					n = files * (line_number / 4);
				
					switch (line_number % 4) {
					case 0: /* Line1: Name and Comment */
						assert(*p1 == '@');
						seqs[n].name   = p1 + 1;
						seqs[n+1].name = p2 + 1;
					
						while (*p1 && !isspace((unsigned char)*p1)) {p1++;p2++;}
						if (*(p1-2) == '/' && isdigit((unsigned char)*(p1-1))) {*(p1-2) = '\0'; *(p2-2) = '\0';}
						if (*p1) {*p1++ = '\0'; *p2++ = '\0';}
						seqs[n].comment = (copy_comment != 0) ? p1 : NULL;
						seqs[n].sam = NULL;
						seqs[n+1].comment = (copy_comment != 0) ? p2 : NULL;
						seqs[n+1].sam = NULL;					
						break;
					case 1: /* Line2: Sequence */
						seqs[n].seq = p1;
						seqs[n].l_seq = q1 - p1;
						seqs[n+1].seq = p2;
						seqs[n+1].l_seq = q2 - p2;

						bases += seqs[n].l_seq;
						bases += seqs[n+1].l_seq;
						break;
					case 2: /* Line3: Ignored */
						assert(*p1 == '+');
						break;
					case 3: /* Line4: Quality */
						seqs[n].qual = p1;
						seqs[n+1].qual = p2;
						break; }
					p1 = ++q1; 
					p2 = ++q2; 
					line_number++; 
				}
				//fprintf(stderr, "rank %d ::: Parse %d lines \n",rank_num, line_number );
			}
		
			aft = MPI_Wtime();
			fprintf(stderr, "%s: parsed sequences (%.02f)\n", __func__, aft - bef);
			total_time_parsing += (aft - bef);
			fprintf(stderr, "rank %d ::: [M::%s] read %zu sequences (%ld bp)...\n",rank_num , __func__, reads, (long)bases);
			/* Datas computation ... */
			bef = MPI_Wtime();
			//fprintf(stderr, "rank %d ::::  Call memseqs with count sequences %zu \n", rank_num, reads);

			mem_process_seqs(opt, indix.bwt, indix.bns, indix.pac, 0, (int)reads, seqs, pes0);

			aft = MPI_Wtime();
			fprintf(stderr, "%s: computed mappings (%.02f)\n", __func__, aft - bef);

			//MPI_Barrier(MPI_COMM_WORLD);
			/* Write results ... */
			bef = MPI_Wtime();
            	
			//pointer to the sam line
			char *currentCarac;
                        char *start_sam_line;
                        int chr, mchr, i;
                        size_t sam_line_size;
                        
                        int nbchr = indix.bns->n_seqs + incrmnt;

                        char *tmp_chr = malloc( MAX_CHR_NAME_SIZE * sizeof(char));
                        tmp_chr[0] = 0;

                        int next;
                        char currentLine[MAX_CHAR_SIZE];

			//first we cont how many sam line we have

			size_t total_sam_line = 0;
                        n = 0;
			int ret_code = 0;
                        if (dofixmate){

                                struct thread_data *td;
                                td = malloc (NUM_THREADS * sizeof(struct thread_data));
                                bef = MPI_Wtime();
                                pthread_attr_init(&attr);
                                pthread_attr_setstacksize(&attr, SMALL_STACK);
                                pthread_attr_setdetachstate(&attr, 0);

                                for( n = 0; n < NUM_THREADS; n++ ){
                                        td[n].total_thread = NUM_THREADS;
                                        td[n].thread_id = n;
                                        assert(seqs);
                                        td[n].seqs_thr = seqs;
                                        td[n].job_rank = rank_num;
                                        td[n].total_reads = reads;
                                        td[n].start_index_seqs = 0;
                                        td[n].final_index_seqs = reads-1;
                                        td[n].indix_thr = &indix;
                                        td[n].total_lines = 0;
                                        ret_code = pthread_create(&threads[n], &attr, call_fixmate, (void *)(&td[n]));
                                        if (ret_code) {
                                                fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                                        }
                                }
                                for(n=0; n<NUM_THREADS; n++)
                                        pthread_join(threads[n], (void *)(&td[n]));
                                for(n=0; n<NUM_THREADS; n++)
                                        total_sam_line += td[n].total_lines;
                                pthread_attr_destroy(&attr);
                                free(td);
                                
                        }
			else{
                                
                                n = 0;                      
                                while ( n < reads){
                                        next = tokenizer(seqs[n].sam,'\n', currentLine);
                                        while (next) {
                                                total_sam_line++;
                                                next = tokenizer(NULL, '\n', currentLine);
                                        }
                                        for(i=0; i<MAX_CHAR_SIZE; i++) currentLine[i]=0;
                                        next = tokenizer(seqs[n+1].sam,'\n', currentLine);
                                        while (next) {
                                                total_sam_line++;
                                                next = tokenizer(NULL, '\n', currentLine);
                                        }
                                        for(i=0; i<MAX_CHAR_SIZE; i++) currentLine[i]=0;
                                        n = n + 2;

                                }

                        }

                        if (dofixmate){
                                 aft = MPI_Wtime();
                                 fprintf(stderr, "rank: %d :: %s: time spend in fixmate (%.02f) \n", rank_num, __func__, aft - bef);
                        }

			int *sam_buff_dest      = calloc ( total_sam_line, sizeof(int) );
                        int *add_in_disc        = calloc ( total_sam_line, sizeof(int) );
			char **start_addr       = malloc ( total_sam_line * sizeof(char*));
                        int *line_size_to_cpy   = malloc ( total_sam_line * sizeof(int));
                        size_t incr_line =0;

                        for(i=0; i< MAX_CHAR_SIZE; i++){
                                currentLine[i]=0;
                        }
			if (dofixmate){
                        	for (n = 0; n < reads; n++) {
                                	int total_sam_line_size = 0;
                                	//seqs[n].l_seq = strlen(seqs[n].sam);
                                	next = tokenizer(seqs[n].sam,'\n', currentLine);
                                	currentCarac = currentLine;
                                	start_sam_line = seqs[n].sam;
					while (next){

        	                                start_sam_line = seqs[n].sam + total_sam_line_size;
                	                        sam_line_size = strlen(currentLine) + 1;
                        	                //we update the offset in the
                                	        //source file
                                        	currentLine[sam_line_size - 1] = '\n';
                                        	currentLine[sam_line_size] = '\0';
                                        	currentCarac = currentLine;
                                        	//WE PASS READ NAME
                                        	currentCarac = strstr (currentCarac, "\t");
                                        	currentCarac++;

                                        	//WE PASS FLAG  
						currentCarac = strstr (currentCarac, "\t");
                                        	currentCarac++;

                                        	if ( currentCarac[0] == '*') chr =  nbchr-1;
                                        	else chr = getChr(currentCarac - 1, files_out_sam_name, nbchr-1, tmp_chr);
      

                                        	if (chr < (nbchr - incrmnt)){
							chr_buff_size[chr]              += sam_line_size;
                                                	sam_buff_dest[incr_line]        = chr;
                                                	start_addr[incr_line]           = start_sam_line;

                                        	}
                                        	else {
							chr_buff_size[nbchr - incrmnt]  += sam_line_size;
                                                	sam_buff_dest[incr_line]        = nbchr-1;
                                                	start_addr[incr_line]           = start_sam_line;
                                        	}

                                        	line_size_to_cpy[incr_line] = sam_line_size;
                                        	total_sam_line_size += sam_line_size;
                                        	incr_line++;
                                        	next = tokenizer(NULL, '\n', currentLine);
                                	}
                        	}
			}
			else{
				 for (n = 0; n < reads; n++) {
                                        int total_sam_line_size = 0;
                                        next = tokenizer(seqs[n].sam,'\n', currentLine);
                                        currentCarac = currentLine;
                                        start_sam_line = seqs[n].sam;

                                        while (next){

                                                start_sam_line = seqs[n].sam + total_sam_line_size;
                                                sam_line_size = strlen(currentLine) + 1;
						currentLine[sam_line_size - 1] = '\n';
                                                currentLine[sam_line_size] = '\0';
                                                currentCarac = currentLine;
                                                //READ NAME
                                                currentCarac = strstr (currentCarac, "\t");
                                                currentCarac++;
                                                // FLAG
                                                currentCarac = strstr (currentCarac, "\t");
                                                currentCarac++;
                                                // CHROMOSOM 
						if ( currentCarac[0] == '*') chr =  nbchr - incrmnt;
                                                else chr = getChr(currentCarac - 1, files_out_sam_name, nbchr - incrmnt, tmp_chr);

                                                // COORD
                                                currentCarac = strstr (currentCarac + 1, "\t");
                                                
                                                // MAPQ
                                                strtoull(currentCarac, &currentCarac, 10);

                                                strtoull(currentCarac, &currentCarac, 10);                                               

						currentCarac = strstr (currentCarac + 1, "\t");
                                                currentCarac++;

                                                if ( currentCarac[0] == '=') mchr =  chr;
                                                else if ( currentCarac[0] == '*') mchr =  nbchr - 2;
                                                else mchr = getChr(currentCarac - 1, files_out_sam_name, nbchr - incrmnt, tmp_chr);

                                                if (chr < (nbchr - incrmnt)){
                                                        chr_buff_size[chr]              += sam_line_size;
                                                        sam_buff_dest[incr_line]        = chr;
                                                        start_addr[incr_line]           = start_sam_line;

                                                }
                                                else {
                                                        chr_buff_size[nbchr - 1] 	+= sam_line_size;
                                                        sam_buff_dest[incr_line]        = nbchr-1;
                                                        start_addr[incr_line]           = start_sam_line;
                                                }

                                                if ( (chr < (nbchr - incrmnt)) && (mchr < (nbchr - incrmnt)) && (chr != mchr) ){
                                                        chr_buff_size[nbchr - 2]        += sam_line_size;
                                                        add_in_disc[incr_line]           = 1;
                                                        
                                                }


                                                line_size_to_cpy[incr_line] = sam_line_size;
                                                total_sam_line_size += sam_line_size;
                                                incr_line++;
						next = tokenizer(NULL, '\n', currentLine);
					}
				}
			}

			free(tmp_chr);

			if ( write_format == 2 ){
				bef = MPI_Wtime();
				//now we fill up the buffer_out_vec
				for (n = 0; n < (indix.bns->n_seqs + incrmnt); n++) {

                                	assert(chr_buff_size[n] <= INT_MAX);

                                	if (chr_buff_size[n]){
                                        	buffer_out_vec[n] = calloc( chr_buff_size[n] + 1, sizeof(char));
                                        	assert(buffer_out_vec[n] != NULL);
                                        	buffer_out_vec[n][chr_buff_size[n]] = '\0';
                                	}
                        	}
                        	size_t *actual_size = calloc(indix.bns->n_seqs + incrmnt, sizeof(size_t));
                        	char *p_temp2;
                        	for (n = 0; n < total_sam_line; n++) {

                                	p_temp2 = buffer_out_vec[sam_buff_dest[n]] + actual_size[sam_buff_dest[n]];
                                	memmove(p_temp2, start_addr[n], line_size_to_cpy[n]);
                                	actual_size[sam_buff_dest[n]] += line_size_to_cpy[n];
					if ( add_in_disc[n] ){
                                	        p_temp2 = buffer_out_vec[nbchr - 2 ] + actual_size[nbchr - 2 ];
                                        	memmove(p_temp2, start_addr[n], line_size_to_cpy[n]);
                                        	actual_size[nbchr -2] += line_size_to_cpy[n];
                                	}
                        	}
			
                        	for (n = 0; n < reads; n++) free(seqs[n].sam);
				free(add_in_disc);
                        	free(sam_buff_dest);
				free(seqs);
                        	free(actual_size);
                        	free(start_addr);
                        	free(line_size_to_cpy);
				for (n = 0; n < (indix.bns->n_seqs + incrmnt); n++) {
					if ( chr_buff_size[n] ) {
						res = MPI_File_write_shared(fh_out[n], buffer_out_vec[n], chr_buff_size[n], MPI_CHAR, &status);
                                		assert(res == MPI_SUCCESS);
                                		res = MPI_Get_count(&status, MPI_CHAR, &count);
                                		assert(res == MPI_SUCCESS);
                                		assert(count == (int)chr_buff_size[n]);
                                		free(buffer_out_vec[n]);
                                		chr_buff_size[n] = 0;
					}
				}

				aft = MPI_Wtime();
                                fprintf(stderr, "rank: %d :: %s: wrote results (%.02f) \n", rank_num, __func__, aft - bef);
                                total_time_writing += (aft - bef);

			}

			if ( write_format == 0 || write_format == 1) {
				                      	      		
                                bef = MPI_Wtime();

                                int ret_code = 0;
                                struct thread_data_compress_by_chr *tdc;
                                tdc = malloc (NUM_THREADS * sizeof(struct thread_data_compress_by_chr));
                                bef = MPI_Wtime();
                                pthread_attr_init(&attr);
                                pthread_attr_setstacksize(&attr, SMALL_STACK);
                                pthread_attr_setdetachstate(&attr, 0);

                                for( n = 0; n < NUM_THREADS; n++ ){
                                	tdc[n].total_thread = NUM_THREADS;
                                        tdc[n].thread_id = n;
					tdc[n].indix_thr = &indix;
					tdc[n].add_in_disc = add_in_disc;
					tdc[n].start_addr = start_addr;
					tdc[n].line_size_to_cpy = line_size_to_cpy;
					tdc[n].sam_buff_dest = sam_buff_dest;	
					tdc[n].total_sam_line = total_sam_line;
					tdc[n].incrmnt = incrmnt;
                                        tdc[n].job_rank = rank_num;
                                        tdc[n].comp_level = compression_level;
                                        tdc[n].fh_out = fh_out;
                                        ret_code = pthread_create(&threads[n], &attr, compress_thread_by_chr, (void *)(&tdc[n]));
                                        if (ret_code) {
                                        	fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                                        }
                                }
                                for(n=0; n<NUM_THREADS; n++)
                                	pthread_join(threads[n], (void *)(&tdc[n]));
				
				pthread_attr_destroy(&attr);
                                free(tdc);
                                total_time_writing += (aft - bef);
                                	                                		
                                for (n = 0; n < reads; n++) free(seqs[n].sam);
                                free(add_in_disc);
                                free(sam_buff_dest);
                                free(seqs);
                                free(start_addr);
                                free(line_size_to_cpy);
						
				aft = MPI_Wtime();
                        	fprintf(stderr, "rank: %d :: %s: wrote results (%.02f) \n", rank_num, __func__, aft - bef);
                        	total_time_writing += (aft - bef);
                                                       	
                        }

                        free(buffer_r1);
                        free(buffer_r2);
                        fprintf(stderr, "rank: %d :: finish for chunck %zu \n", rank_num, u1);

			 //we update u1
			 //get the index_chunk
			if ( proc_num > 1) {	 
				bef = MPI_Wtime();

                        	MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank_target, 0, win);
                        	MPI_Request req;
                        	MPI_Rget(&index_chunk, 1, MPI_UINT64_T, rank_target, 0, 1, MPI_UINT64_T, win, &req);
                        	MPI_Wait(&req, MPI_STATUS_IGNORE);
                        	u1 = (uint64_t)index_chunk;
                        	MPI_Fetch_and_op(&incr,&index_chunk, MPI_UINT64_T, rank_target, 0, MPI_SUM, win);
                        	MPI_Win_flush(rank_target, win);
                        	MPI_Win_unlock(rank_target, win);

                        	aft = MPI_Wtime();
                        	fprintf(stderr, "rank %d ::: update u1 = %zu :: time %.02f \n",rank_num, u1, aft-bef);
			}
			else u1++;
		} //end for (u1 = 0; u1 < chunk_count; u1++){


		MPI_Barrier(MPI_COMM_WORLD);
                if ( (write_format == 1) && (rank_num == 0)) {
                	for (n = 0; n < (indix.bns->n_seqs + incrmnt); n++){
                        	static uint8_t magic[28] =  "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00";
                                res = MPI_File_write_shared(fh_out[n], magic, 28, MPI_BYTE, &status);
                                assert(res == MPI_SUCCESS);
                                res = MPI_Get_count(&status, MPI_BYTE, &count);
                                assert(res == MPI_SUCCESS);
                                assert(count == 28);
                        }
                }


	
		//free data structures and close sam files 
		for (n = 0; n < (indix.bns->n_seqs + incrmnt); n++)  {
			free(files_out_sam_name[n]);
                	free(file_map_by_chr[n]);
                	res = MPI_File_close(&fh_out[n]);
                	assert(res == MPI_SUCCESS);
        	}

	} //end if (file_r2 != NULL && stat_r1.st_size == stat_r2.st_size)
	
	if (file_r1 != NULL && file_r2 != NULL && stat_r1.st_size != stat_r2.st_size) {
	   
	   	/*
		 *	We are in the case the reads are paired and trimmed
		 */

		aft = 0; aft++;
		bef = 0; bef++;

		//global offset contains the starting offset in the fastq for each process
		//TODO?: use only one vector for all the global offsets iot reduce communications
		size_t *goff 	= NULL;
		size_t *goff2 	= NULL;
		goff 	= malloc((proc_num * NUM_THREADS_1 + 1) * sizeof(size_t));
		goff2 	= malloc((proc_num * NUM_THREADS_1 + 1) * sizeof(size_t));

		//this function is used to fill the goff vectors
		 size_t *goff_inter = calloc( (proc_num * NUM_THREADS_1 + 1) , sizeof(size_t));
                size_t *goff_inter_2 = calloc( (proc_num * NUM_THREADS_1 + 1) , sizeof(size_t));

                //this function is used to fill the goff vectors
                bef = MPI_Wtime();
                find_process_starting_offset_mt(goff, stat_r1.st_size, file_r1, proc_num, rank_num, NUM_THREADS_1);
                find_process_starting_offset_mt(goff2, stat_r2.st_size, file_r2, proc_num, rank_num, NUM_THREADS_1);
                aft = MPI_Wtime();
		if (rank_num == 0)
                	fprintf(stderr, "%s: rank %d time spend in finding process start offset = (%.02f) \n", __func__, rank_num, aft - bef);

		//now we exchange the goff buffer between all proc
		//rank 0 gather the vector
		int i12=0;
                for ( i12 = 0; i12 < proc_num * NUM_THREADS_1 + 1; i12++ ){
                        goff_inter[i12] = goff[i12];
                        goff_inter_2[i12] = goff2[i12];
                }

                res = MPI_Allgather(&goff_inter[rank_num*NUM_THREADS_1], NUM_THREADS_1, MPI_LONG_LONG_INT, goff , NUM_THREADS_1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
                assert(res == MPI_SUCCESS);

                res = MPI_Allgather(&goff_inter_2[rank_num*NUM_THREADS_1], NUM_THREADS_1, MPI_LONG_LONG_INT, goff2 , NUM_THREADS_1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
                assert(res == MPI_SUCCESS);

                free(goff_inter);
                free(goff_inter_2);

	
		//we compute the new size according to the shift
		//We calculate the size to read for each process
				
                MPI_Barrier(MPI_COMM_WORLD);


		//opt->chunk_size=100000000;
		//fprintf(stderr,"rank %d ::: chunck size = %zu \n", rank_num, opt->chunk_size);

		///resources needed to find the offsets and size of each read.	
		size_t grand_total_num_reads 	= 0; 
		size_t grand_total_num_reads_2	= 0;
		size_t local_num_reads 		= 0;
		size_t local_num_reads_2	= 0;
		size_t total_num_reads 		= 0;
		size_t total_num_reads_2	= 0;

		//Here I decided to keep separated vectors because otherwise with all the reallocs they would be too big and take too much contiguous memory space
		int *local_read_size    	= NULL;
		int *local_read_size_2  	= NULL;
		size_t *local_read_bytes    	= NULL;
		size_t *local_read_bytes_2    	= NULL;
		size_t *local_read_offsets 	= NULL;
		size_t *local_read_offsets_2 	= NULL;

		///find offsets and sizes for the first file
		bef = MPI_Wtime();

		pthread_attr_t attr;
                pthread_attr_init(&attr);
                pthread_attr_setstacksize(&attr, SMALL_STACK);
                pthread_attr_setdetachstate(&attr, 0);

                pthread_t threads_1[NUM_THREADS_1];

                struct struct_data_thread_1 *td_1;

                size_t *local_num_reads_t            = calloc(NUM_THREADS_1, sizeof(size_t));
                size_t *total_num_reads_t            = calloc(NUM_THREADS_1, sizeof(size_t));
                size_t **local_read_offsets_t        = calloc(NUM_THREADS_1, sizeof(size_t*));
                size_t **local_read_bytes_t          = calloc(NUM_THREADS_1, sizeof(size_t*));
                int **local_read_size_t              = calloc(NUM_THREADS_1, sizeof(int*));

                td_1 = malloc(NUM_THREADS_1 * sizeof(struct struct_data_thread_1));

                int goff_idx = 0;
                for ( n = 0; n < NUM_THREADS_1; n++){

                        goff_idx = (rank_num * NUM_THREADS_1) + n;
                        td_1[n].offset_in_file_mt         = goff[goff_idx];
                        td_1[n].size2read_mt              = goff[goff_idx + 1] - goff[goff_idx];
                        td_1[n].file_r1_mt                = file_r1;
                        td_1[n].local_num_reads_mt        = &local_num_reads_t[n];
                        td_1[n].total_num_reads_mt        = &total_num_reads_t[n];
                        td_1[n].local_read_offsets_mt     = &local_read_offsets_t[n];
                        td_1[n].local_read_size_mt        = &local_read_size_t[n];
                        td_1[n].local_read_bytes_mt       = &local_read_bytes_t[n];
                        td_1[n].proc_num_mt               = proc_num;
                        td_1[n].rank_num_mt               = rank_num;
                        td_1[n].thread_num_mt             = n;
                        td_1[n].previous_read_num         = 0;
                        pthread_create(&threads_1[n], &attr, find_reads_size_and_offsets_mt, (void *)(&td_1[n]));

                }

                total_num_reads = 0;
                for (n = 0; n < NUM_THREADS_1; n++){
                        pthread_join(threads_1[n], (void *)(&td_1[n]));
                        total_num_reads += *(td_1[n].total_num_reads_mt);
                }

		local_read_offsets  = calloc(total_num_reads, sizeof(size_t));
                local_read_size     = calloc(total_num_reads, sizeof(int));
                local_read_bytes    = calloc(total_num_reads, sizeof(size_t));

                assert(local_read_offsets);
                assert(local_read_size);
                assert(local_read_bytes);

                size_t tmp_var = 0;
                for (n = 0; n < NUM_THREADS_1; n++){
                        td_1[n].local_read_offsets     = local_read_offsets;
                        td_1[n].local_read_size        = local_read_size;
                        td_1[n].local_read_bytes       = local_read_bytes;
                        td_1[n].previous_read_num      = tmp_var;
                        tmp_var                        += *(td_1[n].total_num_reads_mt);
			pthread_create(&threads_1[n], &attr, copy_local_read_info_mt, (void *)(&td_1[n]));

                }

                for (n = 0; n < NUM_THREADS_1; n++){
                        pthread_join(threads_1[n], (void *)(&td_1[n]));

                        free(local_read_offsets_t[n]);
                        free(local_read_bytes_t[n]);
                        free(local_read_size_t[n]);
                }

                free(local_num_reads_t);
                free(total_num_reads_t);

                free(local_read_offsets_t);
                free(local_read_bytes_t);
                free(local_read_size_t);

                pthread_attr_destroy(&attr);
                free(td_1);

                aft = MPI_Wtime();
                fprintf(stderr, "%s: rank %d num reads parsed: %zu ::: time spend reading and parsing entire buffer 1 = (%.02f) \n", __func__, rank_num, total_num_reads, aft - bef);

                if (goff) free(goff);



		/*
		find_reads_size_and_offsets(	goff[ind],
				 		siz2read,
						file_r1,
						&local_num_reads,
						&total_num_reads,
						&local_read_offsets,
						&local_read_size,
						&local_read_bytes,
						proc_num,
						rank_num);

		*/
	
		//communications about the first file 
		res = MPI_Reduce(&total_num_reads, &grand_total_num_reads, 1, MPI_LONG_LONG_INT, MPI_SUM, 0,MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);	

		res = MPI_Bcast(&grand_total_num_reads, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);	

		//we dispatch the previous result
		if ( rank_num == 0 || rank_num ==1 )
			fprintf(stderr, "rank %d ::: total_num_reads (1rst file)= %zu \n", rank_num, grand_total_num_reads);

		MPI_Barrier(MPI_COMM_WORLD);

		bef = MPI_Wtime();
		///find offsets and sizes for second file + time the process
		/*
		find_reads_size_and_offsets(	goff2[ind],
						siz2read_2,
						file_r2,
						&local_num_reads_2,
						&total_num_reads_2,
						&local_read_offsets_2,
						&local_read_size_2,
						&local_read_bytes_2,
						proc_num,
						rank_num);
		*/

		pthread_attr_t attr2;
                pthread_attr_init(&attr2);
                pthread_attr_setstacksize(&attr2, SMALL_STACK);
                pthread_attr_setdetachstate(&attr2, 0);

                pthread_t threads_2[NUM_THREADS_1];

                struct struct_data_thread_1 *td_2;

                size_t *local_num_reads_t2            = calloc(NUM_THREADS_1, sizeof(size_t));
                size_t *total_num_reads_t2            = calloc(NUM_THREADS_1, sizeof(size_t));
                size_t **local_read_offsets_t2        = calloc(NUM_THREADS_1, sizeof(size_t*));
                size_t **local_read_bytes_t2          = calloc(NUM_THREADS_1, sizeof(size_t*));
                int **local_read_size_t2              = calloc(NUM_THREADS_1, sizeof(int*));

                td_2 = malloc(NUM_THREADS_1 * sizeof(struct struct_data_thread_1));

                int goff_idx2 = 0;
                for ( n = 0; n < NUM_THREADS_1; n++){

                        goff_idx2 = (rank_num * NUM_THREADS_1) + n;
                        td_2[n].offset_in_file_mt         = goff2[goff_idx2];
                        td_2[n].size2read_mt              = goff2[goff_idx2 + 1] - goff2[goff_idx2];
                        td_2[n].file_r1_mt                = file_r2;
                        td_2[n].local_num_reads_mt        = &local_num_reads_t2[n];
                        td_2[n].total_num_reads_mt        = &total_num_reads_t2[n];
                        td_2[n].local_read_offsets_mt     = &local_read_offsets_t2[n];
                        td_2[n].local_read_size_mt        = &local_read_size_t2[n];
                        td_2[n].local_read_bytes_mt       = &local_read_bytes_t2[n];
                        td_2[n].proc_num_mt               = proc_num;
                        td_2[n].rank_num_mt               = rank_num;
                        td_2[n].thread_num_mt             = n;
                        td_2[n].previous_read_num         = 0;
                        pthread_create(&threads_2[n], &attr2, find_reads_size_and_offsets_mt, (void *)(&td_2[n]));

                }


                total_num_reads_2 = 0;
                for (n = 0; n < NUM_THREADS_1; n++){
                        pthread_join(threads_2[n], (void *)(&td_2[n]));
                        total_num_reads_2 += *(td_2[n].total_num_reads_mt);
                }

		local_read_offsets_2  = calloc(total_num_reads_2, sizeof(size_t));
                local_read_size_2     = calloc(total_num_reads_2, sizeof(int));
                local_read_bytes_2    = calloc(total_num_reads_2, sizeof(size_t));

                assert(local_read_offsets_2);
                assert(local_read_size_2);
                assert(local_read_bytes_2);

                size_t tmp_var2 = 0;
                for (n = 0; n < NUM_THREADS_1; n++){
                        td_2[n].local_read_offsets     = local_read_offsets_2;
                        td_2[n].local_read_size        = local_read_size_2;
                        td_2[n].local_read_bytes       = local_read_bytes_2;
                        td_2[n].previous_read_num      = tmp_var2;
                        tmp_var2                       += *(td_2[n].total_num_reads_mt);

			pthread_create(&threads_2[n], &attr2, copy_local_read_info_mt, (void *)(&td_2[n]));
                }

                for (n = 0; n < NUM_THREADS_1; n++){
                        pthread_join(threads_2[n], (void *)(&td_2[n]));

                        free(local_read_offsets_t2[n]);
                        free(local_read_bytes_t2[n]);
                        free(local_read_size_t2[n]);
                }

                free(local_num_reads_t2);
                free(total_num_reads_t2);

                free(local_read_offsets_t2);
                free(local_read_bytes_t2);
                free(local_read_size_t2);

                pthread_attr_destroy(&attr2);
                free(td_2);

                aft = MPI_Wtime();
                fprintf(stderr, "%s: rank %d num reads parsed: %zu ::: time spend reading and parsing entire buffer 2 = (%.02f) \n", __func__, rank_num, total_num_reads_2, aft - bef);

                if (goff2) free (goff2);

		//commmunications about the second file
		res = MPI_Reduce(&total_num_reads_2, &grand_total_num_reads_2, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);	
	
		res = MPI_Bcast(&grand_total_num_reads_2, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);	

		//we dispatch the previous result
		if ( rank_num == 0 || rank_num ==1 )
			fprintf(stderr, "rank %d ::: total_num_reads (for 2nd file) = %zu \n", rank_num, grand_total_num_reads_2);
		
		assert ( grand_total_num_reads_2 == grand_total_num_reads );
	
		///Now find the required information to create the chunks, such as total number of bases
		local_num_reads 	= total_num_reads;
		local_num_reads_2 	= total_num_reads_2;
	
		//change chunck_size if necessary
		//opt->chunk_size=10000000;

		size_t i;
		//since all the sizes are not the same anymore, you have to loop

		size_t min_num_read = min(local_num_reads, local_num_reads_2);
		size_t bases_tmp 	= 0;
		size_t chunck_num 	= 0;
		size_t chunck_num_2	= 0;

		for(i=0; i < min_num_read; i++)   {

			bases_tmp  += (local_read_size[i] + local_read_size_2[i]);

			if ( bases_tmp > fixed_chunk_size ){

				bases_tmp = 0;
				chunck_num++;
				chunck_num_2++;
			}
		}

		chunck_num 	+= 2;
		chunck_num_2 	+= 2;
	
		fprintf(stderr,"Rank %d :: chunk_num = %zu\n", rank_num, chunck_num);
		fprintf(stderr,"Rank %d :: chunk_num_2 = %zu\n", rank_num, chunck_num_2);

		///assert that the sizes and offsets found earlier are real.
		//changed the comparisons here for the size. they were compared to blen which has become impossible
		size_t h=0;
		/*
		for (h=0; h < total_num_reads; h++)
		{
			assert(local_read_size[h] > 0 );
			assert(local_read_offsets[h] >= 0 );
		}
		for(h=0; h < total_num_reads_2; h++)
		{
			assert(local_read_size_2[h] > 0 );
			assert(local_read_offsets_2[h] >= 0 );
		}
		*/
		///Now that we know how many chunks we need, we create the vectors associated
		///and we update their informations according to the file they are treating

		// we allocate vector for chunks offset
		size_t *begin_offset_chunk 	= calloc(chunck_num, sizeof(size_t));
		size_t *chunk_size 		= calloc(chunck_num, sizeof(size_t));
		size_t *reads_in_chunk 		= calloc(chunck_num, sizeof(size_t));

		assert( begin_offset_chunk != NULL );
		assert( chunk_size != NULL );
		assert( reads_in_chunk != NULL );
	
		size_t chunk_count = 0;

		size_t *begin_offset_chunk_2 	= calloc(chunck_num_2, sizeof(size_t));
		size_t *chunk_size_2 		= calloc(chunck_num_2, sizeof(size_t));
		size_t *reads_in_chunk_2 	= calloc(chunck_num_2, sizeof(size_t));

		assert( begin_offset_chunk_2 != NULL );
		assert( chunk_size_2 != NULL );
		assert( reads_in_chunk_2 != NULL );


		maxsiz = ( fixed_chunk_size ); 
		MPI_Barrier(MPI_COMM_WORLD);
		fprintf(stderr,"rank %d ::: Call find_chunks_info \n", rank_num);
		// the detail of he paramters is at the function definition
		// fprintf(stderr, "rank %d ::: begin_offset_chunk = %zu \n", rank_num, begin_offset_chunk[0]);
		// fprintf(stderr, "rank %d ::: begin_offset_chunk_2 = %zu \n", rank_num, begin_offset_chunk_2[0]);
		// fprintf(stderr, "rank %d ::: chunk_size = %zu \n", rank_num, chunk_size);
		// fprintf(stderr, "rank %d ::: chunk_size_2 = %zu \n", rank_num, chunk_size_2);

		bef = MPI_Wtime();
		find_chunks_info_trim(begin_offset_chunk,
				 begin_offset_chunk_2,
				 chunk_size,
				 chunk_size_2,
				 reads_in_chunk,
				 reads_in_chunk_2,
				 local_read_size,
				 local_read_size_2,
				 local_read_bytes,
				 local_read_bytes_2,
				 local_read_offsets,
				 local_read_offsets_2,
				 rank_num,
				 proc_num,
				 local_num_reads,
				 local_num_reads_2,
				 grand_total_num_reads,
				 grand_total_num_reads_2,
				 maxsiz,
				 &chunk_count);

		aft = MPI_Wtime();
		fprintf(stderr, "%s ::: rank %d ::: evaluating offsets chuncks and sizes (%.02f) found %zu chuncks \n", __func__, rank_num, aft - bef, chunk_count);
	
		MPI_Barrier(MPI_COMM_WORLD);

		if (local_read_offsets) 	free(local_read_offsets);
		if (local_read_size) 		free(local_read_size);
		if (local_read_offsets_2) 	free(local_read_offsets_2);
		if (local_read_size_2)		free(local_read_size_2);

		//verify the number of reads

		size_t num_reads_1 = 0;
		size_t num_reads_2 = 0;
		size_t total_num_reads_v1 = 0;
		size_t total_num_reads_v2 = 0;

		for(h=0; h < chunk_count; h++){
			 assert( reads_in_chunk[h] == reads_in_chunk_2[h]);
			 num_reads_1 += reads_in_chunk[h];
			 num_reads_2 += reads_in_chunk_2[h];
		}

		fprintf(stderr, "rank %d ::: after find_chunks_info ::: num_reads_1 = %zu \n", rank_num, num_reads_1);
		fprintf(stderr, "rank %d ::: after find_chunks_info ::: num_reads_2 = %zu \n", rank_num, num_reads_2);

		res = MPI_Reduce(&num_reads_1, &total_num_reads_v1, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);

		res = MPI_Bcast(&total_num_reads_v1, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);

		if ( rank_num == 0 || rank_num ==1 )
			fprintf(stderr, "rank %d ::: after find_chunks_info ::: total_num_reads forward = %zu \n", rank_num, total_num_reads_v1);

		res = MPI_Reduce(&num_reads_2, &total_num_reads_v2, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);

		res = MPI_Bcast(&total_num_reads_v2, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);

		if ( rank_num == 0 || rank_num ==1 )
			fprintf(stderr, "rank %d ::: after find_chunks_info ::: total_num_reads backward = %zu \n", rank_num, total_num_reads_v2);
		
		assert( total_num_reads_v1 == total_num_reads_v2 );
		assert( grand_total_num_reads == total_num_reads_v1);
		assert( grand_total_num_reads_2 == total_num_reads_v2);
		MPI_Barrier(MPI_COMM_WORLD);

		//display test to check that each vector was well filled with the right information
		//fprintf(stderr,"rank %d ::: begin offset %zu, %zu\n", rank_num, begin_offset_chunk[0], begin_offset_chunk_2[0]);
		//fprintf(stderr,"rank %d ::: chunk size %zu, %zu\n", rank_num, chunk_size[0], chunk_size_2[0]);
		//fprintf(stderr,"rank %d ::: reads in chunk %zu, %zu\n", rank_num, reads_in_chunk[0], reads_in_chunk_2[0]);

		bef = MPI_Wtime();
		size_t total_chunks = 0;
                res = MPI_Reduce(&chunk_count, &total_chunks, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
                assert ( res == MPI_SUCCESS);
                res = MPI_Bcast(&total_chunks, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
                assert ( res == MPI_SUCCESS);

                size_t chunks_per_rank[proc_num];
                res = MPI_Allgather(&chunk_count, 1, MPI_LONG_LONG_INT, chunks_per_rank, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);

                all_begin_offset_chunk = malloc(total_chunks * sizeof(size_t));
                all_begin_offset_chunk[0]=0;
                all_chunk_size = malloc(total_chunks * sizeof(size_t));
                all_chunk_size[0] = 0;
                all_reads_in_chunk =  malloc(total_chunks * sizeof(size_t));
                all_reads_in_chunk[0] = 0;

                all_begin_offset_chunk_2 = malloc(total_chunks * sizeof(size_t));
                all_begin_offset_chunk_2[0]=0;
                all_chunk_size_2 = malloc(total_chunks * sizeof(size_t));
                all_chunk_size_2[0] = 0;
                all_reads_in_chunk_2 =  malloc(total_chunks * sizeof(size_t));
                all_reads_in_chunk_2[0] = 0;

                int displ_chunk[proc_num];
                displ_chunk[0] = 0;

                for (i = 1; i < proc_num; i++) displ_chunk[i] = (displ_chunk[i-1] + chunks_per_rank[i-1]);

                int indx=displ_chunk[rank_num];
                for (i = 0; i <  chunks_per_rank[rank_num]; i++) {
                        all_chunk_size[indx+i]=chunk_size[i];
                        all_begin_offset_chunk[indx+i]=begin_offset_chunk[i];
                        all_reads_in_chunk[indx+i]=reads_in_chunk[i];

                        all_chunk_size_2[indx+i]=chunk_size_2[i];
                        all_begin_offset_chunk_2[indx+i]=begin_offset_chunk_2[i];
                        all_reads_in_chunk_2[indx+i]=reads_in_chunk_2[i];
                }

		if (rank_num > 0){
                        res=MPI_Send(chunk_size, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
                        assert(res == MPI_SUCCESS);
                }
                else{
                        for (i = 1; i < proc_num; i++){
                                res=MPI_Recv(&(all_chunk_size[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                assert(res == MPI_SUCCESS);
                        }
                }

                if (rank_num > 0){
                        res=MPI_Send(chunk_size_2, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
                        assert(res == MPI_SUCCESS);
                }
                else{
                        for (i = 1; i < proc_num; i++){
                                res=MPI_Recv(&(all_chunk_size_2[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                assert(res == MPI_SUCCESS);
                        }
                }

                if (rank_num > 0){
                        res=MPI_Send(begin_offset_chunk, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
                        assert(res == MPI_SUCCESS);
                }
                else{
                        for (i = 1; i < proc_num; i++){
                                res=MPI_Recv(&(all_begin_offset_chunk[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                assert(res == MPI_SUCCESS);
                        }
                }
	
		if (rank_num > 0){
                        res=MPI_Send(begin_offset_chunk_2, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
                        assert(res == MPI_SUCCESS);
                }
                else{
                for (i = 1; i < proc_num; i++){
                        res=MPI_Recv(&(all_begin_offset_chunk_2[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        assert(res == MPI_SUCCESS);
                        }
                }

                if (rank_num > 0){
                        res=MPI_Send(reads_in_chunk, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
                        assert(res == MPI_SUCCESS);
                }
                else{
                        for (i = 1; i < proc_num; i++){
                                res=MPI_Recv(&(all_reads_in_chunk[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                assert(res == MPI_SUCCESS);
                        }
                }

                if (rank_num > 0){
                        res=MPI_Send(reads_in_chunk_2, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
                        assert(res == MPI_SUCCESS);
                }
                else{
                        for (i = 1; i < proc_num; i++){
                                res=MPI_Recv(&(all_reads_in_chunk_2[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                assert(res == MPI_SUCCESS);
                        }
                }

                res=MPI_Bcast(all_chunk_size, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
                res=MPI_Bcast(all_reads_in_chunk, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
                res=MPI_Bcast(all_begin_offset_chunk, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

                res=MPI_Bcast(all_chunk_size_2, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
                res=MPI_Bcast(all_reads_in_chunk_2, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
                res=MPI_Bcast(all_begin_offset_chunk_2, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

                assert(res == MPI_SUCCESS);

		for ( i = 0; i < total_chunks; i++)
                        assert(all_chunk_size[i] != 0);

                free(chunk_size);
                free(begin_offset_chunk);
                free(reads_in_chunk);
                free(chunk_size_2);
                free(begin_offset_chunk_2);
                free(reads_in_chunk_2);

		aft = MPI_Wtime();
                fprintf(stderr, "%s: rank %d time spend gathering chunks info = (%.02f) \n", __func__, rank_num, aft - bef);

		/// Map reference genome indexes in shared memory (by host)
		bef = MPI_Wtime();
		map_indexes(file_map, &count, &indix, &ignore_alt, &win_shr, shared_mem);
		aft = MPI_Wtime();
		fprintf(stderr, "%s ::: rank %d ::: mapped indexes (%.02f)\n", __func__, rank_num, aft - bef);

		///This is a testline to stop the program wherever I want
		//if(0){ MPI_Barrier(MPI_COMM_WORLD); MPI_Finalize(); return 0;}


		//we create a vector with chromosom names 
		int s;
		int incrmnt = 2;
		if (dofixmate) incrmnt = 1;
		
		MPI_File *fh_out        = malloc( (indix.bns->n_seqs + incrmnt) * sizeof(MPI_File));
                char *files_out_sam_name[indix.bns->n_seqs + incrmnt];
                char *file_map_by_chr[indix.bns->n_seqs + incrmnt];

                for (s = 0; s < indix.bns->n_seqs; ++s){
	                files_out_sam_name[s] = malloc( (strlen( indix.bns->anns[s].name) + 1 ) * sizeof(char));
                        files_out_sam_name[s][strlen( indix.bns->anns[s].name)] = 0;
                        char *p0 = files_out_sam_name[s];
        	        memmove(p0, indix.bns->anns[s].name, strlen( indix.bns->anns[s].name));
                }
		if (!dofixmate){
                        char DISCORDANT[] = "discordant";
                        files_out_sam_name[s++] = strdup(DISCORDANT);
                }
                char UNMAPPED[]   = "unmapped";
                files_out_sam_name[s++] = strdup(UNMAPPED);

                //char *file_map_by_chr[(indix.bns->n_seqs + 2)];
                size_t file_name_len = 0;
                
		for (s = 0; s < (indix.bns->n_seqs + incrmnt); ++s){
                        /* Derived file names */
                        file_name_len = strlen(output_path) + strlen(files_out_sam_name[s]) + 6;
                        file_map_by_chr[s] = calloc( file_name_len, sizeof(char));
                        if ( write_format == 2  )
                                sprintf(file_map_by_chr[s], "%s/%s.sam", output_path, files_out_sam_name[s]);
                        if ( write_format == 1  )
                                sprintf(file_map_by_chr[s], "%s/%s.bam", output_path, files_out_sam_name[s]);
                        if ( write_format == 0  )  
                                sprintf(file_map_by_chr[s], "%s/%s.gz", output_path, files_out_sam_name[s]);
    

                }

		if ( write_format == 2  )
                        create_sam_header_by_chr_file(file_map_by_chr, &indix, &count, hdr_line, rg_line, pg_line, rank_num);

                if ( write_format == 1 || write_format == 0)
                        create_bam_header_by_chr_file(file_map_by_chr, &indix, &count, hdr_line, rg_line, pg_line, rank_num, compression_level, dofixmate);

                        
		 ///Create SAM header
		 //TODO: Add line for BWA version
		//create_sam_header_by_chr_file(file_map_by_chr, &indix, &count, hdr_line, rg_line, rank_num);

		//we create a file which contains only the header
		//this file will be used by the sorting to get individual chromosom file
		//create_sam_header(file_out, &indix, &count, hdr_line, rg_line, rank_num);

		//we add the header in the discordant SAM
		//create_sam_header(file_map_by_chr[indix.bns->n_seqs], &indix, &count, hdr_line, rg_line, rank_num);
		
		if (dofixmate){
			for (s = 0; s < (indix.bns->n_seqs + 1); ++s){		
				res = MPI_File_open(MPI_COMM_WORLD, file_map_by_chr[s], MPI_MODE_CREATE|MPI_MODE_WRONLY|MPI_MODE_APPEND, MPI_INFO_NULL, &fh_out[s]);
				assert(res == MPI_SUCCESS);
			}
		}
		else{
			for (s = 0; s < (indix.bns->n_seqs + 2); ++s){
                                res = MPI_File_open(MPI_COMM_WORLD, file_map_by_chr[s], MPI_MODE_CREATE|MPI_MODE_WRONLY|MPI_MODE_APPEND, MPI_INFO_NULL, &fh_out[s]);
                                assert(res == MPI_SUCCESS);
                        }

		}
		//now we open fastq files backward and forward
		if (file_r1 != NULL) {
			//fprintf(stderr, "rank %d ::: open file %s \n",rank_num, file_r1);
			res = MPI_File_open(MPI_COMM_WORLD, file_r1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_r1);
			assert(res == MPI_SUCCESS);
		}
		if (file_r2 != NULL) {
			//fprintf(stderr, "rank %d ::: open file %s \n",rank_num, file_r2);
			res = MPI_File_open(MPI_COMM_WORLD, file_r2, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_r2);
			assert(res == MPI_SUCCESS);
		}

		buffer_r1 = buffer_r2 = NULL; seqs = NULL;
		
		bef = MPI_Wtime();
                if ( proc_num > 1) {
                        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank_target, 0, win);
                        MPI_Request req;
                        MPI_Rget(&index_chunk, 1, MPI_UINT64_T, rank_target, 0, 1, MPI_UINT64_T, win, &req);
                        MPI_Wait(&req, MPI_STATUS_IGNORE);
                        u1 = (uint64_t)index_chunk;
                        MPI_Fetch_and_op(&incr, &index_chunk, MPI_UINT64_T, rank_target, 0, MPI_SUM, win);
                        MPI_Win_flush(rank_target, win);
                        MPI_Win_unlock(rank_target, win);
                }
                
                aft = MPI_Wtime();
                fprintf(stderr, "rank %d ::: initial u1 = %zu :: time %.02f \n",rank_num, u1, aft-bef);

		before_local_mapping = MPI_Wtime();

		// here we loop until there's nothing to read
		//we loop the chunck_count
		
		 //localsize_vec contain the size of the buffer to write in the sam file
		
		int *chr_buff_size  = calloc ( (indix.bns->n_seqs + incrmnt), sizeof(int) );
                //char *buffer_out_vec[indix.bns->n_seqs + incrmnt];
 
		while ( u1 < total_chunks ){

			offset_chunk 		= all_begin_offset_chunk[u1];
			size_chunk   		= all_chunk_size[u1];
			offset_chunk_2 		= all_begin_offset_chunk_2[u1];
			size_chunk_2 		= all_chunk_size_2[u1];

			assert(size_chunk 	!= 0);
			assert(size_chunk_2 != 0);
		
			bef = MPI_Wtime();
			///allocate the buffer sizes to store an entire chunk in it
			buffer_r1 = malloc(size_chunk + 1);
			assert(buffer_r1 != NULL);
			buffer_r1[size_chunk] = 0;
		
			buffer_r2 = malloc(size_chunk_2 + 1);
			assert(buffer_r2 != NULL);
			buffer_r2[size_chunk_2]=0;
			
			struct struct_pread_fastq *td_pread1;
                        td_pread1 = malloc (NUM_THREADS_1 * sizeof(struct struct_pread_fastq));
                        bef = MPI_Wtime();
                        pthread_attr_t attr4;
                        pthread_attr_init(&attr4);
                        pthread_attr_setstacksize(&attr4, BIG_STACK);
                        pthread_attr_setdetachstate(&attr4, 0);

                        for( n = 0; n < NUM_THREADS_1; n++ ){
                                td_pread1[n].total_thread = NUM_THREADS_1;
                                td_pread1[n].thread_id = n;
                                td_pread1[n].job_rank = rank_num;
                                td_pread1[n].offset= offset_chunk;
                                td_pread1[n].size = size_chunk;
                                td_pread1[n].buffer = buffer_r1;
                                td_pread1[n].fd  = fh_r1;
                                int ret_code = pthread_create(&threads[n], &attr4, pread_fastq_chunck, (void *)(&td_pread1[n]));
                                if (ret_code) {
                                        fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                                }
                        }
                        for(n=0; n<NUM_THREADS_1; n++)
                        pthread_join(threads[n], (void *)(&td_pread1[n]));

                        pthread_attr_destroy(&attr4);
                        free(td_pread1);

                        assert(strlen(buffer_r1) == size_chunk);
                        assert(*buffer_r1 == '@');

                        struct struct_pread_fastq *td_pread2;
                        td_pread2 = malloc (NUM_THREADS_1 * sizeof(struct struct_pread_fastq));
                        bef = MPI_Wtime();
                        pthread_attr_t attr5;
                        pthread_attr_init(&attr5);
                        pthread_attr_setstacksize(&attr5, BIG_STACK);
                        pthread_attr_setdetachstate(&attr5, 0);

			for( n = 0; n < NUM_THREADS_1; n++ ){
                                td_pread2[n].total_thread = NUM_THREADS_1;
                                td_pread2[n].thread_id = n;
                                td_pread2[n].job_rank = rank_num;
                                td_pread2[n].offset= offset_chunk_2;
                                td_pread2[n].size = size_chunk_2;
                                td_pread2[n].buffer = buffer_r2;
                                td_pread2[n].fd  = fh_r2;

                                int ret_code = pthread_create(&threads[n], &attr5, pread_fastq_chunck, (void *)(&td_pread2[n]));
                                if (ret_code) {
                                        fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                                }
                        }
                        for(n=0; n<NUM_THREADS_1; n++)
                                pthread_join(threads[n], (void *)(&td_pread2[n]));

                        pthread_attr_destroy(&attr5);
                        free(td_pread2);

                        assert(strlen(buffer_r2) == size_chunk_2);
                        assert(*buffer_r2 == '@');

                        if (u1 < total_chunks){
                                MPI_File_seek(fh_r1, (MPI_Offset)all_begin_offset_chunk[u1 + proc_num], MPI_SEEK_CUR ) ;
                                MPI_File_seek(fh_r2, (MPI_Offset)all_begin_offset_chunk[u1 + proc_num], MPI_SEEK_CUR ) ;
                        }

			///some stats
			aft = MPI_Wtime();
			fprintf(stderr, "%s: num_rank = %d :: read sequences (%.02f)\n", __func__, rank_num, aft - bef);
			total_time_reading_seq += (aft - bef);

			///operations on the number of reads
			bef = MPI_Wtime();
			reads_r1 = all_reads_in_chunk[u1];
			reads_r2 = all_reads_in_chunk_2[u1];
			reads = reads_r1 + reads_r2;
			bases = 0;
			total_local_reads_aligned += reads/2;
			assert(reads <= INT_MAX);

			/* Parse sequences ... */
			seqs = malloc(reads * sizeof(*seqs));
			assert(seqs != NULL);

			if (file_r1 != NULL) {
				p1 = q1 = buffer_r1; e1 = buffer_r1 + size_chunk; line_number = 0;
							
				while (q1 < e1) {
					if (*q1 != '\n') { q1++; continue; }
					/* We have a full line ... process it */
					*q1 = '\0'; 
					//n will take 4 times in row the same value from 0 with a step of 2
					//it is because the 2 operations are using int types, so it is cast as int twice (within and then outside the brackets)
					n = files * (line_number / 4);
				
					switch (line_number % 4) {
					case 0: /* Line1: Name and Comment */
						assert(*p1 == '@');
						seqs[n].name   = p1 + 1;
						while (*p1 && !isspace((unsigned char)*p1)) {p1++;}
						if (*(p1-2) == '/' && isdigit((unsigned char)*(p1-1))) {*(p1-2) = '\0';}
						if (*p1) {*p1++ = '\0';}
						seqs[n].comment = (copy_comment != 0) ? p1 : NULL;
						seqs[n].sam = NULL;
						break;
					case 1: /* Line2: Sequence */
						seqs[n].seq = p1;
						seqs[n].l_seq = q1 - p1;
						bases += seqs[n].l_seq;
						break;
					case 2: /* Line3: Ignored */
						assert(*p1 == '+');
						break;
					case 3: /* Line4: Quality */
						seqs[n].qual = p1;
						break; }
					p1 = ++q1; 
					line_number++; 
				}
			}
			assert( (line_number/4) == reads_r1);
		
			if (file_r2 != NULL) {
				p2 = q2 = buffer_r2; e2 = buffer_r2 + size_chunk_2; line_number = 0;
							
				while (q2 < e2) {
					if (*q2 != '\n') { q2++; continue; }
					/* We have a full line ... process it */
					*q2 = '\0'; 
					//n will take 4 times in row the same value from 0 with a step of 2
					//it is because the 2 operations are using int types, so it is cast as int twice (within and then outside the brackets)
					n = files * (line_number / 4);
				
					switch (line_number % 4) {
					case 0: /* Line1: Name and Comment */
						assert(*p2 == '@');
						seqs[n+1].name   = p2 + 1;
						while (*p2 && !isspace((unsigned char)*p2)) {p2++;}
						if (*(p2-2) == '/' && isdigit((unsigned char)*(p2-1))) {*(p2-2) = '\0';}
						if (*p2) {*p2++ = '\0';}
						seqs[n+1].comment = (copy_comment != 0) ? p2 : NULL;
						seqs[n+1].sam = NULL;
						break;
					case 1: /* Line2: Sequence */
						seqs[n+1].seq = p2;
						seqs[n+1].l_seq = q2 - p2;
						bases += seqs[n+1].l_seq;
						break;
					case 2: /* Line3: Ignored */
						assert(*p2 == '+');
						break;
					case 3: /* Line4: Quality */
						seqs[n+1].qual = p2;
						break; }
					p2 = ++q2; 
					line_number++; 
				}
			}

			assert( (line_number/4) == reads_r2);
			aft = MPI_Wtime();

			fprintf(stderr, "rank %d :: %s: parsed %zu read of sequences (%ld bp) in (%.02f)\n", rank_num, __func__, reads, (long)bases, aft - bef);
			total_time_parsing += (aft - bef);
			//fprintf(stderr, "rank %d ::: [M::%s] read %zu sequences (%ld bp)...\n",rank_num , __func__, reads, (long)bases);
			/* Datas computation ... */
			bef = MPI_Wtime();
			fprintf(stderr, "rank %d ::::  Call memseqs with count sequences %zu \n", rank_num, reads);
			mem_process_seqs(opt, indix.bwt, indix.bns, indix.pac, 0, (int)reads, seqs, pes0);
			aft = MPI_Wtime();
			fprintf(stderr, "rank %d :: %s: computed mappings (%.02f)\n", rank_num, __func__, aft - bef);

			total_time_mapping += (aft - bef);

			/* Write results ... */
			bef = MPI_Wtime();
			
			//int *sam_buff_dest  = calloc ( reads, sizeof(int) );
			//vecto to check if a discordant need to be add in discordant.sam we set 1 if true 0 else
			//int *add_in_disc    =  calloc ( reads, sizeof(int) );
			char *currentCarac;
			char *start_sam_line;
			int chr, mchr;
			size_t sam_line_size;
			int nbchr = indix.bns->n_seqs + incrmnt;
						
			int next;
			char currentLine[MAX_CHAR_SIZE];
			
			//first we cont how many sam line we have
			size_t total_sam_line = 0;
			n = 0;
                       	int ret_code = 0;
                        if (dofixmate){

				struct thread_data *td;
        			td = malloc (NUM_THREADS * sizeof(struct thread_data));
				bef = MPI_Wtime();
				pthread_attr_init(&attr);
        	                pthread_attr_setstacksize(&attr, SMALL_STACK);
                	        pthread_attr_setdetachstate(&attr, 0);
			
				for( n = 0; n < NUM_THREADS; n++ ){
					td[n].total_thread = NUM_THREADS;
					td[n].thread_id = n;
					assert(seqs);
					td[n].seqs_thr = seqs;
					td[n].job_rank = rank_num;
					td[n].total_reads = reads;
					td[n].start_index_seqs = 0;
					td[n].final_index_seqs = reads-1;
					td[n].indix_thr = &indix;
					td[n].total_lines = 0;
					ret_code = pthread_create(&threads[n], &attr, call_fixmate, (void *)(&td[n]));
					if (ret_code) {
                                        	fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                                	}
				}
				for(n=0; n<NUM_THREADS; n++)
					pthread_join(threads[n], (void *)(&td[n]));
				for(n=0; n<NUM_THREADS; n++)
				 	total_sam_line += td[n].total_lines;
				pthread_attr_destroy(&attr);
				free(td);
				
			}
			else{
				
				n = 0;                      
                        	while ( n < reads){
					next = tokenizer(seqs[n].sam,'\n', currentLine);
                               		while (next) {
                                       		total_sam_line++;
                                       		next = tokenizer(NULL, '\n', currentLine);
                               		}
					for(i=0; i<MAX_CHAR_SIZE; i++) currentLine[i]=0;
					next = tokenizer(seqs[n+1].sam,'\n', currentLine);
                               		while (next) {
                                       		total_sam_line++;
                                       		next = tokenizer(NULL, '\n', currentLine);
                        		}
					for(i=0; i<MAX_CHAR_SIZE; i++) currentLine[i]=0;
	                                n = n + 2;
			
				}
				
                        }
		
			if (dofixmate){
				 aft = MPI_Wtime();
                                 fprintf(stderr, "rank: %d :: %s: time spend in fixmate (%.02f) \n", rank_num, __func__, aft - bef);
			}
			
			char tmp_chr[MAX_CHR_NAME_SIZE];
			int *sam_buff_dest  	= calloc ( total_sam_line, sizeof(int) );
			int *add_in_disc	= calloc ( total_sam_line, sizeof(int) );
			char **start_addr    	= malloc ( total_sam_line * sizeof(char*));
			int *line_size_to_cpy   = malloc ( total_sam_line * sizeof(int));	
			size_t incr_line =0;			
			for(i=0; i< MAX_CHAR_SIZE; i++) currentLine[i]=0;
                        bef = MPI_Wtime();

			if (dofixmate){
				for (n = 0; n < reads; n++) {
					int total_sam_line_size = 0;					
					next = tokenizer(seqs[n].sam,'\n', currentLine);
					currentCarac = currentLine;
					start_sam_line = seqs[n].sam;
					//fprintf(stderr, "rank: %d :: sam = %s \n", rank_num, seqs[n].sam);	
					while (next){
				
						start_sam_line = seqs[n].sam + total_sam_line_size;	
						sam_line_size = strlen(currentLine) + 1;
                         			//we update the offset in the
                                        	//source file
                                        	currentLine[sam_line_size - 1] = '\n';
                                        	currentLine[sam_line_size] = '\0';
						currentCarac = currentLine;
						//WE PASS READ NAME
						currentCarac = strstr (currentCarac, "\t");
						currentCarac++;
					 				
						//WE PASS FLAG					
						currentCarac = strstr (currentCarac, "\t");	
						currentCarac++;
					
						if ( currentCarac[0] == '*') chr =  nbchr-1;
						else chr = getChr(currentCarac - 1, files_out_sam_name, nbchr - incrmnt , tmp_chr);

						if (chr < (nbchr - incrmnt)){
                                        		//the read goes in the SAM it belongs 
                                        		chr_buff_size[chr]    		+= sam_line_size;
                                         		sam_buff_dest[incr_line]       	= chr;
							start_addr[incr_line]		= start_sam_line;
						
                                		}
						else {
							//the read goes in unmapped sam
							chr_buff_size[nbchr - incrmnt] 	+= sam_line_size;
                                        		sam_buff_dest[incr_line]        = nbchr - incrmnt ;
							start_addr[incr_line]           = start_sam_line;
						}
					
						line_size_to_cpy[incr_line] = sam_line_size;
						total_sam_line_size += sam_line_size;
						incr_line++;				
						next = tokenizer(NULL, '\n', currentLine);
					}
                
				}
			}
			else{

				for (n = 0; n < reads; n++) {
                                        int total_sam_line_size = 0;
                                        next = tokenizer(seqs[n].sam,'\n', currentLine);
                                        currentCarac = currentLine;
                                        start_sam_line = seqs[n].sam;

                                        while (next){

                                                start_sam_line = seqs[n].sam + total_sam_line_size;
                                                sam_line_size = strlen(currentLine) + 1;
						//we update the offset in the
						//source file
						currentLine[sam_line_size - 1] = '\n';
                                                currentLine[sam_line_size] = '\0';
                                                currentCarac = currentLine;
						//READ NAME
						currentCarac = strstr (currentCarac, "\t");
                                                currentCarac++;
						// FLAG
						currentCarac = strstr (currentCarac, "\t");
                                                currentCarac++;
						// CHROMOSOM	
						if ( currentCarac[0] == '*') chr =  nbchr - incrmnt;
                                                else chr = getChr(currentCarac - 1, files_out_sam_name, nbchr - incrmnt, tmp_chr);
						
						// COORD
						currentCarac = strstr (currentCarac + 1, "\t");
                                                
						// MAPQ
          					strtoull(currentCarac, &currentCarac, 10);
			
						strtoull(currentCarac, &currentCarac, 10);			
						// CIGAR					
						currentCarac = strstr (currentCarac + 1, "\t");
                  				currentCarac++;

						if ( currentCarac[0] == '=') mchr =  chr;
						else if ( currentCarac[0] == '*') mchr =  nbchr - 2;
                                                else mchr = getChr(currentCarac - 1, files_out_sam_name, nbchr - incrmnt, tmp_chr);

						if ( (chr < (nbchr - incrmnt))) {
							chr_buff_size[chr]              += sam_line_size;
                                                        sam_buff_dest[incr_line]        = chr;
                                                        start_addr[incr_line]           = start_sam_line;

                                                }
                                                else {
                                                        chr_buff_size[nbchr - 1] 	+= sam_line_size;
                                                        sam_buff_dest[incr_line]        = nbchr - 1;
                                                        start_addr[incr_line]           = start_sam_line;
                                                }
							
						if ( (chr < (nbchr - incrmnt)) && (mchr < (nbchr - incrmnt)) && (chr != mchr) ){
							chr_buff_size[nbchr - 2] 	+= sam_line_size;
                                                        add_in_disc[incr_line]			= 1;
						}
						
                                                line_size_to_cpy[incr_line] = sam_line_size;
                                                total_sam_line_size += sam_line_size;
                                                incr_line++;
                                                next = tokenizer(NULL, '\n', currentLine);

					}
				}        
			}


			if ( write_format == 2 ){

				bef = MPI_Wtime();
				char *buffer_out_vec[indix.bns->n_seqs + incrmnt];

				//now we fill up the buffer_out_vec
				for (n = 0; n < (indix.bns->n_seqs + incrmnt); n++) {
					assert(chr_buff_size[n] <= INT_MAX);
					if (chr_buff_size[n]){

						buffer_out_vec[n] = (char *)calloc( (chr_buff_size[n] + 1), sizeof(char));
						assert(buffer_out_vec[n] != NULL);
						buffer_out_vec[n][chr_buff_size[n]] = '\0';
					}
				}			
                       		size_t *actual_size = calloc(indix.bns->n_seqs + incrmnt, sizeof(size_t));
				char *p_temp2;
				for (n = 0; n < total_sam_line; n++) {

					p_temp2 = buffer_out_vec[sam_buff_dest[n]] + actual_size[sam_buff_dest[n]];
					memmove(p_temp2, start_addr[n], line_size_to_cpy[n]);
					actual_size[sam_buff_dest[n]] += line_size_to_cpy[n];
					if ( add_in_disc[n] ){
						p_temp2 = buffer_out_vec[nbchr - 2 ] + actual_size[nbchr - 2 ];
						memmove(p_temp2, start_addr[n], line_size_to_cpy[n]);
                                		actual_size[nbchr -2] += line_size_to_cpy[n];
					}
				}
				for (n = 0; n < reads; n++) free(seqs[n].sam); 

				free(add_in_disc);
				free(sam_buff_dest);
				free(seqs);
				free(actual_size);
				free(start_addr);
				free(line_size_to_cpy);

				for (n = 0; n < (indix.bns->n_seqs + incrmnt); n++) {
				
					if (chr_buff_size[n]) {
						res = MPI_File_write_shared(fh_out[n], buffer_out_vec[n], chr_buff_size[n], MPI_CHAR, &status);
						assert(res == MPI_SUCCESS);
						res = MPI_Get_count(&status, MPI_CHAR, &count);
						assert(res == MPI_SUCCESS);
						assert(count == (int)chr_buff_size[n]);
						free(buffer_out_vec[n]);
						chr_buff_size[n] = 0;
					}
				}
			
				aft = MPI_Wtime();
				fprintf(stderr, "rank: %d :: %s: wrote results (%.02f) \n", rank_num, __func__, aft - bef);
				total_time_writing += (aft - bef);
			}
			
			if ( write_format == 0 || write_format == 1) {

                                bef = MPI_Wtime();

                                int ret_code = 0;
                                struct thread_data_compress_by_chr *tdc;
                                tdc = malloc (NUM_THREADS * sizeof(struct thread_data_compress_by_chr));
                                bef = MPI_Wtime();
                                pthread_attr_init(&attr);
                                pthread_attr_setstacksize(&attr, SMALL_STACK);
                                pthread_attr_setdetachstate(&attr, 0);

                                for( n = 0; n < NUM_THREADS; n++ ){
                                        tdc[n].total_thread = NUM_THREADS;
                                        tdc[n].thread_id = n;
                                        tdc[n].indix_thr = &indix;
                                        tdc[n].add_in_disc = add_in_disc;
                                        tdc[n].start_addr = start_addr;
                                        tdc[n].line_size_to_cpy = line_size_to_cpy;
                                        tdc[n].sam_buff_dest = sam_buff_dest;
                                        tdc[n].total_sam_line = total_sam_line;
                                        tdc[n].incrmnt = incrmnt;
                                        tdc[n].job_rank = rank_num;
                                        tdc[n].comp_level = compression_level;
                                        tdc[n].fh_out = fh_out;
                                        ret_code = pthread_create(&threads[n], &attr, compress_thread_by_chr, (void *)(&tdc[n]));
                                        if (ret_code) {
                                                fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                                        }
                                }
                                for(n=0; n<NUM_THREADS; n++)
                                        pthread_join(threads[n], (void *)(&tdc[n]));

                                pthread_attr_destroy(&attr);
                                free(tdc);
                                total_time_writing += (aft - bef);

                                for (n = 0; n < reads; n++) free(seqs[n].sam);
                                free(add_in_disc);
                                free(sam_buff_dest);
                                free(seqs);
                                free(start_addr);
				free(line_size_to_cpy);

                                aft = MPI_Wtime();
                                fprintf(stderr, "rank: %d :: %s: wrote results (%.02f) \n", rank_num, __func__, aft - bef);
                                total_time_writing += (aft - bef);

                        }
			free(buffer_r1);
			free(buffer_r2);
			fprintf(stderr, "rank: %d :: finish for chunck %zu \n", rank_num, u1);


			if ( proc_num > 1) {
		        	bef = MPI_Wtime();
	                       	MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank_target, 0, win);
         		        MPI_Request req;
                        	MPI_Rget(&index_chunk, 1, MPI_UINT64_T, rank_target, 0, 1, MPI_UINT64_T, win, &req);
                        	MPI_Wait(&req, MPI_STATUS_IGNORE);
                        	u1 = (uint64_t)index_chunk;
                        	MPI_Fetch_and_op(&incr,&index_chunk, MPI_UINT64_T, rank_target, 0, MPI_SUM, win);
                        	MPI_Win_flush(rank_target, win);
                        	MPI_Win_unlock(rank_target, win);

                        	aft = MPI_Wtime();
                        	fprintf(stderr, "rank %d ::: update u1 = %zu :: time %.02f \n",rank_num, u1, aft-bef);
			}
			else u1++;

		} //end for loop on chunks
		
		MPI_Barrier(MPI_COMM_WORLD);
                if ( (write_format == 1) && (rank_num == 0)) {
                	for (n = 0; n < (indix.bns->n_seqs + incrmnt); n++){
                        	static uint8_t magic[28] =  "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00";
                                res = MPI_File_write_shared(fh_out[n], magic, 28, MPI_BYTE, &status);
                                assert(res == MPI_SUCCESS);
                                res = MPI_Get_count(&status, MPI_BYTE, &count);
                                assert(res == MPI_SUCCESS);
                                assert(count == 28);
                        }
                }		

		for (n = 0; n < (indix.bns->n_seqs + incrmnt); n++)  {
			free(files_out_sam_name[n]);
			free(file_map_by_chr[n]);
			res = MPI_File_close(&fh_out[n]);
			assert(res == MPI_SUCCESS);
		}
	}// end else case files are trimmed

	if (file_r1 != NULL && file_r2 == NULL){

		/*
		 *
		 * We are in the case reads are single
		 * We do both case trimmed or not 
		 *
		 * We use the same algo as in paired trimmed
		 */
		

		aft = 0; aft++;
		bef = 0; bef++;

		//global offset contains the starting offset in the fastq for each process
		//TODO?: use only one vector for all the global offsets iot reduce communications
		size_t *goff    = malloc((proc_num * NUM_THREADS_1 + 1) * sizeof(size_t));
                size_t *goff_inter = calloc( (proc_num * NUM_THREADS_1 + 1) , sizeof(size_t));

		 bef = MPI_Wtime();
                find_process_starting_offset_mt(goff, stat_r1.st_size, file_r1, proc_num, rank_num, NUM_THREADS_1);
                aft = MPI_Wtime();
                fprintf(stderr, "%s: rank %d time spend in finding process start offset = (%.02f) \n", __func__, rank_num, aft - bef);

                int i12=0;
                for ( i12 = 0; i12 < proc_num * NUM_THREADS_1 + 1; i12++ )  goff_inter[i12] = goff[i12];

		res = MPI_Allgather(&goff_inter[rank_num*NUM_THREADS_1], NUM_THREADS_1, MPI_LONG_LONG_INT, goff , NUM_THREADS_1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
                assert(res == MPI_SUCCESS);

                free(goff_inter);	
			
		//we compute the new size according to the shift
		//We calculate the size to read for each process

		MPI_Barrier(MPI_COMM_WORLD);
	
		///resources needed to find the offsets and size of each read.	
		size_t grand_total_num_reads 	= 0; 
		size_t local_num_reads 		= 0;
		size_t total_num_reads 		= 0;
		
		//Here I decided to keep separated vectors because otherwise with all the reallocs they would be too big and take too much contiguous memory space
		int *local_read_size    	= NULL;
		size_t *local_read_bytes    	= NULL;
		size_t *local_read_offsets	= NULL;
		
		///find offsets and sizes for the first file
		bef = MPI_Wtime();
		/*
		find_reads_size_and_offsets(goff[ind],
					    siz2read,
					    file_r1,
					    &local_num_reads,
					    &total_num_reads,
					    &local_read_offsets,
					    &local_read_size,
					    &local_read_bytes,
					    proc_num,
					    rank_num);
		*/

		pthread_attr_t attr;
                pthread_attr_init(&attr);
                pthread_attr_setstacksize(&attr, SMALL_STACK);
                pthread_attr_setdetachstate(&attr, 0);

                pthread_t threads_1[NUM_THREADS_1];

                struct struct_data_thread_1 *td_1;

                size_t *local_num_reads_t            = calloc(NUM_THREADS_1, sizeof(size_t));
                size_t *total_num_reads_t            = calloc(NUM_THREADS_1, sizeof(size_t));
                size_t **local_read_offsets_t        = calloc(NUM_THREADS_1, sizeof(size_t*));
                size_t **local_read_bytes_t          = calloc(NUM_THREADS_1, sizeof(size_t*));
                int **local_read_size_t              = calloc(NUM_THREADS_1, sizeof(int*));

                td_1 = malloc(NUM_THREADS_1 * sizeof( struct struct_data_thread_1));

                int goff_idx = 0;
                for ( n = 0; n < NUM_THREADS_1; n++){

                        goff_idx = (rank_num * NUM_THREADS_1) + n;
                        td_1[n].offset_in_file_mt         = goff[goff_idx];
                        td_1[n].size2read_mt              = goff[goff_idx + 1] - goff[goff_idx];
                        td_1[n].file_r1_mt                = file_r1;
                        td_1[n].local_num_reads_mt        = &local_num_reads_t[n];
                        td_1[n].total_num_reads_mt        = &total_num_reads_t[n];
                        td_1[n].local_read_offsets_mt     = &local_read_offsets_t[n];
                        td_1[n].local_read_size_mt        = &local_read_size_t[n];
                        td_1[n].local_read_bytes_mt       = &local_read_bytes_t[n];
                        td_1[n].proc_num_mt               = proc_num;
                        td_1[n].rank_num_mt               = rank_num;
                        td_1[n].thread_num_mt             = n;
                        td_1[n].previous_read_num         = 0;
                        pthread_create(&threads_1[n], &attr, find_reads_size_and_offsets_mt, (void *)(&td_1[n]));

                }

                total_num_reads = 0;
                for (n = 0; n < NUM_THREADS_1; n++){
                        pthread_join(threads_1[n], (void *)(&td_1[n]));
                        total_num_reads += *(td_1[n].total_num_reads_mt);
                }

		local_read_offsets  = calloc(total_num_reads, sizeof(size_t));
                local_read_size     = calloc(total_num_reads, sizeof(int));
                local_read_bytes    = calloc(total_num_reads, sizeof(size_t));

                assert(local_read_offsets);
                assert(local_read_size);
                assert(local_read_bytes);

                size_t tmp_var = 0;
                for (n = 0; n < NUM_THREADS_1; n++){
                        td_1[n].local_read_offsets     = local_read_offsets;
                        td_1[n].local_read_size        = local_read_size;
                        td_1[n].local_read_bytes       = local_read_bytes;
                        td_1[n].previous_read_num      = tmp_var;
                        tmp_var                        += *(td_1[n].total_num_reads_mt);
			pthread_create(&threads_1[n], &attr, copy_local_read_info_mt, (void *)(&td_1[n]));

                }

                for (n = 0; n < NUM_THREADS_1; n++){
                        pthread_join(threads_1[n], (void *)(&td_1[n]));

                        free(local_read_offsets_t[n]);
                        free(local_read_bytes_t[n]);
                        free(local_read_size_t[n]);
                }

                free(local_num_reads_t);
                free(total_num_reads_t);

                free(local_read_offsets_t);
                free(local_read_bytes_t);
                free(local_read_size_t);

                pthread_attr_destroy(&attr);
                free(td_1);

		aft = MPI_Wtime();
		fprintf(stderr, "%s: rank %d num reads parsed: %zu ::: time spend reading and parsing entire buffer = (%.02f) \n", __func__, rank_num, total_num_reads, aft - bef);

		if (goff) free(goff);

		//communications about the first file 
		res = MPI_Reduce(&total_num_reads, &grand_total_num_reads, 1, MPI_LONG_LONG_INT, MPI_SUM, 0,MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);	

		res = MPI_Bcast(&grand_total_num_reads, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);	

		//we dispatch the previous result
		if ( rank_num == 0 || rank_num ==1 )
			fprintf(stderr, "rank %d ::: total_num_reads in the file = %zu \n", rank_num, grand_total_num_reads);

		MPI_Barrier(MPI_COMM_WORLD);

			
		///Now find the required information to create the chunks, such as total number of bases
		local_num_reads 	= total_num_reads;
		
		//change chunck_size if necessary
		//opt->chunk_size=10000000;

		size_t i;
		//since all the sizes are not the same anymore, you have to loop

		size_t min_num_read = local_num_reads;
		size_t bases_tmp 	= 0;
		size_t chunck_num 	= 0;

		for(i=0; i < min_num_read; i++)   {

			bases_tmp  += (local_read_size[i]);

			if ( bases_tmp > fixed_chunk_size ){

				bases_tmp = 0;
				chunck_num++;
			}
		}

		chunck_num 	+= 2;
			
		fprintf(stderr,"Rank %d :: chunk_num = %zu\n", rank_num, chunck_num);
		
		///assert that the sizes and offsets found earlier are real.
		//changed the comparisons here for the size. they were compared to blen which has become impossible
		size_t h=0;
		/*
		for (h=0; h < total_num_reads; h++)
		{
			assert(local_read_size[h] > 0 );
			assert(local_read_offsets[h] >= 0 );
		}
		for(h=0; h < total_num_reads_2; h++)
		{
			assert(local_read_size_2[h] > 0 );
			assert(local_read_offsets_2[h] >= 0 );
		}
		*/
		///Now that we know how many chunks we need, we create the vectors associated
		///and we update their informations according to the file they are treating

		// we allocate vector for chunks offset
		size_t *begin_offset_chunk 	= calloc(chunck_num, sizeof(size_t));
		size_t *chunk_size 		= calloc(chunck_num, sizeof(size_t));
		size_t *reads_in_chunk 		= calloc(chunck_num, sizeof(size_t));

		assert( begin_offset_chunk != NULL );
		assert( chunk_size != NULL );
		assert( reads_in_chunk != NULL );
	
		size_t chunk_count = 0;

		maxsiz = fixed_chunk_size; 
		MPI_Barrier(MPI_COMM_WORLD);
		fprintf(stderr,"rank %d ::: Call find_chunks_info \n", rank_num);
		// the detail of he paramters is at the function definition
		// fprintf(stderr, "rank %d ::: begin_offset_chunk = %zu \n", rank_num, begin_offset_chunk[0]);
		// fprintf(stderr, "rank %d ::: begin_offset_chunk_2 = %zu \n", rank_num, begin_offset_chunk_2[0]);
		// fprintf(stderr, "rank %d ::: chunk_size = %zu \n", rank_num, chunk_size);
		// fprintf(stderr, "rank %d ::: chunk_size_2 = %zu \n", rank_num, chunk_size_2);

		bef = MPI_Wtime();
		find_chunks_info(begin_offset_chunk,
				 chunk_size,
				 reads_in_chunk,
				 local_read_size,
				 local_read_bytes,
				 local_read_offsets,
				 rank_num,
				 proc_num,
				 local_num_reads,
				 grand_total_num_reads,
				 maxsiz,
				 &chunk_count);

		aft = MPI_Wtime();
		fprintf(stderr, "%s ::: rank %d ::: evaluating offsets chuncks and sizes (%.02f) found %zu chuncks \n", __func__, rank_num, aft - bef, chunk_count);
	
		MPI_Barrier(MPI_COMM_WORLD);

		if (local_read_offsets) 	free(local_read_offsets);
		if (local_read_size) 		free(local_read_size);
		if (local_read_bytes)		free(local_read_bytes);	
		//verify the number of reads

		size_t num_reads_1 = 0;
		size_t total_num_reads_v1 = 0;
		
		for(h=0; h < chunk_count; h++){
			num_reads_1 += reads_in_chunk[h];
		}

		fprintf(stderr, "rank %d ::: after find_chunks_info ::: num_reads_1 = %zu \n", rank_num, num_reads_1);
		
		res = MPI_Reduce(&num_reads_1, &total_num_reads_v1, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);

		res = MPI_Bcast(&total_num_reads_v1, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);

		if ( rank_num == 0 || rank_num ==1 )
			fprintf(stderr, "rank %d ::: after find_chunks_info ::: total_num_reads forward = %zu \n", rank_num, total_num_reads_v1);

		MPI_Barrier(MPI_COMM_WORLD);

		//display test to check that each vector was well filled with the right information
		//fprintf(stderr,"rank %d ::: begin offset %zu, %zu\n", rank_num, begin_offset_chunk[0], begin_offset_chunk_2[0]);
		//fprintf(stderr,"rank %d ::: chunk size %zu, %zu\n", rank_num, chunk_size[0], chunk_size_2[0]);
		//fprintf(stderr,"rank %d ::: reads in chunk %zu, %zu\n", rank_num, reads_in_chunk[0], reads_in_chunk_2[0]);

		bef = MPI_Wtime();
                size_t total_chunks = 0;

                res = MPI_Reduce(&chunk_count, &total_chunks, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
                assert ( res == MPI_SUCCESS);
                res = MPI_Bcast(&total_chunks, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
                assert ( res == MPI_SUCCESS);

                size_t chunks_per_rank[proc_num];
                res = MPI_Allgather(&chunk_count, 1, MPI_LONG_LONG_INT, chunks_per_rank, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
                assert( res == MPI_SUCCESS);

                all_begin_offset_chunk = malloc(total_chunks * sizeof(size_t));
                all_begin_offset_chunk[0]=0;
                all_chunk_size = malloc(total_chunks * sizeof(size_t));
                all_chunk_size[0] = 0;
                all_reads_in_chunk =  malloc(total_chunks * sizeof(size_t));
                all_reads_in_chunk[0] = 0;

                int displ_chunk[proc_num];
                displ_chunk[0] = 0;

                for (i = 1; i < proc_num; i++) displ_chunk[i] = (displ_chunk[i-1] + chunks_per_rank[i-1]);

                int indx=displ_chunk[rank_num];
                for (i = 0; i <  chunks_per_rank[rank_num]; i++) {
                        all_chunk_size[indx+i]=chunk_size[i];
                        all_begin_offset_chunk[indx+i]=begin_offset_chunk[i];
                        all_reads_in_chunk[indx+i]=reads_in_chunk[i];
                }

                if (rank_num > 0){
                        res=MPI_Send(chunk_size, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
                        assert(res == MPI_SUCCESS);
                }
                else{
                for (i = 1; i < proc_num; i++){
                        res=MPI_Recv(&(all_chunk_size[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        assert(res == MPI_SUCCESS);
                        }
                }

		if (rank_num > 0){
                        res=MPI_Send(begin_offset_chunk, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
                        assert(res == MPI_SUCCESS);
                }
                else{
                for (i = 1; i < proc_num; i++){
                        res=MPI_Recv(&(all_begin_offset_chunk[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        assert(res == MPI_SUCCESS);
                        }
                }

                if (rank_num > 0){
                        res=MPI_Send(reads_in_chunk, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
                        assert(res == MPI_SUCCESS);
                }
                else{
                        for (i = 1; i < proc_num; i++){
                                res=MPI_Recv(&(all_reads_in_chunk[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                assert(res == MPI_SUCCESS);
                        }
                }

                res=MPI_Bcast(all_chunk_size, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
                res=MPI_Bcast(all_reads_in_chunk, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
                res=MPI_Bcast(all_begin_offset_chunk, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

                assert(res == MPI_SUCCESS);

                size_t total = 0;

                for ( i = 0; i < total_chunks; i++){
                        assert(all_chunk_size[i] != 0);
                        total += all_reads_in_chunk[i];
                }

                assert (total_num_reads_v1 == total);

                free(chunk_size);
                free(begin_offset_chunk);
                free(reads_in_chunk);

                aft = MPI_Wtime();
		fprintf(stderr, "%s: rank %d time spend gathering chunks info = (%.02f) \n", __func__, rank_num, aft - bef);
		

		/// Map reference genome indexes in shared memory (by host)
		bef = MPI_Wtime();
		map_indexes(file_map, &count, &indix, &ignore_alt, &win_shr, shared_mem);
		aft = MPI_Wtime();
		fprintf(stderr, "%s ::: rank %d ::: mapped indexes (%.02f)\n", __func__, rank_num, aft - bef);

		//we create a vector with chromosom names 
		int s;
		
		MPI_File *fh_out          = malloc( (indix.bns->n_seqs + 2) * sizeof(MPI_File));
		char *files_out_sam_name[indix.bns->n_seqs + 1];
		char *file_map_by_chr[(indix.bns->n_seqs + 1)];
		for (s = 0; s < indix.bns->n_seqs; ++s){
			files_out_sam_name[s] = malloc( (strlen( indix.bns->anns[s].name) + 1 ) * sizeof(char));
		        files_out_sam_name[s][strlen( indix.bns->anns[s].name)] = 0;
		        char *p0 = files_out_sam_name[s];
		        memmove(p0, indix.bns->anns[s].name, strlen( indix.bns->anns[s].name));
		}
		char UNMAPPED[]   = "unmapped";
		files_out_sam_name[s++] = strdup(UNMAPPED);
		int file_name_len = 0;
	        for (s = 0; s < (indix.bns->n_seqs + 1); ++s){
	        	/* Derived file names */
	        	file_name_len = strlen(output_path) + strlen(files_out_sam_name[s]) + 6;
	        	file_map_by_chr[s] = calloc( file_name_len, sizeof(char));
			if ( write_format == 2  )
                                sprintf(file_map_by_chr[s], "%s/%s.sam", output_path, files_out_sam_name[s]);
                        if ( write_format == 1  )
                                sprintf(file_map_by_chr[s], "%s/%s.bam", output_path, files_out_sam_name[s]);
                        if ( write_format == 0  )
                                sprintf(file_map_by_chr[s], "%s/%s.gz", output_path, files_out_sam_name[s]);
		}
		
		if ( write_format == 2  )
                        create_sam_header_by_chr_file(file_map_by_chr, &indix, &count, hdr_line, rg_line, pg_line, rank_num);

                if ( write_format == 1 || write_format == 0)
                        create_bam_header_by_chr_file(file_map_by_chr, &indix, &count, hdr_line, rg_line, pg_line, rank_num, compression_level, dofixmate);

		//create_sam_header_by_chr_file(file_map_by_chr, &indix, &count, hdr_line, rg_line, rank_num);
		
		//we create a file which contains only the header
		//this file will be used by the sorting to get individual chromosom file
		//create_sam_header(file_out, &indix, &count, hdr_line, rg_line, rank_num);


		///This is a testline to stop the program wherever I want
		//if(0){ MPI_Barrier(MPI_COMM_WORLD); MPI_Finalize(); return 0;}

		for (s = 0; s < (indix.bns->n_seqs + 1); ++s){
                        res = MPI_File_open(MPI_COMM_WORLD, file_map_by_chr[s], MPI_MODE_CREATE|MPI_MODE_WRONLY|MPI_MODE_APPEND, MPI_INFO_NULL, &fh_out[s]);
                        assert(res == MPI_SUCCESS);
                }

		if (file_r1 != NULL) {
			//fprintf(stderr, "rank %d ::: open file %s \n",rank_num, file_r1);
			res = MPI_File_open(MPI_COMM_WORLD, file_r1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_r1);
			assert(res == MPI_SUCCESS);
		}
		
		buffer_r1 = NULL; seqs = NULL;
		//localsize_vec contain the size of the buffer to write in the sam file
		int *chr_buff_size  = calloc ( (indix.bns->n_seqs + 1), sizeof(int) );

		if ( proc_num > 1) {
			bef = MPI_Wtime();

                	MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank_target, 0, win);
                	MPI_Request req;
                	MPI_Rget(&index_chunk, 1, MPI_UINT64_T, rank_target, 0, 1, MPI_UINT64_T, win, &req);
                	MPI_Wait(&req, MPI_STATUS_IGNORE);
                	u1 = (uint64_t)index_chunk;
                	MPI_Fetch_and_op(&incr, &index_chunk, MPI_UINT64_T, rank_target, 0, MPI_SUM, win);
                	MPI_Win_flush(rank_target, win);
                	MPI_Win_unlock(rank_target, win);

                	aft = MPI_Wtime();
                	fprintf(stderr, "rank %d ::: initial u1 = %zu :: time %.02f \n",rank_num, u1, aft-bef);

			before_local_mapping = MPI_Wtime();
		}
		// here we loop until there's nothing to read
		//we loop the chunck_count
	 
		while ( u1 < total_chunks ){

			offset_chunk 		= all_begin_offset_chunk[u1];
			size_chunk   		= all_chunk_size[u1];
			
			assert(size_chunk 	!= 0);
					
			bef = MPI_Wtime();
			///allocate the buffer sizes to store an entire chunk in it
			buffer_r1 = malloc(size_chunk + 1);
			assert(buffer_r1 != NULL);
			buffer_r1[size_chunk] = 0;

			struct struct_pread_fastq *td_pread1;
                        td_pread1 = malloc (NUM_THREADS_1 * sizeof(struct struct_pread_fastq));
                        bef = MPI_Wtime();
                        pthread_attr_t attr4;
                        pthread_attr_init(&attr4);
                        pthread_attr_setstacksize(&attr4, BIG_STACK);
                        pthread_attr_setdetachstate(&attr4, 0);

                        for( n = 0; n < NUM_THREADS_1; n++ ){
                                td_pread1[n].total_thread = NUM_THREADS_1;
                                td_pread1[n].thread_id = n;
                                td_pread1[n].job_rank = rank_num;
                                td_pread1[n].offset= offset_chunk;
                                td_pread1[n].size = size_chunk;
                                td_pread1[n].buffer = buffer_r1;
                                td_pread1[n].fd  = fh_r1;
                                int ret_code = pthread_create(&threads[n], &attr4, pread_fastq_chunck, (void *)(&td_pread1[n]));
                                if (ret_code) {
                                        fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                                }
                        }
                        for(n=0; n < NUM_THREADS_1; n++)
                                pthread_join(threads[n], (void *)(&td_pread1[n]));

                        pthread_attr_destroy(&attr4);
                        free(td_pread1);

                        assert(strlen(buffer_r1) == size_chunk);
                        assert(*buffer_r1 == '@');

			///some stats
			aft = MPI_Wtime();
			fprintf(stderr, "%s: num_rank = %d :: read sequences (%.02f)\n", __func__, rank_num, aft - bef);
			total_time_reading_seq += (aft - bef);

			///operations on the number of reads
			bef = MPI_Wtime();
			reads_r1 = all_reads_in_chunk[u1];
			reads = reads_r1;
			bases = 0;
			total_local_reads_aligned += reads;
			assert(reads <= INT_MAX);

			/* Parse sequences ... */
			seqs = malloc(reads * sizeof(*seqs));
			assert(seqs != NULL);

			if (file_r1 != NULL) {
				p1 = q1 = buffer_r1; e1 = buffer_r1 + size_chunk; line_number = 0;
							
				while (q1 < e1) {
					if (*q1 != '\n') { q1++; continue; }
					/* We have a full line ... process it */
					*q1 = '\0'; 
					//n will take 4 times in row the same value from 0 with a step of 2
					//it is because the 2 operations are using int types, so it is cast as int twice (within and then outside the brackets)
					n = files * (line_number / 4);
				
					switch (line_number % 4) {
					case 0: /* Line1: Name and Comment */
						assert(*p1 == '@');
						seqs[n].name   = p1 + 1;
						while (*p1 && !isspace((unsigned char)*p1)) {p1++;}
						if (*(p1-2) == '/' && isdigit((unsigned char)*(p1-1))) {*(p1-2) = '\0';}
						if (*p1) {*p1++ = '\0';}
						seqs[n].comment = (copy_comment != 0) ? p1 : NULL;
						seqs[n].sam = NULL;
						break;
					case 1: /* Line2: Sequence */
						seqs[n].seq = p1;
						seqs[n].l_seq = q1 - p1;
						bases += seqs[n].l_seq;
						break;
					case 2: /* Line3: Ignored */
						assert(*p1 == '+');
						break;
					case 3: /* Line4: Quality */
						seqs[n].qual = p1;
						break; }
					p1 = ++q1; 
					line_number++; 
				}
			}
			assert( (line_number/4) == reads_r1);
					
			aft = MPI_Wtime();

			fprintf(stderr, "rank %d :: %s: parsed %zu read of sequences (%ld bp) in (%.02f)\n", rank_num, __func__, reads, (long)bases, aft - bef);
			total_time_parsing += (aft - bef);
			//fprintf(stderr, "rank %d ::: [M::%s] read %zu sequences (%ld bp)...\n",rank_num , __func__, reads, (long)bases);
			/* Datas computation ... */
			bef = MPI_Wtime();
			fprintf(stderr, "rank %d ::::  Call memseqs with count sequences %zu \n", rank_num, reads);
			mem_process_seqs(opt, indix.bwt, indix.bns, indix.pac, 0, (int)reads, seqs, pes0);
			aft = MPI_Wtime();
			fprintf(stderr, "rank %d :: %s: computed mappings (%.02f)\n", rank_num, __func__, aft - bef);
			
			total_time_mapping += (aft - bef);

			//first we cont how many sam line we have
			size_t total_sam_line = 0;
			n = 0;
			char currentLine[MAX_CHAR_SIZE];
			int next;

			while ( n < reads){
                        	next = tokenizer(seqs[n].sam,'\n', currentLine);
                                while (next) {
                                	total_sam_line++;
                                        next = tokenizer(NULL, '\n', currentLine);
                                }
				//for(i=0; i<MAX_CHAR_SIZE; i++) current_line[i]=0;
                                n = n + 1;
                        }

			/* Write results ... */
			bef = MPI_Wtime();
                        
                        int *sam_buff_dest  	= calloc ( total_sam_line, sizeof(int) );
			char **start_addr       = malloc ( total_sam_line * sizeof(char*));
			int *line_size_to_cpy   = malloc ( total_sam_line * sizeof(int));

                        char *start_sam_line;
                        char *currentCarac;
                        int chr;
                        size_t sam_line_size;
                        int nbchr = indix.bns->n_seqs + 1;
			size_t incr_line =0;
			char *tmp_chr = malloc( 200 * sizeof(char));
                        tmp_chr[0] = 0;
			
                        for (n = 0; n < reads; n++) {
				int total_sam_line_size = 0;
				next = tokenizer(seqs[n].sam,'\n', currentLine);
                                currentCarac = currentLine;
                                start_sam_line = seqs[n].sam;

				while (next){
				
					start_sam_line = seqs[n].sam + total_sam_line_size;
                                        sam_line_size = strlen(currentLine) + 1;

					currentLine[sam_line_size - 1] = '\n';
                                        currentLine[sam_line_size] = '\0';
                                        currentCarac = currentLine;
					//READ NAME
					currentCarac = strstr (currentCarac, "\t");
					currentCarac++;
					
                                	//GO TO FLAG
					currentCarac = strstr(currentLine, "\t");                               
					currentCarac++;
					//GO TO RNAME (Chr name)
					currentCarac = strstr(currentCarac + 1, "\t");
                                	if ( currentCarac[1] == '*') chr =  nbchr-1;
                                	else chr = getChr(currentCarac, files_out_sam_name, nbchr-1, tmp_chr);

                                	if ((chr < (nbchr - 1))){
						//then we found concordant reads                                
				        	chr_buff_size[chr]    		+= sam_line_size;
				        	sam_buff_dest[incr_line]       	= chr;
						start_addr[incr_line]           = start_sam_line;
					}
					else if ((chr == (nbchr - 1))){
                                        	//we found unmapped reads
                                        	chr_buff_size[nbchr - 1] 	+= sam_line_size;
                                        	sam_buff_dest[incr_line]        = nbchr - 1;
						start_addr[incr_line]           = start_sam_line;
                                	}
                                	else{
                                        	//we never go here
                                        	chr_buff_size[nbchr - 1] 	+= sam_line_size;
                                        	sam_buff_dest[incr_line]        = nbchr - 1;
                                	}

					line_size_to_cpy[incr_line] = sam_line_size;
                                        total_sam_line_size += sam_line_size;
                                        incr_line++;
                                        next = tokenizer(NULL, '\n', currentLine);
				}
			}
                        free(tmp_chr);                                                                                                                                        
			//now we fill up the buffer_out_vec
			if ( write_format == 2 ){

				bef = MPI_Wtime();

                                int ret_code = 0;
                                struct thread_data_compress_by_chr_single *tdc;
                                tdc = malloc (NUM_THREADS * sizeof(struct thread_data_compress_by_chr_single));

                                pthread_attr_init(&attr);
                                pthread_attr_setstacksize(&attr, SMALL_STACK);
                                pthread_attr_setdetachstate(&attr, 0);

                                for( n = 0; n < NUM_THREADS; n++ ){
                                        tdc[n].total_thread = NUM_THREADS;
                                        tdc[n].thread_id = n;
                                        tdc[n].indix_thr = &indix;
                                        tdc[n].start_addr = start_addr;
                                        tdc[n].line_size_to_cpy = line_size_to_cpy;
                                        tdc[n].sam_buff_dest = sam_buff_dest;
                                        tdc[n].total_sam_line = total_sam_line;
                                        tdc[n].incrmnt = 1;
                                        tdc[n].job_rank = rank_num;
                                        tdc[n].comp_level = compression_level;
                                        tdc[n].fh_out = fh_out;
                                        ret_code = pthread_create(&threads[n], &attr, write_sam_thread_by_chr_single, (void *)(&tdc[n]));
                                        if (ret_code) {
                                                fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                                        }
                                }
                                for(n=0; n<NUM_THREADS; n++)
                                        pthread_join(threads[n], (void *)(&tdc[n]));

                                pthread_attr_destroy(&attr);
                                free(tdc);

                                for (n = 0; n < reads; n++) free(seqs[n].sam);

                                free(sam_buff_dest);
                                free(seqs);
                                free(start_addr);
                                free(line_size_to_cpy);
                                aft = MPI_Wtime();
                                fprintf(stderr, "rank: %d :: %s: wrote results (%.02f) \n", rank_num, __func__, aft - bef);
                                total_time_writing += (aft - bef);



			}
			
			if ( write_format == 0 || write_format == 1) {

                                bef = MPI_Wtime();

                                int ret_code = 0;
                                struct thread_data_compress_by_chr_single *tdc;
                                tdc = malloc (NUM_THREADS * sizeof(struct thread_data_compress_by_chr_single));
                                
                                pthread_attr_init(&attr);
                                pthread_attr_setstacksize(&attr, SMALL_STACK);
                                pthread_attr_setdetachstate(&attr, 0);

                                for( n = 0; n < NUM_THREADS; n++ ){
                                        tdc[n].total_thread = NUM_THREADS;
                                        tdc[n].thread_id = n;
                                        tdc[n].indix_thr = &indix;
                                        tdc[n].start_addr = start_addr;
                                        tdc[n].line_size_to_cpy = line_size_to_cpy;
                                        tdc[n].sam_buff_dest = sam_buff_dest;
                                        tdc[n].total_sam_line = total_sam_line;
                                        tdc[n].incrmnt = 1;
                                        tdc[n].job_rank = rank_num;
                                        tdc[n].comp_level = compression_level;
                                        tdc[n].fh_out = fh_out;
                                        ret_code = pthread_create(&threads[n], &attr, compress_thread_by_chr_single, (void *)(&tdc[n]));
                                        if (ret_code) {
                                                fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                                        }
                                }
                                for(n=0; n<NUM_THREADS; n++)
                                        pthread_join(threads[n], (void *)(&tdc[n]));

                                pthread_attr_destroy(&attr);
                                free(tdc);
				                            
                                for (n = 0; n < reads; n++) free(seqs[n].sam);
                                
                                free(sam_buff_dest);
                                free(seqs);
                                free(start_addr);
                                free(line_size_to_cpy);
                                aft = MPI_Wtime();
                                fprintf(stderr, "rank: %d :: %s: wrote results (%.02f) \n", rank_num, __func__, aft - bef);
                                total_time_writing += (aft - bef);

                        }
                        free(buffer_r1);
                        fprintf(stderr, "rank: %d :: finish for chunck %zu \n", rank_num, u1);
			

			if ( proc_num > 1) {
                        	bef = MPI_Wtime();
                        	MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank_target, 0, win);
                        	MPI_Request req;
                        	MPI_Rget(&index_chunk, 1, MPI_UINT64_T, rank_target, 0, 1, MPI_UINT64_T, win, &req);
                        	MPI_Wait(&req, MPI_STATUS_IGNORE);
                        	u1 = (uint64_t)index_chunk;
                        	MPI_Fetch_and_op(&incr,&index_chunk, MPI_UINT64_T, rank_target, 0, MPI_SUM, win);
                        	MPI_Win_flush(rank_target, win);
                        	MPI_Win_unlock(rank_target, win);
			
                        	aft = MPI_Wtime();
                        	fprintf(stderr, "rank %d ::: update u1 = %zu :: time %.02f \n",rank_num, u1, aft-bef);
			}
			else u1++;
                } //end for (u1 = 0; u1 < chunk_count; u1++){

		MPI_Barrier(MPI_COMM_WORLD);
                if ( (write_format == 1) && (rank_num == 0)) {
                	for (n = 0; n < (indix.bns->n_seqs + 1); n++){
                        	static uint8_t magic[28] =  "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00";
                                res = MPI_File_write_shared(fh_out[n], magic, 28, MPI_BYTE, &status);
                                assert(res == MPI_SUCCESS);
                                res = MPI_Get_count(&status, MPI_BYTE, &count);
                                assert(res == MPI_SUCCESS);
                                assert(count == 28);
                         }
                }
	
		//free data structures and close sam files		
		for (n = 0; n < (indix.bns->n_seqs + 1); n++)  {
                	free(files_out_sam_name[n]);
                	free(file_map_by_chr[n]);
               		res = MPI_File_close(&fh_out[n]);
                	assert(res == MPI_SUCCESS);
        	}
	}

	/*
	*
	*   Print some statistics
	*
	*/

	if (hdr_line) free(hdr_line);
        if (rg_line) free(rg_line);
        if (pg_line) free(pg_line);
	
	if ( proc_num > 1) MPI_Win_free(&win);

	free(output_path);
	after_local_mapping	 = MPI_Wtime();
	total_time_local_mapping = after_local_mapping - before_local_mapping;

	res = MPI_Allreduce(&total_time_local_mapping, &total_time_mapping, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	assert ( res == MPI_SUCCESS);

	res = MPI_Allreduce(&total_local_reads_aligned, &total_reads_check, 1,MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
	assert ( res == MPI_SUCCESS);

	res = MPI_Allreduce(&total_time_reading_seq, &grand_total_time_reading_seq, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	assert ( res == MPI_SUCCESS);

	res = MPI_Allreduce(&total_time_parsing, &grand_total_time_parsing, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	assert ( res == MPI_SUCCESS);

	res = MPI_Allreduce(&total_time_writing, &grand_total_time_writing, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	assert ( res == MPI_SUCCESS);
	
	fprintf(stderr, "rank %d :::: total_time_reading_seq = %.02f seconds \n", rank_num, grand_total_time_reading_seq);
	fprintf(stderr, "rank %d :::: total_time_writing = %.02f seconds \n", rank_num, grand_total_time_writing);
	fprintf(stderr, "rank %d :::: total_time_parsing = %.02f  seconds \n", rank_num, grand_total_time_parsing);
	fprintf(stderr, "rank %d :::: total_time_local_mapping = %.02f seconds \n", rank_num, total_time_local_mapping);
	fprintf(stderr, "rank %d :::: total_time_mapping = %.02f seconds \n", rank_num, total_time_mapping);
	fprintf(stderr, "rank %d :::: total_reads_check for all ranks= %zu \n", rank_num, total_reads_check);
	fprintf(stderr, "rank %d :::: total_local_reads_aligned = %zu \n", rank_num, total_local_reads_aligned);
	

	fprintf(stderr, "rank %d :::: finish mappings for reads \n", rank_num);
	if (all_begin_offset_chunk != NULL) free(all_begin_offset_chunk);
	if (all_chunk_size != NULL) free(all_chunk_size);
	if (all_reads_in_chunk != NULL) free(all_reads_in_chunk);

	if (all_begin_offset_chunk_2 != NULL) free(all_begin_offset_chunk_2);
        if (all_chunk_size_2 != NULL) free(all_chunk_size_2);
        if (all_reads_in_chunk_2 != NULL) free(all_reads_in_chunk_2);

	if (opt != NULL) free(opt);

	if (file_r2 != NULL) {
		res = MPI_File_close(&fh_r2);
		assert(res == MPI_SUCCESS);
	}
	if (file_r1 != NULL) {
		res = MPI_File_close(&fh_r1);
		assert(res == MPI_SUCCESS);
	}

	res = MPI_Win_free(&win_shr);
	assert(res == MPI_SUCCESS);
	
	//res = MPI_File_close(&fh_out);
	//assert(res == MPI_SUCCESS);

	res = MPI_Finalize();
	assert(res == MPI_SUCCESS);
	
	return 0;
}

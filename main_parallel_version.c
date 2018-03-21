/*
This file is part of mpiBWA

NGS aligner inspired by BWA

The project was developped by Frederic Jarlier from Institut Curie and Nicolas Joly from Institut Pasteur

Copyright (C) 2016-2017  Institut Curie / Institut Pasteur



This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

#include "bwa.h"
#include "bwamem.h"
#include "utils.h"

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


void init_goff(size_t *goff, MPI_File mpi_filed, size_t fsize,int numproc,int rank){


	char * current_line = NULL;
	MPI_Status status;
	int i = 0;
	int j = 0;
	// TODO problem here when we have only one job!
	size_t lsize = fsize/numproc;
	goff[0]=0;
	for(i=1;i<numproc;i++){goff[i]=lsize*i;}
	goff[numproc]=fsize;
	for(i=1;i<numproc;i++){
		current_line =(char*)calloc(2000,sizeof(char));
		//we read 2K caracters
		MPI_File_read_at(mpi_filed, (MPI_Offset)goff[i], current_line, 2000, MPI_CHAR, &status);
		j=0;
		while (j < 2000 ){
			//we check the \n after and before
			j++;
			if (current_line[j] == '+' && current_line[j+1] == '\n' && current_line[j-1] == '\n') {j++;break;}}
		//then we look for return after quality
		while( j < 2000){j++;if (current_line[j] == '\n') {break;}}
		goff[i]+=(j+1); free(current_line);}

	return;
}

int main(int argc, char *argv[]) {
	double bef, aft;

	char *file_r1 = NULL, *file_r2 = NULL;
	struct stat stat_r1, stat_r2;
	char *buffer_r1, *buffer_r2;
	int fd_in1, fd_in2;
	uint8_t *a, *addr, *addr_map;

	char *file_out = NULL;
	char *buffer_out;

	char *file_ref = NULL;

	char file_map[PATH_MAX], file_tmp[PATH_MAX];

	int proc_num, rank_num, rank_shr;
	size_t localsize;
	size_t n = 0;

	MPI_Aint size_shr;
	MPI_Comm comm_shr;
	MPI_File fh_map, fh_r1, fh_r2, fh_out, fh_tmp;
	MPI_Offset *coff, m, size_map, size_tot;
	MPI_Status status;
	MPI_Win win_shr;

	mem_opt_t *opt, opt0;
	mem_pestat_t pes[4], *pes0 = NULL;
	bwaidx_t indix;
	bseq1_t *seqs;

	int res, count;
	int files, nargs;

	char *p, *q, *s, *e;
	off_t locoff, locsiz, *alloff, *curoff, maxsiz, totsiz, filsiz;

	int c, copy_comment = 0;
	char *rg_line = NULL, *hdr_line = NULL;
	int ignore_alt = 0;
	const char *mode = NULL;

	char *progname = basename(argv[0]);

	if (argc < 2) {
		fprintf(stderr, "Program: MPI version of BWA MEM\n\n"
			"Version: v%s\n\n"
			"Contact 1: Frederic Jarlier (frederic.jarlier@curie.fr)\n\n"
			"usage : mpirun -n TOTAL_PROC %s mem -t 8 -o RESULTS REFERENCE FASTQ_R1 FASTQ_R2\n\n"
			"Requirements : After the creation of reference file with BWA you need to create a referenced\n"
			"	   map file genome with pidx like this  pidx ref.fasta generate a ref.fasta.map.\n"
			"	   The .map file is a copy of memory mapped reference used for shared memory purpose.\n",
			VERSION, progname);
		return 1; }

	/* Validate provided command (first argument) */
	if (strcmp(argv[1], "mem") != 0) {
		fprintf(stderr, "%s: unsupported %s command\n", progname, argv[1]);
		return 1; }

	/* initialize the BWA-MEM parameters to the default values */
	opt = mem_opt_init();
	memset(&opt0, 0, sizeof(opt0));
	while ((c = getopt(argc-1, argv+1, "1paMCSPVYjk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:W:x:G:h:y:K:X:H:o:")) >= 0) {
		if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
		else if (c == '1') ; /* FIXME: unsupported */
		else if (c == 'x') mode = optarg;
		else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;
		else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1;
		else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1;
		else if (c == 'O') {
			opt0.o_del = opt0.o_ins = 1;
			opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->o_ins = strtol(p+1, &p, 10);
		}
		else if (c == 'E') {
			opt0.e_del = opt0.e_ins = 1;
			opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->e_ins = strtol(p+1, &p, 10);
		}
		else if (c == 'L') {
			opt0.pen_clip5 = opt0.pen_clip3 = 1;
			opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->pen_clip3 = strtol(p+1, &p, 10);
		}
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
		else if (c == 'K') ; /* FIXME: unsupported */
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
	res = MPI_Init(&argc, &argv);
	assert(res == MPI_SUCCESS);
	res = MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	assert(res == MPI_SUCCESS);
	res = MPI_Comm_rank(MPI_COMM_WORLD, &rank_num);
	assert(res == MPI_SUCCESS);

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
	if (file_r2 != NULL && stat_r1.st_size != stat_r2.st_size) {
		fprintf(stderr, "file sizes for R1 and R2 do not match\n");
		//res = MPI_Finalize();
		//assert(res == MPI_SUCCESS);
		//exit(2);
	}

	/* Work around build warning in non timing case */
	aft = 0; aft++;
	bef = 0; bef++;
	/*
	 * Rank 0 estimate the size of a read
	 *
	 */
	size_t slen, blen;
	off_t tmp_sz = 1024;
	if (rank_num == 0){
		// 512 Mo are
		int fd_tmp = open(file_r1, O_RDONLY, 0666);
		char *buffer = malloc(tmp_sz  + 1);
		buffer[tmp_sz] = '\0';
		read(fd_tmp, buffer, tmp_sz);
		assert(strlen(buffer) == tmp_sz);
		assert( *buffer == '@');
		/* Estimate sequence size (bases and bytes) */
		/* With this estimation we compute approximately the chuncks size and offset*/

		s = buffer;
		e = buffer + tmp_sz;
		p = q = s;
		while (q < e && *q != '\n') q++; p = ++q;
		while (q < e && *q != '\n') q++; blen = q - p; p = ++q;
		while (q < e && *q != '\n') q++; p = ++q;
		while (q < e && *q != '\n') q++; slen = q - s + 1;

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

	size_t read_number=0;
	assert(fd_in1 != -1);
	size_t *goff = NULL; //global offset contain the start offset in the fastq
	goff = malloc((proc_num + 1) * sizeof(size_t));
	char *current_line = NULL;

	int i = 0;
	int j = 0;
	// TODO problem here when we have only one job!
	size_t lsize = stat_r1.st_size / proc_num;
	goff[0]=0;
	for(i = 1 ; i < proc_num; i++){goff[i] = lsize*i;}
	goff[proc_num] = stat_r1.st_size;
	char *buffer_r0 = malloc( tmp_sz + 1);
	buffer_r0[tmp_sz] = '\0';

	res = MPI_File_read_at(mpi_fd_in1, (MPI_Offset)goff[rank_num], buffer_r0, tmp_sz, MPI_CHAR, &status);
	assert(res == MPI_SUCCESS);	
	p = buffer_r0;
	e = buffer_r0 + tmp_sz;

	while (p < e) {
		if (*p != '@') { p++; continue; }
		if (p != buffer_r0 && *(p-1) != '\n') { p++; continue; }
		q = p + 1;
		while (q < e && *q != '\n') q++; q++;
		while (q < e && *q != '\n') q++; q++;
		if (q < e && *q == '+') break;
		p++;
	}
	assert(*p == '@');
	//we update begining offsets of windows in order to start at the begining of a read
	goff[rank_num] += p - buffer_r0;	
	free(buffer_r0);
	res = MPI_File_close(&mpi_fd_in1);
	assert(res == MPI_SUCCESS);
	//now we exchange the goff buffer between all proc
	size_t goff_inter = goff[rank_num]; //avoid memcpy overlap
	//rank 0 gather the vector
	res = MPI_Allgather(&goff_inter, 1, MPI_LONG_LONG_INT, goff , 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
	assert(res == MPI_SUCCESS);

	//we compute the new size according to the shift
	//We calculate the size to read for each process
	int ind = rank_num;
	size_t siz2read = goff[ind+1]-goff[ind];
	MPI_Barrier(MPI_COMM_WORLD);
	j=0;
	size_t off_in_file = goff[ind]; //Current offset in file sam
	
	size_t total_siz_read = 0;
	size_t read_buffer_sz = 0;
	size_t offset_read = 0;
	
	size_t offset_buffer_r1 = 0;
	size_t offset_tmp  = 0;
	size_t count1;

	MPI_Datatype arraytype;
	MPI_Datatype arraytype0;	
	MPI_Datatype arraytype_r1;
	MPI_Datatype arraytype_r2;

	if ( siz2read < DEFAULT_INBUF_SIZE ) read_buffer_sz = siz2read;
	else read_buffer_sz = DEFAULT_INBUF_SIZE;

	size_t local_num_reads = 0;
	size_t total_num_reads = 0;
	size_t *local_read_offsets = calloc(1 , sizeof(size_t));
	int *local_read_size    = calloc(1, sizeof(int));
	assert( local_read_offsets != NULL);
	assert( local_read_size != NULL);


	bef = MPI_Wtime();
	char *b, *r, *t;
	size_t offset_end_buff;
	size_t pos_in_vect_offset = 0;
	size_t lines = 0;
	size_t total_parsing = 0;
	size_t new_offset = 0;
	size_t new_off_in_file =0;
	size_t g=0;
	size_t w=0;
	size_t total_read = 0;
	size_t total_lines = 0;
	while (1){

		buffer_r1 = malloc(read_buffer_sz + 1);
		assert( buffer_r1 != NULL );
		buffer_r1[read_buffer_sz] = 0;		
		
		// FOR TEST //
		// MPI_Type_contiguous(read_buffer_sz, MPI_CHAR, &arraytype0);
		// MPI_Type_commit(&arraytype0);
		// MPI_File_set_view(mpi_fd_in2, (MPI_Offset)off_in_file, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL ) ; 
		// res = MPI_File_read(mpi_fd_in2, b , 1, arraytype0, &status);

		res = MPI_File_read_at(mpi_fd_in2, (MPI_Offset)off_in_file, buffer_r1, read_buffer_sz, MPI_CHAR, MPI_STATUS_IGNORE);
		assert(res == MPI_SUCCESS);
		assert(*buffer_r1 == '@');	
		//we search the last read 
		b = buffer_r1;
		r  = b + read_buffer_sz;

		if ( read_buffer_sz == DEFAULT_INBUF_SIZE){			
			//go to last read			
			while (r-- != b){if (*r == '\n' && *(r+1) == '+') {r--; break;}}					
			while (r-- != b){if (*r == '\n') break;}
			while (r-- != b){if (*r == '@') break;			}			
			r--; //stop at \n
			offset_end_buff = (r - b);
		}
		else
			offset_end_buff = (r - b);

		
		//we count the number of lines
		t = b ;
		e = buffer_r1 + offset_end_buff;
		m = 0;
		lines = 0;		
		while (t++ < e){if (*t == '\n') lines++;}

		local_num_reads =  lines/4;
		total_lines += lines;		
		total_num_reads += local_num_reads;
			
		local_read_size  	= (int *)realloc(local_read_size, sizeof(int) * (total_num_reads + 1));
		local_read_offsets	= (size_t *)realloc(local_read_offsets, sizeof(size_t) * (total_num_reads + 1));

		assert( local_read_offsets != NULL);
		assert( local_read_size != NULL);
		
		local_read_size[total_num_reads] = 0;
		local_read_offsets[total_num_reads] = 0;
		
		//for the first first read
		local_read_offsets[pos_in_vect_offset] = off_in_file;
		local_read_size[pos_in_vect_offset]    = blen;				
		pos_in_vect_offset++;
		g = 0;
		t = buffer_r1;
		size_t w = 0;
		//we move to the next line because we test *(t-1)
		lines = 0;
		
		while (t < e ){
			if (*t == '\n') lines++;
			if (lines == 4) {
				local_read_size[pos_in_vect_offset] = blen;
				local_read_offsets[pos_in_vect_offset] = off_in_file + g + 1; //+1 for \n caracters
				pos_in_vect_offset++;				
				w++; lines = 0;//t+=blen;g+=blen; 
			}			
			t++;g++;
		}
	
		assert(g == offset_end_buff);
		total_parsing 	+= offset_end_buff; 				
		
		if (total_parsing == siz2read) {free(buffer_r1);break;}

		if ((siz2read - total_parsing) < DEFAULT_INBUF_SIZE)
			read_buffer_sz = siz2read - total_parsing;			
		else read_buffer_sz = DEFAULT_INBUF_SIZE;
					
		off_in_file  += offset_end_buff + 1;
	
		free(buffer_r1);				
	}
	fprintf(stderr, "rank %d : in wait\n", rank_num);
	MPI_Barrier(MPI_COMM_WORLD);
	aft = MPI_Wtime();
	fprintf(stderr, "%s: rank %d num reads parsed: %zu ::: time spend reading and parsing entire buffer = (%.02f) \n", __func__, rank_num, total_num_reads, aft - bef);
	lines = 0;	
	if (goff) free(goff);
	

	MPI_File_close(&mpi_fd_in2);	
	size_t grand_total_num_reads = 0;
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
	size_t chunck_num = (local_num_reads * blen) / (( opt->chunk_size * opt->n_threads) / 2);
	chunck_num += 2; //the last chunk hold the remain bases

	size_t h=0;
	for ( h = 0; h < total_num_reads; h++){
		assert( local_read_size[h] == blen );
		assert( local_read_offsets[h] >= 0 );
	}

	/*
	 * We compute number of chunks for forward reads
	 * and send it to
	 */
	bef = MPI_Wtime();

	size_t bases_previous_chunk      = 0;
	size_t offset_previous_chunk 	 = 0;
	size_t size_previous_chunk 	 = 0;

	if (rank_num > 0){

			//we wait the previous rank to send the final size chunk
			//and offset of the chunk
			MPI_Recv(&bases_previous_chunk,
					1,
					MPI_LONG_LONG_INT,
					rank_num - 1,
					0,
					MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);

			MPI_Recv(&offset_previous_chunk,
					1,
					MPI_LONG_LONG_INT,
					rank_num - 1,
					0,
					MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);

			MPI_Recv(&size_previous_chunk,
					1,
					MPI_LONG_LONG_INT,
					rank_num - 1,
					0,
					MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);
		}

	/*
	 * We update information from previous rank
	 */

	// we allocate vector for chunks offset
	size_t *begin_offset_chunk 	= calloc(chunck_num, sizeof(size_t));
	size_t *chunk_size 		= calloc(chunck_num, sizeof(size_t));
	size_t *reads_in_chunk 		= calloc(chunck_num, sizeof(size_t));

	assert( begin_offset_chunk != NULL );
	assert( chunk_size != NULL );
	assert( reads_in_chunk != NULL );
	
	size_t u1=0;
	
	size_t counter_bases  	 = 0;
	size_t bytes_in_chunk 	 = 0;

	maxsiz = ( opt->chunk_size * opt->n_threads) / 2; //half the normal number of bases for bwa mem
	size_t index_in_chunk = 0;

	size_t begin_offset = 0;
	size_t chunk_count = 0;

	if (rank_num == 0){

		for (u1 = 0; u1 < local_num_reads; u1++){
			//we for offsets multiple of maxsize

			if (counter_bases == 0)  begin_offset = local_read_offsets[u1];

			counter_bases 	+= local_read_size[u1];
			bytes_in_chunk  += local_read_offsets[u1+1] - local_read_offsets[u1];

			if ( counter_bases > maxsiz){
				//then we have the number of beses wanted
				begin_offset_chunk[index_in_chunk] = begin_offset;
				chunk_size[index_in_chunk] 	   = bytes_in_chunk;
				reads_in_chunk[index_in_chunk]	   = counter_bases / blen; 
				chunk_count 			   +=1;
				index_in_chunk 			   +=1;
				//we reset counter of bases
				counter_bases  			   = 0;
				bytes_in_chunk 			   = 0;
			}
		}
	}
	else{

		counter_bases 			= bases_previous_chunk;
		bytes_in_chunk 			= size_previous_chunk;
		begin_offset			= offset_previous_chunk;

		for (u1 = 0; u1 < local_num_reads; u1++){
			//we for offsets multiple of maxsize

			if (counter_bases == 0) begin_offset = local_read_offsets[u1];

			counter_bases 	+= local_read_size[u1];
			bytes_in_chunk  += local_read_offsets[u1+1] - local_read_offsets[u1];

			if ( counter_bases > maxsiz){
				//then we have the number of bases wanted
				begin_offset_chunk[index_in_chunk] = begin_offset;
				chunk_size[index_in_chunk] 	   = bytes_in_chunk;
				reads_in_chunk[index_in_chunk]	   = counter_bases / blen; 
				chunk_count 			   +=1;
				index_in_chunk 			   +=1;
				//we reset counter of bases
				counter_bases  			   = 0;
				bytes_in_chunk 			   = 0;
			}
		}
	}
	if (rank_num == (proc_num - 1)){

		// we complete the last chunk
		// with the last reads

		begin_offset_chunk[index_in_chunk] 	= begin_offset;
		chunk_size[index_in_chunk] 		= bytes_in_chunk;
		reads_in_chunk[index_in_chunk]		= counter_bases / blen;

		index_in_chunk 	+=1;
		chunk_count 	+=1;

	}
	if (rank_num < (proc_num -1)){

		//we send to rank + 1
		MPI_Send(&counter_bases,
				1,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				0,
				MPI_COMM_WORLD);

		MPI_Send(&begin_offset,
				1,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				0,
				MPI_COMM_WORLD);

		MPI_Send(&bytes_in_chunk,
				1,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				0,
				MPI_COMM_WORLD);
	}

	aft = MPI_Wtime();
	fprintf(stderr, "%s: rank %d time spend evaluating chunks = (%.02f) \n", __func__, rank_num, aft - bef);

	free(local_read_offsets);
	free(local_read_size);


	/*
	 * Map reference genome indexes in shared memory (by host)
	 */
	
	MPI_Info win_info;
	MPI_Info_create(&win_info);
	MPI_Info_set(win_info, "alloc_shared_non_contig", "true");

	bef = MPI_Wtime();
	res = MPI_File_open(MPI_COMM_WORLD, file_map, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_map);
	assert(res == MPI_SUCCESS);
	res = MPI_File_get_size(fh_map, &size_map);
	assert(res == MPI_SUCCESS);
	res = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_shr);
	assert(res == MPI_SUCCESS);
	res = MPI_Comm_rank(comm_shr, &rank_shr);
	assert(res == MPI_SUCCESS);
	size_shr = (rank_shr == 0) ? size_map : 0;
	res = MPI_Win_allocate_shared(size_shr, 0, win_info, comm_shr, &addr, &win_shr);
	assert(res == MPI_SUCCESS);
	MPI_Info_free(&win_info);
	res = MPI_Win_shared_query(win_shr, MPI_PROC_NULL, &size_shr, &res, &addr_map);
	assert(res == MPI_SUCCESS);

	m = size_map; a = addr_map; size_tot = 0;
	while (rank_shr == 0) {
		res = MPI_File_read(fh_map, a, INT_MAX/2, MPI_UINT8_T, &status);
		assert(res == MPI_SUCCESS);
		res = MPI_Get_count(&status, MPI_UINT8_T, &count);
		assert(res == MPI_SUCCESS);
		if (count == 0) break;
		m -= count; a += count; size_tot += count; }
	assert(size_tot == 0 || size_tot == size_map);
	res = MPI_Win_fence(0, win_shr);
	assert(res == MPI_SUCCESS);
	bwa_mem2idx(size_map, addr_map, &indix);
	if (ignore_alt)
		for (c = 0; c < indix.bns->n_seqs; ++c)
		indix.bns->anns[c].is_alt = 0;

	res = MPI_File_close(&fh_map);
	assert(res == MPI_SUCCESS);
	aft = MPI_Wtime();
	fprintf(stderr, "%s: mapped indexes (%.02f)\n", __func__, aft - bef);

	/*
	 * Create SAM header
	 * TODO: Add line for BWA version
	 */
	if (rank_num == 0) {
		int s, len;
		char *buff;

		res = MPI_File_delete(file_out, MPI_INFO_NULL);
		assert(res == MPI_SUCCESS || res == MPI_ERR_NO_SUCH_FILE || res == MPI_ERR_IO);
		res = MPI_File_open(MPI_COMM_SELF, file_out, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_out);
		assert(res == MPI_SUCCESS);
		/* Add reference sequence lines */
		for (s = 0; s < indix.bns->n_seqs; ++s) {
		len = asprintf(&buff, "@SQ\tSN:%s\tLN:%d\n", indix.bns->anns[s].name, indix.bns->anns[s].len);
		res = MPI_File_write(fh_out, buff, len, MPI_CHAR, &status);
		assert(res == MPI_SUCCESS);
		res = MPI_Get_count(&status, MPI_CHAR, &count);
		assert(res == MPI_SUCCESS);
		assert(count == len);
		free(buff);
		}
		/* Add header lines */
		if (hdr_line != NULL) {
		len = asprintf(&buff, "%s\n", hdr_line);
		res = MPI_File_write(fh_out, buff, len, MPI_CHAR, &status);
		assert(res == MPI_SUCCESS);
		res = MPI_Get_count(&status, MPI_CHAR, &count);
		assert(res == MPI_SUCCESS);
		assert(count == len);
		free(buff);
		}
		/* Add read group line */
		if (rg_line != NULL) {
		len = asprintf(&buff, "%s\n", rg_line);
		res = MPI_File_write(fh_out, buff, len, MPI_CHAR, &status);
		assert(res == MPI_SUCCESS);
		res = MPI_Get_count(&status, MPI_CHAR, &count);
		assert(res == MPI_SUCCESS);
		assert(count == len);
		free(buff);
		}
		res = MPI_File_close(&fh_out);
		assert(res == MPI_SUCCESS);
	}
	bef = MPI_Wtime();
	res = MPI_Barrier(MPI_COMM_WORLD);
	assert(res == MPI_SUCCESS);
	aft = MPI_Wtime();
	xfprintf(stderr, "%s: synched processes (%.02f)\n", __func__, aft - bef);


	res = MPI_File_open(MPI_COMM_WORLD, file_out, MPI_MODE_WRONLY|MPI_MODE_APPEND, MPI_INFO_NULL, &fh_out);
	assert(res == MPI_SUCCESS);

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


	// here we loop until there's nothing to read
	// in the offset and size file
	//char *p, *q, *e;
	char *p1, *q1, *e1, *p2, *q2, *e2;
	int line_number, line_number2;
	size_t reads_r1, reads_r2, reads;
	int64_t bases;
	size_t offset_chunk;
	size_t size_chunk;
	size_t total_local_reads_aligned= 0;
	size_t total_reads_check 		= 0;
	size_t local_time_spend_mapping = 0;
	double before_local_mapping = 0;
	double after_local_mapping	= 0;
	size_t total_time_local_mapping  = 0;
	size_t total_time_mapping = 0;
	size_t total_time_reading_seq = 0;
	size_t grand_total_time_reading_seq = 0;
	size_t total_time_parsing = 0;
	size_t grand_total_time_parsing = 0;
	size_t total_time_writing = 0;
	size_t grand_total_time_writing = 0;
	before_local_mapping = MPI_Wtime();

	//we loop the chunck_count
	for (u1 = 0; u1 < chunk_count; u1++){

		offset_chunk = begin_offset_chunk[u1];
		size_chunk   = chunk_size[u1];
		assert(size_chunk != 0);
		/*
		 * Read sequence datas ...
		 *
		 */
		bef = MPI_Wtime();

		buffer_r1 = malloc(size_chunk+1);
		assert(buffer_r1 != NULL);
		buffer_r1[size_chunk] = 0;
		
		// FOR TEST //
		//MPI_Type_contiguous(size_chunk, MPI_CHAR, &arraytype_r1);
		//MPI_Type_commit(&arraytype_r1);
		//MPI_File_set_view(fh_r2, (MPI_Offset)offset_chunk, MPI_CHAR, MPI_CHAR, "native", finfo ) ; 
		//res = MPI_File_read(fh_r2, buffer_r1, 1, arraytype_r1, &status);

		res = MPI_File_read_at(fh_r1, offset_chunk, buffer_r1, size_chunk, MPI_CHAR, &status);		
		assert(res == MPI_SUCCESS);
		res = MPI_Get_count(&status, MPI_CHAR, &count);
		assert(res == MPI_SUCCESS);
		assert(count == (int)size_chunk && *buffer_r1 == '@');


		buffer_r2 = malloc(size_chunk+1);
		assert(buffer_r2 != NULL);
		buffer_r2[size_chunk]=0;

		// FOR TEST //
		//MPI_Type_contiguous(size_chunk, MPI_CHAR, &arraytype_r2);
		//MPI_Type_commit(&arraytype_r2);
		//MPI_File_set_view(fh_r1, (MPI_Offset)offset_chunk, MPI_CHAR, MPI_CHAR, "native", finfo ) ; 
		//res = MPI_File_read(fh_r1, buffer_r2, 1, arraytype_r2, &status);		

		res = MPI_File_read_at(fh_r2, offset_chunk, buffer_r2, size_chunk, MPI_CHAR, &status);		
		assert(res == MPI_SUCCESS);
		res = MPI_Get_count(&status, MPI_CHAR, &count);
		assert(res == MPI_SUCCESS);
		assert(count == (int) size_chunk && *buffer_r2 == '@');

		aft = MPI_Wtime();
		fprintf(stderr, "%s: read sequences (%.02f)\n", __func__, aft - bef);
		total_time_reading_seq += (aft - bef);

		bef = MPI_Wtime();
		reads_r1 = reads_r2 = reads_in_chunk[u1];
		reads = reads_r1 + reads_r2; bases = 0;
		total_local_reads_aligned += reads/2;
		assert(reads <= INT_MAX);
		fprintf(stderr, "%s: num_rank = %d :: number of reads = %zu \n", __func__, rank_num, reads);

		/* Parse sequences ... */
		seqs = malloc(reads * sizeof(*seqs));
		assert(seqs != NULL);
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
		localsize = 0;
		for (n = 0; n < reads; n++) {
			/* Reuse .l_seq to store SAM line length to avoid multiple strlen() calls */
			seqs[n].l_seq = strlen(seqs[n].sam);
			localsize += seqs[n].l_seq; }
		assert(localsize <= INT_MAX);
		buffer_out = malloc(localsize);
		assert(buffer_out != NULL);
		p = buffer_out;
		for (n = 0; n < reads; n++) {
			memmove(p, seqs[n].sam, seqs[n].l_seq);
			p += seqs[n].l_seq;
			free(seqs[n].sam); }
		free(seqs);
		res = MPI_File_write_shared(fh_out, buffer_out, localsize, MPI_CHAR, &status);
		assert(res == MPI_SUCCESS);
		res = MPI_Get_count(&status, MPI_CHAR, &count);
		assert(res == MPI_SUCCESS);
		assert(count == (int)localsize);
		free(buffer_out);
		aft = MPI_Wtime();
		fprintf(stderr, "%s: wrote results (%.02f)\n", __func__, aft - bef);
		total_time_writing += (aft - bef);
		free(buffer_r1);
		free(buffer_r2);

	}

	after_local_mapping	 = MPI_Wtime();
	total_time_local_mapping = after_local_mapping - before_local_mapping;

	res = MPI_Allreduce(&total_time_local_mapping, &total_time_mapping, 1,MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
	assert ( res == MPI_SUCCESS);

	res = MPI_Allreduce(&total_local_reads_aligned, &total_reads_check, 1,MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
	assert ( res == MPI_SUCCESS);

	res = MPI_Allreduce(&total_time_reading_seq, &grand_total_time_reading_seq, 1,MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
	assert ( res == MPI_SUCCESS);

	res = MPI_Allreduce(&total_time_parsing, &grand_total_time_parsing, 1,MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
	assert ( res == MPI_SUCCESS);

	res = MPI_Allreduce(&total_time_writing, &grand_total_time_writing, 1,MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
	assert ( res == MPI_SUCCESS);
	
	fprintf(stderr, "rank %d :::: total_time_reading_seq = %zu \n", rank_num, grand_total_time_reading_seq);
	fprintf(stderr, "rank %d :::: total_time_writing = %zu \n", rank_num, grand_total_time_writing);
	fprintf(stderr, "rank %d :::: total_time_parsing = %zu \n", rank_num, grand_total_time_parsing);
	fprintf(stderr, "rank %d :::: total_time_local_mapping = %zu \n", rank_num, total_time_local_mapping);
	fprintf(stderr, "rank %d :::: total_time_mapping = %zu \n", rank_num, total_time_mapping);

	fprintf(stderr, "rank %d :::: local_num_reads = %zu \n", rank_num, local_num_reads);
	fprintf(stderr, "rank %d :::: total_reads_check = %zu \n", rank_num, total_reads_check);
	fprintf(stderr, "rank %d :::: total_local_reads_aligned = %zu \n", rank_num, total_local_reads_aligned);
	fprintf(stderr, "rank %d :::: total_num_reads = %zu \n", rank_num, grand_total_num_reads);

	assert( total_reads_check == grand_total_num_reads );

	fprintf(stderr, "rank %d :::: finish mappings for reads \n", rank_num);
	if (begin_offset_chunk != NULL) free(begin_offset_chunk);
	if (chunk_size != NULL) free(chunk_size);
	if (reads_in_chunk != NULL) free(reads_in_chunk);

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

	res = MPI_File_close(&fh_out);
	assert(res == MPI_SUCCESS);

	bef = MPI_Wtime();
	res = MPI_Finalize();
	assert(res == MPI_SUCCESS);
	aft = MPI_Wtime();
	fprintf(stderr, "%s: finalized synched processes (%.02f)\n", __func__, aft - bef);

	return 0;
}

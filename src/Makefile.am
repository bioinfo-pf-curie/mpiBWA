
AM_CPPFLAGS = -DUSE_MALLOC_WRAPPERS -DHAVE_PTHREAD

bin_PROGRAMS = mpiBWA mpiBWAByChr mpiBWAIdx
noinst_LIBRARIES = libbwa.a
noinst_HEADERS = $(HDR)

libbwa_a_SOURCES = $(SRC)

mpiBWA_SOURCES = mainParallel.c
mpiBWA_LDADD = libbwa.a
mpiBWAByChr_SOURCES = mainParallelByChromosome.c
mpiBWAByChr_LDADD = libbwa.a
mpiBWAIdx_SOURCES = pidx.c
mpiBWAIdx_LDADD = libbwa.a

SRC = utils.c kthread.c kstring.c ksw.c bwt.c bntseq.c bwa.c bwamem.c \
	bwamem_pair.c bwamem_extra.c malloc_wrap.c tokenizer.c fixmate.c parallel_aux.c bgzf.c
HDR = utils.h kstring.h ksw.h bwt.h bntseq.h bwa.h bwamem.h malloc_wrap.h \
	kbtree.h khash.h kseq.h ksort.h kvec.h tokenizer.h fixmate.h parallel_aux.h bgzf.h


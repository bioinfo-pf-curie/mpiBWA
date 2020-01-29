
AM_CPPFLAGS = -DUSE_MALLOC_WRAPPERS -DHAVE_PTHREAD

bin_PROGRAMS = mpiBWA mpiBWAIdx
noinst_LIBRARIES = libbwa.a
noinst_HEADERS = $(HDR)

libbwa_a_SOURCES = $(SRC)

pbwa7_SOURCES = main_parallel_version.c
pbwa7_LDADD = libbwa.a
pidx_SOURCES = pidx.c
pidx_LDADD = libbwa.a

SRC = utils.c kthread.c kstring.c ksw.c bwt.c bntseq.c bwa.c bwamem.c \
	bwamem_pair.c bwamem_extra.c malloc_wrap.c
HDR = utils.h kstring.h ksw.h bwt.h bntseq.h bwa.h bwamem.h malloc_wrap.h \
	kbtree.h khash.h kseq.h ksort.h kvec.h

EXTRA_DIST = INSTALL

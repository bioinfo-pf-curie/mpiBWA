/*
This file is part of <mpiBWA>

pidx create a reference .map a copy of what's in memory.

Copyright (C) 2016  Institut Curie / Institut Pasteur

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

#include <assert.h>
#include <fcntl.h>
#include <libgen.h>
#include <limits.h>
#include <stdio.h>

#include "bwa.h"

extern int bwa_idx2mem(bwaidx_t *idx);
extern int bwa_mem2idx(int64_t l_mem, uint8_t *mem, bwaidx_t *idx);

int main(int argc, char **argv) {
	char *indx, imem[PATH_MAX];
	int fd;
	off_t off;
	ssize_t siz, rem, tot;

	bwaidx_t *idx;

	if (argc != 2) {
		fprintf(stderr, "usage: %s <index>\n", basename(*argv));
		return 1; }

	indx = *(argv+1);

	sprintf(imem, "%s.map", indx);

	idx = bwa_idx_load(indx, BWA_IDX_ALL);
	assert(idx != NULL);
	assert(bwa_idx2mem(idx) == 0);

	fd = open(imem, O_WRONLY|O_CREAT|O_TRUNC, 0666);
	assert(fd != -1);
	off = 0; rem = idx->l_mem; tot = 0;
	while (rem > 0) {
		siz = write(fd, idx->mem + off, (size_t)rem);
		assert(siz != -1);
		off += siz; rem -= siz; tot += siz; }
	assert(tot == (ssize_t)idx->l_mem);
	assert(close(fd) != -1);

	return 0; }

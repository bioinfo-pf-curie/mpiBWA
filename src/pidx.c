/*

This file is part of mpiBWA

The project was developped by Frederic Jarlier from Institut Curie and Nicolas Joly from Institut Pasteur

pidx create a reference .map a copy of what's in memory.

Copyright (C) 2016-2020  Institut Curie / Institut Pasteur

You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
   
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
   
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

*/

#include <assert.h>
#include <fcntl.h>
#include <libgen.h>
#include <limits.h>
#include <stdio.h>
#include <unistd.h>

#include "bwa.h"

int main(int argc, char **argv) {
	char *indx, imem[PATH_MAX];
	int fd;
	off_t off;
	ssize_t siz, rem, tot;

	bwaidx_t *idx;

	if (argc != 2) {
		fprintf(stderr, "
		This file is part of mpiBWA

		Copyright Institut Curie 2020

		This software is a computer program whose purpose is to map the reference in memory and write it in a binary file.

		You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).

		The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 

		The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

		usage: %s <index>\n", basename(*argv));
		
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

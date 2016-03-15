
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#include <assert.h>
#include <fcntl.h>
#include <libgen.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "bwa.h"

#ifdef COMPAT_BWA
int bwa_idx2mem(bwaidx_t *idx);
int bwa_mem2idx(int64_t l_mem, uint8_t *mem, bwaidx_t *idx);

int bwa_idx2mem(bwaidx_t *idx)
{
        int i;
        int64_t k, x, tmp;
        uint8_t *mem;

        // copy idx->bwt
        x = idx->bwt->bwt_size * 4;
        mem = realloc(idx->bwt->bwt, sizeof(bwt_t) + x); idx->bwt->bwt = 0;
        memmove(mem + sizeof(bwt_t), mem, x);
        memcpy(mem, idx->bwt, sizeof(bwt_t)); k = sizeof(bwt_t) + x;
        x = idx->bwt->n_sa * sizeof(bwtint_t); mem = realloc(mem, k + x); memcpy(mem + k, idx->bwt->sa, x); k += x;
        free(idx->bwt->sa);
        free(idx->bwt); idx->bwt = 0;

        // copy idx->bns
        tmp = idx->bns->n_seqs * sizeof(bntann1_t) + idx->bns->n_holes * sizeof(bntamb1_t);
        for (i = 0; i < idx->bns->n_seqs; ++i) // compute the size of heap-allocated memory
                tmp += strlen(idx->bns->anns[i].name) + strlen(idx->bns->anns[i].anno) + 2;
        mem = realloc(mem, k + sizeof(bntseq_t) + tmp);
        x = sizeof(bntseq_t); memcpy(mem + k, idx->bns, x); k += x;
        x = idx->bns->n_holes * sizeof(bntamb1_t); memcpy(mem + k, idx->bns->ambs, x); k += x;
        free(idx->bns->ambs);
        x = idx->bns->n_seqs * sizeof(bntann1_t); memcpy(mem + k, idx->bns->anns, x); k += x;
        for (i = 0; i < idx->bns->n_seqs; ++i) {
                x = strlen(idx->bns->anns[i].name) + 1; memcpy(mem + k, idx->bns->anns[i].name, x); k += x;
                x = strlen(idx->bns->anns[i].anno) + 1; memcpy(mem + k, idx->bns->anns[i].anno, x); k += x;
                free(idx->bns->anns[i].name); free(idx->bns->anns[i].anno);
        }
        free(idx->bns->anns);

        // copy idx->pac
        x = idx->bns->l_pac/4+1;
        mem = realloc(mem, k + x);
        memcpy(mem + k, idx->pac, x); k += x;
        free(idx->bns); idx->bns = 0;
        free(idx->pac); idx->pac = 0;

        return bwa_mem2idx(k, mem, idx);
}

int bwa_mem2idx(int64_t l_mem, uint8_t *mem, bwaidx_t *idx)
{
        int64_t k = 0, x;
        int i;

        // generate idx->bwt
        x = sizeof(bwt_t); idx->bwt = malloc(x); memcpy(idx->bwt, mem + k, x); k += x;
        x = idx->bwt->bwt_size * 4; idx->bwt->bwt = (uint32_t*)(mem + k); k += x;
        x = idx->bwt->n_sa * sizeof(bwtint_t); idx->bwt->sa = (bwtint_t*)(mem + k); k += x;

        // generate idx->bns and idx->pac
        x = sizeof(bntseq_t); idx->bns = malloc(x); memcpy(idx->bns, mem + k, x); k += x;
        x = idx->bns->n_holes * sizeof(bntamb1_t); idx->bns->ambs = (bntamb1_t*)(mem + k); k += x;
        x = idx->bns->n_seqs  * sizeof(bntann1_t); idx->bns->anns = malloc(x); memcpy(idx->bns->anns, mem + k, x); k += x;
        for (i = 0; i < idx->bns->n_seqs; ++i) {
                idx->bns->anns[i].name = (char*)(mem + k); k += strlen(idx->bns->anns[i].name) + 1;
                idx->bns->anns[i].anno = (char*)(mem + k); k += strlen(idx->bns->anns[i].anno) + 1;
        }
        idx->pac = (uint8_t*)(mem + k); k += idx->bns->l_pac/4+1;
        assert(k == l_mem);

        idx->l_mem = k; idx->mem = mem;
        return 0;
}
#endif

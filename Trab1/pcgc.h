#ifndef __PCGC_H__
#define __PCGC_H__

#include <stdlib.h>
#include "utils.h"
#include "sislin.h"

inline real_t** Sempc();
inline real_t** pcJacobi();
real_t** pcGausSeidel();
real_t** pcSSOR();
void resolveGC(real_t *A, real_t *b, real_t **M, real_t *X, int maxit, int e, real_t w, real_t n, real_t k, rtime_t *tempo);
void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k, real_t **M, rtime_t *tempo);

#endif // __PCGC_H__
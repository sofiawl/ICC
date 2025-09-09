#ifndef __GAUSS_SEIDEL_EDO__
#define __GAUSS_SEIDEL_EDO__

#include "utils.h"

// Equação Diferencial Ordinária
typedef struct {
 int n; // número de pontos internos na malha
 real_t a, b; // intervalo
 real_t ya, yb; // condições contorno
 real_t (* p)(real_t), (* q)(real_t), (* r)(real_t);
} EDo;


void prnVetor (real_t *v, unsigned int n);
rtime_t gaussSeidel_EDO (EDo *edoeq, real_t *Y, unsigned int maxiter);
real_t normaL2_EDO (EDo *edoeq, real_t *Y);

#endif
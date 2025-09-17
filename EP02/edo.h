#ifndef __EQDIFF_H__
#define __EQDIFF_H__

typedef double real_t;

#define FORMAT "%23.15e"
#define ERRO 1e-5
#define MAXIT 100


// Sistema linear Tri-diagonal
typedef struct {
  real_t *D, *Di, *Ds, *B;
  int n;
} Tridiag;

// Equação Diferencial Ordinária
typedef struct {
  int n; // número de pontos internos na malha
  real_t a, b; // intervalo
  real_t ya, yb; // condições contorno
  real_t p, q, r1, r2, r3, r4; // coeficientes EDO genérica
} EDo;

// Funções auxiliares

Tridiag *genTridiag (EDo *edoeq);
void prnEDOsl (EDo *edoeq);
void prnVetor (real_t *v, unsigned int n);
real_t gaussSeidel_3Diag (Tridiag *sl, real_t *Y, unsigned int *maxiter);
real_t normaL2_3Diag (Tridiag *sl, real_t *Y);
real_t gaussSeidel_EDO (EDo *edoeq, real_t *Y, unsigned int *maxiter);
real_t normaL2_EDO (EDo *edoeq, real_t *Y);

#endif // __EQDIFF_H__


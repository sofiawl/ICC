/*-------------------------------------------------------------
									SOFIA WAMSER LIMA 
										GRR: 20240495
											20/09/2025
---------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fenv.h>

#include "utils.h"
#include "edo.h"

Tridiag *genTridiag (EDo *edo)
{
  Tridiag *sl;
  real_t x, rx;
  int n = edo->n;
  
  sl = (Tridiag *) malloc (sizeof(Tridiag));
  sl->n = edo->n;

  sl->D = (real_t *) calloc(n, sizeof(real_t));
  sl->Di = (real_t *) calloc(n, sizeof(real_t));
  sl->Ds = (real_t *) calloc(n, sizeof(real_t));
  sl->B = (real_t *) calloc(n, sizeof(real_t));

  real_t h = (edo->b - edo->a)/(n+1);

  for (int i=0; i < n; ++i) {
    x = edo->a + (i+1)*h;
    rx = edo->r1*x + edo->r2*x*x + edo->r3*cos(x) + edo->r4*exp(x);
    
    sl->B[i] = h*h * rx;
    sl->Di[i] = 1 - h * edo->p/2.0;
    sl->D[i] = -2 + h*h * edo->q;
    sl->Ds[i] = 1 + h * edo->p/2.0;
  }

  sl->B[0] -= edo->ya * (1 - h*edo->p/2.0);
  sl->B[n-1] -= edo->yb * (1 + h*edo->p/2.0);
  
  return sl;
}


// Exibe SL na saída padrão
void prnEDOsl (EDo *edoeq)
{
  int n = edoeq->n, i, j;
  real_t x, b, d, di, ds,rx;
  real_t h = (edoeq->b - edoeq->a)/(n+1);

  printf ("%d\n", n);

  for (i=0; i < n; ++i) {
    x = edoeq->a + (i+1)*h;
    rx = edoeq->r1*x + edoeq->r2*x*x + edoeq->r3*cos(x) + edoeq->r4*exp(x);
    
    b = h*h * rx; 
    di = 1 - h * edoeq->p/2.0;
    d = -2 + h*h * edoeq->q;
    ds = 1 + h * edoeq->p/2.0;
    
    for (j=0; j < n; ++j) {
      if (i == j)
	printf (FORMAT,d);
      else if (j == i-1 && i)
	printf (FORMAT,di);
      else if (j == i+1 && i != n-1)
	printf (FORMAT,ds);
      else
	printf(FORMAT, 0.0);
    }
      
    if (i == 0)
      b -= edoeq->ya * (1 - h*edoeq->p/2.0);
    else if (i == n-1)
      b -= edoeq->yb * (1 + h*edoeq->p/2.0);

    printf (FORMAT, b);
      
    printf ("\n");
  }
}


void prnVetor (real_t *v, unsigned int n)
{
  printf (" ");
  for(unsigned int i=0; i < n; ++i)
      printf (FORMAT, v[i]);
  printf ("\n");
}


// Algoritmo  Gauss-Seidel com  vetores das diagonais
// Parâmentros, respectivamente: sistema tridiagonal, vetor de resultados
// ponteiro para valor máximo de iterações, ponteiro para valor da norma L2
rtime_t gaussSeidel_3Diag (Tridiag *sl, real_t *Y, unsigned int *maxiter, real_t *normaL2)
{
  real_t tTotal = timestamp();
  *normaL2 = normaL2_3Diag(sl, Y);

  unsigned int k = 0;
  do { 
    // Primeiro valor do vetor
    Y[0] = (sl->B[0] - sl->Ds[0]*Y[1])/ sl->D[0];

    for(int i = 1; i < sl->n-1; i++)
      Y[i] = (sl->B[i] - sl->Di[i-1]*Y[i-1] - sl->Ds[i]*Y[i+1]) / sl->D[i];

    Y[sl->n-1] = (sl->B[sl->n-1] - sl->Di[sl->n-2]*Y[sl->n-2]) / sl->D[sl->n-1];

    k++;
    *normaL2 = normaL2_3Diag(sl, Y);

  } while(k < *maxiter && *normaL2 > ERRO);

  tTotal = timestamp() - tTotal;

  *maxiter = k;  

  return tTotal;
}


// algoritmo para calcular Norma L2 com  vetores   das  diagonais   e  termos
// independentes do SL
// Parâmetros, respectivamente: sistema tridiagonal, vetor de resultados
real_t normaL2_3Diag (Tridiag *sl, real_t *Y)
{
    int n = sl->n;
    if (n <= 0) return -1.0;

    real_t soma = 0.0;
  
    real_t r0 = sl->B[0] - sl->D[0]*Y[0] - sl->Ds[0]*Y[1];
    soma += r0 * r0;

    for (int i = 1; i < n-1; ++i) {
        r0 = sl->B[i] - sl->Di[i-1]*Y[i-1] - sl->D[i]*Y[i] - sl->Ds[i]*Y[i+1];
        soma += r0 * r0;
    }

    r0 = sl->B[n-1] - sl->Di[n-2]*Y[n-2] - sl->D[n-1]*Y[n-1];
    soma += r0 * r0;

    return sqrt(soma);
}
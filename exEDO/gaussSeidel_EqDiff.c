#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "edo.h"
#include "gaussSeidel_EqDiff.h"

// Gera diagonais e termos independentes a partir da definição do EDO
Tridiag *genTridiag (EDo *edo)
{
  Tridiag *sl;
  real_t x, xnext;
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
    sl->B[i] = h*h * edo->r(x);
    sl->Di[i] = 1 - h * edo->p(x)/2.0;
    sl->D[i] = -2 + h*h * edo->q(x);
    sl->Ds[i] = 1 + h * edo->p(x)/2.0;
  }

  sl->B[0] -= edo->ya * (1 - h*edo->p(edo->a+h)/2.0);
  sl->B[n-1] -= edo->yb * (1 + h*edo->p(edo->b-h)/2.0);
  
  return sl;
}



rtime_t gaussSeidel_3Diag (Tridiag *sl, real_t *Y, unsigned int maxiter)
{
  int n = sl->n;
  
  rtime_t tTotal = timestamp();

  // algoritmo  Gauss-Seidel   com  vetores   das  diagonais   e  termos
  // independentes do SL
  while(maxiter)
  { 
    Y[0] = (sl->B[0] - sl->Ds[0]*Y[1])/ sl->Di[0];

    for(int i = 1; i < sl->n-2; i++)
      Y[i] = (sl->B[i] - sl->Di[i-1]*Y[i-1] - sl->Ds[i]*Y[i+1]) / sl->D[i];

    Y[sl->n-1] = (sl->B[sl->n-1] - sl->Di[sl->n-2]*Y[sl->n-2]) / sl->D[sl->n-1];

    maxiter--;
  }
  
  return timestamp() - tTotal;
}


real_t normaL2_3Diag (Tridiag *sl, real_t *Y)
{
  int n = sl->n;
  real_t normaL2;

  normaL2 = 0.0;

  // algoritmo para calcular Norma L2 com  vetores   das  diagonais   e  termos
  // independentes do SL

  //1: resíduo
  real_t *R = malloc(sizeof(real_t) * sl->n);
  if (!R) return -1;


  R[0] = sl->B[0] - sl->D[0]*Y[0] - sl->Ds[0]*Y[1];

  for(int i=1; i < sl->n-1; ++i) {
    R[i] = sl->B[i] - sl->Di[i]*Y[i-1] - sl->D[i]*Y[i] - sl->Ds[i]*Y[i+1];
  }

  R[sl->n-1] = sl->B[sl->n-1] - sl->Di[sl->n-2]*Y[sl->n-2] - sl->D[sl->n-1]*Y[sl->n-1]; 

  //2: norma L2
  for(int i = 0; i < sl->n; i++){
    normaL2 = R[i]*R[i];
  }
  
  normaL2 = sqrt(normaL2);

  return normaL2;
  
}


rtime_t gaussSeidel_EDO (EDo *edoeq, real_t *Y, unsigned int maxiter)
{
  int n = edoeq->n;
  real_t *x, b, yi, d, di, ds, h;

  rtime_t tTotal = timestamp();

  h = (edoeq->b - edoeq->a) / (n+1);
  ds = 2 + h;
  d = 2*h*h;
  di = 2 - h;
  b = 2*h*h;

  // algoritmo Gauss-Seidel usando parâmetros EDO, sem usar vetores para
  // diagonais e termos independentes do SL

  for (int k=0; k < maxiter; ++k) {
    real_t xi = edoeq->a + h;
    Y[0] = (b* edoeq->r(xi) - ds*x[1]) / d*edoeq->q(xi)-4;
    
    for(int i = 1; i < edoeq->n-2; i++){
      xi += h;
      Y[i] = (b* edoeq->r(xi) - di* edoeq->p(xi)* Y[i-1] - ds*Y[i+1]) / (d*edoeq->q(xi) - 4);
    }

    xi += h; 
    Y[edoeq->n-1] = (d* edoeq->r(xi) - di*Y[edoeq->n-1]) / (d*edoeq->q(xi) - 4);
  }

  return timestamp() - tTotal;
}

real_t normaL2_EDO (EDo *edoeq, real_t *Y)
{
  int n=edoeq->n, i;
  real_t normaL2, res, x, b, d, di, ds, h;

  normaL2 = 0.0;

  h = (edoeq->b-edoeq->a)/(n+1);
  ds = 2 + h;
  d = 2*h*h;
  di = 2 - h;
  b = 2*h*h;

  // algoritmo para calcular Norma L2 usando parâmetros EDO, sem usar vetores para
  // diagonais e termos independentes do SL  

  //resíduo + norma L2
  real_t aux;

  real_t xi = edoeq->a + h;
  aux = b*edoeq->r(xi) - (d*Y[0]) - ds*Y[1];
  res = aux*aux; 

  for(int i=1; i < edoeq->n-1; ++i) {
    xi += h; 
    aux = b*edoeq->r(xi) - di*Y[i-1] - d*Y[i] - ds*Y[i+1];
    res += aux*aux;
  }

  xi += h;
  aux = b*edoeq->r(xi) - di*Y[edoeq->n-2] - d*Y[edoeq->n-1]; 
  res = aux*aux;


  normaL2 = sqrt(res);

  return normaL2;
}

// Exibe um vetor na saída padrão
void prnVetor (real_t *v, unsigned int n)
{
  int i;

  //  printf ("\n");
  printf ("   ");
  for(i=0; i < n; ++i)
      printf ("%20.12e", v[i]);
  printf ("\n");

}

void prnTriDiagonal(Tridiag *sl)
{
    prnVetor(sl->Di,sl->n);
    prnVetor(sl->D,sl->n);
    prnVetor(sl->Ds,sl->n);
    prnVetor(sl->B,sl->n);

}


// Exibe SL na saída padrão
void prnEDOsl (EDo *edoeq, int diagOnly)
{
  int n = edoeq->n, i, j;
  real_t x, b, d, di, ds;
  real_t h = (edoeq->b - edoeq->a)/(n+1);

  printf ("\tn = %d, \t H = %.8g\nSL:\n", n, h);
  if (diagOnly) {
    prnTriDiagonal(genTridiag(edoeq));
  }
  else {
    for (i=0; i < n; ++i) {
      x = edoeq->a + (i+1)*h;
      b = h*h * edoeq->r(x);      
      di = 1 - h * edoeq->p(x)/2.0;
      d = -2 + h*h * edoeq->q(x);
      ds = 1 + h * edoeq->p(x)/2.0;
      
      for (j=0; j < n; ++j) {
	if (i == j)
	  printf ("%20.12e",d);
	else if (j == i-1 && i)
	  printf ("%20.12e",di);
	else if (j == i+1 && i != n-1)
	  printf ("%20.12e",ds);
	else
	  printf("%20.12e", 0.0);
      }
      
      if (i == 0)
	b -= edoeq->ya * (1 - h*edoeq->p(edoeq->a+h)/2.0);
      else if (i == n-1)
	b -= edoeq->yb * (1 + h*edoeq->p(edoeq->b-h)/2.0);
      
      printf (" | %20.12e\n", b);
    }
  }
}



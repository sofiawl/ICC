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
    /*
    ds = 2.0 + h*edoeq->p;
    d = 2.0*h*h*edoeq->q - 4.0;
    di = 2.0 - h*edoeq->p;
    b = 2.0*h*h * rx;
    */
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


real_t gaussSeidel_3Diag (Tridiag *sl, real_t *Y, unsigned int *maxiter)
{
  
  // algoritmo  Gauss-Seidel   com  vetores   das  diagonais   e  termos
  // independentes do SL
  real_t normaL2 = normaL2_3Diag(sl, Y);

  unsigned int k = 1;
  while(k < *maxiter && normaL2 > ERRO)
  { 
    Y[0] = (sl->B[0] - sl->Ds[0]*Y[1])/ sl->Di[0];

    for(int i = 1; i < sl->n-2; i++)
      Y[i] = (sl->B[i] - sl->Di[i-1]*Y[i-1] - sl->Ds[i]*Y[i+1]) / sl->D[i];

    Y[sl->n-1] = (sl->B[sl->n-1] - sl->Di[sl->n-2]*Y[sl->n-2]) / sl->D[sl->n-1];

    k++;
    normaL2 = normaL2_3Diag(sl, Y);
  }

  *maxiter = k;  

  return normaL2;
}


// algoritmo para calcular Norma L2 com  vetores   das  diagonais   e  termos
// independentes do SL
real_t normaL2_3Diag (Tridiag *sl, real_t *Y)
{
  int n = sl->n;
  real_t normaL2;

  normaL2 = 0.0;

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
    normaL2 += R[i]*R[i];
  }
  
  free(R);
  return sqrt(normaL2);
}

// algoritmo Gauss-Seidel usando parâmetros EDO, sem usar vetores para
// diagonais e termos independentes do SL
real_t gaussSeidel_EDO (EDo *edoeq, real_t *Y, unsigned int *maxiter)
{
    int n = edoeq->n;
    real_t b, yi, d, di, ds, h; 

    h = (edoeq->b - edoeq->a) / (n+1);
    fesetround(h);
    /*ds = 2.0 + h*edoeq->p;
    d = 2.0*h*h*edoeq->q - 4.0;
    di = 2.0 - h*edoeq->p;
    b = 2.0*h*h;
    */
    b = h*h; 
    di = 1 - h * edoeq->p/2.0;
    d = -2 + h*h * edoeq->q;
    ds = 1 + h * edoeq->p/2.0;
    
    fesetround(ds); fesetround(d); fesetround(di); fesetround(b);

    real_t normaL2 = normaL2_EDO(edoeq, Y);

    unsigned int k = 1;
    for (; k < *maxiter && normaL2 > ERRO; k++)
    {
        real_t yi = edoeq->a + h;
        real_t r = edoeq->r1*yi + edoeq->r2*yi*yi + edoeq->r3*cos(yi) + edoeq->r4*exp(yi);
        // dúvida da formula di ou ds?
        Y[0] = (b*r - ds*Y[1]) / d;
        //printf("b: %fl r: %fl d: %fl Y0: %fl\n", b, r, d, Y[0]);
        fesetround(Y[0]);

        for(int i = 1; i < edoeq->n-2; i++){
            yi += h;
            r = edoeq->r1*yi + edoeq->r2*yi*yi + edoeq->r3*cos(yi) + edoeq->r4*exp(yi);
            Y[i] = (b*r - di*Y[i-1] - ds*Y[i+1]) / d;
            fesetround(Y[i]);
        }

        yi += h; 
        r = edoeq->r1*yi + edoeq->r2*yi*yi + edoeq->r3*cos(yi) + edoeq->r4*exp(yi);
        // dúvida da formula di ou ds?
        Y[edoeq->n-1] = (d*r - di*Y[edoeq->n-2]) / d;
        fesetround(Y[edoeq->n-1]);

        //prnVetor(Y, edoeq->n);
        //printf("\n");
        normaL2 = normaL2_EDO(edoeq, Y);
        
    }

  *maxiter = k;

  return normaL2;
}


// algoritmo para calcular Norma L2 usando parâmetros EDO, sem usar vetores para
// diagonais e termos independentes do SL    
real_t normaL2_EDO (EDo *edoeq, real_t *Y)
{
  int n=edoeq->n, i;
  real_t normaL2, res, b, d, di, ds, h;

  normaL2 = 0.0;

  h = (edoeq->b-edoeq->a)/(n+1);
  fesetround(h);
  /*ds = 2.0 + h*edoeq->p;
  d = 2.0*h*h*edoeq->q - 4.0;
  di = 2.0 - h*edoeq->p;
  b = 2.0*h*h;
  */
 
  b = h*h; 
  di = 1 - h * edoeq->p/2.0;
  d = -2 + h*h * edoeq->q;
  ds = 1 + h * edoeq->p/2.0;
    
  fesetround(ds); fesetround(d); fesetround(di); fesetround(b);

  //resíduo + norma L2
  real_t aux;

  real_t yi = edoeq->a + h;
  real_t r = edoeq->r1*yi + edoeq->r2*yi*yi + edoeq->r3*cos(yi) + edoeq->r4*exp(yi);
  fesetround(r);

  // dúvida da formula di ou ds?
  aux = b*r - d*Y[0] - ds*Y[1];
  fesetround(aux);
  res = aux*aux; 

  for(int i=1; i < edoeq->n-1; ++i) {
    yi += h;
    r = edoeq->r1*yi + edoeq->r2*yi*yi + edoeq->r3*cos(yi) + edoeq->r4*exp(yi);
    aux = b*r - di*Y[i-1] - d*Y[i] - ds*Y[i+1];
    fesetround(aux);
    res += aux*aux;
  }

  yi += h;
  r = edoeq->r1*yi + edoeq->r2*yi*yi + edoeq->r3*cos(yi) + edoeq->r4*exp(yi);
  // dúvida da formula di ou ds?
  aux = b*r - di*Y[edoeq->n-2] - d*Y[edoeq->n-1]; 
  fesetround(aux);
  res += aux*aux;


  normaL2 = sqrt(res);

  return normaL2;
}


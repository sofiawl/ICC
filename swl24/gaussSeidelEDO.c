#include <stdlib.h>
#include <stdio.h>
#include "utils.h"
#include "sislin.h"
#include "gaussSeidelEDO.h"


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

    // normaL2 debe ser menor que 10^-5
    int k = 1;
    do {
        real_t xi = edoeq->a + h;
        Y[0] = (b* edoeq->r(xi) - ds*x[1]) / d*edoeq->q(xi)-4;
    
        for(int i = 1; i < edoeq->n-2; i++){
            xi += h;
            Y[i] = (b* edoeq->r(xi) - di* edoeq->p(xi)* Y[i-1] - ds*Y[i+1]) / (d*edoeq->q(xi) - 4);
        }

        xi += h; 
        Y[edoeq->n-1] = (d* edoeq->r(xi) - di*Y[edoeq->n-1]) / (d*edoeq->q(xi) - 4);

        k++;

    } while(k < maxiter && normaL2_EDO(edoeq, Y) < 10^-5);

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
  aux = b*edoeq->r(xi) - (d*Y[0]-4) - ds*Y[1];
  res = aux*aux; 

  for(int i=1; i < edoeq->n-1; ++i) {
    xi += h; 
    aux = b*edoeq->r(xi) - di*Y[i-1] - d*Y[i]-4 - ds*Y[i+1];
    res += aux*aux;
  }

  // rever aqui
  xi += h;
  aux = b*edoeq->r(xi) - di*Y[edoeq->n-2] - d*Y[edoeq->n-1]-4; 
  res = aux*aux;


  normaL2 = sqrt(res);

  return normaL2;
}


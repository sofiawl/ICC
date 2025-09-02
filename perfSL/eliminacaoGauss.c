/* Matriz  'normal' (vetor  de  ponteiros (linhas  matriz) para  vetores
   (colunas da matriz), estilo 'Mazieiro/Prog 2'
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fenv.h>

#include "utils.h"
#include "sislin.h"
#include "eliminacaoGauss.h"

// tamanho n, coluna do pivo k
int encontraMax(SistLinear_t *C, int k, int n)
{

  int max = k;
  for(int i = k+1; i < n; i++){
    if(C->A[i][k] > C->A[max][k]) 
      max = i;
  
  }   

  return max;
}

// linhas k e p
static void trocaLinha (SistLinear_t *C, int k, int p, int n)
{
  // Troca coeficientes das linhas
  real_t aux;
  for(int i = 0; i < C->n; i++){
    aux = C->A[k][i];
    C->A[k][i] = C->A[p][i];
    C->A[p][i] = aux;
  }

  // Troca termos independentes
  aux = C->b[k];
  C->b[k] = C->b[p];
  C->b[p] = aux;
}

/* Seja um S.L. de ordem 'n'
   C = A|B em Ax=B
 */
void triangulariza( SistLinear_t *C )
{
  for(int i = 0; i < C->n; i++){
    int pivo = encontraMax(C, i, C->n);
    if(i != pivo){
      trocaLinha(C, i, pivo, C->n);
    }

    for(int k = i+1; k < C->n; k++){
      for(int j = i+1; j < C->n; j++){
          C->A[k][j] = C->A[k][j]* C->A[i][i] - C->A[i][j] * C->A[k][i];
          fesetround(C->A[k][j]);
      }
      C->b[k] = C->b[k] * C->A[i][i] - C->b[i] * C->A[k][i];
      fesetround(C->b[k]);
      C->A[k][i] = 0.0;
    }
  }
}  


void retrosubst( SistLinear_t *C, real_t *X )
{
  for(int i = C->n-1; i >= 0; i--){
    real_t soma = 0.0;
    for(int j = C->n-1; j > i; j--){
      soma += C->A[i][j] * X[j];
      fesetround(soma);
    }
    X[i] = (C->b[i] - soma) / C->A[i][i];
    fesetround(X[i]);
  } 
}

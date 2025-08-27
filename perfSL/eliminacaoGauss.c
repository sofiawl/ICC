/* Matriz  'normal' (vetor  de  ponteiros (linhas  matriz) para  vetores
   (colunas da matriz), estilo 'Mazieiro/Prog 2'
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"
#include "eliminacaoGauss.h"

// linha n coluna k 
real_t encontraMax(SistLinear_t *C, int k, int n)
{
  real_t max = C->A[n][k]; 
  while(n--){
    if(C->A[n][k] > max) 
      max = C->A[n][k];
  }

  return max;
}

// linhas k e p
static void trocaLinha (SistLinear_t *C, int k, int p, int n)
{
  real_t aux;
  for(int i = 1; i <= C->n; i++){
    aux = C->A[k][i];
    C->A[k][i] = C->A[p][i];
    C->A[p][i] = aux;
  }

}

/* Seja um S.L. de ordem 'n'
   C = A|B em Ax=B
 */
void triangulariza( SistLinear_t *C )
{
  for(int i = 1; i <= C->n; i++){
    real_t pivo = 
    encontraMax(C, );
    if(){
      
    }

    for(int j = 1; <= C->n; j++){

    }
  }
  
}

void retrosubst( SistLinear_t *C, real_t *X )
{
 
}

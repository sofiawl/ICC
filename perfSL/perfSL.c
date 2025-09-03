#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "eliminacaoGauss.h"
#include "sislin.h"
#include "utils.h"
#include "gaussSeidel.h"

int main(){

    SistLinear_t *mtx = lerSisLin();

    printf("MATRIZ ORIGINAL\n");
    prnSisLin(mtx);

/*
    triangulariza(mtx);
    printf("MATRIZ TRIANGULARIZADA\n");
    prnSisLin(mtx);

    retrosubst(mtx, mtx->b);
    //printf("MATRIZ DEPOIS DA RETROSUBSTITUIÇÃO\n");
    //prnSisLin(mtx);
    printf("RESULTADO\n");
    prnVetor(mtx->b, mtx->n);
*/

    real_t *X = malloc(sizeof(real_t) * mtx->n);
    if (!X) {return 1;}  
    for(int i = 0; i < mtx->n; i++){
        X[i] = 0.0;
    } 

    real_t norma;
    if (gaussSeidel(mtx, X, DBL_EPSILON, 60, &norma)) {return 1;}

    printf("MATRIZ DEPOIS DE GAUSS-SEIDEL\n");
    prnSisLin(mtx);
    printf("RESULTADO\n");
    prnVetor(X, mtx->n);

    liberaSisLin(mtx);
    return 0;
}
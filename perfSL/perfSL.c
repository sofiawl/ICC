#include <stdio.h>
#include <stdlib.h>
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

    real_t norma = normaL2(mtx->b, mtx->n);

    real_t *X = malloc(sizeof(real_t) * mtx->n);
    if (!X) {return 1;}  
    for(int i = 0; i < mtx->n; i++){
        X[i] = 0.0;
    } 

    gaussSeidel(mtx, X, __FLT_EPSILON__, 60, &norma);

    liberaSisLin(mtx);
    return 0;
}
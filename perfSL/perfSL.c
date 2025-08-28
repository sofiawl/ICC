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

    triangulariza(mtx);
    printf("MATRIZ TRIANGULARIZADA\n");
    //prnSisLin(mtx);

    //liberaSisLin(mtx);
    return 0;
}
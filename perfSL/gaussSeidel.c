#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"
#include "gaussSeidel.h"

int criterioSassenfeld(SistLinear_t *C){

    real_t soma = 0.0;
    real_t bmax = 0.0;

    // Verifica se há convergência pelo método de Sassenfeld
    for(int i = 0; i < C->n; i++){
        real_t b;
        soma = 0.0;
        for(int j = 0; j <= i-1; j++){
            soma += fabs(C->A[i][j])*b;
        }
        for(int j = i+1; j < C->n; j++){
            soma += fabs(C->A[i][j]);
        }
        b = soma / fabs(C->A[i][i]);
        if (b > bmax)
            bmax = b;
    }

    if (bmax < 1.0){
        printf("Matriz fornecida converge pelo critério de Sassenfeld\n");
        return 1;
    }
    else{
        printf("Matriz fornecida não converge pelo critério de Sassenfeld\n");
        return 0;
    }
}

// Método de Gauss-Seidel clássico
int gaussSeidel (SistLinear_t *C, real_t *X, real_t erro, int maxit, real_t *norma)
{
    if (!criterioSassenfeld(C)) {return 1;}

    while(maxit && &&){
        for(int i = 0; i < C->n; i++){
            X[i] = C->b[i];
            for(int j = 0; j <= i-1; j++){
                X[i] -= C->A[i][j]*X[j];
            }
            for(int j = i+1; j < C->n; j++){
                X[i] -= C->A[i][j]*X[j];
            }

            X[i] /= C->A[i][i]; 
        }

        maxit--;
    }

    return 0;
}



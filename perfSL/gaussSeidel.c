#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "utils.h"
#include "sislin.h"
#include "gaussSeidel.h"

#define ZERO __FLT_EPSILON__

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
        printf("Matriz fornecida converge pelo critério de Sassenfeld\n\n");
        return 1;
    }
    else{
        printf("Matriz fornecida não converge pelo critério de Sassenfeld\n\n");
        return 0;
    }
}

// Método de Gauss-Seidel clássico
int gaussSeidel (SistLinear_t *C, real_t *X, real_t erro, int maxit, real_t *norma)
{
    if (!criterioSassenfeld(C)) {return 1;}

    real_t *Xant = malloc(sizeof(real_t)*C->n);
    if (!Xant) {return 1;}

    real_t *R = malloc(sizeof(real_t)*C->n);
    if (!R) {return 1;}

    real_t res;

    do {
        for (int i = 0; i < C->n; i++) {Xant[i] = X[i];}

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
        *norma = normaMax(X, Xant, C->n);

        residuo(C, X, R, C->n);
        res = normaL2(R, C->n);

    } while(maxit && (*norma > erro) && (res > erro));

    free(Xant);
    free(R);
    return 0;
}



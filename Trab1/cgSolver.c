#include <stdlib.h>
#include <stdio.h>
#include "sislin.h"
#include "utils.h"
#include "pcgc.h"

int main(){
    srandom(20252);

    int n, k, maxit;
    real_t w, e;
    rtime_t tempo = timestamp();

    // n k w maxit, erromax
    // se n for menor que 2 falhar
    scanf("%d", &n);
    scanf("%d", &k); 
    scanf("%lf", &w);
    scanf("%d", &maxit);
    //erro max absoluto, considerando norma maxima em x ( max ( |xi - xi-1| ) < ε )
    scanf("%lf", &e);

    real_t *A, *b;
    criaKDiagonal(n, k, &A, &b);
    prnKDiagonal(n, k, A, b);

    real_t *ASP = malloc(n * k * sizeof(real_t));
    if(!ASP){
        fprintf(stderr, "Erro ao alocar memória.\n");
        exit(EXIT_FAILURE);
    }

    real_t *bsp = malloc(n * sizeof(real_t));
    if(!bsp){
        fprintf(stderr, "Erro ao alocar memória.\n");
        exit(EXIT_FAILURE);
    }

    genSimetricaPositiva(A, b, n, k, ASP, bsp, &tempo);
    printf("\n");
    prnKDiagonal(n, k, ASP, bsp);
    //calcResiduoSL(A, b, X, r, n, k, tempo);
    
    destroiKDiagonal(n, k, A, b);
    destroiKDiagonal(n, k, ASP, bsp);

    return 0;
}
#include <stdio.h>
#include <stdlib.h>
#include <likwid.h>
#include <float.h>
#include <sys/time.h>
#include <math.h>

#include "gaussSeidelEDO.h"
#include "sislin.h"
#include "utils.h"

int main(){

    EDo edo;

    //1ª linha: quantidade de pontos da malha da EDO;
    scanf("%d", &edo.n);

    //2ª linha: intervalo a e b onde a EDO é válida;
    scanf("%d", &edo.a);
    scanf("%d", &edo.b);

    //3ª linha: os valores de contorno  y(a) e y(b);
    scanf("%d", &edo.ya);
    scanf("%d", &edo.yb);

    //4ª linha: os coeficientes p e q da EDO genérica;
    scanf("%d", &edo.p);
    scanf("%d", &edo.q);

    real_t *Y = malloc(sizeof(real_t) * edo.n);
    for(int i = 0; i < edo.n; i++) Y[i] = 0.0;

    //5ª linha em diante: uma ou mais linhas contendo os coeficientes r1, r2, r3  e r4 da definição da função r(x),
    //representando diversas  EDO's que diferem apenas no valor de r(x).
    while(){
        
        gaussSeidel_EDO(edo, Y, 100);
        prnVetor(Y);
    }

    free(Y);

    return 0;
}

// r(x) = r1x + r2x² + r3cos(x) + r4ex
real_t r1(real_t x){
    return x;
}

real_t r2(real_t x){
    return x*x;
}

real_t r3(real_t x){
    return cos(x);
}

real_t r4(real_t x){
    return exp(x);
}
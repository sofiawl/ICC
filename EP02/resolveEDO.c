#include <stdio.h>
#include <stdlib.h>
#include <likwid.h>
#include <float.h>
#include <sys/time.h>


// Equação Diferencial Ordinária
typedef struct {
 int n; // número de pontos internos na malha
 real_t a, b; // intervalo
 real_t ya, yb; // condições contorno
 real_t (* p)(real_t), (* q)(real_t), (* r)(real_t);
} EDo;

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

    //5ª linha em diante: uma ou mais linhas contendo os coeficientes r1, r2, r3  e r4 da definição da função r(x),
    //representando diversas  EDO's que diferem apenas no valor de r(x).

    return 0;
}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <likwid.h>
#include <float.h>
#include <sys/time.h>
#include <math.h>
#include "edo.h"
#include "utils.h"


int main(){
    LIKWID_MARKER_INIT;
    EDo edo;
    Tridiag *sl = malloc(sizeof(Tridiag));
    if (!sl) return 1; 

    //1ª linha: quantidade de pontos da malha da EDO;
    scanf("%d", &edo.n);

    //2ª linha: intervalo a e b onde a EDO é válida;
    scanf("%lf", &edo.a);
    scanf("%lf", &edo.b);

    //3ª linha: os valores de contorno  y(a) e y(b);
    scanf("%lf", &edo.ya);
    scanf("%lf", &edo.yb);

    //4ª linha: os coeficientes p e q da EDO genérica;
    scanf("%lf", &edo.p);
    scanf("%lf", &edo.q);

    real_t *Y = malloc(sizeof(real_t) * edo.n);
    if(!Y) return 1;

    real_t normaL2 = 0.0;
    //5ª linha em diante: uma ou mais linhas contendo os coeficientes r1, r2, r3  e r4 da definição da função r(x),
    //representando diversas  EDO's que diferem apenas no valor de r(x).
    while(scanf("%lf %lf %lf %lf", &edo.r1, &edo.r2, &edo.r3, &edo.r4) == 4)
    {
        memset(Y, 0, sizeof(real_t) * edo.n);

        unsigned int iter = MAXIT;
        sl = genTridiag(&edo);
    
        LIKWID_MARKER_START("GS3D");
        rtime_t time = gaussSeidel_3Diag(sl, Y, &iter, &normaL2);
		    LIKWID_MARKER_STOP("GS3D");

        prnEDOsl(&edo);

        printf("\n");
        prnVetor(Y, sl->n);

        printf("%d\n", iter);
        printf("%23.15e\n", normaL2);
        printf("%16.8e\n", time);
    }
	
    free(Y);
    free(sl);
		LIKWID_MARKER_CLOSE;
    return 0;
}

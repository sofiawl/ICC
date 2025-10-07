#include <stdlib.h>
#include "utils.h"
#include "sislin.h"

// Sem pre-condicionador
inline real_t** Sempc(){

}

// Pre-condicionador de Jacobi
inline real_t** pcJacobi(){

}

// Pre-condicionador de Gauss Seidel
real_t** pcGausSeidel(){

}

// Pré-condicionador de SSOR
real_t** pcSSOR(){

}

//Calcula a solução da matriz usando o método do Gradiente Conjugado
void resolveGC(real_t *A, real_t *b, real_t **M, real_t *X, int maxit, int e, real_t w, real_t n, real_t k, rtime_t *tempo){
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

    // Ax = b, onde A é simétrica positiva

}


/**
 * Devolve matriz M⁻¹
 *
 */
void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k,
		 real_t **M, rtime_t *tempo)
{
  *tempo = timestamp();

    if (w == -1.0) {
        // case -1.0, sem pré-condiconador
        //*M = M⁻1;
        // M = genSimetricaPositiva();
        return;
    }
    else if (w == 0.0) {
        // case 0.0, pré-condicionador de Jacobi
        //*M = D;
        // M = genSimetricaPositiva();
    }
    else if (w == 1.0) {
        // case 1.0, pré-condicionador de Gauss Seidel
        // M = (D + ωL)D−1 (D + ωU )
        // M = genSimetricaPositiva();
    }
    else if (w > 1.0 && w < 2.0) { 
        // case 1.0 < w < 2.0, pré-condicionador SSOR
        // M = (D + ωL)D−1 (D + ωU )
        // M = genSimetricaPositiva();
    }
    else {
        // default logic

    }

  *tempo = timestamp() - *tempo;
}
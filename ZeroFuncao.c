#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"
//DOUBLE E ARREDONDAMENTO PARA BAIXO

/*
// Aproximação aceitável como valor zero
#define ZERO DBL_EPSILON

// Parâmetros para teste de convergência
#define MAXIT 500
#define EPS 1.0e-7
#define ULPS 3

typedef struct {
  real_t *p;
  int grau;
} Polinomio;
*/

/*  Retorna tempo em milisegundos

    Forma de uso:
 
    double tempo;
    tempo = timestamp();
    <trecho de programa do qual se deseja medir tempo>
    tempo = timestamp() - tempo;
*/


// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t newtonRaphson (Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz)
{
    switch(criterioParada){
        // Criterio-01: |xk - xk-1| / |xk| <= 10 ^-17
        case 1: 

            break;
        // Criterio-02: |f(xk)| <= DLB_EPSOLON
        case 2:

            break;
        // Criterio-03: ULP's entre xk e xk-1 <= 3
        case 3:

            break;
    }
}


// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao (Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz, int calcPolinomio)
{
    real_t c_old, c_new;
    real_t xa, xm;
    real_t erro;
    
    switch (calcPolinomio)
    {
        case 1: // Calculo Rápido
            c_new = (a + b) / 2;
            calcPolinomio_rapido(p, a, &xa, NULL);
            calcPolinomio_rapido(p, c_new, &xm, NULL);
            if(xa*xm < 0)
                a = c_new;
            else if(xa*xm > 0)
                b = c_new;
            else {
                *raiz = c_new;
                return 0;
            }

            switch(criterioParada){
                // Criterio-01: |xk - xk-1| / |xk| <= 10 ^-17
                case 1: 
                    do{
                        c_old = c_new;
                        c_new = (a + b)/ 2;
                        calcPolinomio_rapido(p, a, &xa, NULL);
                        calcPolinomio_rapido(p, c_new, &xm, NULL);
                        if(xa*xm < 0)
                            a = c_new;
                        else if(xa*xm > 0)
                            b = c_new;
                        else {
                            *raiz = c_new;
                            return 0;
                        }
                        (*it)++;
                        erro = (c_new - c_old) / c_new;
                    }while(erro > EPS);
                    break;
                // Criterio-02: |f(xk)| <= DLB_EPSOLON
                case 2:
                    do{ 
                        c_old = c_new;
                        c_new = (a + b)/ 2;
                        calcPolinomio_rapido(p, a, &xa, NULL);
                        calcPolinomio_rapido(p, c_new, &xm, NULL);
                        if(xa*xm < 0)
                            a = c_new;
                        else if(xa*xm > 0)
                            b = c_new;
                        else {
                            *raiz = c_new;
                            return 0;
                        }
                        (*it)++;
                        erro = xm;
                    }while(xm > ZERO);
                    break;
                // Criterio-03: ULP's entre xk e xk-1 <= 3
                case 3:
                    do {    
                        c_old = c_new;
                        c_new = (a + b)/ 2;
                        calcPolinomio_rapido(p, a, &xa, NULL);
                        calcPolinomio_rapido(p, c_new, &xm, NULL);
                        if(xa*xm < 0)
                            a = c_new;
                        else if(xa*xm > 0)
                            b = c_new;
                        else {
                            *raiz = c_new;
                            return 0;
                        }
                        (*it)++;
                        erro = xm;
                    }while(xm > ZERO);

                    break;
            }
            break;

        case 2: // Cálculo Lento
            c_new = (a + b) / 2;
            calcPolinomio_lento(p, a, &xa, NULL);
            calcPolinomio_lento(p, c_new, &xm, NULL);
            if(xa*xm > 0)
                a = c_new;
            else if(xa*xm < 0)
                b = c_new;
            else {
                *raiz = c_new;
                return 0;
            }

            switch(criterioParada){
                // Criterio-01: |xk - xk-1| / |xk| <= 10 ^-17
                case 1: 

                    break;
                // Criterio-02: |f(xk)| <= DLB_EPSOLON
                case 2:

                    break;
                // Criterio-03: ULP's entre xk e xk-1 <= 3
                case 3:

                    break;
            }
            break;

    }    
    
}


void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx)
{
    real_t b = 0;
    real_t c = 0;
    for(int i = p.grau; i > 0; --i) {
         b = b*x + p.p[i];
         c = c*x + b;
    }
    b = b*x + p.p[0];
    *px = b;
    *dpx = c;
}


void calcPolinomio_lento(Polinomio p, real_t x, real_t *px, real_t *dpx)
{
    real_t *px = 0;
    for(int i = p.grau; i >= 0; --i){
        *px += p.p[i]*pow(x, i);
        *dpx += p.p[i]*i*pow(x, i-1);
    }
}

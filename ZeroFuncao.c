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

void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx)
{
    *px = 0;
    *dpx = 0;
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
    *px = 0;
    *dpx = 0;
    for(int i = p.grau; i >= 0; --i){
        *px += p.p[i]*pow(x, i);
        *dpx += p.p[i]*i*pow(x, i-1);
    }
}


// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t newtonRaphson (Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz, int calcPolinomio)
{
    real_t x_new, dx, fx_old, fx_new;
    real_t erro = 0;
    switch (calcPolinomio)
    {
        case 1: // Cálculo Rápido
            switch(criterioParada){
                // Criterio-01: |xk - xk-1| / |xk| <= 10 ^-17
                case 1: 
                    do{
                        calcPolinomio_rapido(p, x0, &fx_old, &dx);
                        x_new = x0 - (fx_old / dx);
                        calcPolinomio_rapido(p, x_new, &fx_new, NULL);
                        erro = fabs((fx_new - fx_old) / fx_new);
                        x0 = x_new;
                    }while(erro > EPS);
                    break;
                // Criterio-02: |f(xk)| <= DLB_EPSOLON
                case 2:
                    do{
                        calcPolinomio_rapido(p, x0, &fx_old, &dx);
                        x_new = x0 - (fx_old / dx);
                        calcPolinomio_rapido(p, x_new, &fx_new, NULL);
                        erro = fabs(fx_new);
                        x0 = x_new;
                    }while(erro > ZERO);
                    break;
                // Criterio-03: ULP's entre xk e xk-1 <= 3
                case 3:
                    do{
                        calcPolinomio_rapido(p, x0, &fx_old, &dx);
                        x_new = x0 - (fx_old / dx);
                        calcPolinomio_rapido(p, x_new, &fx_new, NULL);
                        erro = fabs((x_new - x0) - 1);
                        x0 = x_new;
                    }while(erro > ULPS);
                    break;
            }

            break;

        case 2: // Cálculo Lento
            switch(criterioParada){
                // Criterio-01: |xk - xk-1| / |xk| <= 10 ^-17
                case 1: 
                    do{
                        calcPolinomio_lento(p, x0, &fx_old, &dx);
                        x_new = x0 - (fx_old / dx);
                        calcPolinomio_lento(p, x_new, &fx_new, NULL);
                        erro = fabs((fx_new - fx_old) / fx_new);
                        x0 = x_new;
                    }while(erro > EPS);
                    break;
                // Criterio-02: |f(xk)| <= DLB_EPSOLON
                case 2:
                    do{
                        calcPolinomio_lento(p, x0, &fx_old, &dx);
                        x_new = x0 - (fx_old / dx);
                        calcPolinomio_lento(p, x_new, &fx_new, NULL);
                        erro = fabs(fx_new);
                        x0 = x_new;
                    }while(erro > ZERO);
                    break;
                // Criterio-03: ULP's entre xk e xk-1 <= 3
                case 3:
                    do{
                        calcPolinomio_lento(p, x0, &fx_old, &dx);
                        x_new = x0 - (fx_old / dx);
                        calcPolinomio_lento(p, x_new, &fx_new, NULL);
                        erro = fabs((x_new - x0) - 1);
                        x0 = x_new;
                    }while(erro > ULPS);
                    break;
            }

            break;        
    }

    return erro;
    
}


// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao (Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz, int calcPolinomio)
{
    real_t c_old, c_new;
    real_t fxa, fxm;
    real_t erro = 0;
    
    switch (calcPolinomio)
    {
        case 1: // Cálculo Rápido
            c_new = (a + b) / 2;
            calcPolinomio_rapido(p, a, &fxa, NULL);
            calcPolinomio_rapido(p, c_new, &fxm, NULL);
            if(fxa*fxm < 0)
                a = c_new;
            else if(fxa*fxm > 0)
                b = c_new;
            else {
                *raiz = c_new;
                return erro;
            }

            switch(criterioParada){
                // Criterio-01: |xk - xk-1| / |xk| <= 10 ^-17
                case 1: 
                    do{
                        c_old = c_new;
                        c_new = (a + b)/ 2;
                        calcPolinomio_rapido(p, a, &fxa, NULL);
                        calcPolinomio_rapido(p, c_new, &fxm, NULL);
                        if(fxa*fxm < 0)
                            a = c_new;
                        else if(fxa*fxm > 0)
                            b = c_new;
                        else {
                            *raiz = c_new;
                            return 0;
                        }
                        (*it)++;
                        erro = fabs((c_new - c_old) / c_new);
                    }while(erro > EPS);
                    break;
                // Criterio-02: |f(xk)| <= DLB_EPSOLON
                case 2:
                    do{ 
                        c_old = c_new;
                        c_new = (a + b)/ 2;
                        calcPolinomio_rapido(p, a, &fxa, NULL);
                        calcPolinomio_rapido(p, c_new, &fxm, NULL);
                        if(fxa*fxm < 0)
                            a = c_new;
                        else if(fxa*fxm > 0)
                            b = c_new;
                        else {
                            *raiz = c_new;
                            return 0;
                        }
                        (*it)++;
                        erro = fabs(fxm);
                    }while(erro > ZERO);
                    break;
                // Criterio-03: ULP's entre xk e xk-1 <= 3
                case 3:
                    do {    
                        c_old = c_new;
                        c_new = (a + b)/ 2;
                        calcPolinomio_rapido(p, a, &fxa, NULL);
                        calcPolinomio_rapido(p, c_new, &fxm, NULL);
                        if(fxa*fxm < 0)
                            a = c_new;
                        else if(fxa*fxm > 0)
                            b = c_new;
                        else {
                            *raiz = c_new;
                            return 0;
                        }
                        (*it)++;
                        erro = fabs((fxa - fxm) - 1);
                    }while(erro > ULPS);

                    break;
            }
            break;

        case 2: // Cálculo Lento
            c_new = (a + b) / 2;
            calcPolinomio_lento(p, a, &fxa, NULL);
            calcPolinomio_lento(p, c_new, &fxm, NULL);
            if(fxa*fxm > 0)
                a = c_new;
            else if(fxa*fxm < 0)
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
                        calcPolinomio_lento(p, a, &fxa, NULL);
                        calcPolinomio_lento(p, c_new, &fxm, NULL);
                        if(fxa*fxm < 0)
                            a = c_new;
                        else if(fxa*fxm > 0)
                            b = c_new;
                        else {
                            *raiz = c_new;
                            return erro;
                        }
                        (*it)++;
                        erro = fabs((c_new - c_old) / c_new);
                    }while(erro > EPS);
                    break;
                // Criterio-02: |f(xk)| <= DLB_EPSOLON
                case 2:
                    do{ 
                        c_old = c_new;
                        c_new = (a + b)/ 2;
                        calcPolinomio_lento(p, a, &fxa, NULL);
                        calcPolinomio_lento(p, c_new, &fxm, NULL);
                        if(fxa*fxm < 0)
                            a = c_new;
                        else if(fxa*fxm > 0)
                            b = c_new;
                        else {
                            *raiz = c_new;
                            return 0;
                        }
                        (*it)++;
                        erro = fabs(fxm);
                    }while(erro > ZERO);
                    break;
                // Criterio-03: ULP's entre xk e xk-1 <= 3
                case 3:
                    do {    
                        c_old = c_new;
                        c_new = (a + b)/ 2;
                        calcPolinomio_lento(p, a, &fxa, NULL);
                        calcPolinomio_lento(p, c_new, &fxm, NULL);
                        if(fxa*fxm < 0)
                            a = c_new;
                        else if(fxa*fxm > 0)
                            b = c_new;
                        else {
                            *raiz = c_new;
                            return 0;
                        }
                        (*it)++;
                        erro = fabs((fxa - fxm) - 1);
                    }while(erro > ULPS);

                    break;
            }
            break;

    }    
    
    return erro;
}


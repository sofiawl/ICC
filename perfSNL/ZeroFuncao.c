#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"


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
    real_t x_new, dpx_old, dpx_new, fx_old, fx_new;
    real_t erro = 0;
    *it = 0;

    switch (calcPolinomio)
    {
        case 1: // Cálculo Rápido
            switch(criterioParada){
                // Criterio-01: |xk - xk-1| / |xk| <= 10 ^-17
                case 1: 
                    do{
                        calcPolinomio_rapido(p, x0, &fx_old, &dpx_old);
                        x_new = x0 - (fx_old / dpx_old);
                        calcPolinomio_rapido(p, x_new, &fx_new, &dpx_new);
                        erro = fabs((x_new - x0) / x_new);
                        x0 = x_new;
                        (*it)++;
                        if (*it >= MAXIT){
                            *raiz = x0;
                            return erro;
                        }
                    }while(erro > EPS);
                    break;
                // Criterio-02: |f(xk)| <= DLB_EPSOLON
                case 2:
                    do{
                        calcPolinomio_rapido(p, x0, &fx_old, &dpx_old);
                        x_new = x0 - (fx_old / dpx_old);
                        calcPolinomio_rapido(p, x_new, &fx_new, &dpx_new);
                        erro = fabs(fx_new);
                        x0 = x_new;
                        (*it)++;
                        if (*it >= MAXIT){
                            *raiz = x0;
                            return erro;
                        }
                    }while(erro > ZERO);
                    break;
                // Criterio-03: ULP's entre xk e xk-1 <= 3
                case 3:
                    do{
                        calcPolinomio_rapido(p, x0, &fx_old, &dpx_old);
                        x_new = x0 - (fx_old / dpx_old);
                        calcPolinomio_rapido(p, x_new, &fx_new, &dpx_new);
                        erro = fabs((fx_new - x0) - 1);
                        x0 = x_new;
                        (*it)++;
                        if (*it >= MAXIT){
                            *raiz = x0;
                            return erro;
                        }
                    }while(erro > ULPS);
                    break;
            }

            break;

        case 2: // Cálculo Lento
            switch(criterioParada){
                // Criterio-01: |xk - xk-1| / |xk| <= 10 ^-17
                case 1: 
                    do{
                        calcPolinomio_lento(p, x0, &fx_old, &dpx_old);
                        x_new = x0 - (fx_old / dpx_old);
                        calcPolinomio_lento(p, x_new, &fx_new, &dpx_new);
                        erro = fabs((x_new - x0) / x_new);
                        x0 = x_new;
                        (*it)++;
                        if (*it >= MAXIT){
                            *raiz = x0;
                            return erro;
                        }
                    }while(erro > EPS);
                    break;
                // Criterio-02: |f(xk)| <= DLB_EPSOLON
                case 2:
                    do{
                        calcPolinomio_lento(p, x0, &fx_old, &dpx_old);
                        x_new = x0 - (fx_old / dpx_old);
                        calcPolinomio_lento(p, x_new, &fx_new, &dpx_new);
                        erro = fabs(fx_new);
                        x0 = x_new;
                        (*it)++;
                        if (*it >= MAXIT){
                            *raiz = x0;
                            return erro;
                        }
                    }while(erro > ZERO);
                    break;
                // Criterio-03: ULP's entre xk e xk-1 <= 3
                case 3:
                    do{
                        calcPolinomio_lento(p, x0, &fx_old, &dpx_old);
                        x_new = x0 - (fx_old / dpx_old);
                        calcPolinomio_lento(p, x_new, &fx_new, &dpx_new);
                        erro = fabs((fx_new - x0) - 1);
                        x0 = x_new;
                        (*it)++;
                        if (*it >= MAXIT){
                            *raiz = x0;
                            return erro;
                        }
                    }while(erro > ULPS);
                    break;
            }

            break;        
    }

    *raiz = x0;
    return erro;
    
}


// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao (Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz, int calcPolinomio)
{
    real_t c_old, c_new, dpx_a, dpx_m;
    real_t fxa, fxm;
    real_t erro = 0;
    *it = 0;

    switch (calcPolinomio)
    {
        case 1: // Cálculo Rápido
            c_new = (a + b) / 2;
            calcPolinomio_rapido(p, a, &fxa, &dpx_a);
            calcPolinomio_rapido(p, c_new, &fxm, &dpx_m);
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
                        calcPolinomio_rapido(p, a, &fxa, &dpx_a);
                        calcPolinomio_rapido(p, c_new, &fxm, &dpx_m);
                        if(fxa*fxm < 0)
                            a = c_new;
                        else if(fxa*fxm > 0)
                            b = c_new;
                        else {
                            *raiz = c_new;
                            return 0;
                        }
                        erro = fabs((c_new - c_old) / c_new);
                        (*it)++;
                        if (*it >= MAXIT){
                            *raiz = c_new;
                            return erro;
                        }
                    }while(erro > EPS);
                    break;
                // Criterio-02: |f(xk)| <= DLB_EPSOLON
                case 2:
                    do{ 
                        c_old = c_new;
                        c_new = (a + b)/ 2;
                        calcPolinomio_rapido(p, a, &fxa, &dpx_a);
                        calcPolinomio_rapido(p, c_new, &fxm, &dpx_m);
                        if(fxa*fxm < 0)
                            a = c_new;
                        else if(fxa*fxm > 0)
                            b = c_new;
                        else {
                            *raiz = c_new;
                            return 0;
                        }
                        erro = fabs(fxm);
                        (*it)++;
                        if (*it >= MAXIT){
                            *raiz = c_new;
                            return erro;
                        }
                    }while(erro > ZERO);
                    break;
                // Criterio-03: ULP's entre xk e xk-1 <= 3
                case 3:
                    do {    
                        c_old = c_new;
                        c_new = (a + b)/ 2;
                        calcPolinomio_rapido(p, a, &fxa, &dpx_a);
                        calcPolinomio_rapido(p, c_new, &fxm, &dpx_m);
                        if(fxa*fxm < 0)
                            a = c_new;
                        else if(fxa*fxm > 0)
                            b = c_new;
                        else {
                            *raiz = c_new;
                            return 0;
                        }
                        erro = fabs((fxa - fxm) - 1);
                        (*it)++;
                        if (*it >= MAXIT){
                            *raiz = c_new;
                            return erro;
                        }
                    }while(erro > ULPS);

                    break;
            }
            break;

        case 2: // Cálculo Lento
            c_new = (a + b) / 2;
            calcPolinomio_lento(p, a, &fxa, &dpx_a);
            calcPolinomio_lento(p, c_new, &fxm, &dpx_m);
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
                        calcPolinomio_lento(p, a, &fxa, &dpx_a);
                        calcPolinomio_lento(p, c_new, &fxm, &dpx_m);
                        if(fxa*fxm < 0)
                            a = c_new;
                        else if(fxa*fxm > 0)
                            b = c_new;
                        else {
                            *raiz = c_new;
                            return erro;
                        }
                        erro = fabs((c_new - c_old) / c_new);
                        (*it)++;
                        if (*it >= MAXIT){
                            *raiz = c_new;
                            return erro;
                        }
                    }while(erro > EPS);
                    break;
                // Criterio-02: |f(xk)| <= DLB_EPSOLON
                case 2:
                    do{ 
                        c_old = c_new;
                        c_new = (a + b)/ 2;
                        calcPolinomio_lento(p, a, &fxa, &dpx_a);
                        calcPolinomio_lento(p, c_new, &fxm, &dpx_m);
                        if(fxa*fxm < 0)
                            a = c_new;
                        else if(fxa*fxm > 0)
                            b = c_new;
                        else {
                            *raiz = c_new;
                            return 0;
                        }
                        erro = fabs(fxm);
                        (*it)++;
                        if (*it >= MAXIT){
                            *raiz = c_new;
                            return erro;
                        }
                    }while(erro > ZERO);
                    break;
                // Criterio-03: ULP's entre xk e xk-1 <= 3
                case 3:
                    do {    
                        c_old = c_new;
                        c_new = (a + b)/ 2;
                        calcPolinomio_lento(p, a, &fxa, &dpx_a);
                        calcPolinomio_lento(p, c_new, &fxm, &dpx_m);
                        if(fxa*fxm < 0)
                            a = c_new;
                        else if(fxa*fxm > 0)
                            b = c_new;
                        else {
                            *raiz = c_new;
                            return 0;
                        }
                        erro = fabs((fxa - fxm) - 1);
                        (*it)++;
                        if (*it >= MAXIT){
                            *raiz = c_new;
                            return erro;
                        }
                    }while(erro > ULPS);

                    break;
            }
            break;

    }    

    *raiz = c_new;
    return erro;
}


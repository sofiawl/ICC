//------------------------------------------------------------

// Instodução a Computação Científica
// EP-O1 - SOLUÇÕES DE EQUAÇÕES NÃO-LINEARES
// Programa feito po Sofia Wamser Lima, 14-19/08/2025

//------------------------------------------------------------
// labZero.h com funções auxiliares: aplicação dos métodos
// Bisseção e Newton-Raphson e cálculos de polinômios ensinados

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <fenv.h>

#include "utils.h"
#include "ZeroFuncao.h"

#define CRITERIOPARADA01 1
#define CRITERIOPARADA02 2
#define CRITERIOPARADA03 3

#define CALCULORAPIDO 1
#define CALCULOLENTO 2

#define BISSECCAO 1
#define NEWTONRAPHSON 2

// Função auxiliar para imprimir dos dados
// Imprime respectivamente: MÉTODO, RAIZ, ERRO, ITERAÇÕES, TEMPO
void imprimir_solucoes(int velocidade, int metodo, real_t raiz, real_t erro, int it, real_t tempo){
    
    switch(metodo)
    {
      case BISSECCAO:
        printf("bissec  %.15e %.15e %4d  %.8e\n", raiz, erro, it, tempo);
        break;

      case NEWTONRAPHSON:
        printf("newton  %.15e %.15e %4d  %.8e\n", raiz, erro, it, tempo);
        break;
    }

}


int main ()
{

  real_t a, b;
  Polinomio pol;

  scanf("%d", &pol.grau);

  pol.p = malloc((pol.grau + 1) * sizeof(real_t));
  if (!pol.p) return 1;

  for (int i = pol.grau; i >=0; --i)
    scanf("%lf", &pol.p[i]);

  scanf("%lf %lf", &a, &b); 

  if (fesetround(FE_DOWNWARD) != 0) 
    printf("Falha ao setar arredondamento.\n");

  real_t raiz, tempo, erro;
  int it = 0;
  real_t x0 = a + ((real_t)rand() / RAND_MAX) * (b - a);

  // A raiz do polinômio é calculada utilizando os métodos: Bisseção e Newton-Raphson
  // São aplicados três critérios de parada para cada método

  //------------------------------------------------------------
  // Execução com o método de cálculo de polinômio *mais* eficiente 

  printf("\nRAPIDO\n\n");
  tempo = timestamp();
  erro = bisseccao(pol, a, b, CRITERIOPARADA01, &it, &raiz, CALCULORAPIDO);
  tempo = timestamp() - tempo;
  imprimir_solucoes(CALCULORAPIDO, BISSECCAO, raiz, erro, it, tempo);

  tempo = timestamp();
  erro = bisseccao(pol, a, b, CRITERIOPARADA02, &it, &raiz, CALCULORAPIDO);
  tempo = timestamp() - tempo;
  imprimir_solucoes(CALCULORAPIDO, BISSECCAO, raiz, erro, it, tempo);

  tempo = timestamp();
  erro = bisseccao(pol, a, b, CRITERIOPARADA03, &it, &raiz, CALCULORAPIDO);
  tempo = timestamp() - tempo;
  imprimir_solucoes(CALCULORAPIDO, BISSECCAO, raiz, erro, it, tempo);
  
  tempo = timestamp(); // nan no exemplo 03
  erro = newtonRaphson(pol, x0, CRITERIOPARADA01, &it, &raiz, CALCULORAPIDO);
  tempo = timestamp() - tempo;
  imprimir_solucoes(CALCULORAPIDO, NEWTONRAPHSON, raiz, erro, it, tempo);
  
  tempo = timestamp();
  erro = newtonRaphson(pol, x0, CRITERIOPARADA02, &it, &raiz, CALCULORAPIDO);
  tempo = timestamp() - tempo;
  imprimir_solucoes(CALCULORAPIDO, NEWTONRAPHSON, raiz, erro, it, tempo);
  
  tempo = timestamp();
  erro = newtonRaphson(pol, x0, CRITERIOPARADA03, &it, &raiz, CALCULORAPIDO);
  tempo = timestamp() - tempo;
  imprimir_solucoes(CALCULORAPIDO, NEWTONRAPHSON, raiz, erro, it, tempo);
  
  //------------------------------------------------------------
  // Execução com o método de cálculo de polinômio *menos* eficiente

  printf("\n\nLENTO\n\n");
  tempo = timestamp();
  erro = bisseccao(pol, a, b, CRITERIOPARADA01, &it, &raiz, CALCULORAPIDO);
  tempo = timestamp() - tempo;
  imprimir_solucoes(CALCULOLENTO, BISSECCAO, raiz, erro, it, tempo);

  tempo = timestamp();
  erro = bisseccao(pol, a, b, CRITERIOPARADA02, &it, &raiz, CALCULORAPIDO);
  tempo = timestamp() - tempo;
  imprimir_solucoes(CALCULOLENTO, BISSECCAO, raiz, erro, it, tempo);
  
  tempo = timestamp();
  erro = bisseccao(pol, a, b, CRITERIOPARADA03, &it, &raiz, CALCULORAPIDO);
  tempo = timestamp() - tempo;
  imprimir_solucoes(CALCULOLENTO, BISSECCAO, raiz, erro, it, tempo);
  
  tempo = timestamp();
  erro = newtonRaphson(pol, x0, CRITERIOPARADA01, &it, &raiz, CALCULORAPIDO);
  tempo = timestamp() - tempo;
  imprimir_solucoes(CALCULOLENTO, NEWTONRAPHSON, raiz, erro, it, tempo);
  
  tempo = timestamp();
  erro = newtonRaphson(pol, x0, CRITERIOPARADA02, &it, &raiz, CALCULORAPIDO);
  tempo = timestamp() - tempo;
  imprimir_solucoes(CALCULOLENTO, NEWTONRAPHSON, raiz, erro, it, tempo);
  
  tempo = timestamp();
  erro = newtonRaphson(pol, x0, CRITERIOPARADA03, &it, &raiz, CALCULORAPIDO);
  tempo = timestamp() - tempo;
  imprimir_solucoes(CALCULOLENTO, NEWTONRAPHSON, raiz, erro, it, tempo);

  free(pol.p);

  return 0;
}


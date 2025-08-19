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

  scanf("%lf %lf", &a, &b); // intervalo onde está uma das raizes.

  if (fesetround(FE_DOWNWARD) != 0) 
    printf("Falha ao setar arredondamento.\n");

  real_t raiz, tempo, erro;
  int it = 0;
  real_t x0 = a + ((real_t)rand() / RAND_MAX) * (b - a);

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
  
  tempo = timestamp();
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


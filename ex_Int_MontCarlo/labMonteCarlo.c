#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"

#define DIFF 0.0

#define NRAND    ((real_t) random() / RAND_MAX)  // drand48() 
#define SRAND(a) srandom(a) // srand48(a)
#define RAND_AB(a, b) ((a) + NRAND * ((b) - (a))) // random entre [a,b]

inline real_t f(real_t *x, int dim){
  real_t resultado = 0.0;
  
  for (int i=0; i<dim; i++){
    resultado += x[i]*x[i]*x[i]*x[i];
    resultado -= 16*x[i]*x[i];
    resultado += 5*x[i];
    resultado /= 2;
  }

  return resultado;
}

// Integral Monte Carlo da função Styblinski-Tang de 2 variáveis
real_t styblinskiTang(real_t a, real_t b, int namostras, int dim)
{
  if (dim != 2 || dim != 4 || dim != 8){
    printf("Dimensões só podem ser: 2, 4 ou 8\n");
    return 0.0;
  }
    
  real_t resultado;
  real_t x, soma = 0.0;
  
  printf("Metodo de Monte Carlo (x, y).\n");
  printf("a = (%f), b = (%f), n = (%d), variaveis = 2\n", a, b, namostras);
  
  rtime_t t_inicial = timestamp();

  real_t areaAB = 1.0;
  real_t *aleats = malloc(sizeof(real_t)*dim);

  
  // 2, 4, ou 8 DIMENSÕES
  for (int i=0; i<dim; i++){
    areaAB *= (b-a);
    aleats[i] = RAND_AB(a, b);
  }

  if (dim =2){
    for (int i=0; i<namostras; i++){
      // x = a + h
      // 'h' é uma porcentagem aleat do tamanho do intervalo [a, b]
      x = a + NRAND * areaAB;
      soma += f(aleats, dim); 
    }
  }

  resultado = (soma / namostras) * areaAB;
  
  rtime_t t_final = timestamp();
  printf("Tempo decorrido: %f seg.\n", t_final - t_inicial);
  
  free(aleats);
  return resultado;
}

// 2 DIMENSÕES
real_t retangulos_xy(real_t a, real_t b, int npontos) {

  real_t k = ceil(sqrt(npontos));
  real_t h = (b-a)/ k;
  real_t resultado;
  real_t soma = 0;
  real_t *x, *y;
  y = malloc(sizeof(real_t)*k);
  x = malloc(sizeof(real_t)*k);
  for (int i=0; i<k; i++){
    y[i] = RAND_AB(a, b);
    x[i] = RAND_AB(a, b);
  }

  printf("Metodo dos Retangulos (x, y).\n");
  printf("a = (%f), b = (%f), n = (%d), h = (%lg)\n", a, b, npontos, h);
  
  rtime_t t_inicial = timestamp();
  
  for (int i=1; i<npontos; i++)
    for (int j=1; j<npontos; j++){
      resultado += f(x, npontos); 
    }
      
  resultado *= h*h;

  rtime_t t_final = timestamp();
  printf("Tempo decorrido: %f seg.\n", t_final - t_inicial);
  
  return resultado;
}


int main(int argc, char **argv) {

  if (argc < 5) {
    printf("Utilização: %s inicial final n_amostras n_variaveis\n", argv[0]);
    return 1;
  }

  // INICIAR VALOR DA SEMENTE
  SRAND(20252);
    
  // CHAMAR FUNÇÕES DE INTEGRAÇÃO E EXIBIR RESULTADOS
  
  return 0;
}


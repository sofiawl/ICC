#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"

#define DIFF 0.0

#define NRAND    ((real_t) random() / RAND_MAX)  // drand48() 
#define SRAND(a) srandom(a) // srand48(a)
#define RAND_AB(a, b) ((a) + NRAND * ((b) - (a))) // random entre [a,b]

/*
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
*/

// Integral Monte Carlo da função Styblinski-Tang de 2 variáveis
real_t styblinskiTang(real_t a, real_t b, int namostras, int dim)
{
  if (dim != 2 || dim != 4 || dim != 8){
    printf("Dimensões só podem ser: 2, 4 ou 8\n");
    return 0.0;
  }
    
  real_t resultado = 0.0;
  real_t x, qx, soma = 0.0;
  
  printf("Metodo de Monte Carlo (x, y).\n");
  printf("a = (%f), b = (%f), n = (%d), variaveis = %d\n", a, b, namostras, dim);
  
  rtime_t t_inicial = timestamp();

  real_t areaAB = 1.0;
  
  // 2, 4, ou 8 DIMENSÕES
  for (int i=0; i<dim; i++){
    areaAB *= (b-a);
  }

  int h = b-a;
  /* 
     Solução genérica para todas dimensões
     Se fossemos fazer uma solução específica poderíamos ter: x, y, qx, qy  
  */
  real_t imc = 0.0;
  for (int i=0; i<namostras; i++){
    // x = a + h
    // 'h' é uma porcentagem aleat do tamanho do intervalo [a, b]
    imc += soma;
    for (int j=0; j<dim; j++){
      x = a + NRAND * h;  
      qx = x*x;
      soma += qx*qx - 16*qx + 5*x; 
    }
    soma /= 2.0;
  }
  

  resultado = (imc / namostras) * areaAB;
  
  rtime_t t_final = timestamp();
  printf("Tempo decorrido: %f seg.\n", t_final - t_inicial);
  
  return resultado;

  /* OTIMIZANDO:  
    Unroll and jam das somas
    Soma somas no final
    resultado = (somaS / namostras) * areaAB
  */
}

// 2 DIMENSÕES
real_t retangulos_xy(real_t a, real_t b, int npontos) {

  real_t k = ceil(sqrt(npontos));
  real_t h = (b-a)/ k;
  real_t resultado;
  real_t soma = 0;
  real_t x, y, qx, qy; 

  printf("Metodo dos Retangulos (x, y).\n");
  printf("a = (%f), b = (%f), n = (%d), h = (%lg)\n", a, b, npontos, h);
  
  rtime_t t_inicial = timestamp();
  
  for (int i=a; i<=b; i+=h){
    qx = x*x;
    for (int j=a; j<=b; j+=h){
      qy = y*y;
      soma += (qx*qx + qy*qy -16*(qx+qy) +5*(x+y))/2; 
    }
  }

  resultado = h*h*soma;

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


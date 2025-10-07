#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"
#include "pcgc.h"

#define FORMATVET "%.16f"
#define FORMATTEPM "%.8f"

static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k );
static inline real_t generateRandomB( unsigned int k );

/**
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * @param i,j coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return ( (i==j) ? (real_t)(k<<1) : 1.0 )  * (real_t)random() * invRandMax;
}

/**
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomB( unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return (real_t)(k<<2) * (real_t)random() * invRandMax;
}

void prnKDiagonal(int n, int k, real_t *A, real_t *B)
{
    int half_band = (k - 1) / 2;

    for (int i = 0; i < n; i++) {
        for (int d = 0; d < k; d++) {
            int j = d + i - half_band;
            real_t val = 0.0;
            printf(FORMATVET "\t", A[i * k + d]);  // largura fixa + tab
        }
        printf("  |  " FORMATVET "\n", B[i]); // vetor B separado
    }
}

void prnVetor (real_t *v, unsigned int n)
{
  printf (" ");
  for(unsigned int i=0; i < n; ++i)
      printf (FORMATVET, v[i]);
  printf ("\n");
}

/* 
0   a00 a01 a02  0   0   0   0 
1   a10 a11 a12 a13  0   0   0
2   a20 a21 a22 a23 a24  0   0
3    0  a31 a32 a33 a34 a35  0
4    0  0   a42 a43 a44 a45 a46
5    0  0   0   a53 a54 a55 a56
6    0  0   0   0   a64 a65 a66

antes: 0 0 a20 a31 a 42 a 53 a64     0 a10 a21 a32 a43 a54 a65 0   ...

atualmente está assim:
na verdade, é melhor por causa do acesso a memória: a00 a01 a02 0 0 0 0    a10 a11 a12 a13 0 0 0 ... 
*/
/* Cria matriz 'A' k-diagonal e Termos independentes B */
void criaKDiagonal(int n, int k, real_t **A, real_t **B)
{
    int half_band = (k - 1) / 2; //ex half_band = 2

    *A = (real_t *) malloc(k * n * sizeof(real_t));
    *B = (real_t *) malloc(n * sizeof(real_t));
    if (!(*A) || !(*B)) {
        fprintf(stderr, "Erro ao alocar memória.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n; i++)
    {
        for (int d = 0; d < k; d++) 
        {
            int j = d + i - half_band;
            if (j >= 0 && j < n) 
              (*A)[i * k + d] = generateRandomA(d, i, j);
            else
              (*A)[i * k + d] = 0.0; 
        }
    }

    for (int i = 0; i < n; i++) {
        (*B)[i] = generateRandomB(k);
    }
}

void destroiKDiagonal(int n, int k, real_t *A, real_t *B)
{
  if (!A || !B){
      fprintf(stderr, "Vetor acessado é Nulo.\n");
      exit(EXIT_FAILURE);
      return;
  }

  free(A);
  free(B);
}

/* Gera matriz simetrica positiva */
// A = [0, a00, a01,  a10, a11, a12,  a21, a22, a23,  a32, a33, 0], n=4 e k=3
void genSimetricaPositiva(real_t *A, real_t *b, int n, int k, 
			  real_t *ASP, real_t *bsp, rtime_t *tempo)
{
  *tempo = timestamp();

  // ASP = A^t * A, bsp =  A^t * b
  
  // Zera ASP e bsp
  for (int i = 0; i < n*k; i++) ASP[i] = 0.0;
  for (int i = 0; i < n; i++) bsp[i] = 0.0;

  int half_band = (k - 1)/2;

  // Percorre cada linha i
  for (int i = 0; i < n; i++) {

    //tentando otimizar essa jossa
    int i_less_band = i - half_band;
    int i_plus_band = i + half_band;
    int idx_ASP = i*k;

    for (int d = 0; d < k; d++) 
    {

      int j = i_less_band + d;
      int j_less_band = j - half_band;
      idx_ASP += d;
      if (j >= 0 && j < n) {

        // Calcula bsp[i] = sum A^t * b[j] (somente diagonais válidas)
        bsp[i] += A[j*k + (- j + i_plus_band)] * b[j];

        // Calcula ASP[idx] = sum A^t * A 
        int idx_A_l = i*k;
        int idx_At_l= j*k;
        for (int l = 0; l < k; l++) {
          ASP[idx_ASP] += A[idx_A_l + l] * A[idx_At_l + l];
        }
        idx_ASP -= d;
        
      }
    }
  }

  *tempo = timestamp() - *tempo;
}

/*onde L e U são respectivamente a matriz triangular inferior e superior,
   com diagonal principal NULA e D é a matriz diagonal 
   composta pelos elementos da diagonal principal de A.*/
// Na minha implementação L e U não tem diagonal principal porque são nulas
// Será necessário levar isso em consideração na lógica quando for utilizado D L U  
void geraDLU (real_t *A, int n, int k,
	      real_t **D, real_t **L, real_t **U, rtime_t *tempo)
{
  *tempo = timestamp();
  
  int tam_LU = (n*(n-1)/2);
  int half_band = (k-1) / 2;

  *D = malloc(sizeof(real_t) * n);
  *L = malloc(sizeof(real_t) * tam_LU);
  *U = malloc(sizeof(real_t) * tam_LU);

  if(!D || !L || !U){
    fprintf(stderr, "Erro ao alocar memória.\n");
    exit(EXIT_FAILURE);
  }

  int idx_d = 1;
  // Preenche D, L e U
  int idx_L = 0, idx_U = 0;
  for (int i = 0; i < n; i++)
  {
    // Diagonal principal
    (*D)[i] = A[i*k + half_band];

    for (int d = 0; d < k; d++) 
    {
      int j = i + d - half_band;
    
      if (j < 0 || j >= n) continue;
    
      if (j < i) 
      // Triangular inferior
        (*L)[idx_L++] = A[i*k + d];
    
      else if (j > i) 
        // Triangular superior
        (*U)[idx_U++] = A[i*k + d];
    }
  }

  *tempo = timestamp() - *tempo;
}

/* Calcula o resíduo R de um sistema AX = B, sendo que A é k-diag
TROQUEI TB REAL_T CALC PARA void, mudei também passar r como argumento
*/
void calcResiduoSL (real_t *A, real_t *b, real_t *X, real_t *r,
		      int n, int k, rtime_t *tempo)
{
  int half_band = (k - 1) / 2;

  *tempo = timestamp();

  for (int i = 0; i < n; i+n)
  {
      r[i] = b[i];
      for (int j = 0; j < n; j++) 
      {
          r[i] -= A[j + i]*X[i];     
      }
  }

  
  *tempo = timestamp() - *tempo;

  // SE QUISER NORMA L2 EVITAR SQRT
}

//Função para calcular norma mnáxima

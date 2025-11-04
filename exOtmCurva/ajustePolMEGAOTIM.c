#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <fenv.h>
#include <math.h>
#include <stdint.h>

#include <likwid.h>

#include "utils.h"

#define FATORUNROLL 4

/////////////////////////////////////////////////////////////////////////////////////
//   AJUSTE DE CURVAS
/////////////////////////////////////////////////////////////////////////////////////
inline double powOTIM(double num, int exp){
  for (int i = 0; i < exp; i++)
    num *= num;
  return num;
}

void montaSL(double ** restrict A, double *  restrict b, int n, long long int p, double * restrict x, double * restrict y) {
  int maxi = n - (n%FATORUNROLL);
  
  // Unroll and Jam
  for (int i = 0; i < maxi; i+=FATORUNROLL) {

    b[i] = b[i+1] = b[i+2] = b[i+3] = 0.0;
    
    for (int j = 0; j < n; ++j) {

      // Cálculo de A
      A[i][j] = A[i+1][j] = A[i+2][j] = A[i+3][j] = 0.0;
      int sumij = i+j;
      for (long long int k = 0; k < p; ++k) {
	      A[i][j] += powOTIM(x[k], sumij);
      }
      sumij++;
      for (long long int k = 0; k < p; ++k) {
	      A[i+1][j] += powOTIM(x[k], sumij);
      }
      sumij++;
      for (long long int k = 0; k < p; ++k) {
	      A[i+2][j] += powOTIM(x[k], sumij);
      }
      sumij++;
      for (long long int k = 0; k < p; ++k) {
	      A[i+3][j] += powOTIM(x[k], sumij);
      }
    
      // Cálculo de b
      b[i] += powOTIM(x[j],i) * y[j];
      b[i+1] += powOTIM(x[j],i+1) * y[j];
      b[i+2] += powOTIM(x[j],i+2) * y[j];
      b[i+3] += powOTIM(x[j],i+3) * y[j]; 
    }
  }

  // Resíduo do laço
  for (int i = n-(n%FATORUNROLL); i < n; ++i){
    b[i] = 0.0;

    // Cálculo do resíduo de A
    for (int j = 0; j < n; ++j){
      A[i][j] = 0.0;
      int sumij = i+j;
      for (long long int k = 0; k < p; ++k) {
	      A[i][j] += powOTIM(x[k], sumij);
      }
    
      // Cálculo do resíduo de B
      b[i] += powOTIM(x[j],i) * y[j];
    }
  }
}

void eliminacaoGauss(double ** restrict A, double * restrict b, int n) {
  for (int i = 0; i < n; ++i) {
    int iMax = i;
    for (int k = i+1; k < n; ++k)
      if (A[k][i] > A[iMax][i])
	      iMax = k;
    if (iMax != i) {
      double *tmp, aux;
      tmp = A[i];
      A[i] = A[iMax];
      A[iMax] = tmp;

      aux = b[i];
      b[i] = b[iMax];
      b[iMax] = aux;
    }

    for (int k = i+1; k < n; ++k) {
      double m = A[k][i] / A[i][i];
      A[k][i]  = 0.0;

      for (int j = i+1; j < n; ++j)
	      A[k][j] -= A[i][j]*m;
      
      b[k] -= b[i]*m;
    }
  }
}

void retrossubs(double ** restrict A, double * restrict b, double * restrict x, int n) {
  for (int i = n-1; i >= 0; --i) {
    x[i] = b[i];
    for (int j = i+1; j < n; ++j)
      x[i] -= A[i][j]*x[j];
    x[i] /= A[i][i];
  }
}

double Pol(double x, int G, double * restrict alpha) {
  double Px = alpha[0];
  for (int i = 1; i <= G; ++i)
    Px += alpha[i]*powOTIM(x,i);
  
  return Px;
}

int main() {

  int G, g; // G -> grau do polinomio
  long long int P, p; // P -> no. de pontos
  string_t marker;

  scanf("%d %lld", &G, &P);
  
  g = G+1; // tamanho do SL (G + 1)
  if (g < FATORUNROLL){
    printf("O grau do polinômio deve ser maior ou igual a %d\n", FATORUNROLL);
    return 1;
  }
  p = P;   // quantidade de pontos

  double *x = (double *) malloc(sizeof(double)*p);
  double *y = (double *) malloc(sizeof(double)*p);

  // ler numeros
  for (long long int i = 0; i < p; ++i)
    scanf("%lf %lf", x+i, y+i);

  double **A = (double **) malloc(sizeof(double *)*g);
  for (int i = 0; i < g; ++i)
    A[i] = (double *) malloc(sizeof(double)*g);
  
  double *b = (double *) malloc(sizeof(double)*g);
  double *alpha = (double *) malloc(sizeof(double)*g); // coeficientes ajuste

  LIKWID_MARKER_INIT;
  
  // (A) Gera SL
  marker = markerName("SL",p);
  LIKWID_MARKER_START (marker);
  double tSL = timestamp();
  montaSL(A, b, g, p, x, y);
  tSL = timestamp() - tSL;
  LIKWID_MARKER_STOP(marker);
  free(marker);

  // (B) Resolve SL
  marker = markerName("EG",p);
  LIKWID_MARKER_START(marker);
  double tEG = timestamp();
  eliminacaoGauss(A, b, g); 
  retrossubs(A, b, alpha, g); 
  tEG = timestamp() - tEG;
  LIKWID_MARKER_STOP(marker);
  free(marker);

  LIKWID_MARKER_CLOSE;

  // Imprime coeficientes
  for (int i = 0; i < g; ++i)
    printf("%1.15e ", alpha[i]);
  puts("");

  // Imprime resíduos
  for (long long int i = 0; i < p; ++i)
    printf("%1.15e ", fabs(y[i] - Pol(x[i],G,alpha)) );
  puts("");

  // Imprime os tempos
  printf("%lld %1.10e %1.10e\n", P, tSL, tEG);

  return 0;
}

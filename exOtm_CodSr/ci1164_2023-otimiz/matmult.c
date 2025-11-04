#include <stdio.h>
#include <stdlib.h>    /* exit, malloc, calloc, etc. */
#include <string.h>
#include <getopt.h>    /* getopt */
#include <time.h>
#include <likwid.h>

#include "matriz.h"

/**
 * Exibe mensagem de erro indicando forma de uso do programa e termina
 * o programa.
 */

static void usage(char *progname)
{
  fprintf(stderr, "Forma de uso: %s [ <ordem> ] \n", progname);
  exit(1);
}



/**
 * Programa principal
 * Forma de uso: matmult [ -n <ordem> ]
 * -n <ordem>: ordem da matriz quadrada e dos vetores
 *
 */

int main (int argc, char *argv[]) 
{
  LIKWID_MARKER_INIT;
  int n=DEF_SIZE;
  
  MatRow mRow_1, mRow_2, resMat, resMat_otm;
  Vetor vet, res, res_otm;
  
  /* =============== TRATAMENTO DE LINHA DE COMANDO =============== */

  if (argc < 2)
    usage(argv[0]);

  n = atoi(argv[1]);
  
  /* ================ FIM DO TRATAMENTO DE LINHA DE COMANDO ========= */
 
  srandom(20232);
      
  res = geraVetor (n, 0); // (real_t *) malloc (n*sizeof(real_t));
  resMat = geraMatRow(n, n, 1);
  res_otm = geraVetor (n, 0); // (real_t *) malloc (n*sizeof(real_t));
  resMat_otm = geraMatRow(n, n, 1);

  mRow_1 = geraMatRow (n, n, 0);
  mRow_2 = geraMatRow (n, n, 0);

  vet = geraVetor (n, 0);

  if (!res || !resMat || !res_otm || !resMat_otm || !mRow_1 || !mRow_2 || !vet) {
    fprintf(stderr, "Falha em alocação de memória !!\n");
    liberaVetor ((void*) mRow_1);
    liberaVetor ((void*) mRow_2);
    liberaVetor ((void*) resMat);
    liberaVetor ((void*) vet);
    liberaVetor ((void*) res);
    exit(2);
  }
    
#ifdef _DEBUG_
    prnMat (mRow_1, n, n);
    prnMat (mRow_2, n, n);
    prnVetor (vet, n);
    printf ("=================================\n\n");
#endif /* _DEBUG_ */

  LIKWID_MARKER_START("MatVet");
  multMatVet (mRow_1, vet, n, n, res);
  LIKWID_MARKER_STOP("MatVet");

  LIKWID_MARKER_START("MatMat");
  multMatMat (mRow_1, mRow_2, n, resMat);
  LIKWID_MARKER_STOP("MatMat");

  LIKWID_MARKER_START("MatVet_Otm");
  multMatVet_otm (mRow_1, vet, n, n, res_otm);
  LIKWID_MARKER_STOP("MatVet_Otm");

  LIKWID_MARKER_START("MatMat_Otm");
  multMatMat_otm (mRow_1, mRow_2, n, resMat_otm);
  LIKWID_MARKER_STOP("MatMat_Otm");
    
#ifdef _DEBUG_
    prnVetor (res, n);
    prnMat (resMat, n, n);
    prnVetor (res_otm, n);
    prnMat (resMat_otm, n, n);
#endif /* _DEBUG_ */

  liberaVetor ((void*) mRow_1);
  liberaVetor ((void*) mRow_2);
  liberaVetor ((void*) resMat);
  liberaVetor ((void*) vet);
  liberaVetor ((void*) res);
  liberaVetor ((void*) resMat_otm);
  liberaVetor ((void*) res_otm);
  
  LIKWID_MARKER_CLOSE;

  return 0;
}


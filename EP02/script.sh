#!/bin/bash
# -------------------------------------------------------------
#                   SOFIA WAMSER LIMA
#                    GRR: 20240495
#                     20/09/2025
# -------------------------------------------------------------

likwid-perfctr -C 0 -g FLOPS_DP -m ./resolveEDO < "$1" > out.dat 2>&1


# Seria para filtrar o relatório do LIKWID para pegar apenas a linha com o contador.
# Mas não está acontecendo isso na impressão. Por isso a saída não ficou exatamente como no exemplo. 
grep FP_ARITH_INST_RETIRED_SCALAR_DOUBLE out.dat | awk -F '|' '{print "FP_ARITH_INST_RETIRED_SCALAR_DOUBLE," $4}' | sed 's/ //g'

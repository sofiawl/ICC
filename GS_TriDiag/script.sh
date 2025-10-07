#!/bin/bash
# -------------------------------------------------------------
#                   SOFIA WAMSER LIMA
#                    GRR: 20240495
#                     20/09/2025
# -------------------------------------------------------------
# Exemplo para rodar o script: ./script.sh arquivo_de_entrada.dat

./resolveEDO < "$1"

echo ""

likwid-perfctr -C 3 -g FLOPS_DP -m ./resolveEDO < "$1" > out.dat 2>&1


# Eu usei grep 'FLOPS'... porque o contador FP_ARITH_INST_RETIRED_SCALAR_DOUBLE
# existe apenas em CPUs Intel. Quando tentei rodar usando apenas o comando da Intel,
# não obtive a saída esperada como no exemplo do professor. Assim consigo capturar os
# eventos de FLOPS em qualquer CPU, inclusive AMD, mantendo o script portável.

# Se quiser ver a saída completa detalhada, pode consultar o arquivo out.dat.
grep 'FLOPS' out.dat | grep '|' | awk -F '|' 'NF>=4 {gsub(/ /,""); print $2 ", " $3 ", " $4}'

#grep FP_ARITH_INST_RETIRED_SCALAR_DOUBLE out.dat | awk -F '|' '{print "FP_ARITH_INST_RETIRED_SCALAR_DOUBLE," $4}' | sed 's/ //g'

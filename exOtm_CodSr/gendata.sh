#!/bin/bash

PROG=matmult
CPU=3

DATA_DIR="Dados/"

mkdir -p ${DATA_DIR}

echo "performance" > /sys/devices/system/cpu/cpufreq/policy${CPU}/scaling_governor

make purge matmult

METRICA="FLOPS_DP L3CACHE ENERGY"
TEMPOS="${DATA_DIR}/Tempos.csv"
TAMANHOS="64 100 128 1024 2000 2048 3000 4096 6000 7000 10000 50000 60000 70000 100000"

for m in ${METRICA}
do
    LIKWID_CSV="${DATA_DIR}/${m}.csv"
    rm -f ${TEMPOS}

    for n in $TAMANHOS
    do
	LIKWID_OUT="${DATA_DIR}/${m}_${n}.txt"
	
	echo "--->>  $m: ./${PROG} $n" >/dev/tty
	# Assume que programa 'matmult' imprime em stdout (via printf) os valores de tempo
	# medidos para cada função com 5 colunas: N, t_matVet, t_matVet_otim, t_matMat, t_matmat_otim
	likwid-perfctr -O -C ${CPU} -g ${m} -o ${LIKWID_OUT} -m ./${PROG} ${n} >>${TEMPOS}
    done

    # Colocar aqui comando(s) que, a partir dos arquivos '.txt', gera (para cada métrica) um arquivo
    # CSV com 5 colunas:
    # N, metrica_matvet,metrica_matvet_otim, metrica_matmat, metrica_matmat_otim
    # 
done

echo "powersave" > /sys/devices/system/cpu/cpufreq/policy${CPU}/scaling_governor 


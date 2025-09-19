#!/bin/bash

likwid_perfctr -C 0 -g FLOPS_DP -m ./resolveEDO < "$1" > out.dat 2>&1 

grep FP_ARITH_INST_RETIRED_SCALAR_DOUBLE | awk -F '|' '{print "FP_ARITH_INST_RETIRED_SCALAR_DOUBLE," $4}' | sed 's/ //g'

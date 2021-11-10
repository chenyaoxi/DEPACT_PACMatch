#!/bin/bash

for n in CCC_01 COC_01 NH2_01 OH_01 PO4_01 ncnc2_01 ncnc3_01
do
#	cp bf_${n}.txt bf_${n}.txt.bak
cp bf_${n}.txt.bak bf_${n}.txt
#	head -n 100 bf_${n}.txt.bak > bf_${n}.txt
done

#!/bin/bash

for n in bf_ncnc3_01.txt bf_ncnc2_01.txt bf_COC_01.txt bf_CCC_01.txt bf_NH2_01.txt bf_OH_01.txt bf_PO4_01.txt
do
	cp ${n}.bak ${n}
	sed -i '/ATP/d' ${n}
done

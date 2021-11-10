#!/bin/bash
cat all-sdf.sdf | while read line
do
pdbid=$(echo $line | awk -F '_' '{print $1}');
ligand=$(echo $line | awk -F '_' '{print $2}');
if [ "$ligand" == "ZN" ]; then
echo $pdbid >> ZN.txt
fi
done

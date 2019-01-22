#!/bin/bash

trap "echo Exited!; exit;" SIGINT SIGTERM

counter=0
threads=1

(
for f in `ls simulated_lethal_data`; do
	for L in 100; do
		if [[ -f fits_proper/`expr "$f" : '\(.*_\)'`L$L`expr "$f" : '.*\(_.*\)'` ]]; then
			((i=i%threads)); ((i++==0)) && wait
			{ echo "$file '$f' not found, fitting now."
			./fits_xyz_only.R simulated_lethal_data/$f $L || true; } &
		else
			echo "file '$f' already fitted, ignoring"
		fi
	done
done
)

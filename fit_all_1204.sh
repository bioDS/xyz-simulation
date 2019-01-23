#!/bin/bash

trap "echo Exited!; exit;" SIGINT SIGTERM

counter=0
threads=12

(
for f in `ls simulated_data`; do
	for L in 100; do
		if [[ -f fits_proper/`expr "$f" : '\(.*_\)'`L$L`expr "$f" : '.*\(_.*\)'` ]]; then
			echo "file '$f' already fitted, ignoring"
		else
			((i=i%threads)); ((i++==0)) && wait
			{ echo "$file '$f' not found, fitting now."
			./fits_xyz_only.R simulated_data/$f $L write || true; } &
		fi
	done
done
)

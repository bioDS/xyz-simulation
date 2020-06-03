#!/bin/bash

trap "echo Exited!; exit;" SIGINT SIGTERM

if [[ $1 ]]; then
	threads=$1
else
	threads=2
fi

counter=0

(
for f in `ls simulated_data | grep SNR5 | grep lethals0`; do
	for L in 10 100 1000; do
		p=`expr "$f" : '.*p\([0-9]*\)_'`
		l=$(bc<<<"scale=1; x=sqrt($p) + 0.5; scale=0; x/1")
		#l_file="fits_proper/`expr "$f" : '\(.*_\)'`L$l`expr "$f" : '.*\(_.*\)'`"
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

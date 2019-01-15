#!/bin/bash

counter=0
threads=2

if [[ $1 ]]; then
        multiplier=$1
else
        multiplier=1
fi

n=$((multiplier * 1000))
p=$((multiplier * 100))

(
for SNR_pre in 2 5 10; do
	SNR=$((multiplier * SNR_pre))
	for bij_pre in 5 20 50 100; do
		bij=$((multiplier * bij_pre))
		for bi_pre in 0 20 50 100; do
			((i=i%threads)); ((i++==0)) && wait
			{ bi=$((multiplier * bi_pre))
			./data_generator.R $n $p $SNR $bi $bij 0 || true; } &
		done
	done
done
)

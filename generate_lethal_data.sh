#!/bin/bash

trap "echo Exited!; exit;" SIGINT SIGTERM

counter=0
threads=24

if [[ $1 ]]; then
        multiplier=$1
else
        multiplier=1
fi

if [[ $2 ]]; then
	repetitions=$2
else
	repetitions=1
fi

n=$((multiplier * 1000))
p=$((multiplier * 100))

(
for repetition in `seq 1 $repetitions`; do
	for SNR in 2 5 10; do
		for bij_pre in 5 20 50 100; do
			bij=$((multiplier * bij_pre))
			for bi_pre in 0 20 50 100; do
				((i=i%threads)); ((i++==0)) && wait
				{ bi=$((multiplier * bi_pre))
				./data_generator.R $n $p $SNR 10 $bi $bij 0 || true; } &
			done
		done
	done
done
)

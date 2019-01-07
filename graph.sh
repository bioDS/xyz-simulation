#!/bin/bash

if [[ $1 ]]; then
        alpha=$1
else
        alpha=0.9
fi


for SNR in 2 5 10; do
	for bij in 5 20 50 100; do
		for bi in 0 20 50 100; do
			echo "1000 100 $SNR $bi $bij 0 $alpha"
			./fits_xyz.R 1000 100 $SNR $bi $bij 0 $alpha
		done
	done
done

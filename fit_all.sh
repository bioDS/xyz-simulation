#!/usr/bin/bash

trap "echo Exited!; exit;" SIGINT SIGTERM

for i in `ls simulated_data`; do
	for L in 10 100 1000; do
		./fits_xyz_only.R $i $L
	done
done

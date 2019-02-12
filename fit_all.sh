#!/usr/bin/env bash

trap "echo Exited!; exit;" SIGINT SIGTERM

for i in `ls simulated_data`; do
	for L in -1; do
		./fits_xyz_only.R simulated_data/$i $L write
	done
done

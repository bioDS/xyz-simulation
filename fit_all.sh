#!/usr/bin/bash

trap "echo Exited!; exit;" SIGINT SIGTERM

for i in `ls simulated_data`; do
	./fits_xyz_only.R $i
done

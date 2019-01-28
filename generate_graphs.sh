#!/bin/bash

trap "echo Exited!; exit;" SIGINT SIGTERM

if [ $1 == "large" ]; then
    echo "generating large graphs"
    prf_size="l"
    mult=10
elif [ $1 == "small" ]; then 
    echo "generating small graphs"
    prf_size="s"
    mult=1
elif [ $1 == "store" ]; then
    d=`date +%s`
    for i in PrecRecF1 FXstrength FXdiff NumObservations l_diff; do
        if [[ -d $i ]]; then
            cd $i
            if ls *.pdf 1> /dev/null 2>&1; then
                mkdir pdfs_$d; mv *.pdf pdfs_$d/
            fi
            cd ..
        fi
    done
    exit 0
else 
    echo "argument not understood, acceptable arguments are 'large', 'small', 'store'"
    exit 0
fi

if [ $2 == "n" ]; then
    read_fits="n"
else
    read_fits="y"
fi

echo "PrecRecF1"
cd PrecRecF1; ./PrecRecF1.R $read_fits $prf_size
if [ $1 == "large" ]; then
    echo "PrecRecF1_lethals"
    ./PrecRecF1_lethals.R $read_fits # actually won't work unless 'y' is used.
fi
cd ..
echo "FXstrength"
cd FXstrength; ./FXstrength.R $read_fits $mult
cd ..
echo "FXdiff"
cd FXdiff; ./FXdiff.R $mult
cd ..
echo "NumObservations"
cd NumObservations; ./NumObservations.R $read_fits $mult n
cd ..
if [ $1 == "large" ]; then
    echo "l_diff"
    cd l_diff; ./l_diff.R $prf_size $read_fits
fi

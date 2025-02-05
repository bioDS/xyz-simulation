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

use_xyz=1
use_glinternet=1

if [ $3 == "x" ]; then
    use_glinternet=0
elif [ $3 == "g" ]; then
    use_xyz=0
fi

if [ $use_xyz == 1 ]; then
    echo "regenerating xyz graphs:"
    echo "PrecRecF1"
    cd PrecRecF1; ./PrecRecF1.R $read_fits $prf_size x
    if [ $1 == "large" ]; then
        echo "PrecRecF1_lethals"
        ./PrecRecF1_lethals.R $read_fits # actually won't work unless 'y' is used.
    fi
    cd ..
    echo "FXstrength"
    cd FXstrength; ./FXstrength.R $read_fits $mult x
    cd ..
    echo "FXdiff"
    cd FXdiff; ./FXdiff.R $mult x
    cd ..
    echo "NumObservations"
    cd NumObservations; ./NumObservations.R $read_fits $mult n x
    cd ..
    if [ $1 == "large" ]; then
        echo "l_diff"
        cd l_diff; ./l_diff.R $prf_size $read_fits
    fi
fi

if [ $use_glinternet == 1 ]; then
    echo "regenerating glinternet graphs:"
    echo "PrecRecF1"
    cd PrecRecF1; ./PrecRecF1.R $read_fits $prf_size g
    cd ..
    echo "FXstrength"
    cd FXstrength; ./FXstrength.R $read_fits $mult g
    cd ..
    echo "FXdiff"
    cd FXdiff; ./FXdiff.R $mult g
    cd ..
    echo "NumObservations"
    cd NumObservations; ./NumObservations.R $read_fits $mult n g
    cd ..
fi

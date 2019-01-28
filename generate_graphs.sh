#!/bin/bash
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
            if ls *.pdfg 1> /dev/null 2>&1; then 
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


cd PrecRecF1; ./PrecRecF1.R $prf_size
./PrecRecF1_lethals.R
cd ..
cd FXstrength; ./FXstrength.R y $mult
cd ..
cd FXdiff; ./FXdiff.R $mult
cd ..
cd NumObservations; ./NumObservations.R y $mult n
cd ..
cd l_diff; ./l_diff.R $prf_size y

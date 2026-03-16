#!/bin/bash

# 
SOURCE_CODE="Terwel_JB_code.c"
EXECUTABLE="a.out"
MATRIX_SIZE=200
THREADS=16
RUNS=3

#
gcc "$SOURCE_CODE" -lm -o "$EXECUTABLE"
#
for ((i=1; i<=RUNS; i++))
do
    echo "---------------------------"
    ./$EXECUTABLE $MATRIX_SIZE $THREADS

    # check if program crashed
    if [ $? -ne 0 ]; then
        echo "Error: Program crashed on run #$i"
        exit 1
    fi
done

echo "---------------------------"
echo "All runs completed successfully"
#!/bin/bash

# 
SOURCE_CODE="lab01.c"
EXECUTABLE="lab01"
MATRIX_SIZE=25000
THREADS=8
RUNS=3

#
gcc -o "$EXECUTABLE" "$SOURCE_CODE" -lm
#
for ((i=1; i<=RUNS; i++))
do
    echo "---------------------------"
    ./$EXECUTABLE $MATRIX_SIZE

    # check if program crashed
    if [ $? -ne 0 ]; then
        echo "Error: Program crashed on run #$i"
        exit 1
    fi
done

echo "---------------------------"
echo "All runs completed successfully"
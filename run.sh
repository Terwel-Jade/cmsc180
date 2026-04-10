#!/bin/bash
SOURCE_CODE="lab01.c"
EXECUTABLE="lab01"
MATRIX_SIZE=25000
RUNS=3

echo "Compiling $SOURCE_CODE..."
gcc -O2 -o "$EXECUTABLE" "$SOURCE_CODE" -lm

if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

echo "Compilation successful."
echo "Running program $RUNS times..."

total_time=0

for ((i=1; i<=RUNS; i++))
do
    echo "---------------------------"
    echo "Run #$i"

    start=$(date +%s%N)
    ./$EXECUTABLE $MATRIX_SIZE
    status=$?
    end=$(date +%s%N)

    runtime=$(awk "BEGIN {printf \"%.6f\", ($end - $start) / 1000000000}")
    echo "Execution time: $runtime seconds"

    if [ $status -ne 0 ]; then
        echo "Error: Program crashed on run #$i"
        exit 1
    fi

    total_time=$(awk "BEGIN {printf \"%.6f\", $total_time + $runtime}")
done

echo "---------------------------"
average=$(awk "BEGIN {printf \"%.6f\", $total_time / $RUNS}")
echo "All runs completed successfully"
echo "Average execution time over $RUNS runs: $average seconds"
#!/bin/bash

DIRECTORY="Desktop/temp"
BINARY="lab04"
SIZE=4000
RUNS=1

for ((i=1; i<=RUNS; i++))
do 
    echo "Iteration $i --------------"

    ssh 10.0.4.238 "cd $DIRECTORY && ./$BINARY $SIZE 5002 1" &

    ssh 10.0.4.69 "cd $DIRECTORY && ./$BINARY $SIZE 5001 1" &
    
    ./$BINARY $SIZE 5000 0

    wait
done

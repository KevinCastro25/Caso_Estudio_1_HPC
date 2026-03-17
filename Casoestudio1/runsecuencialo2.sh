#!/bin/bash

for ((j=1; j<=10; j++)); do
    for i in 800 1000 2000 4000 6000 8000; do
        ./matricesecuencialo2 "$i" >> timesecuencialo2.doc
    done
done

echo Finalizado.
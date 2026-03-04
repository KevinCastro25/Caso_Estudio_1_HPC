#!/bin/bash

for ((j=1; j<=10; j++)); do
    for i in 1000 2000 3000; do
        ./matricesecuencial "$i" >> timesecuencial.doc
    done
done

echo Finalizado.
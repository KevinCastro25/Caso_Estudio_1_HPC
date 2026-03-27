#!/bin/bash

for ((j=1; j<=10; j++)); do
    for k in 7 10 13 16; do
        ./jacobi_Red_black "$k" 1e-6 | grep "TIEMPO\|VERIF" >> resultados_jacobired_black.doc
    done
done

echo "Finalizado."
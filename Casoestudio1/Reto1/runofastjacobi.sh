#!/bin/bash

for ((j=1; j<=10; j++)); do
    for k in 7 10 13 16; do
        ./jacobiofast "$k" 1e-6 | grep "TIEMPO\|VERIF" >> resultados_jacobiofast.doc
    done
done

echo "Finalizado."
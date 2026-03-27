#!/bin/bash
# 2 Procesos
for ((j=1; j<=10; j++)); do
    for k in 7 10 13 16; do
        ./jacobiprocesos "$k" 1e-6 2 | grep "TIEMPO\|VERIF" >> procesos_jacobi2.doc
    done
done

echo "Finalizado."
#4 Procesos
for ((j=1; j<=10; j++)); do
    for k in 7 10 13 16; do
        ./jacobiprocesos "$k" 1e-6 4 | grep "TIEMPO\|VERIF" >> procesos_jacobi4.doc
    done
done

echo "Finalizado."
#8 Procesos
for ((j=1; j<=10; j++)); do
    for k in 7 10 13 16; do
        ./jacobiprocesos "$k" 1e-6 8 | grep "TIEMPO\|VERIF" >> procesos_jacobi8.doc
    done
done

echo "Finalizado."
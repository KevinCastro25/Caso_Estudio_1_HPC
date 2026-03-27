#!/bin/bash
# 2 Hilos
for ((j=1; j<=10; j++)); do
    for k in 7 10 13 16; do
        ./jacobihilos "$k" 1e-6 2 | grep "TIEMPO\|VERIF" >> hilos_jacobi2.doc
    done
done

echo "Finalizado."
#4 Hilos
for ((j=1; j<=10; j++)); do
    for k in 7 10 13 16; do
        ./jacobihilos "$k" 1e-6 4 | grep "TIEMPO\|VERIF" >> hilos_jacobi4.doc
    done
done

echo "Finalizado."
#8 Hilos
for ((j=1; j<=10; j++)); do
    for k in 7 10 13 16; do
        ./jacobihilos "$k" 1e-6 8 | grep "TIEMPO\|VERIF" >> hilos_jacobi8.doc
    done
done

echo "Finalizado."
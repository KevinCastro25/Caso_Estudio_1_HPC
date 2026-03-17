#!/bin/bash

# Pruebas con 12 hilos
echo "Ejecutando pruebas con 12 hilos..."
for ((j=1; j<=10; j++)); do
    for i in 800 1000 2000 4000 6000 8000; do
        ./matriceshilos "$i" 12 >> timehilos12.doc
    done
done

echo "Finalizado."

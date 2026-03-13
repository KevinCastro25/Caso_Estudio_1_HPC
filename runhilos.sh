#!/bin/bash

# Pruebas con 8 hilos
echo "Ejecutando pruebas con 8 hilos..."
for ((j=1; j<=10; j++)); do
    for i in 800 1000 2000 4000 6000 8000; do
        ./matriceshilos "$i" 8 >> timehilos8.doc
    done
done

# Pruebas con 16 hilos
echo "Ejecutando pruebas con 16 hilos..."
for ((j=1; j<=10; j++)); do
    for i in 800 1000 2000 4000 6000 8000; do
        ./matriceshilos "$i" 16 >> timehilos16.doc
    done
done

echo "Finalizado."

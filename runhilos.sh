#!/bin/bash

# Pruebas con 2 hilos
echo "Ejecutando pruebas con 2 hilos..."
for ((j=1; j<=10; j++)); do
    for i in 800 1000 2000 4000 6000 8000; do
        ./matriceshilos "$i" 2 >> timehilos2.doc
    done
done

# Pruebas con 4 hilos
echo "Ejecutando pruebas con 4 hilos..."
for ((j=1; j<=10; j++)); do
    for i in 800 1000 2000 4000 6000 8000; do
        ./matriceshilos "$i" 4 >> timehilos4.doc
    done
done

# Pruebas con 8 hilos
echo "Ejecutando pruebas con 8 hilos..."
for ((j=1; j<=10; j++)); do
    for i in 800 1000 2000 4000 6000 8000; do
        ./matriceshilos "$i" 8 >> timehilos8.doc
    done
done

echo "Finalizado."

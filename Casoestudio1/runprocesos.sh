#!/bin/bash

# Pruebas con 2 procesos
echo "Ejecutando pruebas con 2 procesos..."
for ((j=1; j<=10; j++)); do
    for i in 800 1000 2000 4000 6000 8000; do
        ./matricesprocesos "$i" 2 >> timeprocesos2.doc
    done
done

# Pruebas con 4 procesos
echo "Ejecutando pruebas con 4 procesos..."
for ((j=1; j<=10; j++)); do
    for i in 800 1000 2000 4000 6000 8000; do
        ./matricesprocesos "$i" 4 >> timeprocesos4.doc
    done
done

# Pruebas con 8 procesos
echo "Ejecutando pruebas con 8 procesos..."
for ((j=1; j<=10; j++)); do
    for i in 800 1000 2000 4000 6000 8000; do
        ./matricesprocesos "$i" 8 >> timeprocesos8.doc
    done
done

# Pruebas con 12 procesos
echo "Ejecutando pruebas con 12 procesos..."
for ((j=1; j<=10; j++)); do
    for i in 800 1000 2000 4000 6000 8000; do
        ./matricesprocesos "$i" 12 >> timeprocesos12.doc
    done
done

echo "Finalizado."
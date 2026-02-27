#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <windows.h>   // for GetProcessTimes
#include <stdint.h>   // for uint64_t

// Función para crear una matriz dinámica nxn
int** crearMatriz(int n) {
    int** matriz = (int**) malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++) {
        matriz[i] = (int*) malloc(n * sizeof(int));
    }
    return matriz;
}

// Llenar matriz con números aleatorios
void llenarMatriz(int** matriz, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matriz[i][j] = rand() % 101; 
        }
    }
}

// Función para liberar memoria
void liberarMatriz(int** matriz, int n) {
    for (int i = 0; i < n; i++) {
        free(matriz[i]);
    }
    free(matriz);
}

// Multiplicación de matrices
void multiplicarMatrices(int** A, int** B, int** C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = 0;
            for (int k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        printf("Uso: %s <tamano de matriz>\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    if (n <= 0) {
        printf("El tamano de la matriz debe ser un numero positivo.\n");
        return 1;
    }

    // Inicializar semilla para números aleatorios
    srand(time(NULL));

    // Crear matrices
    int** A = crearMatriz(n);
    int** B = crearMatriz(n);
    int** C = crearMatriz(n);

    // Llenar matrices A y B
    llenarMatriz(A, n);
    llenarMatriz(B, n);

    // Medir tiempo de CPU mediante GetProcessTimes (usuario + kernel)
    FILETIME ftCreation, ftExit, ftKernel, ftUser;
    GetProcessTimes(GetCurrentProcess(), &ftCreation, &ftExit, &ftKernel, &ftUser);
    uint64_t start = ((uint64_t)ftKernel.dwHighDateTime << 32 | ftKernel.dwLowDateTime) +
                     ((uint64_t)ftUser.dwHighDateTime << 32 | ftUser.dwLowDateTime);

    multiplicarMatrices(A, B, C, n);

    GetProcessTimes(GetCurrentProcess(), &ftCreation, &ftExit, &ftKernel, &ftUser);
    uint64_t end = ((uint64_t)ftKernel.dwHighDateTime << 32 | ftKernel.dwLowDateTime) +
                   ((uint64_t)ftUser.dwHighDateTime << 32 | ftUser.dwLowDateTime);

    double cpu_time = (end - start) * 1e-7; // FILETIME units are 100 ns

    // Liberar memoria
    liberarMatriz(A, n);
    liberarMatriz(B, n);
    liberarMatriz(C, n);

    // Salida única: tiempo de CPU en segundos
    printf("%f\n", cpu_time);

    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

    // Medir tiempo de wall clock
    struct timespec start_wall, end_wall;
    clock_gettime(CLOCK_MONOTONIC, &start_wall);

    multiplicarMatrices(A, B, C, n);

    clock_gettime(CLOCK_MONOTONIC, &end_wall);
    // Calcular tiempo de wall clock en segundos
    double wall_time = (end_wall.tv_sec - start_wall.tv_sec) + 
                      (end_wall.tv_nsec - start_wall.tv_nsec) / 1e9;

    // Liberar memoria
    liberarMatriz(A, n);
    liberarMatriz(B, n);
    liberarMatriz(C, n);

    // Salida: tiempo de wall clock en segundos
    printf("%.6f\n", wall_time);

    return 0;
}

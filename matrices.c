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

// Función para liberar memoria
void liberarMatriz(int** matriz, int n) {
    for (int i = 0; i < n; i++) {
        free(matriz[i]);
    }
    free(matriz);
}

// Llenar matriz con números aleatorios
void llenarMatriz(int** matriz, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matriz[i][j] = rand() % 10; // Números entre 0 y 9
        }
    }
}

// Imprimir matriz
void imprimirMatriz(int** matriz, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%4d", matriz[i][j]);
        }
        printf("\n");
    }
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

int main() {
    int n;

    printf("Ingrese el tamano de las matrices cuadradas: ");
    scanf("%d", &n);

    // Inicializar semilla para números aleatorios
    srand(time(NULL));

    // Crear matrices
    int** A = crearMatriz(n);
    int** B = crearMatriz(n);
    int** C = crearMatriz(n);

    // Llenar matrices A y B
    llenarMatriz(A, n);
    llenarMatriz(B, n);

    // Medir tiempo de CPU solo para la multiplicación
    clock_t start = clock();
    multiplicarMatrices(A, B, C, n);
    clock_t end = clock();
    double cpu_time = (double)(end - start) / CLOCKS_PER_SEC;

    // Liberar memoria
    liberarMatriz(A, n);
    liberarMatriz(B, n);
    liberarMatriz(C, n);

    // Salida única: tiempo de CPU en segundos
    printf("%f\n", cpu_time);

    return 0;
}

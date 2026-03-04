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

// Multiplicación de matrices optimizada (fila x fila)
// Asume que B ya está traspuesta virtualmente (acceso por B[j][k])
void multiplicarMatrices(int** A, int** B, int** C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = 0;
            for (int k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[j][k];
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
    int** B = crearMatriz(n);  // B se asume ya traspuesta
    int** C = crearMatriz(n);

    // Llenar matrices A y B
    llenarMatriz(A, n);
    llenarMatriz(B, n);

    // Medir tiempo de CPU 
    struct timespec start_cpu, end_cpu;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_cpu);

    multiplicarMatrices(A, B, C, n);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_cpu);
    
    // Calcular tiempo de CPU en segundos
    double cpu_time = (end_cpu.tv_sec - start_cpu.tv_sec) + 
                      (end_cpu.tv_nsec - start_cpu.tv_nsec) / 1e9;

    // Liberar memoria
    liberarMatriz(A, n);
    liberarMatriz(B, n);
    liberarMatriz(C, n);

    // Salida: tiempo de CPU en segundos
    printf("%.6f\n", cpu_time);

    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>

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

// Estructura para pasar argumentos a los threads
typedef struct {
    int** A;
    int** B;
    int** C;
    int n;
    int start_row;
    int end_row;
} ThreadArgs;

// Función que ejecuta cada hilo para multiplicar filas
void* multiplicarFilas(void* arg) {
    ThreadArgs* args = (ThreadArgs*) arg;
    
    for (int i = args->start_row; i < args->end_row; i++) {
        for (int j = 0; j < args->n; j++) {
            args->C[i][j] = 0;
            for (int k = 0; k < args->n; k++) {
                args->C[i][j] += args->A[i][k] * args->B[k][j];
            }
        }
    }
    
    free(args);
    pthread_exit(NULL);
}

// Multiplicación de matrices con POSIX threads
void multiplicarMatricesHilos(int** A, int** B, int** C, int n, int num_threads) {
    pthread_t* threads = (pthread_t*) malloc(num_threads * sizeof(pthread_t));
    int rows_per_thread = n / num_threads;
    int remaining_rows = n % num_threads;
    
    int start_row = 0;
    for (int t = 0; t < num_threads; t++) {
        ThreadArgs* args = (ThreadArgs*) malloc(sizeof(ThreadArgs));
        args->A = A;
        args->B = B;
        args->C = C;
        args->n = n;
        args->start_row = start_row;
        // Distribuir filas extras entre los primeros threads
        args->end_row = start_row + rows_per_thread + (t < remaining_rows ? 1 : 0);
        
        pthread_create(&threads[t], NULL, multiplicarFilas, (void*) args);
        start_row = args->end_row;
    }
    
    // Esperar a que terminen todos los threads
    for (int t = 0; t < num_threads; t++) {
        pthread_join(threads[t], NULL);
    }
    
    free(threads);
}

int main(int argc, char* argv[]) {
    if (argc < 2 || argc > 3) {
        printf("Uso: %s <tamano de matriz> [numero de hilos]\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    if (n <= 0) {
        printf("El tamano de la matriz debe ser un numero positivo.\n");
        return 1;
    }

    // Determinar número de hilos
    int num_threads = (argc == 3) ? atoi(argv[2]) : sysconf(_SC_NPROCESSORS_ONLN);
    if (num_threads <= 0) {
        num_threads = 1;
    }
    if (num_threads > n) {
        num_threads = n;  // No tiene sentido más hilos que filas
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

    multiplicarMatricesHilos(A, B, C, n, num_threads);

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

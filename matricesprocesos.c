/**
 * =============================================================================
 * CONCEPTOS IMPORTANTES USADOS EN ESTE PROGRAMA
 * =============================================================================
 *
 * 1. FORK
 *    -----
 *    fork() crea una copia exacta del proceso actual (proceso hijo).
 *    El proceso original se llama "padre" y la copia se llama "hijo".
 *    Después del fork(), padre e hijo corren al mismo tiempo (en paralelo).
 *
 *    Ejemplo visual:
 *
 *        Proceso Padre
 *             |
 *           fork()
 *           /    \
 *    Hijo 0   Hijo 1   <- trabajan al mismo tiempo
 *       |         |
 *    filas 0-2  filas 3-4
 *       |         |
 *    exit(0)   exit(0)
 *           \    /
 *        waitpid()     <- padre espera a que terminen
 *             |
 *        [fin del programa]
 *
 * 2. FORK-JOIN
 *    ----------
 *    Es el patrón donde el padre:
 *      a) Hace FORK: lanza varios hijos a trabajar en paralelo
 *      b) Hace JOIN: espera con waitpid() a que todos terminen
 *    Solo cuando todos los hijos terminaron, el padre continúa.
 *
 * 3. MEMORIA COMPARTIDA (mmap)
 *    --------------------------
 *    Normalmente cuando se hace fork(), el hijo recibe una COPIA de la
 *    memoria del padre. Si el hijo modifica algo, el padre NO lo ve.
 *
 *    Para que los procesos compartan la misma memoria (y puedan escribir
 *    resultados que el padre pueda leer), usamos mmap() con MAP_SHARED.
 *
 *    Analogía: Es como tener una pizarra común donde todos escriben,
 *    en lugar de que cada uno tenga su propio cuaderno.
 *
 * 4. SEMÁFORO (exclusión mutua)
 *    ---------------------------
 *    Cuando varios procesos pueden escribir en la misma memoria al mismo
 *    tiempo, puede haber conflictos (condición de carrera). Un semáforo
 *    actúa como un "semáforo de tráfico": solo deja pasar a un proceso
 *    a la vez a la zona de escritura (región crítica).
 *
 *    sem_wait() → ROJO:  "espera tu turno, otro proceso está escribiendo"
 *    sem_post() → VERDE: "ya terminé, el siguiente puede pasar"
 *
 * 5. REGIÓN CRÍTICA
 *    ---------------
 *    Es el bloque de código donde se accede a datos compartidos.
 *    Solo un proceso puede ejecutarlo a la vez (protegido por el semáforo).
 *
 * =============================================================================
 * ESTRUCTURA DE LA MEMORIA COMPARTIDA
 * =============================================================================
 *
 *   [SharedMemory]
 *   +---------------------------------------------------------+
 *   | n          -> tamaño de la matriz                       |
 *   | num_procs  -> cuántos procesos hay                      |
 *   | mutex      -> semáforo para exclusión mutua             |
 *   | log_count  -> contador de escrituras (dato compartido)  |
 *   | data[]  ->  bloque contiguo con las 3 matrices:         |
 *   |   [ A(0,0) A(0,1)...A(n-1,n-1) ]  <- matriz A aplanada |
 *   |   [ B(0,0) B(0,1)...B(n-1,n-1) ]  <- matriz B aplanada |
 *   |   [ C(0,0) C(0,1)...C(n-1,n-1) ]  <- matriz C aplanada |
 *   +---------------------------------------------------------+
 *
 *   Las matrices se guardan en un solo bloque de memoria contigua
 *   (en lugar de punteros dispersos) para mejor rendimiento de caché.
 *
 * =============================================================================
 */

#include <stdio.h>      // printf, perror
#include <stdlib.h>     // malloc, free, exit, atoi, rand, srand
#include <time.h>       // clock_gettime, time
#include <unistd.h>     // fork, sysconf, _SC_NPROCESSORS_ONLN
#include <sys/wait.h>   // waitpid
#include <sys/mman.h>   // mmap, munmap, MAP_SHARED, MAP_ANONYMOUS
#include <semaphore.h>  // sem_t, sem_init, sem_wait, sem_post, sem_destroy
#include <string.h>     // (incluido por compatibilidad)

/* =============================================================================
 * ESTRUCTURA: SharedMemory
 * =============================================================================
 * Contiene todos los datos que los procesos comparten: las matrices,
 * el semáforo y variables de control. Todo vive en un bloque de
 * memoria creado con mmap(), visible por padre e hijos.
 */
typedef struct {
    int   n;          // Tamaño de la matriz (n x n)
    int   num_procs;  // Número total de procesos hijos lanzados
    sem_t mutex;      // Semáforo binario para exclusión mutua en región crítica
    int   log_count;  // Contador de celdas escritas en C (protegido por mutex)
    /*
     * data[] es un "flexible array member": un array sin tamaño fijo al final
     * de la estructura. Permite que A, B y C queden en memoria contigua.
     *
     * Distribución dentro de data[]:
     *   Posiciones [0        .. n*n-1    ] -> Matriz A
     *   Posiciones [n*n      .. 2*n*n-1  ] -> Matriz B
     *   Posiciones [2*n*n    .. 3*n*n-1  ] -> Matriz C (resultado)
     */
    int data[];
} SharedMemory;

/* =============================================================================
 * MACROS DE ACCESO A LAS MATRICES
 * =============================================================================
 * Permiten usar A(shm,i,j), B(shm,i,j) y C(shm,i,j) como si fueran
 * matrices 2D normales, aunque internamente están aplanadas en data[].
 *
 * Fórmula para convertir (fila, columna) a índice en array 1D:
 *   índice = fila * n + columna
 *
 * Ejemplo para n=3, elemento [1][2]:
 *   A(shm,1,2) = data[1*3 + 2] = data[5]
 */
#define A(shm, i, j) (shm->data[(i) * shm->n + (j)])
#define B(shm, i, j) (shm->data[shm->n * shm->n + (i) * shm->n + (j)])
#define C(shm, i, j) (shm->data[2 * shm->n * shm->n + (i) * shm->n + (j)])

/* =============================================================================
 * FUNCIÓN: llenarMatrices
 * =============================================================================
 * Rellena las matrices A y B con números enteros aleatorios entre 0 y 100.
 * Se llama ANTES de hacer fork(), así todos los procesos ven los mismos datos.
 *
 * Parámetros:
 *   shm -> puntero a la memoria compartida donde están A y B
 */
void llenarMatrices(SharedMemory* shm) {
    int n = shm->n;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A(shm, i, j) = rand() % 101;  // valor entre 0 y 100
            B(shm, i, j) = rand() % 101;
        }
    }
}

/* =============================================================================
 * FUNCIÓN: trabajoProceso
 * =============================================================================
 * Esta es la función que ejecuta CADA PROCESO HIJO.
 * Calcula las filas de C que le fueron asignadas: C = A x B
 *
 * La multiplicación de matrices funciona así:
 *   C[i][j] = suma de A[i][k] * B[k][j]  para k = 0 hasta n-1
 *
 * Ejemplo para C[0][0] con n=3:
 *   C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0]
 *
 * Parámetros:
 *   shm       -> memoria compartida (contiene A, B y C)
 *   start_row -> primera fila que debe calcular este proceso (inclusive)
 *   end_row   -> última fila que debe calcular (exclusive)
 *   proc_id   -> identificador del proceso (0, 1, 2...) para debug
 *
 * NOTA SOBRE LA REGIÓN CRÍTICA:
 *   Cada proceso escribe en filas distintas de C, por lo que en este
 *   caso específico no habría conflicto real. Sin embargo, usamos el
 *   semáforo para demostrar el mecanismo y para proteger log_count,
 *   que SÍ es un dato compartido que todos los procesos modifican.
 */
void trabajoProceso(SharedMemory* shm, int start_row, int end_row, int proc_id) {
    int n = shm->n;

    for (int i = start_row; i < end_row; i++) {
        for (int j = 0; j < n; j++) {

            /* Calcular el valor de C[i][j]
             * Usamos long long para evitar desbordamiento en matrices grandes
             * (100 * 100 * n puede superar el límite de int para n > 214748) */
            long long sum = 0;
            for (int k = 0; k < n; k++) {
                sum += A(shm, i, k) * B(shm, k, j);
            }

            /* -----------------------------------------------------------------
             * REGIÓN CRÍTICA — inicio
             * -----------------------------------------------------------------
             * sem_wait() intenta "tomar" el semáforo:
             *   - Si está libre (valor=1): lo toma (valor->0) y entra
             *   - Si está ocupado (valor=0): se BLOQUEA hasta que se libere
             *
             * Solo UN proceso puede estar dentro de este bloque a la vez.
             */
            sem_wait(&shm->mutex);

                C(shm, i, j) = (int) sum;  // escribir resultado en C compartida
                shm->log_count++;           // incrementar contador compartido

            /* -----------------------------------------------------------------
             * REGIÓN CRÍTICA — fin
             * -----------------------------------------------------------------
             * sem_post() libera el semáforo (valor->1), permitiendo que
             * el siguiente proceso en espera pueda entrar.
             */
            sem_post(&shm->mutex);
        }
    }
}

/* =============================================================================
 * FUNCIÓN: main
 * =============================================================================
 * Punto de entrada del programa. Orquesta todo el flujo:
 *   1. Validar argumentos
 *   2. Crear memoria compartida (mmap)
 *   3. Inicializar semáforo
 *   4. Llenar matrices A y B
 *   5. FORK -> lanzar procesos hijos
 *   6. JOIN -> esperar que terminen (waitpid)
 *   7. Medir e imprimir tiempo
 *   8. Liberar recursos
 */
int main(int argc, char* argv[]) {

    /* -------------------------------------------------------------------------
     * PASO 1: Validar y leer argumentos de línea de comandos
     * -------------------------------------------------------------------------
     * argc = número de argumentos (incluyendo el nombre del programa)
     * argv[0] = nombre del programa
     * argv[1] = tamaño de la matriz (obligatorio)
     * argv[2] = número de procesos  (opcional)
     */
    if (argc < 2 || argc > 3) {
        printf("Uso: %s <tamaño de matriz> [numero de procesos]\n", argv[0]);
        printf("Ejemplo: %s 1000 4\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);  // convertir string "1000" a entero 1000
    if (n <= 0) {
        printf("Error: el tamaño de la matriz debe ser un número positivo.\n");
        return 1;
    }

    /*
     * Si no se especifica número de procesos, usar todos los núcleos
     * disponibles del CPU (detección automática con sysconf).
     */
    int num_procs = (argc == 3) ? atoi(argv[2]) : sysconf(_SC_NPROCESSORS_ONLN);
    if (num_procs <= 0) num_procs = 1;
    if (num_procs > n)  num_procs = n;  // no tiene sentido más procesos que filas

    /* -------------------------------------------------------------------------
     * PASO 2: Crear memoria compartida con mmap()
     * -------------------------------------------------------------------------
     * mmap() reserva un bloque de memoria y lo mapea en el espacio de
     * direcciones del proceso. Con MAP_SHARED + MAP_ANONYMOUS:
     *   - MAP_SHARED:    los cambios son visibles por todos los procesos
     *   - MAP_ANONYMOUS: no está asociado a ningún archivo, es RAM pura
     *
     * Tamaño del bloque:
     *   sizeof(SharedMemory)    -> cabecera con n, num_procs, mutex, log_count
     *   3 * n * n * sizeof(int) -> espacio para matrices A, B y C
     */
    size_t shm_size = sizeof(SharedMemory) + 3 * n * n * sizeof(int);

    SharedMemory* shm = (SharedMemory*) mmap(
        NULL,                       // dejar que el SO elija la dirección
        shm_size,                   // tamaño total del bloque
        PROT_READ | PROT_WRITE,     // permisos: lectura y escritura
        MAP_SHARED | MAP_ANONYMOUS, // compartido, sin archivo
        -1,                         // sin descriptor de archivo (MAP_ANONYMOUS)
        0                           // offset 0
    );

    if (shm == MAP_FAILED) {
        perror("Error: mmap falló");
        return 1;
    }

    // Inicializar campos de control
    shm->n         = n;
    shm->num_procs = num_procs;
    shm->log_count = 0;

    /* -------------------------------------------------------------------------
     * PASO 3: Inicializar el semáforo
     * -------------------------------------------------------------------------
     * sem_init(semaforo, pshared, valor_inicial)
     *
     *   pshared = 0 -> compartido entre hilos del mismo proceso (pthread)
     *   pshared = 1 -> compartido entre procesos distintos (fork) <- este
     *
     *   valor_inicial = 1 -> semáforo binario (mutex):
     *     1 = libre   (nadie está en la región crítica)
     *     0 = ocupado (alguien está en la región crítica)
     */
    if (sem_init(&shm->mutex, 1, 1) != 0) {
        perror("Error: sem_init falló");
        munmap(shm, shm_size);
        return 1;
    }

    /* -------------------------------------------------------------------------
     * PASO 4: Llenar matrices A y B con datos aleatorios
     * -------------------------------------------------------------------------
     * IMPORTANTE: esto se hace ANTES del fork() para que todos los
     * procesos hijos vean las mismas matrices desde el inicio.
     */
    srand(time(NULL));       // semilla basada en tiempo actual para variedad
    llenarMatrices(shm);

    /* -------------------------------------------------------------------------
     * PASO 5: Calcular distribución de filas entre procesos
     * -------------------------------------------------------------------------
     * Se divide n filas entre num_procs procesos lo más equitativamente
     * posible. Si no divide exacto, los primeros procesos reciben 1 fila extra.
     *
     * Ejemplo: n=10, num_procs=3
     *   rows_per_proc = 10/3 = 3
     *   remaining     = 10%3 = 1
     *   Proceso 0 -> filas 0,1,2,3  (4 filas, recibe la extra)
     *   Proceso 1 -> filas 4,5,6    (3 filas)
     *   Proceso 2 -> filas 7,8,9    (3 filas)
     */
    int rows_per_proc  = n / num_procs;
    int remaining_rows = n % num_procs;

    // Iniciar medición de tiempo ANTES de lanzar los procesos
    struct timespec start_wall, end_wall;
    clock_gettime(CLOCK_MONOTONIC, &start_wall);

    // Array para guardar los PIDs (identificadores) de los procesos hijos
    pid_t* pids   = (pid_t*) malloc(num_procs * sizeof(pid_t));
    int start_row = 0;

    /* -------------------------------------------------------------------------
     * PASO 6: FORK — lanzar procesos hijos
     * -------------------------------------------------------------------------
     * Por cada proceso p:
     *   - Se calcula qué filas le corresponden
     *   - Se llama fork() que duplica el proceso
     *   - El HIJO (pid==0) ejecuta su porción y termina con exit(0)
     *   - El PADRE (pid>0) guarda el PID del hijo y continúa el bucle
     *
     * Después del bucle, el padre ha lanzado todos los hijos y todos
     * están trabajando en paralelo.
     */
    for (int p = 0; p < num_procs; p++) {
        int end_row = start_row + rows_per_proc + (p < remaining_rows ? 1 : 0);

        pids[p] = fork();  // <- aquí se crean los procesos hijos

        if (pids[p] < 0) {
            // fork() devuelve -1 si hubo error (memoria insuficiente, etc.)
            perror("Error: fork falló");
            exit(1);
        }

        if (pids[p] == 0) {
            /* =================================================================
             * CÓDIGO DEL PROCESO HIJO
             * =================================================================
             * fork() devuelve 0 en el proceso hijo.
             * El hijo solo debe hacer su trabajo y terminar.
             * NO debe seguir el bucle del padre (por eso exit(0) al final).
             */
            trabajoProceso(shm, start_row, end_row, p);
            exit(0);  // hijo termina limpiamente
        }

        /* Si llegamos aquí, somos el PADRE.
         * fork() devuelve el PID del hijo al padre (número > 0).
         * Guardamos el PID para luego poder esperar a ese hijo específico. */
        start_row = end_row;  // avanzar al siguiente bloque de filas
    }

    /* -------------------------------------------------------------------------
     * PASO 7: FORK-JOIN — padre espera a que todos los hijos terminen
     * -------------------------------------------------------------------------
     * waitpid(pid, status, opciones) espera a que el proceso con ese PID
     * termine. Es el equivalente a pthread_join() pero para procesos.
     *
     * Sin este paso, el padre podría leer C antes de que los hijos
     * terminen de escribirla -> resultado incorrecto.
     */
    for (int p = 0; p < num_procs; p++) {
        waitpid(pids[p], NULL, 0);  // esperar al hijo p específicamente
    }

    // Detener medición de tiempo DESPUÉS de que todos los hijos terminaron
    clock_gettime(CLOCK_MONOTONIC, &end_wall);

    double wall_time = (end_wall.tv_sec  - start_wall.tv_sec) +
                      (end_wall.tv_nsec - start_wall.tv_nsec) / 1e9;

    /* -------------------------------------------------------------------------
     * PASO 8: Mostrar resultado
     * -------------------------------------------------------------------------
     * Se imprime ÚNICAMENTE el tiempo en segundos con 6 decimales,
     * igual que la versión con hilos. Esto permite redirigir la salida
     * a un archivo con >> y procesarla con scripts de bash.
     *
     * Ejemplo de uso en script:
     *   ./matriz_procesos 1000 4 >> tiempo_procesos4.doc
     */
    printf("%.6f\n", wall_time);

    /* -------------------------------------------------------------------------
     * PASO 9: Liberar todos los recursos
     * -------------------------------------------------------------------------
     * Es importante limpiar en este orden:
     *   1. Destruir el semáforo (sem_destroy)
     *   2. Liberar la memoria compartida (munmap)
     *   3. Liberar el array de PIDs (free)
     */
    sem_destroy(&shm->mutex);   // destruir el semáforo
    munmap(shm, shm_size);      // liberar la memoria compartida
    free(pids);                 // liberar array de PIDs

    return 0;
}
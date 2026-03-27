/*
 * jacobi_poisson_1d_hilos.c
 *
 * Solución iterativa de Jacobi para la ecuación de Poisson en 1D.
 * Versión paralela con POSIX Threads (pthreads).
 *
 * Problema:
 *   -u''(x) = f(x),  x en [0, 1]
 *   u(0) = u(1) = 0
 *   f(x) = -x*(x+3)*exp(x)
 *   u_exacta(x) = x*(x-1)*exp(x)
 *
 * Argumentos:
 *   k         Índice de malla. Define n = 2^k + 1 puntos.
 *   tol       Tolerancia para la norma RMS del residuo (ej: 1e-6).
 *   num_hilos Número de hilos a usar (ej: 2, 4, 8).
 *
 * Uso:
 *   ./jacobi_hilos <k> <tol> <num_hilos>
 *   Ejemplo: ./jacobi_hilos 10 1e-6 4
 *
 * Compilación:
 *   gcc -O2 -o jacobi_hilos jacobi_poisson_1d_hilos.c -lm -lpthread
 *
 * =============================================================================
 * ESTRATEGIA DE PARALELISMO
 * =============================================================================
 * La iteración de Jacobi actualiza cada nodo interior j con:
 *   u_new[j] = (f[j]*h² + u_old[j-1] + u_old[j+1]) / 2
 *
 * Cada nodo es INDEPENDIENTE de los demás en la misma iteración
 * (usa u_old, no u_new), lo que hace que Jacobi sea naturalmente
 * paralelizable: se divide el rango de nodos interiores entre los
 * hilos y cada uno actualiza su porción sin conflictos de escritura.
 *
 * SINCRONIZACIÓN CON BARRERAS:
 * Al final de cada iteración todos los hilos deben haber terminado
 * antes de iniciar la siguiente (la iteración k+1 necesita los valores
 * completos de u_old calculados en la iteración k). Para esto se usan
 * dos barreras por iteración:
 *
 *   BARRERA 1: todos copiaron u → u_old antes de que alguien actualice u
 *   BARRERA 2: todos actualizaron u y calcularon su residuo parcial
 *              → el hilo 0 acumula los parciales y decide si continuar
 *   BARRERA 3: todos esperan la decisión del hilo 0 antes de continuar
 *
 * La variable 'seguir' es escrita por el hilo 0 y leída por todos,
 * protegida por un mutex para evitar condiciones de carrera.
 * =============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

/* ─── Límite de seguridad de iteraciones ─────────────────────────────────── */
#define MAX_ITER 5000000

/* ─── Prototipos ─────────────────────────────────────────────────────────── */
double exact(double x);
double force(double x);
void   solve_direct(int n, double h, double *f, double *ud);

/* =============================================================================
 * CONTEXTO GLOBAL COMPARTIDO ENTRE TODOS LOS HILOS
 * =============================================================================
 * Un único bloque de datos compartido evita tener que pasar referencias
 * cruzadas entre estructuras. Los hilos leen de aquí los parámetros del
 * problema y escriben sus resultados parciales en sus ranuras del array.
 */
typedef struct {
    /* Parámetros del problema (solo lectura) */
    int     n;
    double  h2;
    double  inv_h2;
    double *f;
    double *u;
    double *u_old;
    double  tol;
    int     num_hilos;

    /* Distribución de trabajo por hilo */
    int    *start;          /* start[t] = primer nodo del hilo t           */
    int    *end;            /* end[t]   = último nodo + 1 del hilo t       */

    /* Resultados parciales por hilo (cada hilo escribe solo en su ranura) */
    double *res_sq_parcial;
    double *change_sq_parcial;

    /* Control de iteración (hilo 0 escribe, todos leen) */
    int     seguir;
    int     it_num;

    /* Sincronización */
    pthread_barrier_t barrera;
    pthread_mutex_t   mutex_ctrl;
} ContextoJacobi;

/* =============================================================================
 * FUNCIÓN DEL HILO
 * =============================================================================
 * Cada hilo recibe su id y un puntero al contexto global compartido.
 * Ejecuta iteraciones de Jacobi sobre su rango de nodos hasta que
 * el hilo 0 señale convergencia o límite de iteraciones.
 */
typedef struct {
    int              id;
    ContextoJacobi  *ctx;
} HiloArgs;

void *jacobi_worker(void *arg)
{
    HiloArgs       *a   = (HiloArgs *) arg;
    ContextoJacobi *ctx = a->ctx;
    int             id  = a->id;
    int             s   = ctx->start[id];
    int             e   = ctx->end[id];

    while (1) {
        /* ── PASO 1: Copiar u → u_old en la porción de este hilo ── */
        for (int j = s; j < e; j++)
            ctx->u_old[j] = ctx->u[j];

        /* El hilo 0 copia también las fronteras (necesarias para el residuo) */
        if (id == 0) {
            ctx->u_old[0]          = ctx->u[0];
            ctx->u_old[ctx->n - 1] = ctx->u[ctx->n - 1];
        }

        /* ── BARRERA 1: todos han terminado de copiar u_old ──────────────
         * Ningún hilo puede empezar a actualizar u hasta que todos hayan
         * copiado u_old, porque la actualización de u[j] usa u_old[j±1]
         * que puede pertenecer al rango de otro hilo.
         */
        pthread_barrier_wait(&ctx->barrera);

        /* ── PASO 2: Actualizar u[j] con la fórmula de Jacobi ── */
        for (int j = s; j < e; j++) {
            ctx->u[j] = (ctx->f[j] * ctx->h2
                         + ctx->u_old[j-1]
                         + ctx->u_old[j+1]) * 0.5;
        }

        /* ── PASO 3: Calcular residuo y cambio parciales ── */
        double res_sq    = 0.0;
        double change_sq = 0.0;

        for (int j = s; j < e; j++) {
            double r_j = (2.0 * ctx->u_old[j]
                          - ctx->u_old[j-1]
                          - ctx->u_old[j+1]) * ctx->inv_h2
                         - ctx->f[j];
            res_sq    += r_j * r_j;
            double c_j = ctx->u[j] - ctx->u_old[j];
            change_sq += c_j * c_j;
        }

        /* Escribir en la ranura propia — sin conflicto con otros hilos */
        ctx->res_sq_parcial[id]    = res_sq;
        ctx->change_sq_parcial[id] = change_sq;

        /* ── BARRERA 2: todos han calculado su residuo parcial ──────────
         * El hilo 0 necesita que TODOS hayan escrito su parcial antes
         * de poder acumularlos.
         */
        pthread_barrier_wait(&ctx->barrera);

        /* ── PASO 4: Solo el hilo 0 acumula, calcula RMS y decide ── */
        if (id == 0) {
            double total_res    = 0.0;
            double total_change = 0.0;

            for (int t = 0; t < ctx->num_hilos; t++) {
                total_res    += ctx->res_sq_parcial[t];
                total_change += ctx->change_sq_parcial[t];
            }

            double res_rms = sqrt(total_res / ctx->n);

            ctx->it_num++;

            /* REGIÓN CRÍTICA: escribir 'seguir' */
            pthread_mutex_lock(&ctx->mutex_ctrl);
            if (res_rms <= ctx->tol || ctx->it_num >= MAX_ITER)
                ctx->seguir = 0;
            pthread_mutex_unlock(&ctx->mutex_ctrl);
        }

        /* ── BARRERA 3: todos esperan la decisión del hilo 0 ────────────
         * Ningún hilo puede leer 'seguir' antes de que el hilo 0 lo
         * haya actualizado.
         */
        pthread_barrier_wait(&ctx->barrera);

        /* ── PASO 5: Leer 'seguir' y decidir si continuar ── */
        pthread_mutex_lock(&ctx->mutex_ctrl);
        int continuar = ctx->seguir;
        pthread_mutex_unlock(&ctx->mutex_ctrl);

        if (!continuar) break;
    }

    pthread_exit(NULL);
}

/* ════════════════════════════════════════════════════════════════════════════
 * jacobi_paralelo
 *
 * Orquesta la iteración de Jacobi paralela:
 *   1. Inicializa el contexto y la sincronización
 *   2. Distribuye nodos interiores entre hilos
 *   3. FORK: lanza los hilos con pthread_create
 *   4. JOIN: espera que terminen con pthread_join
 *   5. Libera recursos y retorna el número de iteraciones
 * ════════════════════════════════════════════════════════════════════════════ */
int jacobi_paralelo(int n, double h, double *f, double *u, double tol, int num_hilos)
{
    /* ── Reservar u_old compartido ── */
    double *u_old = (double *) malloc(n * sizeof(double));
    if (!u_old) { fprintf(stderr, "Error: malloc u_old\n"); return -1; }
    for (int j = 0; j < n; j++) u_old[j] = 0.0;

    /* ── Arrays de distribución y parciales ── */
    int    *start            = (int *)    malloc(num_hilos * sizeof(int));
    int    *end              = (int *)    malloc(num_hilos * sizeof(int));
    double *res_sq_parcial   = (double *) malloc(num_hilos * sizeof(double));
    double *change_sq_parcial= (double *) malloc(num_hilos * sizeof(double));

    if (!start || !end || !res_sq_parcial || !change_sq_parcial) {
        fprintf(stderr, "Error: malloc arrays de control\n");
        free(u_old); free(start); free(end);
        free(res_sq_parcial); free(change_sq_parcial);
        return -1;
    }

    /* ── Distribuir nodos interiores entre hilos ── */
    int interior     = n - 2;
    int nodos_x_hilo = interior / num_hilos;
    int restantes    = interior % num_hilos;
    int inicio       = 1;

    for (int t = 0; t < num_hilos; t++) {
        start[t] = inicio;
        end[t]   = inicio + nodos_x_hilo + (t < restantes ? 1 : 0);
        inicio   = end[t];
        res_sq_parcial[t]    = 0.0;
        change_sq_parcial[t] = 0.0;
    }

    /* ── Inicializar contexto compartido ── */
    ContextoJacobi ctx;
    ctx.n                 = n;
    ctx.h2                = h * h;
    ctx.inv_h2            = 1.0 / (h * h);
    ctx.f                 = f;
    ctx.u                 = u;
    ctx.u_old             = u_old;
    ctx.tol               = tol;
    ctx.num_hilos         = num_hilos;
    ctx.start             = start;
    ctx.end               = end;
    ctx.res_sq_parcial    = res_sq_parcial;
    ctx.change_sq_parcial = change_sq_parcial;
    ctx.seguir            = 1;
    ctx.it_num            = 0;

    pthread_barrier_init(&ctx.barrera,    NULL, num_hilos);
    pthread_mutex_init (&ctx.mutex_ctrl,  NULL);

    /* ── Crear argumentos por hilo ── */
    HiloArgs  *hilo_args = (HiloArgs *)  malloc(num_hilos * sizeof(HiloArgs));
    pthread_t *hilos     = (pthread_t *) malloc(num_hilos * sizeof(pthread_t));

    if (!hilo_args || !hilos) {
        fprintf(stderr, "Error: malloc hilos\n");
        free(u_old); free(start); free(end);
        free(res_sq_parcial); free(change_sq_parcial);
        free(hilo_args); free(hilos);
        return -1;
    }

    /* ── FORK: lanzar hilos ── */
    for (int t = 0; t < num_hilos; t++) {
        hilo_args[t].id  = t;
        hilo_args[t].ctx = &ctx;
        pthread_create(&hilos[t], NULL, jacobi_worker, (void *) &hilo_args[t]);
    }

    /* ── JOIN: esperar a todos los hilos ── */
    for (int t = 0; t < num_hilos; t++)
        pthread_join(hilos[t], NULL);

    int resultado = ctx.it_num;

    /* ── Liberar recursos ── */
    pthread_barrier_destroy(&ctx.barrera);
    pthread_mutex_destroy  (&ctx.mutex_ctrl);
    free(u_old); free(start); free(end);
    free(res_sq_parcial); free(change_sq_parcial);
    free(hilo_args); free(hilos);

    return resultado;
}

/* ════════════════════════════════════════════════════════════════════════════
 * MAIN
 * ════════════════════════════════════════════════════════════════════════════ */
int main(int argc, char *argv[])
{
    if (argc != 4) {
        fprintf(stderr, "\nUso: %s <k> <tol> <num_hilos>\n", argv[0]);
        fprintf(stderr, "  k          Indice de malla  (entero >= 0)\n");
        fprintf(stderr, "  tol        Tolerancia RMS   (real > 0, ej: 1e-6)\n");
        fprintf(stderr, "  num_hilos  Numero de hilos  (entero >= 1, ej: 4)\n");
        fprintf(stderr, "\nEjemplo: %s 10 1e-6 4\n\n", argv[0]);
        return 1;
    }

    int    k         = atoi(argv[1]);
    double tol       = atof(argv[2]);
    int    num_hilos = atoi(argv[3]);

    if (k < 0)         { fprintf(stderr, "Error: k debe ser >= 0.\n");          return 1; }
    if (tol <= 0.0)    { fprintf(stderr, "Error: tol debe ser > 0.\n");          return 1; }
    if (num_hilos <= 0){ fprintf(stderr, "Error: num_hilos debe ser >= 1.\n");   return 1; }

    double a  = 0.0, b  = 1.0;
    double ua = 0.0, ub = 0.0;

    int    nk = (1 << k) + 1;
    double hk = (b - a) / (nk - 1);

    if (num_hilos > nk - 2) {
        num_hilos = nk - 2;
        fprintf(stderr, "Aviso: num_hilos reducido a %d (nodos interiores).\n", num_hilos);
    }

    double *xk  = (double *) malloc(nk * sizeof(double));
    double *fk  = (double *) malloc(nk * sizeof(double));
    double *uek = (double *) malloc(nk * sizeof(double));
    double *udk = (double *) malloc(nk * sizeof(double));
    double *ujk = (double *) malloc(nk * sizeof(double));

    if (!xk || !fk || !uek || !udk || !ujk) {
        fprintf(stderr, "Error: memoria insuficiente para n = %d nodos.\n", nk);
        free(xk); free(fk); free(uek); free(udk); free(ujk);
        return 1;
    }

    for (int j = 0; j < nk; j++) {
        xk[j]  = a + j * hk;
        fk[j]  = force(xk[j]);
        uek[j] = exact(xk[j]);
        ujk[j] = 0.0;
    }
    fk[0] = ua; fk[nk-1] = ub;

    for (int j = 0; j < nk; j++) udk[j] = fk[j];
    solve_direct(nk, hk, fk, udk);

    printf("\n");
    printf("========================================\n");
    printf(" JACOBI_POISSON_1D  (C - HILOS)\n");
    printf("========================================\n");
    printf(" Indice de malla  k  = %d\n",   k);
    printf(" Numero de nodos  nk = %d\n",   nk);
    printf(" Espaciado        hk = %.6e\n", hk);
    printf(" Tolerancia RMS      = %.6e\n", tol);
    printf(" Numero de hilos     = %d\n",   num_hilos);
    printf(" Limite iteraciones  = %d\n",   MAX_ITER);
    printf("========================================\n\n");

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int it_num = jacobi_paralelo(nk, hk, fk, ujk, tol, num_hilos);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double seg = (t1.tv_sec  - t0.tv_sec)
               + (t1.tv_nsec - t0.tv_nsec) * 1e-9;

    /* VERIF 2: Residuo RMS final */
    double res_final_sq = 0.0;
    double inv_h2_main  = 1.0 / (hk * hk);
    for (int j = 1; j < nk - 1; j++) {
        double r_j = (2.0*ujk[j] - ujk[j-1] - ujk[j+1]) * inv_h2_main - fk[j];
        res_final_sq += r_j * r_j;
    }
    double res_final_rms = sqrt(res_final_sq / nk);

    /* VERIF 3: Error RMS vs exacta */
    double rms_disc = 0.0;
    for (int j = 0; j < nk; j++) {
        double diff = uek[j] - ujk[j];
        rms_disc += diff * diff;
    }
    rms_disc = sqrt(rms_disc / nk);

    /* Salida idéntica al programa original */
    printf("\n");
    printf("========================================\n");
    printf(" RESULTADO DEL EXPERIMENTO\n");
    printf("========================================\n");
    printf(" k                              : %d\n",     k);
    printf(" n (nodos)                      : %d\n",     nk);
    printf(" Tolerancia solicitada          : %.6e\n",   tol);
    printf("----------------------------------------\n");
    printf(" [TIEMPO]  Wall clock time      : %.9f s\n", seg);
    printf("----------------------------------------\n");
    printf(" [VERIF 1] Iteraciones          : %d\n",     it_num);
    printf(" [VERIF 2] Residuo RMS final    : %.10e\n",  res_final_rms);
    printf(" [VERIF 3] Error RMS vs exacta  : %.10e\n",  rms_disc);
    if (it_num >= MAX_ITER)
        printf(" ADVERTENCIA: limite de iteraciones alcanzado.\n");
    printf("========================================\n\n");
    printf("JACOBI_POISSON_1D: fin normal de ejecucion.\n\n");

    free(xk); free(fk); free(uek); free(udk); free(ujk);
    return 0;
}

/* ════════════════════════════════════════════════════════════════════════════
 * exact  -  u(x) = x*(x-1)*exp(x)
 * ════════════════════════════════════════════════════════════════════════════ */
double exact(double x) { return x * (x - 1.0) * exp(x); }

/* ════════════════════════════════════════════════════════════════════════════
 * force  -  f(x) = -x*(x+3)*exp(x)
 * ════════════════════════════════════════════════════════════════════════════ */
double force(double x) { return -x * (x + 3.0) * exp(x); }

/* ════════════════════════════════════════════════════════════════════════════
 * solve_direct  -  algoritmo de Thomas para sistema tridiagonal
 * ════════════════════════════════════════════════════════════════════════════ */
void solve_direct(int n, double h, double *f, double *ud)
{
    double *c_prime = (double *) malloc(n * sizeof(double));
    double *d_prime = (double *) malloc(n * sizeof(double));
    if (!c_prime || !d_prime) {
        fprintf(stderr, "Error de memoria en solve_direct.\n");
        free(c_prime); free(d_prime); return;
    }

    double h2 = h * h;
    c_prime[0] = 0.0;
    d_prime[0] = f[0];

    for (int j = 1; j < n - 1; j++) {
        double diag_j = 2.0 / h2;
        double sub_j  = -1.0 / h2;
        double sup_j  = -1.0 / h2;
        double denom  = diag_j - sub_j * c_prime[j-1];
        c_prime[j]    = sup_j  / denom;
        d_prime[j]    = (f[j] - sub_j * d_prime[j-1]) / denom;
    }

    c_prime[n-1] = 0.0;
    d_prime[n-1] = f[n-1];

    ud[n-1] = d_prime[n-1];
    for (int j = n-2; j >= 0; j--)
        ud[j] = d_prime[j] - c_prime[j] * ud[j+1];

    free(c_prime); free(d_prime);
}

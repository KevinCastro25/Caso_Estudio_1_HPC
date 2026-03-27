/*
 * jacobi_poisson_1d_procesos.c
 *
 * Solución iterativa de Jacobi para la ecuación de Poisson en 1D.
 * Versión paralela con procesos (fork/waitpid) y memoria compartida (mmap).
 *
 * Uso:
 *   ./jacobi_procesos <k> <tol> <num_procs>
 *   Ejemplo: ./jacobi_procesos 10 1e-6 4
 *
 * Compilación:
 *   gcc -O2 -o jacobi_procesos jacobi_poisson_1d_procesos.c -lm
 *
 * =============================================================================
 * SINCRONIZACIÓN: BARRERA CON DOS SEMÁFOROS (fase entrada + fase salida)
 * =============================================================================
 * Una barrera simple de un semáforo tiene una condición de carrera:
 * el último proceso que llega despierta a los demás, pero él mismo puede
 * avanzar antes de que los demás hayan salido de la barrera, corrompiendo
 * el contador para la siguiente barrera.
 *
 * La solución es una barrera de DOS FASES:
 *
 *   FASE 1 — Entrada: todos decrementan el contador y esperan en sem_entrada.
 *             El último (contador==0) resetea el contador y abre la salida.
 *
 *   FASE 2 — Salida: todos decrementan el contador2 y esperan en sem_salida.
 *             El último en salir resetea contador2 y abre la entrada de nuevo.
 *
 * Esto garantiza que ningún proceso puede entrar a la SIGUIENTE barrera
 * hasta que TODOS hayan salido de la ACTUAL.
 * =============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/mman.h>
#include <semaphore.h>

#define MAX_ITER 5000000

double exact(double x);
double force(double x);
void   solve_direct(int n, double h, double *f, double *ud);

/* =============================================================================
 * MEMORIA COMPARTIDA
 * =============================================================================
 * Layout de data[] al final de la estructura:
 *   [0        .. n-1    ]  → u
 *   [n        .. 2n-1   ]  → u_old
 *   [2n       .. 3n-1   ]  → f
 *   [3n       .. 3n+p-1 ]  → res_sq_parcial   (p doubles)
 *   [3n+p     .. 3n+2p-1]  → change_sq_parcial(p doubles)
 * Seguido de 2*p ints para start[] y end[].
 */
typedef struct {
    int    n;
    int    num_procs;
    double h2;
    double inv_h2;
    double tol;
    int    seguir;
    int    it_num;

    /* Barrera de dos fases */
    sem_t  sem_mutex;       /* protege it_num y seguir                     */
    sem_t  sem_entrada;     /* fase 1: espera a que todos lleguen          */
    sem_t  sem_salida;      /* fase 2: espera a que todos salgan           */
    sem_t  sem_cnt;         /* protege contador1 y contador2               */
    int    contador1;       /* cuántos faltan por llegar (fase entrada)    */
    int    contador2;       /* cuántos faltan por salir  (fase salida)     */

    double data[];
} SharedCtx;

/* =============================================================================
 * BARRERA DE DOS FASES
 * =============================================================================
 * Garantiza que ningún proceso pueda pasar a la siguiente iteración
 * hasta que TODOS hayan completado la actual.
 */
static void barrera(SharedCtx *shm, int np)
{
    /* ── FASE 1: todos llegan ── */
    sem_wait(&shm->sem_cnt);
    shm->contador1--;
    int ultimo = (shm->contador1 == 0);
    sem_post(&shm->sem_cnt);

    if (ultimo) {
        /* Resetear contador2 para la fase salida */
        sem_wait(&shm->sem_cnt);
        shm->contador2 = np;
        sem_post(&shm->sem_cnt);
        /* Abrir la puerta de entrada para todos */
        for (int i = 0; i < np; i++)
            sem_post(&shm->sem_entrada);
    }
    sem_wait(&shm->sem_entrada);  /* esperar a que el último abra */

    /* ── FASE 2: todos salen ── */
    sem_wait(&shm->sem_cnt);
    shm->contador2--;
    int ultimo2 = (shm->contador2 == 0);
    sem_post(&shm->sem_cnt);

    if (ultimo2) {
        /* Resetear contador1 para la próxima barrera */
        sem_wait(&shm->sem_cnt);
        shm->contador1 = np;
        sem_post(&shm->sem_cnt);
        /* Abrir la puerta de salida para todos */
        for (int i = 0; i < np; i++)
            sem_post(&shm->sem_salida);
    }
    sem_wait(&shm->sem_salida);   /* esperar a que el último abra */
}

/* =============================================================================
 * FUNCIÓN DEL PROCESO HIJO
 * =============================================================================
 */
static void proceso_worker(int id, SharedCtx *shm)
{
    int     n    = shm->n;
    int     np   = shm->num_procs;
    double *u    = &shm->data[0];
    double *u_old= &shm->data[n];
    double *f    = &shm->data[2*n];
    double *res_p= &shm->data[3*n];
    int    *start= (int *)&shm->data[3*n + 2*np];
    int    *end  = start + np;

    int s = start[id];
    int e = end[id];

    while (1) {
        /* PASO 1: Copiar u → u_old */
        for (int j = s; j < e; j++)
            u_old[j] = u[j];
        if (id == 0) {
            u_old[0]     = u[0];
            u_old[n - 1] = u[n - 1];
        }

        /* BARRERA 1: todos copiaron u_old */
        barrera(shm, np);

        /* PASO 2: Actualizar u con Jacobi */
        for (int j = s; j < e; j++)
            u[j] = (f[j] * shm->h2 + u_old[j-1] + u_old[j+1]) * 0.5;

        /* PASO 3: Calcular residuo parcial con u_old */
        double res_sq = 0.0;
        for (int j = s; j < e; j++) {
            double r_j = (2.0 * u_old[j] - u_old[j-1] - u_old[j+1])
                         * shm->inv_h2 - f[j];
            res_sq += r_j * r_j;
        }

        /* REGIÓN CRÍTICA: escribir parcial propio */
        sem_wait(&shm->sem_mutex);
        res_p[id] = res_sq;
        sem_post(&shm->sem_mutex);

        /* BARRERA 2: todos calcularon parciales */
        barrera(shm, np);

        /* PASO 4: proceso 0 acumula y decide */
        if (id == 0) {
            double total = 0.0;
            for (int t = 0; t < np; t++)
                total += res_p[t];

            double res_rms = sqrt(total / n);

            sem_wait(&shm->sem_mutex);
            shm->it_num++;
            if (res_rms <= shm->tol || shm->it_num >= MAX_ITER)
                shm->seguir = 0;
            sem_post(&shm->sem_mutex);
        }

        /* BARRERA 3: todos esperan decisión del proceso 0 */
        barrera(shm, np);

        /* PASO 5: leer 'seguir' */
        sem_wait(&shm->sem_mutex);
        int continuar = shm->seguir;
        sem_post(&shm->sem_mutex);

        if (!continuar) break;
    }
}

/* =============================================================================
 * jacobi_procesos
 * =============================================================================
 */
int jacobi_procesos(int n, double h, double *f_in, double *u_out,
                    double tol, int num_procs)
{
    size_t sz = sizeof(SharedCtx)
              + (3*n + 2*num_procs) * sizeof(double)
              + (2*num_procs)       * sizeof(int);

    SharedCtx *shm = (SharedCtx *) mmap(
        NULL, sz, PROT_READ | PROT_WRITE,
        MAP_SHARED | MAP_ANONYMOUS, -1, 0);

    if (shm == MAP_FAILED) { perror("mmap"); return -1; }

    shm->n          = n;
    shm->num_procs  = num_procs;
    shm->h2         = h * h;
    shm->inv_h2     = 1.0 / (h * h);
    shm->tol        = tol;
    shm->seguir     = 1;
    shm->it_num     = 0;
    shm->contador1  = num_procs;
    shm->contador2  = num_procs;

    /* pshared=1 → semáforos compartidos entre procesos */
    sem_init(&shm->sem_mutex,   1, 1);
    sem_init(&shm->sem_entrada, 1, 0);
    sem_init(&shm->sem_salida,  1, 0);
    sem_init(&shm->sem_cnt,     1, 1);

    /* Punteros dentro de data[] */
    double *u_shm = &shm->data[0];
    double *u_old = &shm->data[n];
    double *f_shm = &shm->data[2*n];
    double *res_p = &shm->data[3*n];
    int    *start = (int *)&shm->data[3*n + 2*num_procs];
    int    *end   = start + num_procs;

    for (int j = 0; j < n; j++) {
        u_shm[j] = u_out[j];
        u_old[j] = 0.0;
        f_shm[j] = f_in[j];
    }
    for (int t = 0; t < num_procs; t++) res_p[t] = 0.0;

    /* Distribuir nodos interiores */
    int interior = n - 2, nxp = interior / num_procs,
        rest = interior % num_procs, ini = 1;
    for (int t = 0; t < num_procs; t++) {
        start[t] = ini;
        end[t]   = ini + nxp + (t < rest ? 1 : 0);
        ini      = end[t];
    }

    /* FORK */
    pid_t *pids = malloc(num_procs * sizeof(pid_t));
    for (int t = 0; t < num_procs; t++) {
        pids[t] = fork();
        if (pids[t] < 0) { perror("fork"); exit(1); }
        if (pids[t] == 0) { proceso_worker(t, shm); exit(0); }
    }

    /* JOIN */
    for (int t = 0; t < num_procs; t++)
        waitpid(pids[t], NULL, 0);

    int resultado = shm->it_num;
    for (int j = 0; j < n; j++) u_out[j] = u_shm[j];

    sem_destroy(&shm->sem_mutex);
    sem_destroy(&shm->sem_entrada);
    sem_destroy(&shm->sem_salida);
    sem_destroy(&shm->sem_cnt);
    munmap(shm, sz);
    free(pids);

    return resultado;
}

/* ════════════════════════════════════════════════════════════════════════════
 * MAIN
 * ════════════════════════════════════════════════════════════════════════════ */
int main(int argc, char *argv[])
{
    if (argc != 4) {
        fprintf(stderr, "\nUso: %s <k> <tol> <num_procs>\n", argv[0]);
        fprintf(stderr, "  k          Indice de malla  (entero >= 0)\n");
        fprintf(stderr, "  tol        Tolerancia RMS   (real > 0, ej: 1e-6)\n");
        fprintf(stderr, "  num_procs  Numero de procesos (entero >= 1)\n");
        fprintf(stderr, "\nEjemplo: %s 10 1e-6 4\n\n", argv[0]);
        return 1;
    }

    int k = atoi(argv[1]); double tol = atof(argv[2]); int np = atoi(argv[3]);
    if (k < 0)   { fprintf(stderr,"Error: k >= 0\n");      return 1; }
    if (tol<=0)  { fprintf(stderr,"Error: tol > 0\n");     return 1; }
    if (np <= 0) { fprintf(stderr,"Error: procs >= 1\n");  return 1; }

    int nk = (1<<k)+1; double hk = 1.0/(nk-1);
    if (np > nk-2) { np = nk-2;
        fprintf(stderr,"Aviso: num_procs reducido a %d\n",np); }

    double *xk  = malloc(nk*sizeof(double));
    double *fk  = malloc(nk*sizeof(double));
    double *uek = malloc(nk*sizeof(double));
    double *udk = malloc(nk*sizeof(double));
    double *ujk = malloc(nk*sizeof(double));

    for (int j=0;j<nk;j++){
        xk[j]=j*hk; fk[j]=force(xk[j]); uek[j]=exact(xk[j]); ujk[j]=0.0;
    }
    fk[0]=0.0; fk[nk-1]=0.0;
    for (int j=0;j<nk;j++) udk[j]=fk[j];
    solve_direct(nk,hk,fk,udk);

    printf("\n========================================\n");
    printf(" JACOBI_POISSON_1D  (C - PROCESOS)\n");
    printf("========================================\n");
    printf(" Indice de malla  k  = %d\n",k);
    printf(" Numero de nodos  nk = %d\n",nk);
    printf(" Espaciado        hk = %.6e\n",hk);
    printf(" Tolerancia RMS      = %.6e\n",tol);
    printf(" Numero de procesos  = %d\n",np);
    printf(" Limite iteraciones  = %d\n",MAX_ITER);
    printf("========================================\n\n");

    struct timespec t0,t1;
    clock_gettime(CLOCK_MONOTONIC,&t0);
    int it_num = jacobi_procesos(nk,hk,fk,ujk,tol,np);
    clock_gettime(CLOCK_MONOTONIC,&t1);
    double seg=(t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;

    double res_sq=0.0, inv_h2=1.0/(hk*hk);
    for(int j=1;j<nk-1;j++){
        double r=(2.0*ujk[j]-ujk[j-1]-ujk[j+1])*inv_h2-fk[j];
        res_sq+=r*r;
    }
    double rms=0.0;
    for(int j=0;j<nk;j++){double d=uek[j]-ujk[j];rms+=d*d;}

    printf("\n========================================\n");
    printf(" RESULTADO DEL EXPERIMENTO\n");
    printf("========================================\n");
    printf(" k                              : %d\n",k);
    printf(" n (nodos)                      : %d\n",nk);
    printf(" Tolerancia solicitada          : %.6e\n",tol);
    printf("----------------------------------------\n");
    printf(" [TIEMPO]  Wall clock time      : %.9f s\n",seg);
    printf("----------------------------------------\n");
    printf(" [VERIF 1] Iteraciones          : %d\n",it_num);
    printf(" [VERIF 2] Residuo RMS final    : %.10e\n",sqrt(res_sq/nk));
    printf(" [VERIF 3] Error RMS vs exacta  : %.10e\n",sqrt(rms/nk));
    if(it_num>=MAX_ITER) printf(" ADVERTENCIA: limite alcanzado.\n");
    printf("========================================\n\n");
    printf("JACOBI_POISSON_1D: fin normal de ejecucion.\n\n");

    free(xk);free(fk);free(uek);free(udk);free(ujk);
    return 0;
}

double exact(double x) { return x*(x-1.0)*exp(x); }
double force(double x) { return -x*(x+3.0)*exp(x); }

void solve_direct(int n, double h, double *f, double *ud)
{
    double *c=malloc(n*sizeof(double)), *d=malloc(n*sizeof(double));
    double h2=h*h;
    c[0]=0.0; d[0]=f[0];
    for(int j=1;j<n-1;j++){
        double dg=2.0/h2, sb=-1.0/h2, sp=-1.0/h2;
        double dn=dg-sb*c[j-1];
        c[j]=sp/dn; d[j]=(f[j]-sb*d[j-1])/dn;
    }
    c[n-1]=0.0; d[n-1]=f[n-1];
    ud[n-1]=d[n-1];
    for(int j=n-2;j>=0;j--) ud[j]=d[j]-c[j]*ud[j+1];
    free(c); free(d);
}
/*
 * jacobi_poisson_1d.c
 *
 * Solución iterativa de Jacobi para la ecuación de Poisson en 1D.
 * Traducción y adaptación del código MATLAB de John Burkardt.
 *
 * Problema:
 *   -u''(x) = f(x),  x en [0, 1]
 *   u(0) = u(1) = 0
 *   f(x) = -x*(x+3)*exp(x)
 *   u_exacta(x) = x*(x-1)*exp(x)
 *
 * Argumentos:
 *   k    Índice de malla. Define n = 2^k + 1 puntos.
 *   tol  Tolerancia para la norma RMS del residuo (ej: 1e-6).
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* ─── Límite de seguridad de iteraciones ─────────────────────────────────── */
#define MAX_ITER 5000000

/* ─── Prototipos ─────────────────────────────────────────────────────────── */
double exact(double x);
double force(double x);
void   solve_direct(int n, double h, double *f, double *ud);
int    jacobi(int n, double h, double *f, double *u, double tol);

/* ════════════════════════════════════════════════════════════════════════════
 * MAIN
 * ════════════════════════════════════════════════════════════════════════════ */
int main(int argc, char *argv[])
{
    /* ── Lectura de argumentos por consola ── */
    if (argc != 3) {
        fprintf(stderr, "\nUso: %s <k> <tol>\n", argv[0]);
        fprintf(stderr, "  k    Indice de malla  (entero >= 0, recomendado: 5, 10, 14, 17)\n");
        fprintf(stderr, "  tol  Tolerancia RMS   (real > 0, recomendado: 1e-6)\n");
        fprintf(stderr, "\nEjemplo: %s 5 1e-6\n\n", argv[0]);
        return 1;
    }

    int    k   = atoi(argv[1]);
    double tol = atof(argv[2]);

    if (k < 0) {
        fprintf(stderr, "Error: k debe ser >= 0.\n");
        return 1;
    }
    if (tol <= 0.0) {
        fprintf(stderr, "Error: tol debe ser > 0.\n");
        return 1;
    }

    /* ── Parámetros del dominio (fijos según el documento) ── */
    double a  = 0.0;
    double b  = 1.0;
    double ua = 0.0;   /* u(0) = 0 */
    double ub = 0.0;   /* u(1) = 0 */

    /* ── Geometría de la malla ── */
    int    nk = (1 << k) + 1;          /* nk = 2^k + 1                    */
    double hk = (b - a) / (nk - 1);    /* espaciado uniforme               */

    /* ── Reserva de memoria ── */
    double *xk  = (double *)malloc(nk * sizeof(double)); /* nodos           */
    double *fk  = (double *)malloc(nk * sizeof(double)); /* lado derecho    */
    double *uek = (double *)malloc(nk * sizeof(double)); /* solución exacta */
    double *udk = (double *)malloc(nk * sizeof(double)); /* solución directa*/
    double *ujk = (double *)malloc(nk * sizeof(double)); /* solución Jacobi */

    if (!xk || !fk || !uek || !udk || !ujk) {
        fprintf(stderr, "Error: no hay memoria suficiente para n = %d nodos.\n", nk);
        free(xk); free(fk); free(uek); free(udk); free(ujk);
        return 1;
    }

    /* ── Construcción de la malla y vectores ── */
    for (int j = 0; j < nk; j++) {
        xk[j]  = a + j * hk;
        fk[j]  = force(xk[j]);
        uek[j] = exact(xk[j]);
        ujk[j] = 0.0;              /* aproximación inicial = cero          */
    }

    /* Condiciones de frontera en el lado derecho */
    fk[0]      = ua;
    fk[nk - 1] = ub;

    /* ── Solución directa (sustitución hacia atrás, sistema tridiagonal) ── */
    for (int j = 0; j < nk; j++) udk[j] = fk[j];
    solve_direct(nk, hk, fk, udk);

    /* ── Encabezado ── */
    printf("\n");
    printf("========================================\n");
    printf(" JACOBI_POISSON_1D  (C)\n");
    printf("========================================\n");
    printf(" Indice de malla  k  = %d\n",   k);
    printf(" Numero de nodos  nk = %d\n",   nk);
    printf(" Espaciado        hk = %.6e\n", hk);
    printf(" Tolerancia RMS      = %.6e\n", tol);
    printf(" Limite iteraciones  = %d\n",   MAX_ITER);
    printf("========================================\n\n");

    /* ── Iteración de Jacobi ── */
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int it_num = jacobi(nk, hk, fk, ujk, tol);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double seg = (t1.tv_sec  - t0.tv_sec)
               + (t1.tv_nsec - t0.tv_nsec) * 1e-9;

    /* ── VERIF 2: Residuo RMS final  r = A*u - f  sobre nodos interiores ──
     * Mide que tan bien satisface u la ecuacion discreta.
     * Es independiente de cualquier solucion de referencia externa.         */
    double res_final_sq = 0.0;
    double inv_h2_main  = 1.0 / (hk * hk);
    for (int j = 1; j < nk - 1; j++) {
        double r_j = (2.0*ujk[j] - ujk[j-1] - ujk[j+1]) * inv_h2_main - fk[j];
        res_final_sq += r_j * r_j;
    }
    double res_final_rms = sqrt(res_final_sq / nk);

    /* ── VERIF 3: Error RMS contra la solucion analitica cerrada ──
     * u_exacta(x) = x*(x-1)*exp(x) es una formula fija, independiente
     * de cualquier metodo numerico. Si una version paralela produce un
     * vector u distinto, este valor cambia y delata el error.               */
    double rms_disc = 0.0;
    for (int j = 0; j < nk; j++) {
        double diff = uek[j] - ujk[j];
        rms_disc += diff * diff;
    }
    rms_disc = sqrt(rms_disc / nk);

    /* ══════════════════════════════════════════════════════════════════════
     * SALIDA DEL EXPERIMENTO
     *
     * Formato pensado para ser parseado con:
     *   grep "TIEMPO\|VERIF" salida.txt
     *
     * Regla de verificacion:
     *   Dos implementaciones son equivalentes si y solo si producen
     *   los mismos [VERIF 1], [VERIF 2] y [VERIF 3] para igual k y tol.
     *   El [TIEMPO] puede (y debe) variar entre versiones optimizadas.
     * ══════════════════════════════════════════════════════════════════════ */
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
 * exact  –  solución analítica: u(x) = x*(x-1)*exp(x)
 * ════════════════════════════════════════════════════════════════════════════ */
double exact(double x)
{
    return x * (x - 1.0) * exp(x);
}

/* ════════════════════════════════════════════════════════════════════════════
 * force  –  término forzante: f(x) = -x*(x+3)*exp(x)
 * ════════════════════════════════════════════════════════════════════════════ */
double force(double x)
{
    return -x * (x + 3.0) * exp(x);
}

/* ════════════════════════════════════════════════════════════════════════════
 * solve_direct
 *
 * Resuelve el sistema tridiagonal A*u = f con el algoritmo de Thomas
 * (eliminación gaussiana para matrices tridiagonales). O(n) en tiempo
 * y memoria, equivalente al operador "\" de MATLAB para este sistema.
 *
 * La matriz (sin escalar por h²) tiene la forma:
 *   fila 0    :  [  1,  0,  0, ... ]   (frontera izquierda)
 *   fila j    :  [ -1,  2, -1, ... ]   dividido por h²
 *   fila nk-1 :  [  0,  0,  0,  1 ]   (frontera derecha)
 *
 * Entradas:
 *   n   número de nodos
 *   h   espaciado de malla
 *   f   lado derecho (NO se modifica; se copia internamente)
 *   ud  vector de salida (se sobreescribe con la solución)
 * ════════════════════════════════════════════════════════════════════════════ */
void solve_direct(int n, double h, double *f, double *ud)
{
    double *c_prime = (double *)malloc(n * sizeof(double));
    double *d_prime = (double *)malloc(n * sizeof(double));

    if (!c_prime || !d_prime) {
        fprintf(stderr, "Error de memoria en solve_direct.\n");
        free(c_prime); free(d_prime);
        return;
    }

    double h2 = h * h;

    /* Coeficientes del sistema tridiagonal escalado por h²:
       sub = -1,  diag = 2,  sup = -1   para filas interiores
       diag =  1                         para filas de frontera    */

    /* Barrido hacia adelante */
    c_prime[0] = 0.0;
    d_prime[0] = f[0];          /* u[0] = ua */

    for (int j = 1; j < n - 1; j++) {
        /* diag efectivo después de la eliminación */
        double diag_j = 2.0 / h2;
        double sub_j  = -1.0 / h2;
        double sup_j  = -1.0 / h2;

        double denom  = diag_j - sub_j * c_prime[j - 1];
        c_prime[j]    = sup_j  / denom;
        d_prime[j]    = (f[j] - sub_j * d_prime[j - 1]) / denom;
    }

    c_prime[n - 1] = 0.0;
    d_prime[n - 1] = f[n - 1];  /* u[n-1] = ub */

    /* Sustitución hacia atrás */
    ud[n - 1] = d_prime[n - 1];
    for (int j = n - 2; j >= 0; j--)
        ud[j] = d_prime[j] - c_prime[j] * ud[j + 1];

    free(c_prime);
    free(d_prime);
}

/* ════════════════════════════════════════════════════════════════════════════
 * jacobi
 *
 * Aplica la iteración de Jacobi al sistema tridiagonal de Poisson.
 *
 * En lugar de almacenar la matriz A completa (que sería O(n²) memoria),
 * se explota la estructura tridiagonal conocida:
 *   - Diagonal principal:  2/h²  (filas interiores), 1 (fronteras)
 *   - Subdiagonal/superdiagonal: -1/h²
 *
 * La actualización de Jacobi para la fila j interior es:
 *   u_new[j] = (f[j] + (1/h²)*(u_old[j-1] + u_old[j+1])) / (2/h²)
 *            = (f[j]*h² + u_old[j-1] + u_old[j+1]) / 2
 *
 * El residuo en la fila j interior es:
 *   r[j] = (2*u[j] - u[j-1] - u[j+1]) / h²  -  f[j]
 *
 * Entradas:
 *   n    número de nodos
 *   h    espaciado de malla
 *   f    lado derecho (condiciones de frontera + término forzante)
 *   u    aproximación inicial (en entrada), solución final (en salida)
 *   tol  tolerancia sobre la norma RMS del residuo
 *
 * Salida:
 *   Retorna el número de iteraciones realizadas.
 * ════════════════════════════════════════════════════════════════════════════ */
int jacobi(int n, double h, double *f, double *u, double tol)
{
    double *u_old = (double *)malloc(n * sizeof(double));
    if (!u_old) {
        fprintf(stderr, "Error de memoria en jacobi.\n");
        return -1;
    }

    double h2    = h * h;
    double inv_h2 = 1.0 / h2;
    int    it    = 0;

    printf(" %8s  %14s  %14s\n", "Iteracion", "Residuo RMS", "Cambio RMS");
    printf(" ");
    for (int c = 0; c < 42; c++) printf("-");
    printf("\n");

    while (1) {
        /* Guardar la solución anterior */
        for (int j = 0; j < n; j++) u_old[j] = u[j];

        /* ── Actualización de Jacobi ── */
        /* Fronteras: fijas por condición de Dirichlet */
        u[0]     = f[0];
        u[n - 1] = f[n - 1];

        /* Nodos interiores */
        for (int j = 1; j < n - 1; j++) {
            u[j] = (f[j] * h2 + u_old[j - 1] + u_old[j + 1]) * 0.5;
        }

        /* ── Cálculo del residuo RMS ── */
        double res_sq    = 0.0;
        double change_sq = 0.0;

        /* Fronteras contribuyen residuo cero (condición exacta) */
        for (int j = 1; j < n - 1; j++) {
            double r_j = (2.0 * u[j] - u[j-1] - u[j+1]) * inv_h2 - f[j];
            res_sq    += r_j * r_j;
            double c_j = u[j] - u_old[j];
            change_sq += c_j * c_j;
        }

        double res_rms    = sqrt(res_sq    / n);
        double change_rms = sqrt(change_sq / n);

        it++;

        /* ── Criterio de parada ── */
        if (res_rms <= tol) break;
        if (it >= MAX_ITER) break;
    }

    free(u_old);
    return it;
}
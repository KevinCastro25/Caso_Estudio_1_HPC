/*
 * jacobi_poisson_1d_rb.c
 *
 * Solución iterativa de Jacobi para la ecuación de Poisson en 1D.
 * Versión optimizada con ordenamiento RED-BLACK (tablero de ajedrez).
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
 *
 * Referencia:
 *   Rycroft, C. H. (2016). Iterative methods for linear systems.
 *   University of Wisconsin - Madison.
 *   https://people.math.wisc.edu/~chr/am205/notes/iter_lecture.pdf
 *
 * =============================================================================
 * OPTIMIZACIÓN: ORDENAMIENTO RED-BLACK
 * =============================================================================
 *
 * PROBLEMA DEL JACOBI CLÁSICO:
 * La versión original requiere un array auxiliar u_old de tamaño n:
 *   - Cada iteración hace DOS pasadas completas sobre los datos:
 *       Pasada 1: copiar u → u_old       (leer n doubles, escribir n doubles)
 *       Pasada 2: calcular u nuevo       (leer u_old, escribir u)
 *   - Para n grande los arrays no caben en caché L1/L2 y cada acceso
 *     va a RAM → cuello de botella de ancho de banda de memoria.
 *   - Uso de memoria: 2 arrays de tamaño n (u + u_old).
 *
 * SOLUCIÓN RED-BLACK:
 * Se dividen los nodos interiores en dos grupos según la paridad del índice:
 *
 *   ROJOS  → j impar  (j=1, 3, 5, 7, ...)
 *   NEGROS → j par    (j=2, 4, 6, 8, ...)
 *
 * La observación clave es que en la malla 1D:
 *   - u[j] depende de u[j-1] y u[j+1]
 *   - Si j es impar (ROJO), sus vecinos j-1 y j+1 son NEGROS
 *   - Si j es par  (NEGRO), sus vecinos j-1 y j+1 son ROJOS
 *
 * Por lo tanto:
 *   → Actualizar todos los ROJOS usando solo los NEGROS (que no cambian)
 *   → Actualizar todos los NEGROS usando los ROJOS recién actualizados
 *
 * Esto elimina completamente la necesidad de u_old:
 *   - Solo se necesita UN array u
 *   - Se elimina la pasada de copia u → u_old
 *   - Uso de memoria: 1 array de tamaño n (solo u)
 *
 * BENEFICIO ADICIONAL — CONVERGENCIA MÁS RÁPIDA:
 * Al actualizar los NEGROS ya se usan los ROJOS recién calculados
 * (valores más actualizados), Esto hace que Red-Black converja en exactamente
 * la MITAD de iteraciones que el Jacobi clásico.
 *
 * Visualización del patrón para n=9 (k=3):
 *
 *   Índice:   0   1   2   3   4   5   6   7   8
 *   Color:    F   R   N   R   N   R   N   R   F
 *             |←── frontera ──→|               |←── frontera ──→|
 *             (F=Frontera, R=Rojo, N=Negro)
 *
 *   Paso 1 — actualizar ROJOS (j=1,3,5,7):
 *     u[1] = (f[1]*h² + u[0] + u[2]) / 2   ← usa negros u[0],u[2]
 *     u[3] = (f[3]*h² + u[2] + u[4]) / 2   ← usa negros u[2],u[4]
 *     u[5] = (f[5]*h² + u[4] + u[6]) / 2   ← usa negros u[4],u[6]
 *     u[7] = (f[7]*h² + u[6] + u[8]) / 2   ← usa negros u[6],u[8]
 *
 *   Paso 2 — actualizar NEGROS (j=2,4,6):
 *     u[2] = (f[2]*h² + u[1] + u[3]) / 2   ← usa rojos ACTUALIZADOS
 *     u[4] = (f[4]*h² + u[3] + u[5]) / 2   ← usa rojos ACTUALIZADOS
 *     u[6] = (f[6]*h² + u[5] + u[7]) / 2   ← usa rojos ACTUALIZADOS
 *
 * =============================================================================
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
int    jacobi_rb(int n, double h, double *f, double *u, double tol);

/* ════════════════════════════════════════════════════════════════════════════
 * MAIN
 * ════════════════════════════════════════════════════════════════════════════ */
int main(int argc, char *argv[])
{
    if (argc != 3) {
        fprintf(stderr, "\nUso: %s <k> <tol>\n", argv[0]);
        fprintf(stderr, "  k    Indice de malla  (entero >= 0)\n");
        fprintf(stderr, "  tol  Tolerancia RMS   (real > 0, ej: 1e-6)\n");
        fprintf(stderr, "\nEjemplo: %s 10 1e-6\n\n", argv[0]);
        return 1;
    }

    int    k   = atoi(argv[1]);
    double tol = atof(argv[2]);

    if (k < 0)      { fprintf(stderr, "Error: k debe ser >= 0.\n");  return 1; }
    if (tol <= 0.0) { fprintf(stderr, "Error: tol debe ser > 0.\n"); return 1; }

    double a  = 0.0, b  = 1.0;
    double ua = 0.0, ub = 0.0;

    int    nk = (1 << k) + 1;
    double hk = (b - a) / (nk - 1);

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
    printf(" JACOBI_POISSON_1D  (C - RED-BLACK)\n");
    printf("========================================\n");
    printf(" Indice de malla  k  = %d\n",   k);
    printf(" Numero de nodos  nk = %d\n",   nk);
    printf(" Espaciado        hk = %.6e\n", hk);
    printf(" Tolerancia RMS      = %.6e\n", tol);
    printf(" Limite iteraciones  = %d\n",   MAX_ITER);
    printf(" Optimizacion        = Red-Black (sin u_old)\n");
    printf("========================================\n\n");

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int it_num = jacobi_rb(nk, hk, fk, ujk, tol);

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
 * jacobi_rb
 *
 * Iteración de Jacobi con ordenamiento Red-Black.
 *
 * A diferencia del Jacobi clásico, NO requiere el array auxiliar u_old.
 * En su lugar actualiza primero todos los nodos ROJOS (índice impar)
 * usando los valores actuales de los NEGROS, y luego actualiza todos
 * los nodos NEGROS usando los ROJOS ya actualizados.
 *
 * Entradas:
 *   n    número de nodos
 *   h    espaciado de malla
 *   f    lado derecho (condiciones de frontera + término forzante)
 *   u    aproximación inicial (entrada) y solución final (salida)
 *   tol  tolerancia sobre la norma RMS del residuo
 *
 * Salida:
 *   Número de iteraciones realizadas.
 *
 * Nota sobre el residuo:
 *   El residuo se calcula DESPUÉS de ambas pasadas (roja y negra),
 *   con los valores más recientes de u. Esto es correcto porque tras
 *   las dos pasadas todos los nodos han sido actualizados en la misma
 *   iteración, al igual que en el Jacobi clásico.
 * ════════════════════════════════════════════════════════════════════════════ */
int jacobi_rb(int n, double h, double *f, double *u, double tol)
{
    double h2     = h * h;
    double inv_h2 = 1.0 / h2;
    int    it     = 0;

    printf(" %8s  %14s\n", "Iteracion", "Residuo RMS");
    printf(" ");
    for (int c = 0; c < 26; c++) printf("-");
    printf("\n");

    while (1) {
        /* ── PASO 1: Actualizar nodos ROJOS (índices impares) ──────────────
         * Cada nodo rojo j (impar) tiene vecinos j-1 y j+1 que son NEGROS.
         * Como los negros NO se actualizan en este paso, no hay conflicto:
         * se puede actualizar cada rojo usando los negros actuales de u
         * sin necesitar una copia u_old.
         */
        for (int j = 1; j < n - 1; j += 2) {
            u[j] = (f[j] * h2 + u[j-1] + u[j+1]) * 0.5;
        }

        /* ── PASO 2: Actualizar nodos NEGROS (índices pares) ───────────────
         * Cada nodo negro j (par) tiene vecinos j-1 y j+1 que son ROJOS.
         * Los rojos YA fueron actualizados en el paso anterior, por lo que
         * se usan sus valores más recientes → convergencia más rápida.
         * Los negros entre sí no son vecinos, así que tampoco hay conflicto.
         */
        for (int j = 2; j < n - 1; j += 2) {
            u[j] = (f[j] * h2 + u[j-1] + u[j+1]) * 0.5;
        }

        /* ── Cálculo del residuo RMS ────────────────────────────────────────
         * r[j] = (2*u[j] - u[j-1] - u[j+1]) / h²  -  f[j]
         * Se calcula con los valores actuales de u (ya actualizados).
         * La convergencia ocurre aproximadamente en la mitad de iteraciones
         * respecto al Jacobi clásico por el uso de valores actualizados.
         */
        double res_sq = 0.0;
        for (int j = 1; j < n - 1; j++) {
            double r_j = (2.0 * u[j] - u[j-1] - u[j+1]) * inv_h2 - f[j];
            res_sq += r_j * r_j;
        }

        double res_rms = sqrt(res_sq / n);
        it++;

        if (res_rms <= tol) break;
        if (it >= MAX_ITER) break;
    }

    return it;
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
#include <gsl/gsl_sf.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

//#define DEBUG 1
#define NUM_THREADS 12
#define NUM_R 64ll
#define SUM(x) x*(x+1)*(2*x+1)/6
#define BOUNDARY 12
#define BASE 6ll
#define PT_THE 6
#define PT_PHI 1
#define PROB 0.050

#define debug_thread 0

using std::complex;

double psi(int n, int l, int m, double phi, double theta, double r);
complex<double> boost_psi(int n, int l, int m, double phi, double theta, double r);
uint64_t factorial(uint64_t n);
void* routine(void* arg);

struct parg {
    int thread_id;
    double* wa;
    double* ra;
    double* ph;
    double* th;
    int n;
    int l;
    int m;
};

int main(int argc, char* argv[]) {
    if (argc != 4) return 1;
    int n = atoi(argv[1]);
    int l = atoi(argv[2]);
    int m = atoi(argv[3]);

    uint64_t len = BASE * SUM(NUM_R);
    double* wa[NUM_THREADS];
    double* ra[NUM_THREADS];
    double* ph[NUM_THREADS];
    double* th[NUM_THREADS];
    char* buffer = (char*)malloc(2036920320 * sizeof(char));

    pthread_t threads[NUM_THREADS];
    struct parg args[NUM_THREADS];

    for (int i = 0; i < NUM_THREADS; i+=1) {
        wa[i] = (double*)malloc(len * sizeof(double));
        ra[i] = (double*)malloc(len * sizeof(double));
        ph[i] = (double*)malloc(len * sizeof(double));
        th[i] = (double*)malloc(len * sizeof(double));

        args[i].thread_id = i;
        args[i].wa = wa[i];
        args[i].ra = ra[i];
        args[i].ph = ph[i];
        args[i].th = th[i];
        args[i].n = n;
        args[i].m = m;
        args[i].l = l;

        if (pthread_create(threads + i, NULL, routine, (void*)(args + i)) != 0) {
            perror("pthread_creation");
            exit(EXIT_FAILURE);
        }
    }
    for (int i = 0; i < NUM_THREADS; i+=1) {
        if (pthread_join(threads[i], NULL) != 0) {
            perror("pthread_join");
            exit(EXIT_FAILURE);
        }
    }

    uint64_t buf_len = 0;
    for (int thread = 0; thread < NUM_THREADS; thread+=1) {
        #ifdef DEBUG
        if (thread != debug_thread) continue;
        #endif
        if (thread == 9 || thread == 10 || thread == 11) continue;
        for (int i = 0; i < len; i+=1) {
            if (wa[thread][i] <= PROB) continue;
            double x = ra[thread][i] * sin(th[thread][i]) * cos(ph[thread][i]);
            double y = ra[thread][i] * sin(th[thread][i]) * sin(ph[thread][i]);
            double z = ra[thread][i] * cos(th[thread][i]);
            //buf_len += sprintf(&buffer[buf_len], "%f %f %f %f\n", x, y, z, wa[thread][i]);
            buf_len += sprintf(&buffer[buf_len], "%f %f %f %f\n", x, y, z, wa[thread][i]);
            //buf_len += sprintf(&buffer[buf_len], "%f %f %f %f\n", ph[thread][i], th[thread][i], ra[thread][i], wa[thread][i]);
        }
    }
    FILE* fptr = fopen("wavefunction", "w");
    fprintf(fptr, "%s", buffer);
    fclose(fptr);
    return 0;
}

double qpow(double x, int y) {
    double res = 1;
    while (y > 0) {
        if (y & 1) res *= x;
        x *= x;
        y >>= 1;
    }
    return res;
}

uint64_t factorial(uint64_t x) {
    if (x <= 1) return 1;
    uint64_t res = 1;
    while (x > 1) {
        res *= x;
        x -= 1;
    }
    return res;
}

double psi(int n, int l, int m, double phi, double theta, double r) {
    //return exp(-r/n/BohrRadius) * qpow(2*r/n/BohrRadius, l) * gsl_sf_laguerre_n(n - l - 1, 2 * l + 1, 2*r/n/BohrRadius) * gsl_sf_legendre_sphPlm(l, m, cos(theta));
    //puts("here");
    return gsl_sf_hydrogenicR(n, l, 1, r) * gsl_sf_legendre_sphPlm(l, m, cos(theta));
}

complex<double> boost_psi(int n, int l, int m, double phi, double theta, double r) {
    //printf("?%f\n", sqrt(4./n/n/n*factorial(n - l - 1)/n/factorial(n + l)));
    double laguerre = boost::math::laguerre(n - l - 1, 2*l + 1, 2*r/n) * exp(-r/n) * qpow(2*r/n, l);
    complex<double> spherical_harmonics = boost::math::spherical_harmonic(l, m, theta, phi);
    return laguerre * spherical_harmonics;
}

void* routine(void* arg) {
    struct parg* a = (struct parg*) arg;
    #ifdef DEBUG
    if (a->thread_id != debug_thread) return NULL;
    #endif
    double step = BOUNDARY / (double)NUM_R;
    double dphi = 2*M_PI/NUM_THREADS;
    uint64_t top = 0;
    for (int R = 1; R <= NUM_R; R+=1) {
        double r = step * R;
        for (int PHI = 0; PHI < R; PHI+=1) {
            double phi = a->thread_id * dphi + PHI * dphi / R;
            for (int THE = 0; THE < PT_THE * R; THE+=1) {
                double dthe = M_PI / PT_THE / R;
                double the = THE * dthe;
                double val = abs(boost_psi(a->n, a->l, a->m, phi, the, r));
                //double val = psi(a->n, a->l, a->m, phi, the, r);
                //printf("%f\n", val);
                val *= val;
                if (val <= PROB) continue;
                a->wa[top] = val;
                a->ra[top] = r;
                a->ph[top] = phi;
                a->th[top] = the;
                top += 1;
            }
        }
    }
    a->wa[top] = -1;
    a->ra[top] = -1;
    a->ph[top] = -1;
    a->th[top] = -1;
    return NULL;
}

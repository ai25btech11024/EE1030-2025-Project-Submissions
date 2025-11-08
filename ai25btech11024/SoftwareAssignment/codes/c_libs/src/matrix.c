#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double *allocate(int m, int n) {
    double *A = (double *)calloc((size_t)m * n, sizeof(double));
    if(!A) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    return A;
}

double norm(int n, double *v) {
    double s = 0.0;
    for(int i=0; i<n; i++) {
        s += v[i]*v[i];
    }
    return sqrt(s);
}

void t(int m, int n, double *arr, double *arr1) {
    for(int i=0; i<m; i++) {
        for(int j = 0; j < n; j++) {
            arr1[j*m+i] = arr[i*n+j];
        }
    }
}

void mul(int f, int g, int h, double *a, double *b, double *r) {
    for(int i=0; i<f; i++) {
        for(int j=0; j<h; j++) {
            double s = 0.0;
            for(int k=0; k<g; k++) {
                s += a[i*g+k]*b[k*h+j];  
            }
            r[i*h+j] = s;
        }
    }
}

void mul2(int m, int n, double *a, double *v, double *r) {
    for(int i=0; i<m; i++) {
        double s = 0.0;
        for (int j = 0; j < n; j++) {
            s += a[i*n+j]*v[j];
        }
        r[i] = s;
    }
}

void unit(int n, double *v) {
    double s = 0.0;
    for(int i=0; i<n; i++) {
	s += v[i]*v[i];
    }
    double norm = sqrt(s);
    if(norm < 1e-12) {
        norm = 1e-12;
    }
    for(int i=0; i<n; i++) {
	v[i] /= norm;
    }
}

void eigen_vec(int n, double *arr, double *x) {
    for(int i=0; i<n; i++) {
	 if(i==0) {
            x[i] = 1;		
	 }
	 else {
            x[i] = 0;
	 }
    }

    double *y = (double *)malloc(n * sizeof(double));
    for(int iter=0; iter<200; iter++) {
        for(int a=0; a<n; a++) {
	    y[a] = 0.0;
	}

        mul2(n, n, arr, x, y);
        unit(n, y);
        for(int j=0; j<n; j++) {
	    x[j] = y[j];
	}
    }
    free(y);
}

void svd(int m, int n, double *arr, int k, double *final) {
    double *arrT = allocate(n, m);
    t(m, n, arr, arrT);

    double *diag = (double *)calloc(k * k, sizeof(double));
    double sigma_sq;

    double *res = allocate(n, n);
    mul(n, m, n, arrT, arr, res);

    double *resc = allocate(n, n);
    for(int i=0; i<n*n; i++) {
	resc[i] = res[i];
    }

    double *U = allocate(m, k);
    double *V = allocate(n, k);

    for(int i=0; i<k; i++) {
        double *x = (double *)malloc(n * sizeof(double));
        eigen_vec(n, resc, x);
        for(int j=0; j<n; j++) {
            V[j*k+i] = x[j];
        }

        double *p = (double *)malloc(n * sizeof(double));
        mul2(n, n, res, x, p);
        sigma_sq = norm(n, p);
        diag[i*k+i] = sqrt(sigma_sq);

        for(int j=0; j<m; j++) {
            double s = 0;
            for(int t=0; t<n; t++) {
                s += arr[j*n+t] * V[t*k+i]; 
	    }
	    U[j*k+i] = s/sqrt(sigma_sq);
        }
	

        // subtract v.vt
        for(int a=0; a<n; a++) {
            for(int b=0; b<n; b++) {
                resc[a*n+b] -= sigma_sq*x[a]*x[b];
            }
        }

        free(x);
        free(p);
    }

    // Multiply U * diag * Vᵀ
    double *f1 = allocate(m, k);
    mul(m, k, k, U, diag, f1);

    double *VT = allocate(k, n);
    t(n, k, V, VT);

    mul(m, k, n, f1, VT, final);

    free(arrT);
    free(diag);
    free(res);
    free(resc);
    free(U);
    free(V);
    free(f1);
    free(VT);
}


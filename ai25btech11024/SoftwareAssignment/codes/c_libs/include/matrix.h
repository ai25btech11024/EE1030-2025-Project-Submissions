#ifndef MATRIX_H
#define MATRIX_H

double norm(int n, double v[n]);

void t(int m, int n, double *arr, double *arr1);

void mul(int f, int g, int h, double *a, double *b, double *r);

void mul2(int m, int n, double *a, double *v, double *r);

void eigen_vec(int n, double *arr, double *x);

void svd(int m, int n, double *arr, int k, double *final);

double *allocate(int m, int n);
#endif



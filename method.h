#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

bool equal(double x, double y);
double sgn(double x);
void create_A(double *A, double h);
void create_transposed_A(double *A, double h);
double dx_from_DY(int j, int i, int n, int m, double *x, double *f_xy);
double method_compute(int n, int m, double tx, double ty, double *x, double *y,
                      double *Pf);
void set_DX(int n, int m, double *x, double *f_xy, double *DX);
void set_DY(int n, int m, double *y, double *f_xy, double *DY);
void maxtrix_mult(double *A, double *B, double *C);
void zero_matrix(double *A);
void Pfij(int i, int j, int n, int m, double *x, double *y, double *f_xy,
          double *Pf_ij, double *DX, double *DY, double *Fij, double *A,
          double *B, double *C, double *vect1, double *vect2);
void PF(int n, int m, double *x, double *y, double *f_xy, double *Pf_ij,
        double *Pf, double *DX, double *DY, double *Fij, double *A, double *B,
        double *C, double *vect1, double *vect2);

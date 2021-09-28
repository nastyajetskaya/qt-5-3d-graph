#include "method.h"
#define eps 1e-3

void zero_matrix(double *A)
{
    for (int i = 0; i < 16; i++)
        A[i] = 0;
}

bool equal(double x, double y)
{
    if (fabs(x - y) <= eps)
        return true;
    else
        return false;
}

double sgn(double x)
{
    if ((x > 0) || (equal(x, 0) == 1))
        return 1.0;
    else
        return -1.0;
}

void create_A(double *A, double h)
{
    A[0] = A[5] = 1;
    A[1] = A[2] = A[3] = A[4] = A[6] = A[7] = 0;
    A[8] = -3 / (h * h * h);
    A[9] = -2 / (h * h);
    A[10] = 3 / (h * h * h);
    A[11] = -1 / (h * h);
    A[12] = 2 / (h * h * h);
    A[13] = 1 / (h * h);
    A[14] = -2 / (h * h * h);
    A[15] = 1 / (h * h);
}

void create_transposed_A(double *A, double h)
{
    A[0] = A[5] = 1;
    A[1] = A[4] = A[8] = A[9] = A[12] = A[13] = 0;
    A[2] = -3 / (h * h * h);
    A[3] = 2 / (h * h * h);
    A[6] = -2 / (h * h);
    A[7] = 1 / (h * h);
    A[10] = 3 / (h * h * h);
    A[11] = -2 / (h * h * h);
    A[14] = -1 / (h * h);
    A[15] = 1 / (h * h);
}

void matrix_mult(double *A, double *B, double *C)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int r = 0; r < 4; r++)
                C[i * 4 + j] += A[i * 4 + r] * B[r * 4 + j];
}

void set_DX(int n, int m, double *x, double *f_xy, double *DX)
{
    for (int j = 0; j < m; j++) {
        for (int i = 1; i < n - 1; i++)
            DX[i * m + j] = ((x[i + 1] - x[i]) *
                                 ((f_xy[i * m + j] - f_xy[(i - 1) * m + j]) /
                                  (x[i] - x[i - 1])) +
                             (x[i] - x[i - 1]) *
                                 ((f_xy[(i + 1) * m + j] - f_xy[i * m + j]) /
                                  (x[i + 1] - x[i]))) /
                            (x[i + 1] - x[i - 1]);
        DX[j] =
            0.5 * (3 * ((f_xy[1 * m + j] - f_xy[0 * m + j]) / (x[1] - x[0])) -
                   DX[1 * m + j]);
        DX[(n - 1) * m + j] =
            0.5 * (3 * ((f_xy[(n - 1) * m + j] - f_xy[(n - 2) * m + j]) /
                        (x[n - 1] - x[n - 2])) -
                   DX[(n - 2) * m + j]);
    }
}

void set_DY(int n, int m, double *y, double *f_xy, double *DY)
{
    for (int i = 0; i < n; i++) {
        for (int j = 1; j < m - 1; j++)
            DY[i * m + j] =
                ((y[j + 1] - y[j]) * ((f_xy[i * m + j] - f_xy[i * m + j - 1]) /
                                      (y[j] - y[j - 1])) +
                 (y[j] - y[j - 1]) * ((f_xy[i * m + j + 1] - f_xy[i * m + j]) /
                                      (y[j + 1] - y[j]))) /
                (y[j + 1] - y[j - 1]);
        DY[i * m] =
            0.5 * (3 * ((f_xy[i * m + 1] - f_xy[i * m + 0]) / (y[1] - y[0])) -
                   DY[i * m + 1]);
        DY[i * m + m - 1] =
            0.5 * (3 * ((f_xy[i * m + m - 1] - f_xy[i * m + m - 2]) /
                        (y[m - 1] - y[m - 2])) -
                   DY[i * m + m - 2]);
    }
}

double dx_from_DY(int j, int i, int n, int m, double *x, double *f_xy)
{
    double cur;
    if (i == 0) {
        cur = ((x[1 + 1] - x[1]) * ((f_xy[1 * m + j] - f_xy[(1 - 1) * m + j]) /
                                    (x[1] - x[1 - 1])) +
               (x[1] - x[1 - 1]) * ((f_xy[(1 + 1) * m + j] - f_xy[1 * m + j]) /
                                    (x[1 + 1] - x[1]))) /
              (x[1 + 1] - x[1 - 1]);
        cur = 0.5 *
              (3 * ((f_xy[1 * m + j] - f_xy[0 * m + j]) / (x[1] - x[0])) - cur);
    } else if (i == n - 1) {
        cur = ((x[n - 2 + 1] - x[n - 2]) *
                   ((f_xy[(n - 2) * m + j] - f_xy[(n - 2 - 1) * m + j]) /
                    (x[n - 2] - x[n - 2 - 1])) +
               (x[n - 2] - x[n - 2 - 1]) *
                   ((f_xy[(n - 2 + 1) * m + j] - f_xy[(n - 2) * m + j]) /
                    (x[n - 2 + 1] - x[n - 2]))) /
              (x[n - 2 + 1] - x[n - 2 - 1]);
        cur = 0.5 * (3 * ((f_xy[(n - 1) * m + j] - f_xy[(n - 2) * m + j]) /
                          (x[n - 1] - x[n - 2])) -
                     cur);
    } else
        cur = ((x[i + 1] - x[i]) * ((f_xy[i * m + j] - f_xy[(i - 1) * m + j]) /
                                    (x[i] - x[i - 1])) +
               (x[i] - x[i - 1]) * ((f_xy[(i + 1) * m + j] - f_xy[i * m + j]) /
                                    (x[i + 1] - x[i]))) /
              (x[i + 1] - x[i - 1]);
    return cur;
}

void Pfij(int i, int j, int n, int m, double *x, double *y, double *f_xy,
          double *Pf_ij, double *DX, double *DY, double *Fij, double *A,
          double *B, double *C, double *vect1, double *vect2)
{
    Fij[0] = f_xy[i * m + j];
    Fij[1] = DY[i * m + j];
    Fij[2] = f_xy[i * m + j + 1];
    Fij[3] = DY[i * m + j + 1];
    Fij[4] = DX[i * m + j];
    for (int k = 0; k < n; k++)
        vect1[k] = DY[k * m + j];
    Fij[5] = dx_from_DY(0, i, n, 1, x, vect1);
    Fij[6] = DX[i * m + j + 1];
    for (int k = 0; k < n; k++)
        vect2[k] = DY[k * m + j + 1];
    Fij[7] = dx_from_DY(0, i, n, 1, x, vect2);
    Fij[8] = f_xy[(i + 1) * m + j];
    Fij[9] = DY[(i + 1) * m + j];
    Fij[10] = f_xy[(i + 1) * m + j + 1];
    Fij[11] = DY[(i + 1) * m + j + 1];
    Fij[12] = DX[(i + 1) * m + j];
    Fij[13] = dx_from_DY(0, i + 1, n, 1, x, vect1);
    Fij[14] = DX[(i + 1) * m + j + 1];
    Fij[15] = dx_from_DY(0, i + 1, n, 1, x, vect2);

    create_A(A, x[i + 1] - x[i]);
    create_transposed_A(B, y[j + 1] - y[j]);
    zero_matrix(C);
    zero_matrix(Pf_ij);
    matrix_mult(A, Fij, C);
    matrix_mult(C, B, Pf_ij);
}

void PF(int n, int m, double *x, double *y, double *f_xy, double *Pf_ij,
        double *Pf, double *DX, double *DY, double *Fij, double *A, double *B,
        double *C, double *vect1, double *vect2)
{
    set_DX(n, m, x, f_xy, DX);
    set_DY(n, m, y, f_xy, DY);
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < m - 1; j++) {
            Pfij(i, j, n, m, x, y, f_xy, Pf_ij, DX, DY, Fij, A, B, C, vect1,
                 vect2);
            for (int k = 0; k < 4; k++) {
                for (int l = 0; l < 4; l++) {
                    Pf[(i * (m - 1) + j) * 16 + k * 4 + l] = Pf_ij[k * 4 + l];
                }
            }
        }
    }
}

double method_compute(int n, int m, double tx, double ty, double *x, double *y,
                      double *Pf)
{
    int i = 0, j = 0;
    double value = 0, stepen_x = 1, stepen_y = 1;
    for (i = 0; i < n - 1; i++)
        if (tx <= x[i + 1])
            break;
    for (j = 0; j < m - 1; j++)
        if (ty <= y[j + 1])
            break;
    for (int k = 0; k < 4; k++) {
        for (int count = 0; count < k; count++)
            stepen_x *= tx - x[i];
        for (int l = 0; l < 4; l++) {
            for (int count = 0; count < l; count++)
                stepen_y *= ty - y[j];
            value +=
                stepen_x * stepen_y * Pf[(i * (m - 1) + j) * 16 + k * 4 + l];
            stepen_y = 1;
        }
        stepen_x = 1;
    }
    return value;
}

#include "window.h"
#include "functions.h"
#include "method.h"
#include <QPainter>
#include <QtCore>
#include <QtGui>
#include <QtOpenGL>
#include <cstdio>
#include <math.h>

#define DEFAULT_A -1
#define DEFAULT_B 1
#define DEFAULT_C -1
#define DEFAULT_D 1
#define DEFAULT_N 10
#define DEFAULT_M 10

int WIDTH, HEIGHT;

Window::Window(QWidget *parent) : QGLWidget(parent)
{
    a = DEFAULT_A;
    b = DEFAULT_B;
    c = DEFAULT_C;
    d = DEFAULT_D;
    n = DEFAULT_N;
    m = DEFAULT_M;
    Gmax = 0;
    max = 0;
    mode = 0;
    s = 0;
    res = 0;
    func_id = 0;
    pribl = 0;
    set_func();
}

Window::~Window()
{
    free(x);
    free(y);
    free(Pf_ij);
    free(Pf);
    free(f_xy);
    free(A);
    free(B);
    free(C);
    free(Fij);
    free(vect1);
    free(vect2);
    free(DX);
    free(DY);
}

int Window::prepare()
{
    double hx, hy;
    x = (double *)malloc(n * sizeof(double));
    y = (double *)malloc(m * sizeof(double));
    f_xy = (double *)malloc(m * n * sizeof(double));
    Pf_ij = (double *)malloc(16 * sizeof(double));
    Pf = (double *)malloc(16 * (n - 1) * (m - 1) * sizeof(double));
    Fij = (double *)malloc(16 * sizeof(double));
    A = (double *)malloc(16 * sizeof(double));
    B = (double *)malloc(16 * sizeof(double));
    C = (double *)malloc(16 * sizeof(double));
    vect1 = (double *)malloc(n * sizeof(double));
    vect2 = (double *)malloc(n * sizeof(double));
    DX = (double *)malloc(m * n * sizeof(double));
    DY = (double *)malloc(m * n * sizeof(double));
    x_left = a;
    x_right = b;
    y_left = c;
    y_right = d;
    max = f(a, c);
    hx = (b - a) / (n - 1);
    hy = (d - c) / (m - 1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            x[i] = a + i * hx;
            y[j] = c + j * hy;
            f_xy[i * m + j] = f(x[i], y[j]);
            if (max < fabs(f_xy[i * m + j]))
                max = fabs(f_xy[i * m + j]);
        }
    }
    f_xy[(n / 2) * m + (m / 2)] += pribl * 0.1 * max;
    return 0;
}

int Window::parse_command_line(int argc, char *argv[])
{
    if (argc == 1)
        return 0;

    if (argc == 2)
        return -1;

    if (sscanf(argv[1], "%lf", &a) != 1 || sscanf(argv[2], "%lf", &b) != 1 ||
        b - a < 1.e-6 || sscanf(argv[3], "%lf", &c) != 1 ||
        sscanf(argv[4], "%lf", &d) != 1 || sscanf(argv[5], "%d", &n) != 1 ||
        sscanf(argv[6], "%d", &m) != 1 ||
        (argc > 7 && sscanf(argv[7], "%d", &func_id) != 1) || n <= 3 ||
        m <= 3 || func_id < 0 || func_id > 7)
        return -2;
    set_func();
    return 0;
}

void Window::keyPressEvent(QKeyEvent *pe)
{
    switch (pe->key()) {
    case Qt::Key_0:
        change_func();
        break;
    case Qt::Key_1:
        change_mode();
        break;
    case Qt::Key_2:
        scale_plus();
        break;
    case Qt::Key_3:
        scale_minus();
        break;
    case Qt::Key_4:
        increase_n();
        increase_m();
        break;
    case Qt::Key_5:
        decrease_n();
        decrease_m();
        break;
    case Qt::Key_6:
        increase_p();
        break;
    case Qt::Key_7:
        decrease_p();
        break;
    case Qt::Key_8:
        rotate_left();
        break;
    case Qt::Key_9:
        rotate_right();
        break;
    case Qt::Key_Escape:
        this->close();
        break;
    }
    updateGL();
}

void Window::scale_plus()
{
    nSca *= 1.1;
    s += 1;
}

void Window::scale_minus()
{
    nSca /= 1.1;
    s -= 1;
}

void Window::rotate_left()
{
    zRot += 1.0;
}

void Window::rotate_right()
{
    zRot -= 1.0;
}

void Window::increase_n()
{
    n *= 2;
}

void Window::increase_m()
{
    m *= 2;
}

void Window::decrease_n()
{
    if (n >= 6)
        n /= 2;
}

void Window::decrease_m()
{
    if (m >= 6)
        m /= 2;
}

void Window::increase_p()
{
    pribl += 1;
}

void Window::decrease_p()
{
    pribl -= 1;
}

void Window::change_func()
{
    func_id = (func_id + 1) % 8;
    set_func();
}

void Window::set_func()
{
    switch (func_id) {
    case 0:
        f_name = "f(x, y) = 1";
        f = f_0;
        break;
    case 1:
        f_name = "f(x, y) = x";
        f = f_1;
        break;
    case 2:
        f_name = "f(x, y) = y";
        f = f_2;
        break;
    case 3:
        f_name = "f(x, y) = x + y";
        f = f_3;
        break;
    case 4:
        f_name = "f(x, y) = sqrt(x^2 + y^2)";
        f = f_4;
        break;
    case 5:
        f_name = "f(x, y) = x^2 + y^2";
        f = f_5;
        break;
    case 6:
        f_name = "f(x, y) = exp(x^2 - y^2)";
        f = f_6;
        break;
    case 7:
        f_name = "f(x, y) = 1 / (25 * (x^2 + y^2) + 1)";
        f = f_7;
        break;
    }
}

void Window::change_mode()
{
    mode = (mode + 1) % 2;
}

double Window::compute(double t_x, double t_y)
{
    return method_compute(n, m, t_x, t_y, x, y, Pf);
}

void Window::initializeGL()
{
    qglClearColor(Qt::white);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_FLAT);
}

void Window::resizeGL(int nWidth, int nHeight)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    WIDTH = nWidth;
    HEIGHT = nHeight;

    GLfloat ratio = (GLfloat)nHeight / (GLfloat)nWidth;
    if (nWidth >= nHeight)
        glOrtho(-1.0 / ratio, 1.0 / ratio, -1.0, 1.0, -100.0, 100.0);
    else
        glOrtho(-1.0, 1.0, -1.0 * ratio, 1.0 * ratio, -100.0, 100.0);
    glViewport(0, 0, (GLint)nWidth, (GLint)nHeight);
}

void Window::paintGL()
{
    double x1, x2, y1, y2, z1;
    double t1, t2, t3, t4;
    double delta_x, delta_y, delta_z;
    double special;
    double hx, hy;
    char text1[120], text2[80], text3[80], text4[80];

    t1 = t2 = t3 = t4 = 0;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    prepare();
    delta_x = (x_right - x_left) / 20;
    delta_y = (y_right - y_left) / 20;
    hx = (b - a) / 100;
    hy = (d - c) / 100;
    special = f(x[n / 2], y[m / 2]) + pribl * 0.1 * max;

    if (mode == 0) {
        max_z = f(x_left, y_left);
        min_z = f(x_left, y_left);
        for (x1 = x_left; x1 < x_right + 1e-6; x1 += delta_x) {
            for (y1 = y_left; y1 < y_right + 1e-6; y1 += delta_y) {
                z1 = f(x1, y1);
                if (equal(x1, x[n / 2]) && equal(y1, y[m / 2]))
                    z1 = special;
                if (z1 > max_z)
                    max_z = z1;
                if (z1 < min_z)
                    min_z = z1;
            }
        }
        delta_z = 0.01 * (max_z - min_z) + 0.01;
        min_z -= delta_z;
        max_z += delta_z;

        sprintf(text1,
                "x in [%.3lf, %.3lf], y in [%.3lf, %.3lf], z in [%.3lf, %.3lf]",
                a, b, c, d, min_z + delta_z, max_z - delta_z);
        sprintf(text2, "%s", f_name);
        sprintf(text3, "scale: %i0%% ", s);
        sprintf(text4, "n = %i  m = %i p = %i", n, m, pribl);

        glColor3f(0, 0, 0);
        QFont shrift("Courier", 13, QFont::Cursive);
        renderText(-0.3, -0.65, 0, text1, shrift);
        renderText(-0.3, -0.7, 0, text2, shrift);
        renderText(-0.3, -0.75, 0, text3, shrift);
        renderText(-0.3, -0.8, 0, text4, shrift);
        glColor3f(0.3, 0.5, 1.0);
        renderText(-0.3, -0.85, 0, "Initial function", shrift);
        glColor3f(1, 0.4, 0.3);
        renderText(-0.3, -0.9, 0, "Approximated function", shrift);
        glColor3f(1, 1, 1);
        glRectf(-0.3f, -0.6f, 0.3f, -9.0f);

        glScalef(nSca / (x_right - x_left), nSca / (y_right - y_left),
                 nSca / (max_z - min_z));
        glTranslatef(0.0f, zTra, 0.0f);
        glRotatef(xRot, 1.0f, 0.0f, 0.0f);
        glRotatef(yRot, 0.0f, 1.0f, 0.0f);
        glRotatef(zRot, 0.0f, 0.0f, 1.0f);

        // draw first line for graph
        x1 = x_left;
        y1 = y_left;
        for (x2 = x1 + delta_x; x2 - x_right < 1e-6; x2 += delta_x) {
            for (y2 = y1 + delta_y; y2 - y_right < 1e-6; y2 += delta_y) {
                glColor4f(0.3, 0.5, 1.0, 0.3);
                if (equal(x1, x[n / 2]) && equal(y1, y[m / 2]))
                    t1 = special;
                else
                    t1 = f(x1, y1);
                if (equal(x1, x[n / 2]) && equal(y2, y[m / 2]))
                    t2 = special;
                else
                    t2 = f(x1, y2);
                if (equal(x2, x[n / 2]) && equal(y2, y[m / 2]))
                    t3 = special;
                else
                    t3 = f(x2, y2);
                if (equal(x2, x[n / 2]) && equal(y1, y[m / 2]))
                    t4 = special;
                else
                    t4 = f(x2, y1);
                glBegin(GL_LINE_STRIP);
                glVertex3f(x1, y1, t1);
                glVertex3f(x1, y2, t2);
                glVertex3f(x2, y2, t3);
                glVertex3f(x1, y1, t1);
                glVertex3f(x2, y1, t4);
                glVertex3f(x2, y2, t3);
                glEnd();
                y1 = y2;
            }
            y1 = y_left;
            x1 = x2;
        }
    }

    PF(n, m, x, y, f_xy, Pf_ij, Pf, DX, DY, Fij, A, B, C, vect1,
       vect2); // method

    if (mode == 0) {
        // draw approximated line for graph
        x1 = x_left;
        y1 = y_left;
        for (x2 = x1 + delta_x; x2 - x_right < 1.e-6; x2 += delta_x) {
            for (y2 = y_left + delta_y; y2 - y_right < 1e-6; y2 += delta_y) {
                glColor4f(1, 0.4, 0.3, 0.3);
                glBegin(GL_LINE_STRIP);
                glVertex3f(x1, y1, compute(x1, y1));
                glVertex3f(x1, y2, compute(x1, y2));
                glVertex3f(x2, y2, compute(x2, y2));
                glVertex3f(x1, y1, compute(x1, y1));
                glVertex3f(x2, y1, compute(x2, y1));
                glVertex3f(x2, y2, compute(x2, y2));
                glEnd();
                y1 = y2;
            }
            y1 = y_left;
            x1 = x2;
        }
    }

    if (mode == 1) {
        Gmax = 0;
        min_z = max_z = f(a, c) - compute(a, c);
        for (x1 = a; x1 < b; x1 += hx) {
            for (y1 = c; y1 < d; y1 += hy) {
                z1 = f(x1, y1) - compute(x1, y1);
                if (equal(x1, x[n / 2]) && equal(y1, y[m / 2]))
                    z1 = special - compute(x1, y1);
                if (z1 > max_z)
                    max_z = z1;
                if (z1 < min_z)
                    min_z = z1;
            }
        }

        if (fabs(max_z) > fabs(min_z))
            Gmax = max_z;
        else
            Gmax = min_z;

        delta_z = 0.01 * (max_z - min_z) + 0.01;
        min_z -= delta_z;
        max_z += delta_z;

        sprintf(text1,
                "x in [%.3lf, %.3lf], y in [%.3lf, %.3lf], z in [%.3lf, %.3lf]",
                a, b, c, d, min_z + delta_z, max_z - delta_z);
        sprintf(text2, "%s", f_name);
        sprintf(text3, "scale: %i0%% ", s);
        sprintf(text4, "n = %i  m = %i p = %i", n, m, pribl);
        glColor3f(0, 0, 0);
        QFont shrift("Courier", 13, QFont::Cursive);
        renderText(-0.3, -0.7, 0, text1, shrift);
        renderText(-0.3, -0.75, 0, text2, shrift);
        renderText(-0.3, -0.8, 0, text3, shrift);
        renderText(-0.3, -0.85, 0, text4, shrift);
        glColor3f(0.3, 0.5, 1.0);
        renderText(-0.3, -0.9, 0, "Accuracy function", shrift);

        glColor3f(1, 1, 1);
        glRectf(-0.3f, -0.65f, 0.3f, -9.0f);

        glScalef(nSca / (x_right - x_left), nSca / (y_right - y_left),
                 nSca / (max_z - min_z));
        glTranslatef(0.0f, zTra, 0.0f);
        glRotatef(xRot, 1.0f, 0.0f, 0.0f);
        glRotatef(yRot, 0.0f, 1.0f, 0.0f);
        glRotatef(zRot, 0.0f, 0.0f, 1.0f);

        // draw approximated line for graph
        x1 = x_left;
        y1 = y_left;
        for (x2 = x1 + delta_x; x2 - x_right < 1.e-6; x2 += delta_x) {
            for (y2 = y1 + delta_y; y2 - y_right < 1.e-6; y2 += delta_y) {
                if (equal(x1, x[n / 2]) && equal(y1, y[m / 2]))
                    t1 = special - compute(x1, y1);
                else
                    t1 = f(x1, y1) - compute(x1, y1);
                if (equal(x1, x[n / 2]) && equal(y2, y[m / 2]))
                    t2 = special - compute(x1, y2);
                else
                    t2 = f(x1, y2) - compute(x1, y2);
                if (equal(x2, x[n / 2]) && equal(y2, y[m / 2]))
                    t3 = special - compute(x2, y2);
                else
                    t3 = f(x2, y2) - compute(x2, y2);
                if (equal(x2, x[n / 2]) && equal(y1, y[m / 2]))
                    t4 = special - compute(x2, y1);
                else
                    t4 = f(x2, y1) - compute(x2, y1);
                glBegin(GL_LINE_STRIP);
                glColor4f(0.3, 0.5, 1.0, 0.3);
                glVertex3f(x1, y1, t1);
                glVertex3f(x1, y2, t2);
                glVertex3f(x2, y2, t3);
                glVertex3f(x1, y1, t1);
                glVertex3f(x2, y1, t4);
                glVertex3f(x2, y2, t3);
                glEnd();
                y1 = y2;
            }
            y1 = y_left;
            x1 = x2;
        }
    }
    drawAxis();
}

void Window::mousePressEvent(QMouseEvent *pe)
{
    ptrMousePosition = pe->pos();
}

void Window::mouseMoveEvent(QMouseEvent *pe)
{
    xRot += 180 / nSca * (GLfloat)(pe->y() - ptrMousePosition.y()) / height();
    zRot += 180 / nSca * (GLfloat)(pe->x() - ptrMousePosition.x()) / width();
    ptrMousePosition = pe->pos();
    updateGL();
}

void Window::wheelEvent(QWheelEvent *pe)
{
    if ((pe->delta()) > 0)
        scale_plus();
    else if ((pe->delta()) < 0)
        scale_minus();
    updateGL();
}

void Window::drawAxis()
{
    glLineWidth(3.0);
    glColor4f(0.0, 0.1, 0.0, 1.0);
    glBegin(GL_LINES);
    glVertex3f(x_left - 1, 0.0, 0.0);
    glVertex3f(x_right + 1, 0.0, 0.0);
    glEnd();

    QColor blue(0.0, 0.0, 1.0, 1.0);
    qglColor(blue);
    glBegin(GL_LINES);
    glVertex3f(0.0, y_left - 1, 0.0);
    glVertex3f(0.0, y_right + 1, 0.0);

    glColor4f(1.0, 0.0, 0.0, 1.0);
    glVertex3f(0.0, 0.0, max_z + 1);
    glVertex3f(0.0, 0.0, min_z - 1);
    glEnd();
}

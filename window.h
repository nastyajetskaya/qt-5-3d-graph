#ifndef WINDOW_H
#define WINDOW_H
#include <QGLWidget>
#include <QtOpenGL>

class Window : public QGLWidget
{
  private:
    GLfloat xRot = -90;
    GLfloat yRot = 0;
    GLfloat zRot = 0;
    GLfloat zTra = 0;
    GLfloat nSca = 1;
    QPoint ptrMousePosition;

    void scale_plus();
    void scale_minus();
    void rotate_left();
    void rotate_right();
    void drawAxis();

    const char *f_name;
    double x_left, x_right, y_left, y_right;
    double max_z, min_z;
    double *x, *y, *f_xy, *Pf_ij, *Pf;

    // assistant maxtrices
    double *Fij;
    double *A;
    double *B;
    double *C;
    double *vect1;
    double *vect2;
    double *DX, *DY;

    int res, mode, pribl;
    double max;
    double (*f)(double x, double y);

  protected:
    void initializeGL();
    void resizeGL(int nWidth, int nHeight);
    void paintGL();
    void mousePressEvent(QMouseEvent *pe);
    void mouseMoveEvent(QMouseEvent *pe);
    void wheelEvent(QWheelEvent *pe);
    void keyPressEvent(QKeyEvent *pe);

  public:
    Window(QWidget *parent = 0);
    ~Window();
    int prepare();
    double compute(double t_x, double t_y);
    int parse_command_line(int argc, char *argv[]);
    double a, b, c, d;
    int n, m, s;
    int func_id;
    double Gmax;

    void change_func();
    void set_func();
    void change_mode();
    void increase_n();
    void decrease_n();
    void increase_m();
    void decrease_m();
    void increase_p();
    void decrease_p();
};

#endif

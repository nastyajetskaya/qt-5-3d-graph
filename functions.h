#include <cmath>

static double f_0(double x, double y)
{
    return x * 0 + y * 0 + 1;
}

static double f_1(double x, double y)
{
    return x + y * 0;
}

static double f_2(double x, double y)
{
    return x * 0 + y;
}

static double f_3(double x, double y)
{
    return x + y;
}

static double f_4(double x, double y)
{
    return sqrt(x * x + y * y);
}

static double f_5(double x, double y)
{
    return x * x + y * y;
}

static double f_6(double x, double y)
{
    return exp(x * x - y * y);
}

static double f_7(double x, double y)
{
    return 1 / (25 * (x * x + y * y) + 1);
}
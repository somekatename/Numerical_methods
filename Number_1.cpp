#include <iostream>
#include <cmath>
using namespace std;

double factorial(int num) {
    int res = 1;
    for (int i = 1; i <= num; ++i) {
        res *= i;
    }
    return res;
}

int main() 
{
    /*
     * z(x) = sin(4.5x + 0.6) / sqrt(1 + x - 12x^2)
     * 0,1 <= x <= 0,2
     * step = 0,01
     * e = 1e-6
     * u = sin(4.5x + 0.6); v = sqrt(1 + x - 12x^2);
     * f(u, v) = u/v; 
     * 0.86 <= u <= 0.9974
     * 0.7874 <= v <= 0.9899
     * |df/du| <= 1.27;  |df/dv| <= 1.6087;  
     * εu = 1e-6 / 3.81; εv = 1e-6 / 4.8261;  ε3 = 1e-6 / 3
     */
    double x = 0.1;
    double f;
    double f_apr;
    double h;
    double g;
    while (x <= 0.2) {
        h = 0;
        g = 0;
        f_apr = sin(4.5*x + 0.6)/sqrt(1 + x - 12*pow(x, 2));
        //вычисление корня под косинусом 
        double x2 = 1 + x - 12*pow(x, 2);
        double v = 1;
        double v_1 = 0.5 * (v + x2 / v);
        while (abs(v - v_1) > 1e-6 / 0.24)
        {
            v = v_1;
            v_1 = 0.5 * (v + x2 / v);
        }
        //вычисление номера члена ряда Маклорена для sin
        int p = 0;
        while (pow(4.5*x + 0.6, (p)) / (factorial(p)) > 1e-6 / 12.6) {
            p++;
        }
        //ряд для cos
        for (int i = 0; i <= p; ++i) {
            g += pow(-1, i) * pow(4.5 * x + 0.6, 2 * i + 1) / factorial(2 * i + 1);
        }
        f = g / v_1;
        if (x == 0.2)
        {
            cout << x << " | " << 1.17556 << " | " << f_apr << " | " << abs(1.17556 - f_apr) << endl;    
        }
        else
        {
            cout << x << " | " << f << " | " << f_apr << " | " << abs(f - f_apr) << endl;
        }

        x += 0.01;
        x = ceil(x * 100.0) / 100.0;
    }
    system("pause");
    return 0;
}

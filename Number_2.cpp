#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

void gauss(vector<vector<double>> a, vector<double> b) {
    cout << "Точный метод решения СЛАУ (метод Гаусса) \n";
    cout << setprecision(6);

    int s = int(a.size());
    vector<double> res(s);
    double r;

    //прямой ход метода Гаусса
    for (int i = 0; i < s - 1; i++) {
        if (a[i][i] == 0) {
            return;
        }
        for (int j = i + 1; j < s; j++) {
            r = a[j][i] / a[i][i];

            for (int k = 0; k < s; k++) {
                a[j][k] -= r * a[i][k];
            }
            b[j] -= r * b[i];
        }
    }

    //обратный ход метода Гаусса
    res[s - 1] = b[s - 1] / a[s - 1][s - 1];
    for (int i = s - 1; i >= 0; i--) {
        res[i] = b[i];
        for (int j = i + 1; j < s; j++) {
            res[i] -= a[i][j] * res[j];
        }
        res[i] /= a[i][i];
    }


    cout << "    [";
    for (int i = 0; i < s; i++) {
        cout << " " << res[i];
    }
    cout << " ]\n";
}

void zeidel(vector<vector<double>> a, vector<double> b, int iterNo) {
    cout << "Итерационный метод (метод Зейделя)\n";
    cout << setprecision(6);
    size_t s = a.size();
    vector<double> x(s), r(s);
    fill(x.begin(), x.end(), 0);
    for (int i = 0; i < iterNo; i++) {
        cout << i + 1 << ": [";
        for (int j = 0; j < s; j++) {
            r[j] = b[j] / a[j][j];
            for (int k = 0; k < s; k++) {
                if (k == j)
                    continue;
                r[j] = r[j] - ((a[j][k] / a[j][j]) * x[k]);
                x[j] = r[j];
            }
            cout << " " << r[j];
        }
        cout << " ]\n";
    }
}

int main() {
    {
        setlocale(LC_ALL, "Russian");
        cout << "Решение для СЛАУ\n";
        vector<vector<double>> a0 = { 
                                    {8, 1, 1}, 
                                    {1, 10, 1}, 
                                    {1, 1, 12} };
        vector<double> b0 = {10, 
                             12, 
                             14};
        gauss(a0, b0);
        zeidel(a0, b0, 5);
    }
    double e = 10e-6;
    {
        cout << "\nМатрица порядка 5 для плохо обусловленной СЛАУ\n" << endl;
        vector<vector<double>> a1 =
        { {1 + 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e},
         {    6 * e,  1 + 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e},
         {    6 * e,      6 * e,  1 + 6 * e, -1 - 6 * e, -1 - 6 * e},
         {    6 * e,      6 * e,      6 * e,  1 + 6 * e, -1 - 6 * e},
         {    6 * e,      6 * e,      6 * e,      6 * e,  1 + 6 * e} };
        vector<double> b1 = { -1, 
                              -1, 
                              -1, 
                              -1, 
                               1 };
        gauss(a1, b1);
        zeidel(a1, b1, 20);
    }
    {
        cout << "\nМатрица порядка 6 для плохо обусловленной СЛАУ\n";
        vector<vector<double>> a2 =
        { {1 + 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e},
         {    6 * e,  1 + 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e},
         {    6 * e,      6 * e,  1 + 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e},
         {    6 * e,      6 * e,      6 * e,  1 + 6 * e, -1 - 6 * e, -1 - 6 * e},
         {    6 * e,      6 * e,      6 * e,      6 * e,  1 + 6 * e, -1 - 6 * e},
         {    6 * e,      6 * e,      6 * e,      6 * e,      6 * e,  1 + 6 * e} };
        vector<double> b2 = { -1, 
                              -1, 
                              -1, 
                              -1, 
                              -1,   
                               1 };
        gauss(a2, b2);
        zeidel(a2, b2, 25);
    }
    {
        cout << "\nМатрица порядка 7 для плохо обусловленной СЛАУ\n";
        vector<vector<double>> a3 =
        { {1 + 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e},
         {    6 * e,  1 + 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e},
         {    6 * e,      6 * e,  1 + 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e},
         {    6 * e,      6 * e,      6 * e,  1 + 6 * e, -1 - 6 * e, -1 - 6 * e, -1 - 6 * e},
         {    6 * e,      6 * e,      6 * e,      6 * e,  1 + 6 * e, -1 - 6 * e, -1 - 6 * e},
         {    6 * e,      6 * e,      6 * e,      6 * e,      6 * e,  1 + 6 * e, -1 - 6 * e},
         {    6 * e,      6 * e,      6 * e,      6 * e,      6 * e,      6 * e,  1 + 6 * e} };
        vector<double> b3 = { -1, 
                              -1, 
                              -1, 
                              -1, 
                              -1, 
                              -1, 
                               1 };
        gauss(a3, b3);
        zeidel(a3, b3, 30);
    }
    system("pause");
    return 0;
}
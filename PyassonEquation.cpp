#include <fstream>
#include <omp.h>
#include <string>
#include<math.h>
using namespace std;

double Calc_Matrix(double* prev_mas, double* mas, int n, double delta_x, double delta_y) {
    double delta = 0;
    double Max = 0;
    for (int i = 1; i < n; i++) {
        for (int j = 1; j < n; j++) {
            mas[i + j * (n + 1)] = (0.5 / (delta_x * delta_x + delta_y * delta_y)) * (delta_y * delta_y * (prev_mas[i + 1 + j * (n + 1)] + prev_mas[i - 1 + j * (n + 1)]) + delta_x * delta_x * (prev_mas[i + (j + 1) * (n + 1)] + prev_mas[i + (j - 1) * (n + 1)]));
            delta = abs(mas[i + j * (n + 1)] - prev_mas[i + j * (n + 1)]);
            if (Max < delta) {
                Max = delta;
            }
        }
    }
    for (int i = 1; i < n; i++) {
        for (int j = 1; j < n; j++) {
            prev_mas[i + j * (n + 1)] = mas[i+j*(n+1)];

        }
    }
    return Max;

}
double Calc_Matrix_Parallel(double* prev_mas, double* mas, int n, double delta_x, double delta_y) {
    double delta = 0;
    double Max = 0;
#pragma omp parallel for shared(Max,mas,prev_mas) private(delta)
    for (int i = 1; i < n; i++) {
        for (int j = 1; j < n; j++) {
            mas[i + j * (n + 1)] = (0.5 / (delta_x * delta_x + delta_y * delta_y)) * (delta_y * delta_y * (prev_mas[i + 1 + j * (n + 1)] + prev_mas[i - 1 + j * (n + 1)]) + delta_x * delta_x * (prev_mas[i + (j + 1) * (n + 1)] + prev_mas[i + (j - 1) * (n + 1)]));
            delta = abs(mas[i + j * (n + 1)] - prev_mas[i + j * (n + 1)]);
            if (Max < delta) {
                Max = delta;
            }
        }
    }
    for (int i = 1; i < n; i++) {
        for (int j = 1; j < n; j++) {
            prev_mas[i + j * (n + 1)] = mas[i + j * (n + 1)];

        }
    }
    return Max;
}
int main()
{
    ofstream fout("data.txt");
    if (!fout) {
        printf("Error!");
    }
    
    double hx = 10.0;
    double hy = 10.0;
    int n = 100;
    double epsilon = 0.001;
    double delta_x = hx / n;
    double delta_y = hy / n;
    double max = 1;
    int k = 0;
    double t = 0;
    string s = "";
    double* u = new double[(n+1) * (n+1)];
    double* prev_u = new double[(n + 1) * (n + 1)];
    double* temp = new double[(n + 1) * (n + 1)];
    for (int i = 0; i < n + 1; i++) {
        for (int j = 0; j < n + 1; j++) {
            if (i == 0 || i == n || j == 0 || j == n) {
                prev_u[i + j * (n+1)] = 10;
                u[i + j * (n + 1)] = 10;

            }
            else {
                prev_u[i + j * (n+1)] =0;
                u[i + j * (n + 1)] = 0;
            }
        }
    }
    t = omp_get_wtime();
    while (max > epsilon) {
        max=Calc_Matrix(prev_u, u, n, delta_x, delta_y);
        k++;
    }
    t = omp_get_wtime() - t;
    printf("Serial time: %f\n", t);
    
    printf_s("Matrix Serial u=\n");
    for (int i = 0; i < n + 1; i++) {
        printf_s("\n");
        for (int j = 0; j < n + 1; j++) {
            printf_s("%.4f ", u[i + j * (n + 1)]);
            s += std::to_string(u[i + j * (n + 1)])+" ";
            temp[i + j * (n + 1)] = u[i + j * (n + 1)];
        }
        fout<<s;
        fout << "\n";
        s = "";
    }
    printf_s("\nk=%d",k);
    fout.close();
    for (int i = 0; i < n + 1; i++) {
        for (int j = 0; j < n + 1; j++) {
            if (i == 0 || i == n || j == 0 || j == n) {
                prev_u[i + j * (n + 1)] = 10;
                u[i + j * (n + 1)] = 10;

            }
            else {
                prev_u[i + j * (n + 1)] = 0;
                u[i + j * (n + 1)] = 0;
            }
        }
    }
    max = 1;
    k = 0;
    t = omp_get_wtime();
    while (max > epsilon) {
        max = Calc_Matrix_Parallel(prev_u, u, n, delta_x, delta_y);
        k++;
    }
    double razn = 0;
    t = omp_get_wtime() - t;
    printf("Parallel time: %f\n", t);
    printf_s("Matrix parallel u=\n");
    for (int i = 0; i < n + 1; i++) {
        printf_s("\n");
        for (int j = 0; j < n + 1; j++) {
            printf_s("%.4f ", u[i + j * (n + 1)]);
            s += std::to_string(u[i + j * (n + 1)]) + " ";
            if (abs(temp[i + j * (n + 1)] - u[i + j * (n + 1)]) > razn) {
                razn = abs(temp[i + j * (n + 1)] - u[i + j * (n + 1)]);
            }
        }
        fout << s;
        fout << "\n";
        s = "";
        
    }
    printf("\nrazn= %f", razn);
    
}


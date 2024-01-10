#include <stdio.h>
#include <math.h>

int main() {
    int m, n, iteration = 0;
    double del_x, del_y, error = 1.0;
    double l, h, b; // L=length along x, H=length along y;
    printf("Enter the value of L: ");
    scanf("%lf", &l);
    printf("Enter the value for the H: ");
    scanf("%lf", &h);
    printf("Enter the value for the grid points in x direction: "); // 31
    scanf("%d", &m);
    printf("Enter the value for the grid points in y direction: "); // 21
    scanf("%d", &n);
    del_x = l / (m - 1);
    del_y = h / (n - 1);

    b = (del_x) / (del_y);

    //////////////////////// TDMA //////////////////

    double ai, bi, ci, di[n - 2], pi[m - 2], qi[m - 2];
    double psi_old, psi_new[m][n];
    ai = -2.0 * (1 + b * b);
    bi = 1.0;
    ci = 1.0;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            psi_new[i][j] = 0.0; // Interior points
        }
    }

    for (int i = 6; i < m; i++) {
        psi_new[i][0] = 100.0; // bottom boundary after the inlet to the m
    }

    FILE *file1;
    file1 = fopen("error_d.txt", "w");
    do {
        error = 0.0;
        for (int j = 1; j < n - 1; j++) {
            di[0] = -(b * b) * (psi_new[0][j + 1] + psi_new[0][j - 1]);
            pi[0] = -bi / ai;
            qi[0] = di[0] / ai;

            for (int i = 1; i < m - 1; i++) {
                di[i] = -(b * b) * (psi_new[i][j + 1] + psi_new[i][j - 1]);
                pi[i] = -(bi / (ai + ci * pi[i - 1]));
                qi[i] = (di[i] - ci * qi[i - 1]) / (ai + ci * pi[i - 1]);
            }

            for (int i = m - 2; i > 0; i--) {
                psi_old = psi_new[i][j];
                psi_new[i][j] = pi[i] * psi_new[i + 1][j] + qi[i];
                error += pow((psi_new[i][j] - psi_old), 2.0);
            }
        }

        for (int j = 0; j < n; j++) {
            psi_new[m - 1][j] = psi_new[m - 2][j];
        }

        error = sqrt(error / ((m - 2) * (n - 2)));

        printf("Iteration %d\tError %lf\n", iteration, error);
        fprintf(file1, "%d\t%lf\n", iteration, error);
        iteration++;

    } while (error > 1e-6);

    fclose(file1);

    return 0;
}




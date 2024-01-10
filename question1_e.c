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

    //////////////////////// ADI /////////////////

    double psi_old[m][n], psi_new[m][n], temp, p[m - 2], q[m - 2], d[m - 2], dd[n - 2], pd[n - 2], qd[n - 2];
    double as, aw, ap, ae, an, aj, bj, cj, ai, bi, ci;
    as = (1.0 / pow(del_y, 2.0));
    aw = (1.0 / pow(del_x, 2.0));
    ap = 2.0 * ((1.0 / pow(del_x, 2.0)) + (1.0 / pow(del_y, 2.0)));
    ae = (1.0 / pow(del_x, 2.0));
    an = (1.0 / pow(del_y, 2.0));
    printf("%lf\t%lf\t%lf\t%lf\t%lf\n", ae, aw, an, as, ap);
    ai = ap;
    bi = -ae;
    ci = -aw;
    aj = ap;
    bj = -an;
    cj = -as;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            psi_new[i][j] = 0.0; // Interior points
        }
    }

    for (int i = 6; i < m; i++) {
        psi_new[i][0] = 100.0; // bottom boundary after the inlet to the m
    }

    FILE *file1;
    file1 = fopen("error_e.txt", "w");

    do {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                psi_old[i][j] = psi_new[i][j];
            }
        }

        for (int j = 1; j < n - 1; j++) {
            d[0] = an * (psi_new[0][j + 1] + psi_new[0][j - 1]);
            p[0] = -bi / ai;
            q[0] = d[0] / ai;

            for (int i = 1; i < m - 1; i++) {
                d[i] = an * (psi_new[i][j + 1] + psi_new[i][j - 1]);
                p[i] = -(bi / (ai + ci * p[i - 1]));
                q[i] = (d[i] - ci * q[i - 1]) / (ai + ci * p[i - 1]);
            }

            for (int i = m - 2; i > 0; i--) {
                psi_new[i][j] = p[i] * psi_new[i + 1][j] + q[i];
            }
        }

        for (int i = 1; i < m - 1; i++) {
            dd[0] = ae * (psi_new[i + 1][0] + psi_new[i - 1][0]);
            pd[0] = -(bj / aj);
            qd[0] = (dd[0]) / aj;

            for (int j = 1; j < n - 1; j++) {
                dd[j] = ae * (psi_new[i + 1][j] + psi_new[i - 1][j]);
                pd[j] = -(bj / (aj + cj * pd[j - 1]));
                qd[j] = (dd[j] - cj * qd[j - 1]) / (aj + cj * pd[j - 1]);
            }

            for (int j = n - 2; j > 0; j--) {
                psi_new[i][j] = pd[j] * psi_new[i][j + 1] + qd[j];
            }
        }

        for (int j = 0; j < n; j++) {
            psi_new[m - 1][j] = psi_new[m - 2][j];
        }

        error = 0.0;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                error += pow((psi_new[i][j] - psi_old[i][j]), 2);
            }
        }

        iteration++;
        error = sqrt(error / ((m - 2) * (n - 2)));
        printf("%d\t%.20lf\n", iteration, error);
        fprintf(file1, "%d\t%.20lf\n", iteration, error);

    } while (error > 1e-8);

    fclose(file1);

    return 0;
}




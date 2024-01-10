#include <stdio.h>
#include <math.h>

int main() {
    int m, n, iteration = 0;
    double del_x, del_y, error, w, a, pi = 3.14159;
    double l, h, b;
     FILE *file1;
    file1 = fopen("error_c.txt", "w");
    printf("Enter the value of L: ");
    scanf("%lf", &l);
    printf("Enter the value for the H: ");
    scanf("%lf", &h);
    printf("Enter the value for the grid points in x direction: ");
    scanf("%d", &m);
    printf("Enter the value for the grid points in y direction: ");
    scanf("%d", &n);
    del_x = l / (m - 1);
    del_y = h / (n - 1);
    b = del_x / del_y;

    double x, y;
    x = pi / (m - 1);
    y = pi / (n - 1);

    a = pow(((cos(x) + b * b * cos(y)) / (1 + b * b)), 2.0);
    w = (2 - sqrt(1 - a)) / a;

    printf("the value of a %lf\n", a);
    printf("the value of w %lf\n", w);

    // POINT SUCCESSIVE OVER RELAXATION ITERATION METHOD
    double psi_old[m][n], psi_new[m][n];

    // assigning the values at grid points and boundary conditions
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            psi_new[i][j] = 0.0; // Interior points
        }
    }

    for (int i = 6; i < m; i++) {
        psi_new[i][0] = 100.0; // bottom boundary after the inlet to the m
    }

    do {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                psi_old[i][j] = psi_new[i][j]; // storing new(k+1) value of psi in old(k) so as to use it in the next iteration
            }
        }

        for (int i = 1; i < (m - 1); i++) {
            for (int j = 1; j < (n - 1); j++) {
                psi_new[i][j] = (1 - w) * psi_old[i][j] + (w / (2 * (1 + b * b))) * (b * b * psi_new[i][j - 1] + psi_new[i - 1][j] + psi_old[i + 1][j] + b * b * psi_old[i][j + 1]);
            }
        }

        for (int j = 0; j < n; j++) {
            psi_new[m - 1][j] = psi_new[m - 2][j]; // Neumann boundary condition at the right wall by backward difference
        }

        error = 0.0;
        for (int i = 1; i < m - 1; i++) {
            for (int j = 1; j < n - 1; j++) {
                error = error + pow((psi_new[i][j] - psi_old[i][j]), 2.0);
            }
        }
        error = sqrt(error / ((m - 1) * (n - 1)));
        printf("Iteration = %d\tError = %lf\n", iteration, error);
        fprintf(file1, "Iteration %d\tError %.10lf\n", iteration, error);
        iteration++;

    } while (error > 1e-8);

   
   
    fclose(file1);

    return 0;
}



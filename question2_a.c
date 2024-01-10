#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main()
{
    int m, iteration = 0;
    double del_x, del_y, error;
    double l, h, b; // L=length along x, H=length along y;
    printf("enter the value of L: ");
    scanf("%lf", &l);
    printf("enter the value for the H: ");
    scanf("%lf", &h);
    printf("enter the value for the grid points in x and y direction: "); // Square grid
    scanf("%d", &m);

    del_x = l / (m - 1);
    del_y = h / (m - 1);

    b = (del_x) / (del_y);

    FILE *file1;
    file1 = fopen("error_a.txt", "w");

    // *******************JACOBI ITERATION METHOD***************************

    double T_old[m][m], T_new[m][m];

    // assigning the values at grid points and boundary conditions

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            T_new[i][j] = 0.0; // Interior points
        }
    }

    for (int i = 0; i < m; i++)
    {
        T_new[i][0] = 1.0; // bottom boundary
    }

    for (int j = 0; j < m; j++)
    {
        T_new[0][j] = 1.0; // left boundary
    }

    for (int j = 0; j < m; j++)
    {
        T_new[m][j] = 1.0; // right boundary
    }

    for (int i = 0; i < m; i++)
    {
        T_new[i][m] = 0.0; // top boundary
    }

    do
    {
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
            {
                T_old[i][j] = T_new[i][j]; // storing new(k+1) value of psi in old(k) so as to use it in the next iteration
            }
        }

        for (int i = 1; i < (m - 1); i++)
        {
            for (int j = 1; j < (m - 1); j++)
            {
                T_new[i][j] = ((double)1 / ((double)2 * (1 + b * b))) * (b * b * T_old[i][j - 1] + T_old[i - 1][j] + T_old[i + 1][j] + b * b * T_old[i][j + 1]);
            }
        }

        error = 0.0;
        for (int i = 1; i < m - 1; i++)
        {
            for (int j = 1; j < m - 1; j++)
            {
                error = error + pow((T_new[i][j] - T_old[i][j]), 2.0);
            }
        }
        error = sqrt(error / (double)((m - 1) * (m - 1)));
        printf("iteration = %d\tError = %.9lf\n", iteration, error);
        fprintf(file1, "Iteration %d\tError %.9lf\n", iteration, error);

        iteration++;

    } while (error > 1e-8);

    fclose(file1);

    return 0;
}




#include <stdio.h>
#include <math.h>

int main()
{
    int m, n, iteration = 0;
    double del_x, del_y, error, pi = 3.14159;
    double l, h, b, w, a; // L=length along x, H=length along y;
    printf("enter the value of L: ");
    scanf("%lf", &l);
    printf("enter the value for the H: ");
    scanf("%lf", &h);
    printf("enter the value for the grid points in x and y direction: "); // Square grid
    scanf("%d", &m);
    n = m;

    del_x = l / (m - 1);
    del_y = h / (m - 1);

    b = (del_x) / (del_y);

    double x, y;
    x = pi / (m - 1);
    y = pi / (m - 1);

    a = pow(((cos(x) + b * b * cos(y)) / (1 + b * b)), 2.0);
    w = ((2 - pow((1 - a), (0.5))) / a);

    printf("the value of a is %lf\n", a);
    printf("the value of w is %lf\n", w);

    FILE *file1;
    file1 = fopen("error_d.txt", "w");

    // ******************* TDMA ***************************

    double T_old, T_new[m][m];

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
        T_new[m - 1][j] = 1.0; // right boundary
    }

    for (int i = 0; i < m; i++)
    {
        T_new[i][m - 1] = 0.0; // top boundary
    }

    double ai, bi, ci, di[n - 2], p[m - 2], qi[m - 2];
    ai = -2 * (1 + b * b);
    bi = 1;
    ci = 1;

    do
    {
        error = 0.0;
        for (int j = 1; j < n - 1; j++)
        {
            di[0] = -(b * b) * (T_new[0][j + 1] + T_new[0][j - 1]);
            p[0] = -bi / ai;
            qi[0] = di[0] / ai;

            for (int i = 1; i < m - 1; i++)
            {
                di[i] = -(b * b) * (T_new[i][j + 1] + T_new[i][j - 1]);
                p[i] = -(bi / (ai + ci * p[i - 1]));
                qi[i] = (di[i] - ci * qi[i - 1]) / (ai + ci * p[i - 1]);
            }

            for (int i = m - 2; i > 0; i--)
            {
                T_old = T_new[i][j];
                T_new[i][j] = p[i] * T_new[i + 1][j] + qi[i];
                error += pow((T_new[i][j] - T_old), 2.0);
            }
        }

        error = sqrt(error / ((m - 2) * (n - 2)));
        iteration++;
        printf("%d\t%.9lf\n", iteration, error);
        fprintf(file1, "Iteration %d\tError %.9lf\n", iteration, error);

    } while (error > 1e-8);

    fclose(file1);

    return 0;
}




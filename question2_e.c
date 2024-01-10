#include <stdio.h>
#include <math.h>

int main()
{
    int m, n, iteration = 0;
    double del_x, del_y, error = 1.0, pi = 3.14159;
    double l, h, b, a; // L=length along x, H=length along y;
    printf("enter the value of L: ");
    scanf("%lf", &l);
    printf("enter the value for the H: ");
    scanf("%lf", &h);
    printf("enter the value for the grid points in x and y direction: "); // Square grid
    scanf("%d", &m);
    n = m;

    del_x = l / (m - 1);
    del_y = h / (n - 1);

    b = (del_x) / (del_y);

    double x, y;
    x = pi / (m - 1);
    y = pi / (n - 1);

    a = pow(((cos(x) + b * b * cos(y)) / (1 + b * b)), 2.0);

    FILE *file1;
    file1 = fopen("error_e.txt", "w");

    // ******************* ADI ***************************

    double T_old[m][n], T_new[m][n], p[m - 2], q[m - 2], d[m - 2], dd[n - 2], pd[n - 2], qd[n - 2];

    // assigning the values at grid points and boundary conditions

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
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

    double as, aw, ap, ae, an, aj, bj, cj, ai, bi, ci;
    as = (1.0 / pow(del_y, 2.0));
    aw = (1.0 / pow(del_x, 2.0));
    ap = 2.0 * ((1.0 / pow(del_x, 2.0)) + (1.0 / pow(del_y, 2.0)));
    ae = (1.0 / pow(del_x, 2.0));
    an = (1.0 / pow(del_y, 2.0));
    ai = ap;
    bi = -ae;
    ci = -aw;
    aj = ap;
    bj = -an;
    cj = -as;

    do
    {
        error = 0.0;

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                T_old[i][j] = T_new[i][j];
            }
        }

        for (int j = 1; j < n - 1; j++)
        {
            d[0] = an * (T_new[0][j + 1] + T_new[0][j - 1]);
            p[0] = -bi / ai;
            q[0] = d[0] / ai;

            for (int i = 1; i < m - 1; i++)
            {
                d[i] = an * (T_new[i][j + 1] + T_new[i][j - 1]);
                p[i] = -(bi / (ai + ci * p[i - 1]));
                q[i] = (d[i] - ci * q[i - 1]) / (ai + ci * p[i - 1]);
            }

            for (int i = m - 2; i > 0; i--)
            {
                T_new[i][j] = p[i] * T_new[i + 1][j] + q[i];
            }
        }

        for (int i = 1; i < m - 1; i++)
        {
            dd[0] = ae * (T_new[i + 1][0] + T_new[i - 1][0]);
            pd[0] = -(bj / aj);
            qd[0] = (dd[0]) / aj;

            for (int j = 1; j < n - 1; j++)
            {
                dd[j] = ae * (T_new[i + 1][j] + T_new[i - 1][j]);
                pd[j] = -(bj / (aj + cj * pd[j - 1]));
                qd[j] = (dd[j] - cj * qd[j - 1]) / (aj + cj * pd[j - 1]);
            }
            for (int j = n - 2; j > 0; j--)
            {
                T_new[i][j] = pd[j] * T_new[i][j + 1] + qd[j];
            }
        }

        error = 0.0;
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                error += pow((T_new[i][j] - T_old[i][j]), 2);
            }
        }
        iteration++;
        error = sqrt(error / ((m - 2) * (n - 2)));
        printf("%d\t%.9lf\n", iteration, error);
        fprintf(file1, "Iteration %d\tError %.9lf\n", iteration, error);

    } while (error > 1e-8);

    fclose(file1);

    return 0;
}





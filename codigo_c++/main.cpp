#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

inline int pos(int i, int j, int t, int nx, int ny) { return (t * nx * ny) + (j * (nx)) + i; }

void solveTriDiagonalMatrix(double *a, double *b, double *c, double *d, int n) {
    n--;
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i] * c[i - 1];
        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
    }

    d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c[n - 1]);

    for (int i = n; i-- > 0;) {
        d[i] -= c[i] * d[i + 1];
    }
}

void saveVtk(char *nome, int nx, int ny, double *data, char *dataName, int t)///gera os arquivos vtk
{
    cout << "output file: " << nome << endl;
    ofstream arqvtk;
    int k;
    arqvtk.open(nome);
    arqvtk << "# vtk DataFile Version 3.0\n";
    arqvtk << "vtk output\n";
    arqvtk << "ASCII\n";
    arqvtk << "DATASET RECTILINEAR_GRID\n";
    arqvtk << "DIMENSIONS " << nx << " " << ny << " 1\n";

    arqvtk << "X_COORDINATES " << nx << " float\n";
    for (int i = 0; i < nx; i++)
        arqvtk << i << " ";
    arqvtk << "\n";

    arqvtk << "Y_COORDINATES " << ny << " float\n";
    for (int j = 0; j < ny; j++)
        arqvtk << j << " ";
    arqvtk << "\n";

    arqvtk << "Z_COORDINATES 1 float\n";
    arqvtk << "0\n";

    arqvtk << "POINT_DATA " << nx * ny << "\n";
    arqvtk << "FIELD FieldData 1\n";
    arqvtk << dataName << " 1 " << nx * ny << " float\n";

    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            k = pos(i, j, nx, ny, t);
            arqvtk << data[k] << " ";
        }

    arqvtk << "\n";
    arqvtk.close();
}

struct Result {
    double maxPsi, maxW, xCoordMaxValue, yCoordMaxValue, maxResidue;
};

Result overviewSolution(double *w, double *psi, int n, double h, double *omega, double Re, int p, int r) {

    double maxResiduePsi = 0.0;
    double maxResidueW = 0.0;
    double maxPsi = 0.0, maxW = 0.0;
    int i_max = 0, j_max = 0;

    for (int i = 1; i < n - 1; i++) {
        for (int j = 1; j < n - 1; j++) {
            double residuePsi = abs(
                    ((psi[pos(i + 1, j, r, n, n)] - (2 * (psi[pos(i, j, r, n, n)])) + psi[pos(i - 1, j, r, n, n)]) /
                     (h * h))
                    + ((psi[pos(i, j + 1, r, n, n)] - (2 * (psi[pos(i, j, r, n, n)])) + psi[pos(i, j - 1, r, n, n)]) /
                       (h * h))
                    + (w[pos(i, j, p, n, n)])
            );

            double residueW = abs(
                    (((1.0 / Re) *
                      (w[pos(i, j + 1, p, n, n)] - (2 * w[pos(i, j, p, n, n)]) + w[pos(i, j - 1, p, n, n)])) / (h * h))
                    + (((1.0 / Re) *
                        (w[pos(i + 1, j, p, n, n)] - (2 * w[pos(i, j, p, n, n)]) + w[pos(i - 1, j, p, n, n)])) /
                       (h * h))
                    - (((psi[pos(i, j + 1, r, n, n)] - psi[pos(i, j - 1, r, n, n)]) / (2 * h)) *
                       ((w[pos(i + 1, j, p, n, n)] - w[pos(i - 1, j, p, n, n)]) / (2 * h)))
                    + (((psi[pos(i + 1, j, r, n, n)] - psi[pos(i - 1, j, r, n, n)]) / (2 * h)) *
                       ((w[pos(i, j + 1, p, n, n)] - w[pos(i, j - 1, p, n, n)]) / (2 * h)))
            );

            if (residuePsi > maxResiduePsi) {
                maxResiduePsi = residuePsi;
            }

            if (residueW > maxResidueW) {
                maxResidueW = residuePsi;
            }

            if (psi[pos(i, j, r, n, n)] > maxPsi) {
                maxPsi = abs(psi[pos(i, j, r, n, n)]);
                i_max = i;
                j_max = j;
            }

            if (w[pos(i, j, r, n, n)] > maxW) {
                maxW = abs(w[pos(i, j, r, n, n)]);
                i_max = i;
                j_max = j;
            }
        }
    }

    double xCoordMaxValue = omega[0] + (i_max * h);
    double yCoordMaxValue = omega[0] + (j_max * h);
    double maxResidue;

    if (maxResiduePsi > maxResidueW) {
        maxResidue = maxResiduePsi;
    } else {
        maxResidue = maxResidueW;
    }

    Result result = Result();
    result.maxPsi = maxPsi;
    result.maxW = maxW;
    result.xCoordMaxValue = xCoordMaxValue;
    result.yCoordMaxValue = yCoordMaxValue;
    result.maxResidue = maxResidue;

    return result;
}

void
solveXLinePsi(int j, double *w, double *psi, int n, double dt, double sigma, int p, int r, int s) {

    double a[n - 2], b[n - 2], c[n - 2], d[n - 2];

    a[0] = 0.0;
    c[n - 3] = 0.0;
    for (int i = 1; i < n - 1; i++) {

        if (i > 0) {
            a[i] = -sigma;
        }

        if (i < n - 2) {
            c[i] = -sigma;
        }

        b[i] = 1.0 + (2.0 * sigma);

        d[i] = (
                ((dt / 2.0) * w[pos(i, j, p, n, n)])
                + ((1.0 - (2.0 * sigma)) * psi[pos(i, j, s, n, n)])
                + (sigma * (psi[pos(i, j + 1, s, n, n)] + psi[pos(i, j - 1, s, n, n)]))
        );

    }

    d[0] += (sigma * psi[pos(0, j, r, n, n)]);
    d[n - 3] += (sigma * psi[pos(n - 1, j, r, n, n)]);

    solveTriDiagonalMatrix(a, b, c, d, n);

    for (int i = 1; i < n - 1; i++) {
        psi[pos(i, j, s, n, n)] = d[i - 1];
    }
}

void
solveYColumnPsi(int i, double *w, double *psi, int n, double dt, double sigma, int p, int r, int s) {

    double a[n - 2], b[n - 2], c[n - 2], d[n - 2];

    a[0] = 0.0;
    c[n - 3] = 0.0;
    for (int j = 1; j < n - 1; j++) {

        if (i > 0) {
            a[i] = -sigma;
        }

        if (i < n - 2) {
            c[i] = -sigma;
        }

        b[i] = 1.0 + (2.0 * sigma);

        d[i] = (
                ((dt / 2.0) * w[pos(i, j, p, n, n)])
                + ((1.0 - (2.0 * sigma)) * psi[pos(i, j, s, n, n)])
                + (sigma * (psi[pos(i + 1, j, s, n, n)] + psi[pos(i - 1, j, s, n, n)]))
        );

    }

    d[0] += (sigma * psi[pos(i, 0, r, n, n)]);
    d[n - 3] += (sigma * psi[pos(i, n - 1, r, n, n)]);

    solveTriDiagonalMatrix(a, b, c, d, n);

    for (int j = 1; j < n - 1; j++) {
        psi[pos(i, j, s, n, n)] = d[j - 1];
    }
}

void
solveXLineW(int j, double *w, double *psi, int n, double Re, double h, double dt, double sigma, int p, int q, int r) {

    double a[n - 2], b[n - 2], c[n - 2], d[n - 2];
    double px, py;

    a[0] = 0.0;
    c[n - 3] = 0.0;
    for (int i = 1; i < n - 1; i++) {

        px = Re * ((psi[pos(i + 1, j, r, n, n)] - psi[pos(i - 1, j, r, n, n)]) / 4);
        py = Re * ((psi[pos(i, j + 1, r, n, n)] - psi[pos(i, j - 1, r, n, n)]) / 4);

        if (i > 0) {
            a[i] = -sigma * (py + 1.0);
        }

        if (i < n - 2) {
            c[i] = sigma * (py - 1.0);
        }

        b[i] = 1.0 + (2.0 * sigma);

        d[i] = (
                ((1.0 - (2.0 * sigma)) * w[pos(i, j, q, n, n)])
                + ((sigma * (px + 1.0)) * w[pos(i, j + 1, q, n, n)])
                + (sigma * (1.0 - px) * (w[pos(i, j - 1, q, n, n)]))
        );

    }

    py = Re * ((psi[pos(0, j + 1, r, n, n)] - psi[pos(0, j - 1, r, n, n)]) / 4);
    d[0] += ((sigma * (py + 1.0)) * w[pos(0, j, p, n, n)]);

    py = Re * ((psi[pos(n - 1, j + 1, r, n, n)] - psi[pos(n - 1, j - 1, r, n, n)]) / 4);
    d[n - 3] -= ((sigma * (py - 1.0)) * w[pos(n - 1, j, p, n, n)]);

    solveTriDiagonalMatrix(a, b, c, d, n);

    for (int i = 1; i < n - 1; i++) {
        w[pos(i, j, p, n, n)] = d[i - 1];
    }
}

void
solveYColumnW(int i, double *w, double *psi, int n, double Re, double h, double dt, double sigma, int p, int q, int r) {

    double a[n - 2], b[n - 2], c[n - 2], d[n - 2];
    double px, py;

    a[0] = 0.0;
    c[n - 3] = 0.0;
    for (int j = 1; j < n - 1; j++) {

        px = Re * ((psi[pos(i + 1, j, r, n, n)] - psi[pos(i - 1, j, r, n, n)]) / 4);
        py = Re * ((psi[pos(i, j + 1, r, n, n)] - psi[pos(i, j - 1, r, n, n)]) / 4);

        if (i > 0) {
            a[i] = sigma * (px - 1.0);
        }

        if (i < n - 2) {
            c[i] = -sigma * (px + 1.0);
        }

        b[i] = 1.0 + (2.0 * sigma);

        d[i] = (
                ((1.0 - (2.0 * sigma)) * w[pos(i, j, q, n, n)])
                + ((sigma * (1.0 - py)) * w[pos(i + 1, j, q, n, n)])
                + (sigma * (1.0 + py) * (w[pos(i - 1, j, q, n, n)]))
        );

    }

    px = Re * ((psi[pos(i + 1, 0, r, n, n)] - psi[pos(i - 1, 0, r, n, n)]) / 4);
    d[0] -= ((sigma * (px - 1.0)) * w[pos(i, 0, q, n, n)]);

    px = Re * ((psi[pos(i + 1, n - 1, r, n, n)] - psi[pos(i - 1, n - 1, r, n, n)]) / 4);
    d[n - 3] += ((sigma * (px + 1.0)) * w[pos(i, n - 1, q, n, n)]);

    solveTriDiagonalMatrix(a, b, c, d, n);

    for (int j = 1; j < n - 1; j++) {
        w[pos(i, j, p, n, n)] = d[i - 1];
    }
}

void solver(double *omega, int n, double Re, double tol = 1e-6, int generateDataInterval = 1, int maxIt = 100) {

    double h = (omega[1] - omega[0]) / (n - 1);

    double dt = (h * h) / 2;

    double sigma = dt / (2 * (h * h));

    double *psi = new double[2 * n * n];
    double *w = new double[2 * n * n];

    int p = 0, q = 1, r = 0, s = 1, it = 0, aux;
    double error = 999;

    while (error > tol && it <= maxIt) {

        aux = p;
        p = q;
        q = aux;
        for (int j = 1; j < n - 1; j++) {
            solveXLinePsi(j, w, psi, n, dt, sigma, p, r, s);
        }

        aux = p;
        p = q;
        q = aux;
        for (int i = 1; i < n - 1; i++) {
            solveYColumnPsi(i, w, psi, n, dt, sigma, p, r, s);
        }

        //contour
        for (int k = 0; k < n; k++) {

            //contour (i, 0):
            w[pos(k, 0, p, n, n)] = -2 * ((psi[pos(k, 1, r, n, n)]) / (h * h));
            w[pos(k, 0, q, n, n)] = -2 * ((psi[pos(k, 1, r, n, n)]) / (h * h));

            //contour (I, j):
            w[pos(n - 1, k, p, n, n)] = -2 * ((psi[pos(n - 2, k, r, n, n)]) / (h * h));
            w[pos(n - 1, k, q, n, n)] = -2 * ((psi[pos(n - 2, k, r, n, n)]) / (h * h));

            //contour (0, j):
            w[pos(0, k, p, n, n)] = -2 * ((psi[pos(1, k, r, n, n)]) / (h * h));
            w[pos(0, k, q, n, n)] = -2 * ((psi[pos(1, k, r, n, n)]) / (h * h));

            //contour (i, 0):
            w[pos(k, n - 1, p, n, n)] = -2 * ((psi[pos(k, n - 2, r, n, n)]) / (h * h)) - (2 / h);
            w[pos(k, n - 1, q, n, n)] = -2 * ((psi[pos(k, n - 2, r, n, n)]) / (h * h)) - (2 / h);

        }

        aux = r;
        r = s;
        s = aux;
        for (int j = 1; j < n - 1; j++) {
            solveXLineW(j, w, psi, n, Re, h, dt, sigma, p, q, r);
        }

        aux = r;
        r = s;
        s = aux;
        for (int i = 1; i < n - 1; i++) {
            solveYColumnW(i, w, psi, n, Re, h, dt, sigma, p, q, r);
        }

        if (it % generateDataInterval == 0) {
            char *folderName = "data", fileNamePsi[100], fileNameW[100], *dataNamePsi = "psi", *dataNameW = "w";
            sprintf(fileNamePsi, "%s/%s%d.vtk", folderName, dataNamePsi, it);
            sprintf(fileNameW, "%s/%s%d.vtk", folderName, dataNameW, it);

            Result result = overviewSolution(w, psi, n, h, omega, Re, p, r);
            printf("---------------INTERAÇÃO %d---------------", it);
            printf("\n\nresiduo: %f", result.maxResidue);
            printf("\nmax psi: %f", result.maxPsi);
            printf("\nmax w: %f", result.maxW);
            printf("\nx_max: %f", result.xCoordMaxValue);
            printf("\ny_max: %f\n", result.yCoordMaxValue);
            saveVtk(fileNamePsi, n, n, psi, dataNamePsi, r);
            saveVtk(fileNameW, n, n, psi, dataNameW, r);
            printf("\n---------------------------------------------\n\n\n\n\n");
        }

        it++;

    }

}

int main() {
    double omega[2] = {0.0, 1.0};

    int n = 32;

    double Re = 1000.0;

    double tol = 1e-6;

    int generateDataInterval = 25;

    int maxIt = 5000;

    solver(omega, n, Re, tol, generateDataInterval, maxIt);

    return 0;
}

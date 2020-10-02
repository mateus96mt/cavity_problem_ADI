#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sys/stat.h>

#define checkNaN

using namespace std;

struct Result {
    double maxPsi, maxW, xCoordMaxValue, yCoordMaxValue, maxResidue;
};

inline int pos(int i, int j, int t, int nx, int ny) { return (t * nx * ny) + (j * (nx)) + i; }

void solveTriDiagonalMatrix(double *a, double *b, double *c, double *d, int n) {

//    double *p;
//    printf("\nVETORES RECEBIDOS PARA RESOLVER MATRIZ TRI-DIAGONAL:");
//    for (int j = 0; j < 4; j++) {
//        switch (j) {
//            case 0:
//                printf("\na: ");
//                p = a;
//                break;
//            case 1:
//                printf("\nb: ");
//                p = b;
//                break;
//            case 2:
//                printf("\nc: ");
//                p = c;
//                break;
//            case 3:
//                printf("\nd: ");
//                p = d;
//                break;
//            default:
//                break;
//        }
//        for (int i = 0; i < n; i++) {
//            printf("%f, ", p[i]);
//        }
//    }
//    printf("\n");

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

void saveVtk(char *nome, int nx, int ny, double *data, char *dataName, int t, bool showLog = false) {
    if (showLog) {
        cout << "output file: " << nome << endl;
    }

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

            if (abs(psi[pos(i, j, r, n, n)]) > maxPsi) {
                maxPsi = abs(psi[pos(i, j, r, n, n)]);
                i_max = i;
                j_max = j;
            }

            if (abs(w[pos(i, j, p, n, n)]) > maxW) {
                maxW = abs(w[pos(i, j, p, n, n)]);
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

void saveOutPut(double *w, double *psi, int n, double h, double *omega, double Re, int p, int r, int it, double tol,
                double error, char *nome) {
    ofstream saida;
    saida.open(nome);

    Result result = overviewSolution(w, psi, n, h, omega, Re, p, r);
    saida << "CONVERGIU EM " << it << " INTERAÇÕES" << endl;
    saida << "TOLERANCIA DE RESIDUO PARA CONVERGENCIA: " << tol << endl;
    saida << "RESIDUO FINAL OBTIDO APOS CONVERGIR:     " << error << endl;
    saida << "\n\nresiduo: " << result.maxResidue;
    saida << "\nmax psi: " << result.maxPsi;
    saida << "\nmax w: " << result.maxW;
    saida << "\nx_max: " << result.xCoordMaxValue;
    saida << "\ny_max: " << result.yCoordMaxValue << endl;
    saida.close();
}

void
calculatePsiX(int j, double *w, double *psi, int n, double dt, double sigma, int p, int r, int s) {

//    printf("\nchamou calculatePsiX");

    double a[n - 2], b[n - 2], c[n - 2], d[n - 2];

    a[0] = 0.0;
    c[n - 3] = 0.0;
    for (int i = 1; i < n - 1; i++) {

        if (i > 1) {
            a[i - 1] = -sigma;
        }

        if (i < n - 2) {
            c[i - 1] = -sigma;
        }

        b[i - 1] = 1.0 + (2.0 * sigma);

        d[i - 1] = (
                ((dt / 2.0) * w[pos(i, j, p, n, n)])
                + ((1.0 - (2.0 * sigma)) * psi[pos(i, j, s, n, n)])
                + (sigma * (psi[pos(i, j + 1, s, n, n)] + psi[pos(i, j - 1, s, n, n)]))
        );

#ifdef checkNaN
        if (isnan(a[i - 1]) || isnan(b[i - 1]) || isnan(c[i - 1]) || isnan(d[i - 1])) {
            fprintf(stderr, "\nERRO NaN NA ETAPA 1 PSI: ");
            if (isnan(a[i - 1])) {
                fprintf(stderr, "a[%d],", i - 1);
            }
            if (isnan(b[i - 1])) {
                fprintf(stderr, "b[%d],", i - 1);
            }
            if (isnan(c[i - 1])) {
                fprintf(stderr, "c[%d],", i - 1);
            }
            if (isnan(d[i - 1])) {
                fprintf(stderr, "d[%d],", i - 1);
            }
            fprintf(stderr, "\n");
            exit(1);
        }
#endif

    }

//    d[0] += (sigma * psi[pos(0, j, r, n, n)]);
//    d[n - 3] += (sigma * psi[pos(n - 1, j, r, n, n)]);

    solveTriDiagonalMatrix(a, b, c, d, n - 2);

    for (int i = 1; i < n - 1; i++) {
#ifdef checkNaN
        if (isnan(d[i - 1])) {
            fprintf(stderr, "\nDEU ERRO DE NaN em d  ETAPA 1 psi ao resolver matriz tri-diagonal!: d[%d], ", i - 1);
            exit(1);
        }
#endif
        psi[pos(i, j, s, n, n)] = d[i - 1];
    }
}

void
calculatePsiY(int i, double *w, double *psi, int n, double dt, double sigma, int p, int r, int s) {

//    printf("\nchamou calculatePsiY");

    double a[n - 2], b[n - 2], c[n - 2], d[n - 2];

    a[0] = 0.0;
    c[n - 3] = 0.0;
    for (int j = 1; j < n - 1; j++) {

        if (j > 1) {
            a[j - 1] = -sigma;
        }

        if (j < n - 2) {
            c[j - 1] = -sigma;
        }

        b[j - 1] = 1.0 + (2.0 * sigma);

        d[j - 1] = (
                ((dt / 2.0) * w[pos(i, j, p, n, n)])
                + ((1.0 - (2.0 * sigma)) * psi[pos(i, j, s, n, n)])
                + (sigma * (psi[pos(i + 1, j, s, n, n)] + psi[pos(i - 1, j, s, n, n)]))
        );

#ifdef checkNaN
        if (isnan(a[j - 1]) || isnan(b[j - 1]) || isnan(c[j - 1]) || isnan(d[j - 1])) {
            fprintf(stderr, "\nERRO NaN NA ETAPA 2 PSI: ");
            if (isnan(a[j - 1])) {
                fprintf(stderr, "a[%d],", j - 1);
            }
            if (isnan(b[j - 1])) {
                fprintf(stderr, "b[%d],", j - 1);
            }
            if (isnan(c[j - 1])) {
                fprintf(stderr, "c[%d],", j - 1);
            }
            if (isnan(d[j - 1])) {
                fprintf(stderr, "d[%d],", j - 1);
            }
            fprintf(stderr, "\n");
            exit(1);
        }
#endif


    }

//    d[0] += (sigma * psi[pos(i, 0, r, n, n)]);
//    d[n - 3] += (sigma * psi[pos(i, n - 1, r, n, n)]);

    solveTriDiagonalMatrix(a, b, c, d, n - 2);

    for (int j = 1; j < n - 1; j++) {
#ifdef checkNaN
        if (isnan(d[j - 1])) {
            fprintf(stderr, "\nDEU ERRO DE NaN em d  ETAPA 1 psi ao resolver matriz tri-diagonal!: d[%d], ", j - 1);
            exit(1);
        }
#endif
        psi[pos(i, j, s, n, n)] = d[j - 1];
    }
}

void
calculateWX(int j, double *w, double *psi, int n, double Re, double h, double dt, double sigma, int p, int q, int r) {

//    printf("\nchamou calculateWX");

    double a[n - 2], b[n - 2], c[n - 2], d[n - 2];
    double px, py;

    a[0] = 0.0;
    c[n - 3] = 0.0;
    for (int i = 1; i < n - 1; i++) {

        px = Re * ((psi[pos(i + 1, j, r, n, n)] - psi[pos(i - 1, j, r, n, n)]) / 4);
        py = Re * ((psi[pos(i, j + 1, r, n, n)] - psi[pos(i, j - 1, r, n, n)]) / 4);

        if (i > 1) {
            a[i - 1] = -sigma * (py + 1.0);
        }

        if (i < n - 2) {
            c[i - 1] = sigma * (py - 1.0);
        }

        b[i - 1] = 1.0 + (2.0 * sigma);

        d[i - 1] = (
                ((1.0 - (2.0 * sigma)) * w[pos(i, j, q, n, n)])
                + ((sigma * (px + 1.0)) * w[pos(i, j + 1, q, n, n)])
                + (sigma * (1.0 - px) * (w[pos(i, j - 1, q, n, n)]))
        );

#ifdef checkNaN
        if (isnan(a[i - 1]) || isnan(b[i - 1]) || isnan(c[i - 1]) || isnan(d[i - 1])) {
            fprintf(stderr, "\nERRO NaN NA ETAPA 1 W: ");
            if (isnan(a[i - 1])) {
                fprintf(stderr, "a[%d],", i - 1);
            }
            if (isnan(b[i - 1])) {
                fprintf(stderr, "b[%d],", i - 1);
            }
            if (isnan(c[i - 1])) {
                fprintf(stderr, "c[%d],", i - 1);
            }
            if (isnan(d[i - 1])) {
                fprintf(stderr, "d[%d],", i - 1);
            }
            fprintf(stderr, "\n");
            exit(1);
        }
#endif


    }

    py = Re * ((psi[pos(0, j + 1, r, n, n)] - psi[pos(0, j - 1, r, n, n)]) / 4);
    d[0] += ((sigma * (py + 1.0)) * w[pos(0, j, p, n, n)]);

    py = Re * ((psi[pos(n - 1, j + 1, r, n, n)] - psi[pos(n - 1, j - 1, r, n, n)]) / 4);
    d[n - 3] -= ((sigma * (py - 1.0)) * w[pos(n - 1, j, p, n, n)]);

    solveTriDiagonalMatrix(a, b, c, d, n - 2);

    for (int i = 1; i < n - 1; i++) {
#ifdef checkNaN
        if (isnan(d[i - 1])) {
            fprintf(stderr, "\nDEU ERRO DE NaN em d  ETAPA 1 W ao resolver matriz tri-diagonal!: d[%d], ", i - 1);
            exit(1);
        }
#endif
        w[pos(i, j, p, n, n)] = d[i - 1];
    }
}

void
calculateWY(int i, double *w, double *psi, int n, double Re, double h, double dt, double sigma, int p, int q, int r) {

//    printf("\nchamou calculateWY");

    double a[n - 2], b[n - 2], c[n - 2], d[n - 2];
    double px, py;

    a[0] = 0.0;
    c[n - 3] = 0.0;
    for (int j = 1; j < n - 1; j++) {

        px = Re * ((psi[pos(i + 1, j, r, n, n)] - psi[pos(i - 1, j, r, n, n)]) / 4);
        py = Re * ((psi[pos(i, j + 1, r, n, n)] - psi[pos(i, j - 1, r, n, n)]) / 4);

        if (j > 1) {
            a[j - 1] = sigma * (px - 1.0);
        }

        if (j < n - 2) {
            c[j - 1] = -sigma * (px + 1.0);
        }

        b[j - 1] = 1.0 + (2.0 * sigma);

        d[j - 1] = (
                ((1.0 - (2.0 * sigma)) * w[pos(i, j, q, n, n)])
                + ((sigma * (1.0 - py)) * w[pos(i + 1, j, q, n, n)])
                + (sigma * (1.0 + py) * (w[pos(i - 1, j, q, n, n)]))
        );

#ifdef checkNaN
        if (isnan(a[j - 1]) || isnan(b[j - 1]) || isnan(c[j - 1]) || isnan(d[j - 1])) {
            fprintf(stderr, "\nERRO NaN NA ETAPA 2 W: ");
            if (isnan(a[j - 1])) {
                fprintf(stderr, "a[%d],", j - 1);
            }
            if (isnan(b[j - 1])) {
                fprintf(stderr, "b[%d],", j - 1);
            }
            if (isnan(c[j - 1])) {
                fprintf(stderr, "c[%d],", j - 1);
            }
            if (isnan(d[j - 1])) {
                fprintf(stderr, "d[%d],", j - 1);
            }
            fprintf(stderr, "\n");
            exit(1);
        }
#endif

    }

    px = Re * ((psi[pos(i + 1, 0, r, n, n)] - psi[pos(i - 1, 0, r, n, n)]) / 4);
    d[0] -= ((sigma * (px - 1.0)) * w[pos(i, 0, q, n, n)]);

    px = Re * ((psi[pos(i + 1, n - 1, r, n, n)] - psi[pos(i - 1, n - 1, r, n, n)]) / 4);
    d[n - 3] += ((sigma * (px + 1.0)) * w[pos(i, n - 1, q, n, n)]);

    solveTriDiagonalMatrix(a, b, c, d, n - 2);

    for (int j = 1; j < n - 1; j++) {
#ifdef checkNaN
        if (isnan(d[j - 1])) {
            fprintf(stderr, "\nDEU ERRO DE NaN em d  ETAPA 2 W ao resolver matriz tri-diagonal!: d[%d], ", j - 1);
            exit(1);
        }
#endif

        w[pos(i, j, p, n, n)] = d[j - 1];
    }
}

void solver(double *omega, int n, double Re, double h, double dt, double tol = 1e-6, int generateDataInterval = 1,
            int maxIt = 100, const char *folderName = "data", const char *saidaNome = "saida%dx%d.txt") {

    double sigma = dt / (2 * (h * h));

    double *psi = new double[2 * n * n];
    double *w = new double[2 * n * n];

    //initial condition
    for (int i = 0; i < 2 * n * n; i++) {
        psi[i] = 0.0;
        w[i] = 0.0;
    }

    int p = 0, q = 1, r = 0, s = 1, it = 1, aux;
    double error = 999;

    char fileNamePsi[100], fileNameW[100], *dataNamePsi = "psi", *dataNameW = "w";

    sprintf(fileNamePsi, "%s/%s%d.vtk", folderName, dataNamePsi, 0);
    sprintf(fileNameW, "%s/%s%d.vtk", folderName, dataNameW, 0);

    saveVtk(fileNamePsi, n, n, psi, dataNamePsi, r);
    saveVtk(fileNameW, n, n, w, dataNameW, r);

    while (error > tol && it <= maxIt) {

//        printf("\n---------------INTERAÇÃO %d---------------", it);

        aux = r;
        r = s;
        s = aux;
        for (int j = 1; j < n - 1; j++) {
            calculatePsiX(j, w, psi, n, dt, sigma, p, r, s);
        }

        aux = r;
        r = s;
        s = aux;
        for (int i = 1; i < n - 1; i++) {
            calculatePsiY(i, w, psi, n, dt, sigma, p, r, s);
        }

        //contour (i, 0) LOW:
        for (int k = 0; k < n; k++) {
            w[pos(k, 0, p, n, n)] = -2 * ((psi[pos(k, 1, r, n, n)]) / (h * h));
            w[pos(k, 0, q, n, n)] = -2 * ((psi[pos(k, 1, r, n, n)]) / (h * h));
        }
        //contour (I, j) RIGHT:
        for (int k = 0; k < n; k++) {
            w[pos(n - 1, k, p, n, n)] = -2 * ((psi[pos(n - 2, k, r, n, n)]) / (h * h));
            w[pos(n - 1, k, q, n, n)] = -2 * ((psi[pos(n - 2, k, r, n, n)]) / (h * h));
        }
        //contour (0, j) LEFT:
        for (int k = 0; k < n; k++) {
            w[pos(0, k, p, n, n)] = -2 * ((psi[pos(1, k, r, n, n)]) / (h * h));
            w[pos(0, k, q, n, n)] = -2 * ((psi[pos(1, k, r, n, n)]) / (h * h));
        }
        //contour (i, 0) UP:
        for (int k = 0; k < n; k++) {
            w[pos(k, n - 1, p, n, n)] = -2 * ((psi[pos(k, n - 2, r, n, n)]) / (h * h)) - (2 / h);
            w[pos(k, n - 1, q, n, n)] = -2 * ((psi[pos(k, n - 2, r, n, n)]) / (h * h)) - (2 / h);
        }

        aux = p;
        p = q;
        q = aux;
        for (int j = 1; j < n - 1; j++) {
            calculateWX(j, w, psi, n, Re, h, dt, sigma, p, q, r);
        }

        aux = p;
        p = q;
        q = aux;
        for (int i = 1; i < n - 1; i++) {
            calculateWY(i, w, psi, n, Re, h, dt, sigma, p, q, r);
        }

        if (it % generateDataInterval == 0) {

            sprintf(fileNamePsi, "%s/%s%d.vtk", folderName, dataNamePsi, it);
            sprintf(fileNameW, "%s/%s%d.vtk", folderName, dataNameW, it);

            Result result = overviewSolution(w, psi, n, h, omega, Re, p, r);
            printf("---------------INTERAÇÃO %d---------------", it);
            printf("\n\nresiduo: %.20f", result.maxResidue);
            printf("\nmax psi: %.20f", result.maxPsi);
            printf("\nmax w: %.20f", result.maxW);
            printf("\nx_max: %.20f", result.xCoordMaxValue);
            printf("\ny_max: %.20f\n", result.yCoordMaxValue);
            saveVtk(fileNamePsi, n, n, psi, dataNamePsi, r);
            saveVtk(fileNameW, n, n, w, dataNameW, r);
            printf("\n---------------------------------------------\n\n\n\n\n");

            error = result.maxResidue;
        }

//        printf("\n---------------------------------------------\n\n\n\n\n");
        it++;

    }

    if (it < maxIt) {
        printf("CONVERGIU EM %d INTERAÇÕES\n", it - 1);
        printf("TOLERANCIA DE RESIDUO PARA CONVERGENCIA: %.20f\n", tol);
        printf("RESIDUO FINAL OBTIDO APOS CONVERGIR:     %.20f\n", error);
        char nome[100];

        sprintf(nome, saidaNome, n, n);
        saveOutPut(w, psi, n, h, omega, Re, p, r, it - 1, tol, error, nome);
    }

    printf("\n\nn: %d", n);
    printf("\n\nRe: %f", Re);
    printf("\n\ntol: %f", tol);
    printf("\n\ninterval: %d", generateDataInterval);
    printf("\n\nmaxIt: %d", maxIt);
    printf("\n\nh: %f", h);
    printf("\n\ndt: %f", dt);
    printf("\n\nsigma: %f", sigma);
}

int main() {
    double omega[2] = {0.0, 1.0};

    int n = 32;

    double h = (omega[1] - omega[0]) / (n - 1);

    double dt = (h * h) / 10;

    double Re = 1000.0;

    double tol = 1e-6;

    int generateDataInterval = 100;

    int maxIt = 1000000;

    char outPutFolder[100], saidaNome[100];

    sprintf(outPutFolder, "dados_Re_%d_%dx%d", (int) Re, n, n);

    sprintf(saidaNome, "saida_Re_%d_%dx%d.txt", (int) Re, n, n);

    mkdir(outPutFolder, 0777);

    solver(omega, n, Re, h, dt, tol, generateDataInterval, maxIt, outPutFolder, saidaNome);

    return 0;
}

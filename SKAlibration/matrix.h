#ifndef MATRIX
#define MATRIX
#include <pthread.h>
#include "complex.h"


typedef struct Matrix_i
{
    int *data;
    unsigned int nbRow;
    unsigned int nbColumn;
}Matrix_i;

typedef struct Matrix_f
{
    Complex_f *data;
    unsigned int nbRow;
    unsigned int nbColumn;
}Matrix_f;


int calibration (int argc, char **argv);

/*****************************************************
 *
 *        Prototypes to handle int matrix
 *
 * **************************************************/

void printMatrixI(Matrix_i *a);

Matrix_i * matrixAllocI(unsigned int rows, unsigned int column);

void freeMatrixI(Matrix_i *A);

int matrixGetI(Matrix_i *A, unsigned int i, unsigned int j);

void matrixSetI(Matrix_i *A, unsigned int i, unsigned int j, const int c);

void create_grid(Matrix_i *A, Matrix_i *B);

Matrix_i * make1DvectorFromMatrix(Matrix_i *A);

void compute_LBL_Matrix(Matrix_i *LBL, Matrix_i *A0, Matrix_i *A1);


/*****************************************************
 *
 *    Prototypes to handle simple precision matrix
 *
 * **************************************************/

void printMatrixF(Matrix_f *a);

Matrix_f * matrixAllocF(unsigned int rows, unsigned int column);

Matrix_f * createRandomMatrixF(unsigned int rows, unsigned int column);

void createZeroMatrixF(Matrix_f *J);

Matrix_f *createTestVectorF(unsigned int rows, unsigned int column);

void createIdentityF(Matrix_f *);

void freeMatrixF(Matrix_f *A);

Complex_f matrixGetF(Matrix_f *A, unsigned int i, unsigned int j);

void matrixSetF(Matrix_f *A, unsigned int i, unsigned int j, const Complex_f c);

void scaleMatrixF(Matrix_f *a, const Complex_f c);

void addMatrixF(Matrix_f *a, Matrix_f *b);

void subMatrixF(Matrix_f *a, Matrix_f *b);

void matNormF(Matrix_f *A);

void scaleLineF(Matrix_f *a, unsigned int i, Complex_f f);

int subXLinesF(Matrix_f *a, unsigned int i, unsigned int i2, Complex_f f);

/*****************************************************
 *
 *   Adapted functions to handle our new Jacobian
 *
 * **************************************************/

void giveJfF(Matrix_f *A, Matrix_f * sky_c , Matrix_i * lbl, Matrix_f *j, int nt, int na);

void compute_data_vectorF(Matrix_f *j, Matrix_f *A,  Matrix_f *res, int na, int nt);


#endif /* MATRIX*/


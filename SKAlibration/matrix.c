#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include "complex.h"
#include "matrix.h"
#include "thread_utils.h"


#ifdef __k1__
    #include <mppa_power.h>
    #include <mppa_routing.h>
    #include <mppa_remote.h>
    #include <mppa_rpc.h>
    #include <mppa_async.h>
    #include <assert.h>
    #include <utask.h>// for pthread_barrier_t
    #include <HAL/hal/board/boot_args.h>
#endif



#define NB_POLA (4)


long int MyMem = 0;


/*****************************************************
     *
     *        Prototypes to handle int matrix
     *
     * **************************************************/


/**
     * @brief printMatrixF display a matrix (usefull for testing)
     * @param a the matrix to display
     */
void printMatrixI(Matrix_i *a)
{
    if(a)
    {
        unsigned int i, j;
        for (i=0; i<a->nbRow; i++)
        {
            for(j=0; j<a->nbColumn; j++)
            {
                printf("%d ",matrixGetI(a, i, j));
                printf("\t");
            }
            printf("\n");
        }
    }
}


/**
     * @brief matrixAllocI allocate memory for a matrix and initialize
     * @param rows number of rows
     * @param column number of columns
     * @return the allocated matrix if no errors occured, NULL if not
     */
Matrix_i * matrixAllocI(unsigned int rows, unsigned int column)
{
    Matrix_i *A = (Matrix_i *)malloc(sizeof(Matrix_i));
    int *c=(int*)calloc(rows*column,sizeof(int));
    if (!A || !c)
    {
        fprintf(stderr, " Matrix Error \n in file: %s \n function %s\n line : %d\n An error occured while allocating memory for matrices \n", __FILE__, __FUNCTION__,__LINE__);
        return NULL;
    }
    A->nbColumn = column;
    A->nbRow = rows;
    A->data = c;
    MyMem += sizeof(Matrix_i) + rows*column*sizeof(int);
    return A;
}


/**
     * @brief freeMatrixI free a int matrix
     * @param A the matrix to free
     */
void freeMatrixI(Matrix_i *A)
{
    if(A)
    {
        free (A->data);
        free (A);
    }
}

/**
     * @brief matrixGetI access a specified data of a matrix
     * @param A the matrix
     * @param i the index of rows
     * @param j the index on columns
     * @return the data A(i,j)
     */
int matrixGetI(Matrix_i *A, unsigned int i, unsigned int j)
{
    if(A)
        return A->data[i*A->nbColumn+j];
    return 0;
}

/**
     * @brief matrixSetI set the value of a specified data
     * @param A the matrix
     * @param i the index on rows
     * @param j the index on columns
     * @param c the complex to set at A(i,j)
     */
void matrixSetI(Matrix_i *A, unsigned int i, unsigned int j, const int c)
{
    if(A)
        A->data[i*A->nbColumn+j] = c;
}


/**
     * @brief create_grid taking values from 0 to n (n*n being the size of the TWO matrices)
     * @param A will be the grid with respect to the lines
     * @param B will be the grid with respect to the columns
     */
void create_grid(Matrix_i *A, Matrix_i *B)
{
    if(A&&B)
    {
        unsigned int i, j;
        for(i=0; i< A->nbRow; i++)
        {
            for(j=0; j<A->nbColumn; j++)
            {
                matrixSetI(B, i, j, i);
                matrixSetI(A, i, j, j);
            }
        }
    }
}


/**
     * @brief make1DvectorFromMatrix make an array of the values present in a matrix
     * @param A the 2D atrix we want to set in 1D
     * @return a 1D matrix containing the A matrix value
     * @note the A parameter will be freed
     */
Matrix_i * make1DvectorFromMatrix(Matrix_i *A)
{
    if(A)
    {
        Matrix_i *res = matrixAllocI(1,A->nbRow* A->nbColumn);
        double tmp;
        unsigned int i, j;
        for(i=0; i< A->nbRow; i++)
        {
            for(j=0; j<A->nbColumn; j++)
            {
                tmp = matrixGetI(A,i,j);
                matrixSetI(res, 0, i*A->nbColumn + j, tmp);
            }
        }
        return res;
    }
    return NULL;
}


/**
     * @brief compute_LBL_Matrix compute a matrix that lists all the baselines to compute the jacobian
     * @param LBL an initialized array that will be filled after function call
     * @param A0 a 1D vector containing values from 0 to na-1
     * @param A1 a 1D vector containing values from 0 to na-1 (with a different order)
     */
void compute_LBL_Matrix(Matrix_i *LBL, Matrix_i *A0, Matrix_i *A1)
{
    if( LBL && A0 && A1)
    {
        unsigned int i, countLBL=0;
        int a, b;
        for (i=0; i<A0->nbColumn; i++)
        {
            b = matrixGetI(A1, 0, i);
            a = matrixGetI(A0, 0, i);

            if(a!=b&&countLBL<LBL->nbRow)
            {
                matrixSetI(LBL,countLBL,0, b);
                matrixSetI(LBL,countLBL,1, a);
                ++countLBL;
            }
        }
    }
}


/*****************************************************
     *
     *   Functions to handle a simple precision matrix
     *
     * **************************************************/

static Complex_f getGainComponent(Matrix_f * Gain, int na, int npola)
{
    return matrixGetF(Gain, npola+NB_POLA*na,0);
}


/**
     * @brief printMatrixF display a matrix (usefull for testing)
     * @param a the matrix to display
     */
void printMatrixF(Matrix_f *a)
{
    if(a)
    {
        unsigned int i, j;
        printf("\n\n\n");
        for (i=0; i<a->nbRow; i++)
        {
            for(j=0; j<a->nbColumn; j++)
            {
                printComplexF(matrixGetF(a, i, j));
                printf("\t");
            }
            printf("\n");
        }
    }
}



/**
     * @brief matrixAllocF allocate memory for a matrix and initialize it to zero
     * @param rows number of rows
     * @param column number of columns
     * @return the allocated matrix if no errors occured, NULL if not
     */
Matrix_f * matrixAllocF(unsigned int rows, unsigned int column)
{
    Matrix_f *A = (Matrix_f *)calloc(1,sizeof(Matrix_f));
    Complex_f *c=(Complex_f *)calloc(rows*column, sizeof(Complex_f));
    if (!A || !c)
    {
        fprintf(stderr,"Matrix Error \n in file: %s \n function %s\n line : %d\n Unable to allocate memory \n", __FILE__,__FUNCTION__, __LINE__);
        return NULL;
    }
    A->nbColumn = column;
    A->nbRow = rows;
    A->data = c;
    createZeroMatrixF(A);
    MyMem += sizeof(Matrix_f) + rows*column*sizeof(Complex_f);
    return A;
}


/**
     * @brief createRandomMatrix initialize a random Matrix
     * @param rows number of rows
     * @param column number of columns
     * @return the random matrix
     */
Matrix_f * createRandomMatrixF(unsigned int rows, unsigned int column)
{
    Matrix_f *res = matrixAllocF(rows, column);
    unsigned int i,j;
    for (i=0; i<res->nbRow; i++)
    {
        for(j=0; j<res->nbColumn; j++)
        {
            matrixSetF(res, i, j, createRandomComplexF());
        }
    }
    return res;
}


/**
     * @brief createRandomMatrix initialize a random Matrix
     * @param rows number of rows
     * @param column number of columns
     * @return the random matrix
     */
void createZeroMatrixF(Matrix_f *J)
{
    if(J)
    {
        unsigned int i,j;
        for (i=0; i<J->nbRow; i++)
        {
            for(j=0; j<J->nbColumn; j++)
            {
                matrixSetF(J, i, j, createNulComplexF());
            }
        }
    }
}


/**
     * @brief createTestVector initialize a known Matrix
     * @param rows number of rows
     * @param column number of columns
     * @return the random matrix
     */
Matrix_f *createTestVectorF(unsigned int rows, unsigned int column)
{

    Matrix_f *res = matrixAllocF(rows, column);
    unsigned int i,j;
    Complex_f tmp;
    tmp.im = tmp.re = 0.5f;
    for (i=0; i<rows; i++)
    {
        for(j=0; j<column; j++)
        {
            matrixSetF(res, i, j, tmp);
        }
    }
    return res;
}


/**
     * @brief createIdentityF compute the identoty matrix
     * @param I A pointer to store the result
     */
void createIdentityF(Matrix_f *I)
{
    if(I && I->nbColumn==I->nbRow)
    {
        unsigned int i;
        Complex_f a;
        a.re =1;
        a.im =0;
        for (i=0; i<I->nbRow; i++)
        {
            matrixSetF(I, i, i, a);
        }
    }
}


/**
     * @brief freeMatrixF free a complex float matrix
     * @param A the matrix to free
     */
void freeMatrixF(Matrix_f *A)
{
    if(A)
    {
        free (A->data);
        free (A);
    }
}



/**
     * @brief matrixGetF access a specified data of a matrix
     * @param A the matrix
     * @param i the index of rows
     * @param j the index on columns
     * @return the data A(i,j)
     */
Complex_f matrixGetF(Matrix_f *A, unsigned int i, unsigned int j)
{
    if(A)
        return A->data[i*A->nbColumn+j];
    return createNulComplexF();
}



/**
     * @brief matrixSetF set the value of a specified data
     * @param A the matrix
     * @param i the index on rows
     * @param j the index on columns
     * @param c the complex to set at A(i,j)
     */
void matrixSetF(Matrix_f *A, unsigned int i, unsigned int j, const Complex_f c)
{
    if(A)
        A->data[i*A->nbColumn+j] = c;
}



/**
     * @brief scaleMatrixF multiply each data by a constant (complex)
     * @param a the matrix to scale
     * @param c the complex
     * @todo Parallel programming
     */
void scaleMatrixF(Matrix_f *a, const Complex_f c)
{
    if(a)
    {
        unsigned int i, j;
        Complex_f tmp;
        for (i=0; i<a->nbRow; i++)
        {
            for(j=0; j<a->nbColumn; j++)
            {
                tmp = matrixGetF(a, i, j);
                matrixSetF(a, i, j, multiplyF(tmp, c));
            }
        }
    }
}


/**
     * @brief addMatrixF compute A = A+B
     * @param a complex matrix
     * @param b complex matrix
     * @param c a matrix to store the result
     * @todo Parallel programming
     */
void addMatrixF(Matrix_f *a, Matrix_f *b)
{
    if ((a->nbColumn==b->nbColumn) && (a->nbRow==b->nbRow))
    {
        unsigned int i, j;
        Complex_f tmpa, tmpb;
        for (i=0; i<a->nbRow; i++)
        {
            for(j=0; j<a->nbColumn; j++)
            {
                tmpa = matrixGetF(a, i,j);
                tmpb = matrixGetF(b, i, j);
                matrixSetF(a, i, j, addF(tmpa, tmpb));
            }
        }
    }
}


/**
     * @brief addMatrixF compute A = A+BT
     * @param a complex matrix
     * @param b complex matrix
     * @param c a matrix to store the result
     * @todo Parallel programming
     */
void addMatrixTF(Matrix_f *a, Matrix_f *b)
{
    if ((a->nbColumn==b->nbRow) && (a->nbRow==b->nbColumn))
    {
        unsigned int i, j;
        Complex_f tmpa, tmpb;
        for (i=0; i<a->nbRow; i++)
        {
            for(j=0; j<a->nbColumn; j++)
            {
                tmpa = matrixGetF(a, i,j);
                tmpb = matrixGetF(b, j, i);
                matrixSetF(a, i, j, addF(tmpa, tmpb));
            }
        }
    }
}
/**
     * @brief subMatrixF compute C = A-B
     * @param a complex matrix
     * @param b complex matrix
     * @param c a matrix to store the result
     * @todo Parallel programming
     */
void subMatrixF(Matrix_f *a, Matrix_f *b)
{
    if (a && b)
    {
        if ((a->nbColumn==b->nbColumn) && (a->nbRow==b->nbRow))
        {
            unsigned int i, j;
            Complex_f tmpa, tmpb;
            for (i=0; i<a->nbRow; i++)
            {
                for(j=0; j<a->nbColumn; j++)
                {
                    tmpa = matrixGetF(a, i,j);
                    tmpb = matrixGetF(b, i, j);
                    matrixSetF(a, i, j, subF(tmpa, tmpb));
                }
            }
        }
    }
}


/**
     * @brief matNormF take the first matrix A value as a reference (ie: gain = 1 and no imaginary part)
     * @param A the matrix that we want to "norm"
     */
void matNormF(Matrix_f *A)
{
    if (A)
    {
        Complex_f c = matrixGetF(A, 0, 0);
        Complex_f cconj = conjugateF(c);
        c=scaleF(cconj, 1/gainF(c));
        scaleMatrixF(A, c);
    }
}




/**
     *                          THE DATA VECTOR IS THE VISIBILITY VECTOR
     *
     * @brief compute_data_vector  uses the property to y=J_x.x, with J_x being the jacobian estimated at x
     * @param A the 1D matrix that will produce the data vector
     * @return the data Vector
     *
     * we use our knowledge on the shape of the jacobian to spare another matrix product
     */
void compute_data_vectorF(Matrix_f *jac, Matrix_f *A,  Matrix_f *res, int na, int nt)
{
    if(A&&res&&jac)
    {
        int i, k, l;
        Complex_f tmp, tmpa;
        int nb = (na-1)*nt*NB_POLA;
        for(i=0; i<na; i++)
        {
            for(k=0; k<nb; k++)
            {
                tmpa.re = tmpa.im = 0;
                for(l=0; l<NB_POLA; l++)
                {
                    tmp = multiplyF(matrixGetF(jac,i*nb+k,l), matrixGetF(A,i*NB_POLA+l,0));
                    tmpa = addF(tmpa,tmp);
                }
                matrixSetF(res,i*nb+k, 0,tmpa);
            }
        }
    }
}


/**
     * @brief invert_diag_matrixF as we know that the hessian is a diagonal matrix we just have to invert it's elements
     * @param a the matrix to invert
     * @param i an empty matrix to store the inverse
     */
void invert_diag_matrixF(Matrix_f *a, Matrix_f *i)
{
    if(a && i)
    {
        unsigned int j;
        Complex_f c;
        createIdentityF(i);
        for(j=0; j<a->nbRow; j++)
        {
            c=matrixGetF(a,j,j);
            matrixSetF(i,j, j, invF(c));
        }
    }
}



/**
     * @brief scaleLineF multiply a line with a complex
     * @param a the matrix
     * @param i the index of the line we want to scale
     * @param f the complex that will scale the line
     */
void scaleLineF(Matrix_f *a, unsigned int i, Complex_f f)
{
    if(a)
    {
        if(i<a->nbRow)
        {
            Complex_f c;
            unsigned int j;
            for(j=0; j<a->nbColumn; j++)
            {
                c = multiplyF(f, matrixGetF(a, i, j));
                matrixSetF(a, i, j, c);
            }
        }
    }
}


/**
     * @brief subXLinesF computes line i - f*line i2 and store it into line i
     * @param a the matrix
     * @param i the index of the line to update
     * @param i2 the index of the "reference" line
     * @param f a complex to scale the matrix
     */
int subXLinesF(Matrix_f *a, unsigned int i, unsigned int i2, Complex_f f)
{
    int ret = 1;
    if(a)
    {
        if(i<a->nbRow)
        {
            Complex_f c,tmp, tmp1;
            unsigned int j;
            for(j=0; j<a->nbColumn; j++)
            {
                tmp = matrixGetF(a, i, j);
                tmp1 = matrixGetF(a, i2, j);
                tmp1 =multiplyF(f, tmp1);
                c = subF(tmp, tmp1);
                matrixSetF(a, i, j, c);
                if(c.re >1e-6 || c.im >1e-6)
                {
                    ret = 0;
                }
            }
        }
    }
    return ret;
}



/**
     * @brief giveJfF function that computes the half Jacobian of a complex float matrix
     * @param A the 1D matrix on which we want to know the jacobian
     * @param j the destination of the jacobian
     * @return the half jacobian matrix
     *
     * again, we use the properties of this system to reduce the comlexity
     */

/*
void giveJfF(Matrix_f *A, Matrix_f * sky_c , Matrix_i * lbl, Matrix_f *jac, int nt, int na)
{
    if(A && jac)
    {
        int i,j,k,l,m;
        int b;
        Complex_f insert, tmp;
        for(i=0; i<na; i++)
        {
            for(j=0; j<(na-1); j++)
            {
                b = (int) matrixGetI(lbl, i, 1);
                for(k=0; k<nt; k++)
                {
                    for(l=0; l<NB_POLA; l++)
                    {
                        for (m=0; m<NB_POLA;  m++)
                        {
                            if((l+m)%2 == 0)
                            {
                                insert = conjugateF(getGainComponent(A,b,m)); // get the gain corresponding to Q at pola n
                                tmp = matrixGetF(sky_c,(i*(na-1)+j)*nt+k, l); //get the sky coherency for the curent baseline
                                insert = multiplyF(insert, tmp);
                                matrixSetF(jac,i*(na-1)*nt*NB_POLA + j*nt*NB_POLA +k*NB_POLA + l, m, insert);
                            }
                        }
                    }
                }
            }
        }
    }
}*/

static void matrix_product(Matrix_f *gain, Complex_f *sky_coherency, Matrix_f *c, int baseline)
{
    int i, j, k;
    Complex_f tmp, tmpa;
    for(i=0; i<2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            tmp.re=tmp.im=0;
            for(k=0; k<2; k++)
            {
                tmpa = conjugateF(getGainComponent(gain, baseline, j*2+k));
                tmpa = multiplyF(tmpa, sky_coherency[i*2+k]);
                tmp  = addF(tmp, tmpa);
            }
            matrixSetF(c, i, j, tmp);
        }
    }
}

void giveJfF(Matrix_f *A, Matrix_f * sky_c , Matrix_i * lbl, Matrix_f *jac, int nt, int na)
{
    if(A && jac)
    {
        int i,j,k,l, count_tab;
        int b;
        int tab[] = {0, 5, 8, 13, 2, 7, 10, 15}; 
        Complex_f insert;
        Matrix_f *tmpmat = matrixAllocF(2,2);
        for(i=0; i<na; i++)
        {
            for(j=0; j<(na-1); j++)
            {
                b = (int) matrixGetI(lbl, i, 1);
                for(k=0; k<nt; k++)
                {
                    matrix_product(A, &(sky_c->data[(i*(na-1)+j)*nt+k]), tmpmat, b);
                    count_tab =0;
                    for(l=0; l<NB_POLA; l++)
                    {
                        insert = matrixGetF(tmpmat, 0, l);
                        matrixSetF(jac,i*(na-1)*nt*NB_POLA + j*nt*NB_POLA +k*NB_POLA,tab[count_tab++], insert);
                        matrixSetF(jac,i*(na-1)*nt*NB_POLA + j*nt*NB_POLA +k*NB_POLA, tab[count_tab++], insert);

                    }
                }
            }
        }
        freeMatrixF(tmpmat);
    }
}
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "complex.h"
#include "matrix.h"


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
    tmp.im = tmp.re = 1.0f;
    for (i=0; i<rows; i++)
    {
        for(j=0; j<column; j++)
        {
        	tmp.re /= (i*0.3+j*0.2+1);
        	tmp.im /=(i*0.5+j*0.3+1);
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
void freeMatrixF(Matrix_f **A)
{
    if(*A)
    {
		if((*A)->data) {
			free((*A)->data);
			(*A)->data = NULL;
		}
        free(*A);
		*A = NULL;
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
       float gain  =0.;
        scaleMatrixF(A, cconj);
        c = matrixGetF(A, 0, 0);
        gain = gainF(c);
        gain = 1/gain;
        c.re = gain;
        c.im = 0;
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
void compute_data_vectorF(Matrix_f *jac, Matrix_f *A,  Matrix_f *res, int na)
{
    int i,j;
    Complex_f tmp, tmpa;
    for(i=0; i<na; i++)
    {
    	tmpa.re=tmpa.im = 0;
   	 	for (j=0; j<na-1; j++)
   	 	{
			tmp =  multiplyF(matrixGetF(jac,i*(na-1)+j,0),matrixGetF(A,i,0)); 
			tmpa = addF(tmpa,tmp);
			matrixSetF(res, i*(na-1)+j, 0, tmpa);
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
*	@brief function aimed to give the Jacobian matrix
*	@param A the matrix from which we want the Jacobian
*	@param lbl a matrix that lists all the baselines
*	@param jac a pointer to the INITIALIZED Jacobian matrix
*	@param na the number of antenna in the calibration problem
*/
void giveJfF(Matrix_f *A, Matrix_i * lbl, Matrix_f *jac, int na)
{
    if(A && jac)
    {
        int i,j;
        int b;
        Complex_f insert;
        for(i=0; i<na; i++)
        {
            for(j=0; j<(na-1); j++)
            {
                b = (int) matrixGetI(lbl, 2*i+j, 1);
                insert = conjugateF(matrixGetF(A,b,0));
                matrixSetF(jac,i*(na-1) + j,0, insert);
            }
        }
    }
    //printMatrixF(jac);
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex.h"

/**
 * @file complex.cpp to handle basic complexes functions
 * @author Nicolas Sourbier
 * @date 22/06/2017
 **/


int equals(float a, float b, float precis)
{
  if (fabs(a-b) < precis) return 1;
  return 0;
}

/**********************************************************************
 *          ---- Useful functions for float complexes ----
 * *******************************************************************/

/**
 * @brief conjugateF computes the conjugate of a complex
 * @param c the complex we want to know the conjugate
 * @return a new complex containing the conjugate of the one given in param
 */
Complex_f conjugateF(Complex_f c)
{
    Complex_f res = c;
    res.im *= -1;
    return res;
}

/**
 * @brief addF add two complexes
 * @param a the first operand
 * @param b the complex to add to the first one
 * @return  the sum of the two complexes
 */
Complex_f addF(Complex_f a, Complex_f b)
{
    Complex_f res;
    res.re = a.re + b.re;
    res.im = a.im + b.im;
    return res;
}

/**
 * @brief subF substract two complexes (a-b)
 * @param a the fist operand
 * @param b the complex we want to substract
 * @return the difference of the two
 */
Complex_f subF(Complex_f a, Complex_f b)
{
    Complex_f res;
    res.re = a.re - b.re;
    res.im = a.im - b.im;
    return res;
}

/**
 * @brief scaleF multiply a complex by a constant
 * @return the product of a and c
 */
Complex_f scaleF(Complex_f a, float c)
{
    Complex_f res;
    res.re = a.re*c;
    res.im = a.im*c;
    return res;
}


/**
 * @brief multiplyF multiply two complexes
 * @return the product of a and b
 */
Complex_f multiplyF(Complex_f a, Complex_f b)
{
    Complex_f res;
    res.re = a.re*b.re -a.im*b.im;
    res.im = a.re*b.im +a.im*b.re;
    return res;
}

/**
 * @brief divF divide two complexes (a/b)
 * @return the quotient of the two
 */
Complex_f divF(Complex_f a, Complex_f b)
{
    Complex_f res, tmp;
    if(!equals(b.im, 0.0f, 1e-9))
    {
        tmp = conjugateF(b);
        a = multiplyF(a, tmp);
        b = multiplyF(b, tmp);
        res.re = a.re/b.re;
        res.im = a.im/b.re;
    }
    else{
        res.re = a.re/b.re;
        res.im = a.im/b.re;
    }
    return res;
}

/**
 * @brief invF compute the inverse of a complex
 * @param a the complex to invert
 * @return the inverted complex
 */
Complex_f invF(Complex_f a)
{
    Complex_f tmp;
    tmp.re = 1.0f;
    tmp.im = 0.0f;
    return divF(tmp, a);
}

/**
 * @brief multiplyByConjF multiply a complex by it's conjugate
 * @param a the complex to multiply
 * @return the product of the complex and it's conjugate
 */
Complex_f multiplyByConjF(Complex_f a)
{
    Complex_f res;
    res.im=0;
    res.re=2*a.re + 2*a.im;
    return res;
}

/**
 * @brief gainF compute the gain of a complex
 * @param a the complex we want to know the gain
 * @return the gain of the complex
 */
float gainF(Complex_f a)
{
    return sqrt(pow(a.im,2)+pow(a.re,2));
}

/**
 * @brief createNulComplexF create a complex 0 + 0*i
 * @return a nul complex
 */
Complex_f createNulComplexF(void)
{
    Complex_f res;
    res.re = res.im =0;
    return res;
}

/**
 * @brief createUnitaryComplexF create a complex 1 + 1*i
 * @return a complex with a real an imaginary part of 1
 */
Complex_f createUnitaryComplexF(void)
{
    Complex_f res;
    res.re = res.im =1;
    return res;
}


/**
 * @brief createRandomComplexF create a random complex
 * @return a random complex
 */
Complex_f createRandomComplexF(void)
{
    Complex_f res;
    res.re = (float) rand()/RAND_MAX;
    res.im = (float) rand()/RAND_MAX;
    return res;
}


Complex_f createComplexF(float re, float im)
{
    Complex_f res;
    res.re = re;
    res.im = im;
    return res;
}


/**
 * @brief printComplexF pint a complex float
 * @param a the complex to print
 */
void printComplexF(Complex_f a)
{
    if(a.re >=0 && a.im >=0)
        printf(" %f +%f*j", a.re, a.im);
    else if (a.im>0 && a.re <0)
        printf("%f +%f*j", a.re, a.im);
    else if (a.im<0 && a.re >0)
        printf(" %f %f*j", a.re, a.im);
    else printf("%f %f*j", a.re, a.im);
}

#ifndef COMPLEX_H
#define COMPLEX_H



typedef struct Complex_f
{
    float re;
    float im;
}Complex_f;

int equals(float a, float b, float precis);

/*****************************************************
 *
 *  Prototypes to handle simple precision complexes
 *
 * **************************************************/

Complex_f conjugateF(Complex_f c);

Complex_f addF(Complex_f a, Complex_f b);

Complex_f subF(Complex_f a, Complex_f b);

Complex_f scaleF(Complex_f a, float c);

Complex_f multiplyF(Complex_f a, Complex_f b);

Complex_f divF(Complex_f a, Complex_f b);

Complex_f invF(Complex_f a);

Complex_f multiplyByConjF(Complex_f a);

float gainF(Complex_f a);

Complex_f createNulComplexF(void);

Complex_f createUnitaryComplexF(void);

Complex_f createRandomComplexF(void);

Complex_f createComplexF(float re, float im);

void printComplexF(Complex_f a);


#endif /* COMPLEX_H*/


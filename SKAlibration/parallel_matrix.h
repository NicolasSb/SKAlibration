#ifndef PARALLEL_MATRIX
#define PARALLEL_MATRIX

#include "matrix.h"

#define NB_POLA (4)

#ifndef NB_PE_CLUSTER
#define NB_PE_CLUSTER 1
#endif


typedef struct Pe_Data
{
	Matrix_f ** mat;
	float scaleArg;
	int nb_time_freq;
	int nb_antenna;
	int nb_pola;
	int id;
	int job;
	int running_threads;
}Pe_Data;


void startPE();

void stopPE();


enum JOBS {DO_NOTHING = 0, HESSIAN, INVERT, TMP_STEP, DOT_PRODUCT, SCALE_ADD};

void update_step(void * arg);

void * pe_main(void * args);

void computeHessianF(Pe_Data *pe_data);

void giveHessianPar(Pe_Data *pe_data);

void matInvF(Pe_Data *pe_data);

void tmpStep(Pe_Data *pe_data); /*it only is the product of J0 transposed and conjugate with the visibility vector*/

void dotProductF(Pe_Data *pe_data);
 
void scaleMatAddT(Pe_Data *pe_data);


#endif // PARALLEL_MATRIX
#ifndef PARALLEL_MATRIX
#define PARALLEL_MATRIX

#include "matrix.h"


#ifndef NB_PE_CLUSTER
#define NB_PE_CLUSTER 1
#endif



void startPE();

void stopPE();


void * calibrate_data(void *args);

void update_step(void * arg);

void giveHessian(int index);

void computeJv(int index);


#endif // PARALLEL_MATRIX

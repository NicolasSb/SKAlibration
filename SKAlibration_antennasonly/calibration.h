#ifndef CALIB_H
#define CALIB_H

#include "matrix.h"


#define DATA_MATRIX_SEG_ID (600)
#define JAC_MATRIX_SEG_ID (601)
#define RES_MATRIX_SEG_ID (602)
#define FLAG_SEG_ID (610)


#ifndef NB_CLUSTER
#define NB_CLUSTER 1
#endif

#ifndef NB_PE_CLUSTER
#define NB_PE_CLUSTER 1
#endif

#ifndef TIME_DISPLAY
#define TIME_DISPLAY 0 
#endif

#ifndef TEST_MODE
#define TEST_MODE 0
#endif

#define NB_PE (NB_CLUSTER*NB_PE_CLUSTER)

#include "parallel_matrix.h"


	typedef struct IterationThread
	{
		int currentAntenna;
	}IterationThread;

	typedef struct Cluster_data
	{
		int nb_antenna, current_antenna;
		int nb_pe_required;
	    Complex_f  *dG, *jac;
	    Matrix_f *data;
	    Matrix_f *dg0;
	}Cluster_data;

int calibration ();

void testAllF();

/*****************************************************
 *
 *       Prototypes to handle multi-threading
 *
 * **************************************************/

void * cluster_main(void *);



#endif

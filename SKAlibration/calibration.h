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

	typedef struct InitThread
	{
		int nb_antenna;
		int nb_time_freq;
		int nb_polarization;
	}InitThread;

	typedef struct IterationThread
	{
		int currentAntenna;
	}IterationThread;

	typedef struct Cluster_data
	{
		int nb_antenna, nb_time_freq, nb_polarization, current_antenna;
	    Complex_f  *dG, *jac;
	    Matrix_f *data;
	    Matrix_f *H, *Hinv, *tmpMat, *dg0;
	    Pe_Data *p;
	}Cluster_data;

int calibration (int argc, char **argv);

void testAllF(int nt);

/*****************************************************
 *
 *       Prototypes to handle multi-threading
 *
 * **************************************************/

void * cluster_main(void *);



#endif
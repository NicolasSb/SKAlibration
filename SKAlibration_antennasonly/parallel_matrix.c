#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <utask.h>
#include "complex.h"
#include "matrix.h"
#include "parallel_matrix.h"
#include "calibration.h" 

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



#ifndef NB_PE_CLUSTER
#define NB_PE_CLUSTER 1
#endif



#ifdef __k1__

int running = 1;
int nb_running_thread = 1;

int nb_antenna;
Complex_f hessian[NB_PE_CLUSTER];
Complex_f temporaryScal[NB_PE_CLUSTER];
Matrix_f *J = NULL;
Matrix_f *visibilities = NULL;
Matrix_f *residualGain = NULL;
Complex_f dg0[NB_PE_CLUSTER];
int busy[NB_PE_CLUSTER] = {0};

pthread_mutex_t step_hessian;
pthread_mutex_t step_jv;
pthread_mutex_t step_scale;
pthread_barrier_t barrier;
#endif


/**
*	\brief update running flag of the PE
*/
void startPE()
{
	__builtin_k1_swu(&running, 1);
}



/**
*	\brief update running flag of the PE
*/
void stopPE()
{
	__builtin_k1_swu(&running, 0);
}

/**
*	@brief used to flush the cache memory on a cluster
*
*/
static void memBarrier()
{
	__builtin_k1_wpurge();
	__builtin_k1_fence();
}


/**
*	@brief threaad runnning on one PE that calibrate the data of one antenna
*/
void * calibrate_data(void * args __attribute__((unused)))
{
	int index = __k1_get_cpu_id();
	while(1)
	{
		while(!__builtin_k1_ldu(&busy[index]));
		giveHessian(index);
		hessian[index] = invF(hessian[index]);	
		computeJv(index);
		Complex_f res = multiplyF(hessian[index], temporaryScal[index]);
		res = scaleF(res, 0.33f);
		pthread_mutex_lock(&step_scale);
		dg0[index] = addF(scaleF(matrixGetF(residualGain, index, 0), 0.67f),res);	
		pthread_mutex_unlock(&step_scale);		
		pthread_barrier_wait(&barrier);
		__builtin_k1_sdu(&busy[index], 0);
		if(__builtin_k1_lwu(&running) == 0)
			break;
	}
	return NULL;
}

/**
*	@brief calibrate an antenna on the main PE of a cluster
*
*/
static void calibrate_data_main()
{
	int index = __k1_get_cpu_id();
	giveHessian(index);
	hessian[index] = invF(hessian[index]);
	computeJv(index);	
	Complex_f res = multiplyF(hessian[index], temporaryScal[index]);
	res = scaleF(res, 0.33f);

	pthread_mutex_lock(&step_scale);
	dg0[index] = addF(scaleF(matrixGetF(residualGain, index, 0), 0.67f),res);
	pthread_mutex_unlock(&step_scale);
}

/**
 * @brief run_test the calibration algorithms on the main PE of a cluster 
 *	that computes for an antenna
 * @param arg a structure containing matrices and indexes
 */
void update_step(void * arg)
{
	/* Init */
	Cluster_data * cluster_data = (Cluster_data *)arg;			//get the structure from the main core
	
	nb_antenna = cluster_data->nb_antenna;

	int started_pe = cluster_data->nb_pe_required;
	
	Complex_f *dG = cluster_data->dG, *J0 = cluster_data->jac;  //get the complex arrays

	Matrix_f *jac = (Matrix_f *)malloc(sizeof(Matrix_f));		//create structure to store the jacobian
	Matrix_f *residual = (Matrix_f *)malloc(sizeof(Matrix_f));	//create structure to store the residual matrix

	residual->data = dG;										//
	residual->nbRow = started_pe;
	residual->nbColumn = 1;

	jac->data = J0;
	jac->nbRow = started_pe*(nb_antenna-1);
	jac->nbColumn = 1;
	
	J = jac;	
	residualGain = residual;

	visibilities = cluster_data->data;	

	pthread_mutex_init(&step_hessian, NULL);
	pthread_mutex_init(&step_jv, NULL);
	pthread_mutex_init(&step_scale, NULL);
	pthread_barrier_init(&barrier,NULL, started_pe);
	int i;

	for(i = 0; i<started_pe; i++)
	{
		__builtin_k1_sdu(&busy[i], 1);
	}

	calibrate_data_main(NULL);	
	
	pthread_barrier_wait(&barrier);

	for(i = 0; i<started_pe; i++)
	{
		matrixSetF(residual, i, 0, dg0[i]);
	}

	//free allocated memory
	jac->data = NULL;
	J = NULL;	
	freeMatrixF(&jac);	
	residual->data =NULL;
	freeMatrixF(&residual);
	pthread_mutex_destroy(&step_hessian);
	pthread_mutex_destroy(&step_jv);
	pthread_mutex_destroy(&step_scale);
	pthread_barrier_destroy(&barrier);
	memBarrier();
}


/**
*	@brief compute the hessian term for an antenna
*/
void giveHessian(int index)
{
	int i;
	Complex_f tmp, tmpa;
	int start = index*(nb_antenna-1);
	tmpa.re = tmpa.im =0.f;
	hessian[index].re = hessian[index].im = 0;
	pthread_mutex_lock(&step_hessian);
	for(i= start; i<start+nb_antenna-1; i++)
	{
		tmp = multiplyF(matrixGetF(J, i, 0),conjugateF(matrixGetF(J,i,0)));
		tmpa = addF(tmpa, tmp);
	}
	hessian[index] = tmpa;	
	pthread_mutex_unlock(&step_hessian);
}

/**
*	@brief compute the product of the Jacobian and the visibilities  for ne antenna
*
*/
void computeJv(int index)
{
	int i;
	Complex_f tmp, tmpa;
	int start = index*(nb_antenna-1);
	tmpa.re = tmpa.im =0.f;
	temporaryScal[index].re = temporaryScal[index].im = 0;
	pthread_mutex_lock(&step_jv);
	for(i=start; i<start+nb_antenna-1; i++)
	{	
		tmp = multiplyF(conjugateF(matrixGetF(J,i,0)),matrixGetF(visibilities, i, 0));
		tmpa = addF(tmpa, tmp);
	}
	temporaryScal[index] = tmpa;
	pthread_mutex_unlock(&step_jv);
}



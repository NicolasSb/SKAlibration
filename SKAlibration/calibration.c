	#include <stdlib.h>
	#include <stdio.h>
	#include <math.h>
	#include <pthread.h>
	#include "complex.h"
	#include "matrix.h"
	#include "parallel_matrix.h"
	#include "calibration.h" 
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

	#define CHIP_FREQ ((float)__bsp_frequency/1000) // frequence en kHz du MPPA

	#define NB_POLA (4)


	

Matrix_i *lbl = NULL;

long long mymem=0;

/**
* 	verbose printf : printf("[Cluster %d PE %d]",__k1_get_cluster_id(),__k1_get_cpu_id());
*/

/**
 * @brief calibration a main function to launch the calibration
 * @param argc number of arguments
 * @param argv the value of the arguments
 * @return 0 if everything has been fine
 */
int calibration (int argc, char **argv)
{	
	srand((int)time(NULL));
	if(argc && argv[1])
	{
		testAllF(atoi(argv[1]));
		return 0;
	}
	else
	{
		testAllF(10);
		return 0;
	}
}



	/**
	*	\brief  function sent to our threads 
	*
	*	it will send the data to treat to each clusters.
	*
	*
	**/
void * cluster_main(void * arg __attribute__((unused)))
{
	// init thread

	// Init segments 
	mppa_async_segment_t paramQueue;

	mppa_async_segment_t dataMatrixSegment;
	mppa_async_segment_t jacobianMatrixSegment;
	mppa_async_segment_t residualMatrixSegment;

	mppa_async_segment_t flagSegment;

	// create THREADS
	pthread_t t[NB_PE_CLUSTER-1];

	//printf("[Cluster %d PE %d] Init segment done\n",__k1_get_cluster_id(),__k1_get_cpu_id());
	
	// create queue reception data
	void * bufferQueue = malloc(sizeof(InitThread));
	mymem +=2*sizeof(InitThread);
	if (bufferQueue == NULL){
		printf("[Cluster %d PE %d] Alloc failed\n",__k1_get_cluster_id(),__k1_get_cpu_id());
		while(1);
	}
	//printf("[Cluster %d PE %d] creating queue : %d \n",__k1_get_cluster_id(),__k1_get_cpu_id(), __k1_get_cluster_id());
	if (mppa_async_segment_create(&paramQueue,__k1_get_cluster_id(), bufferQueue, sizeof(InitThread), MPPA_ASYNC_SEGMENT_FLAG_QUEUE0, 0, NULL) != 0) //note MPPA...QUEUE0 ID problem
	{
		printf("[Cluster %d PE %d] segment %d failed to create\n", __k1_get_cluster_id(), __k1_get_cpu_id(), __k1_get_cluster_id());
		while(1);
	}
	//printf("[Cluster %d PE %d] SUCCESS on creating queue : %d \n",__k1_get_cluster_id(),__k1_get_cpu_id(),__k1_get_cluster_id());


	// clone matrix RDMA segment
	
	mppa_async_segment_clone(&dataMatrixSegment,DATA_MATRIX_SEG_ID,NULL, 0, NULL);

	mppa_async_segment_clone(&jacobianMatrixSegment,JAC_MATRIX_SEG_ID,NULL, 0, NULL);
	mppa_async_segment_clone(&residualMatrixSegment,RES_MATRIX_SEG_ID,NULL, 0, NULL);
	

	// clone flag RDMA segment 
	mppa_async_segment_clone(&flagSegment,FLAG_SEG_ID,NULL, 0, NULL); 

	// initialization
	// dequeue
	void *temp = NULL;
	mppa_async_dequeue(&paramQueue,sizeof(InitThread), &temp, NULL);
	InitThread *initValue = (InitThread*) temp;
	int nb_antenna, nb_time_freq, nb_polarization;
	nb_antenna = initValue->nb_antenna;
	nb_time_freq = initValue->nb_time_freq;
	nb_polarization = initValue->nb_polarization;

	// Allocate space for data
	Complex_f * jac = calloc(nb_polarization*nb_time_freq*(nb_antenna-1)*nb_polarization, sizeof(Complex_f));
	Matrix_f * data = matrixAllocF(nb_polarization*(nb_antenna-1)*nb_time_freq,1);
	Complex_f * res = calloc(nb_polarization, sizeof(Complex_f));

	//allocate space for temp matrices
	Matrix_f *dg0 = matrixAllocF(nb_polarization,1);
	Matrix_f *H = matrixAllocF(nb_polarization, nb_polarization);
	Matrix_f *Hinv = matrixAllocF(nb_polarization, nb_polarization);
	Matrix_f *tmpMat = matrixAllocF(nb_polarization,1);
	

	//allocate pe data pointer structure
	Pe_Data * p = (Pe_Data *) malloc(sizeof(Pe_Data));
	p->mat = (Matrix_f **) calloc(3, sizeof(Matrix_f *));

	mymem +=nb_polarization*nb_time_freq*(nb_antenna-1)*nb_polarization*sizeof(Complex_f);				//jac mem
	mymem +=nb_polarization*(nb_antenna-1)*nb_time_freq*sizeof(Complex_f)+sizeof(Matrix_f);	//data mem
	mymem +=nb_polarization*sizeof(Complex_f);															//res mem
	mymem +=nb_polarization*sizeof(Complex_f)+sizeof(Matrix_f);								//dg0 mem
	mymem +=nb_polarization*nb_polarization*sizeof(Complex_f)+sizeof(Matrix_f);							//H mem
	mymem +=nb_polarization*nb_polarization*sizeof(Complex_f)+sizeof(Matrix_f);							//Hinv mem
	mymem +=nb_polarization*sizeof(Complex_f)+sizeof(Matrix_f);											//tmpMat mem
	mymem+=3*sizeof(Matrix_f*)+sizeof(Pe_Data);
	

	IterationThread *iterThread = NULL;

	//creating data structure
	int currentAntenna=0;
	Cluster_data cData;
	
	cData.data = data;
	cData.H = H;
	cData.Hinv = Hinv;
	cData.tmpMat = tmpMat;
	cData.dg0 = dg0;
	cData.nb_antenna = nb_antenna;
	cData.nb_time_freq = nb_time_freq;
	cData.nb_polarization = nb_polarization;
	cData.p = p;

	//Start the PE ant initialize a stop flag
	int i;
	for(i=1; i<NB_PE_CLUSTER; i++)
    {
    	pthread_create(&t[i-1], NULL, pe_main, (void*) p);
    }

    startPE();
	while(1) {

		// dequeue (pop)
		mppa_async_dequeue(&paramQueue,sizeof(IterationThread), &temp, NULL);

		iterThread = (IterationThread*) temp;
		
		currentAntenna = iterThread->currentAntenna;

		cData.current_antenna = currentAntenna;

		if (currentAntenna == -1)
		{
			//printf("[Cluster %d PE %d] exited with success\n",__k1_get_cluster_id(),__k1_get_cpu_id());

			break;
		} 

		// RMDA_get
		// get the Data matrix
		off64_t offset = currentAntenna*nb_time_freq*NB_POLA*(nb_antenna-1)*sizeof(Complex_f);
	
		if(mppa_async_get(data->data, &dataMatrixSegment, offset, nb_polarization*(nb_antenna-1)*nb_time_freq*sizeof(Complex_f), NULL) != 0)
		{
			printf("[Cluster %d] failed to retrieve buffer address %p size %d\n", __k1_get_cluster_id(), data->data,
										nb_polarization*(nb_antenna-1)*nb_time_freq*sizeof(Complex_f));
			while(1);
		}


		// RDMA_get data  
		// the offsets and sizes are not sure for now

		offset = currentAntenna*nb_polarization*(nb_antenna-1)*nb_time_freq*nb_polarization*sizeof(Complex_f);
		// printf("[Cluster %d PE %d] offset jac %lld\n",__k1_get_cluster_id(),__k1_get_cpu_id(),offset);

		if(mppa_async_get(jac, &jacobianMatrixSegment,offset, nb_polarization*(nb_antenna-1)*nb_time_freq*nb_polarization*sizeof(Complex_f), NULL) != 0)
		{
			printf("[Cluster %d] failed to retrieve buffer address %p size %d\n", __k1_get_cluster_id(), jac, 
									nb_polarization*nb_antenna*(nb_antenna-1)*nb_time_freq*sizeof(Complex_f));
			while(1);
		}
		//printf("[c %d] accessing %d data with offset : %d \n", __k1_get_cluster_id(),nb_polarization*(nb_antenna-1)*nb_time_freq*nb_polarization, (int)offset);

		offset = currentAntenna*nb_polarization*sizeof(Complex_f);
		// printf("[Cluster %d PE %d] offset res %lld\n",__k1_get_cluster_id(),__k1_get_cpu_id(),offset);

		if(mppa_async_get(res, &residualMatrixSegment, offset, nb_polarization*sizeof(Complex_f), NULL) != 0)
		{
			printf("[Cluster %d] failed to retrieve buffer address %p size %d\n", __k1_get_cluster_id(), res, nb_polarization*sizeof(Complex_f));
			while(1);
		}


		// prepare structure for the computing
		cData.dG = res;
		cData.jac = jac;


		//printf("[Cluster %d PE %d] Work start on antenna %d \n",__k1_get_cluster_id(),__k1_get_cpu_id(), currentAntenna);
		update_step((void *)&cData);
		//printf("[Cluster %d PE %d] Work done\n",__k1_get_cluster_id(),__k1_get_cpu_id());

		// RDMA_put data
		offset = currentAntenna*nb_polarization*sizeof(Complex_f);
		//printf("[c %d] puting %d data with offset : %d \n", __k1_get_cluster_id(),nb_polarization, currentAntenna*nb_polarization);
		mppa_async_put(cData.dG, &residualMatrixSegment,offset,nb_polarization*sizeof(Complex_f), NULL);

		// (update flag) ->one more antenna has been treated
		mppa_async_poke(&flagSegment, __k1_get_cluster_id()*sizeof(long long), 0);
		mppa_async_postadd(&flagSegment, NB_CLUSTER*sizeof(long long), 1);
	}
	//printf("Total memory allocated in one cluster : %lld o\n\n\n",mymem);


	stopPE();

	for(i=1; i<NB_PE_CLUSTER; i++)
    {
    	pthread_join(t[i-1], NULL);
    }


	// fin thread
	freeMatrixF(data);
	free(jac);
	free(res);

	
	freeMatrixF(H);
	freeMatrixF(Hinv);
	freeMatrixF(tmpMat);
	freeMatrixF(dg0);

	// destroy queue
	mppa_async_segment_destroy(&paramQueue);

	//free the buffer
	free(bufferQueue);

	return 0;
}


	/**
	 * @brief testAll small implementation to see if the computed matrix converges on the true one
	 * @param na the number of antennas taken in account (should be at least 3)
	 */
void testAllF(int nb_time_freq)
{
	#if TEST_MODE == 0
	/*Create random calibration gain vector (the True gain vector to be found)*/
    unsigned int nb_antenna = 64;
	int nb_polarization = 4;
    Matrix_f * observation = createRandomMatrixF(nb_antenna*nb_polarization,1); /*for each time frequency point we have a vector of na * 4 elements containing
	                                                             for each antennas it's four polarizations */
    Matrix_f *dG = createRandomMatrixF(nb_antenna*nb_polarization,1);
    int nbIter = 10;
    #else
    /*Create random calibration gain vector (the True gain vector to be found)*/
    unsigned int nb_antenna = 3;
    int nb_polarization= 4;
    nb_time_freq =1;
    Matrix_f * observation = createTestVectorF(nb_antenna*nb_polarization,1);
    Matrix_f *dG = createTestVectorF(nb_antenna*nb_polarization,1);
    int nbIter = 10;
    #endif
    #if TIME_DISPLAY == 1
    uint64_t time_spent = 0;
    uint64_t program_start = __k1_read_dsu_timestamp();
    uint64_t time_per_iteration = 0;
    #endif
	//Useful matrices declaration
	Matrix_f * sky_coherency = createRandomMatrixF(nb_time_freq*nb_antenna*(nb_antenna-1),nb_polarization); // used ? 
	Matrix_f * data = NULL;
	Matrix_f * J0 = NULL;
    int i; /*nd is useful when we use more than one directions*/
	Matrix_i *a0 = matrixAllocI(nb_antenna, nb_antenna), *a1 = matrixAllocI(nb_antenna, nb_antenna);
	Matrix_i *b0 = NULL, *b1 = NULL;



	unsigned int j;

	//declaration of MPPA segments
	mppa_async_segment_t paramQueue[NB_CLUSTER];
	mppa_async_segment_t dataMatrixSegment;
	mppa_async_segment_t jacobianMatrixSegment;
	mppa_async_segment_t residualMatrixSegment;
	mppa_async_segment_t flagSegment;



	lbl = matrixAllocI(nb_antenna*nb_antenna-nb_antenna,2);
	create_grid(a0, a1);
	b0 = make1DvectorFromMatrix(a0);

	b1 = make1DvectorFromMatrix(a1);

	freeMatrixI(a0);
	freeMatrixI(a1);

	compute_LBL_Matrix(lbl, b0, b1);
	freeMatrixI(b0);
	freeMatrixI(b1);

	data = matrixAllocF(nb_polarization*nb_antenna*(nb_antenna-1)*nb_time_freq, 1); // Jacobian * observation 
	J0 = matrixAllocF(nb_antenna*(nb_antenna-1)*nb_time_freq*nb_polarization, nb_polarization);

    /* we take the first antenna as a reference*/
	matNormF(observation);

	/*Create a random Guess */
	addMatrixF(dG, observation);

	// Declaration of MPPA buffers 
	long long flag[NB_CLUSTER+1] = {0};
	

	giveJfF(observation, sky_coherency, lbl, J0, nb_time_freq, nb_antenna);

	    /*compute the data corresponding to the visibilities*/
	compute_data_vectorF(J0, observation, data, nb_antenna, nb_time_freq);

	    /*Norm the guessed vector*/
	matNormF(dG);

	/*#ifdef TEST_MODE
	printf(" \nThe true matrix we want to have in final\n");
	printMatrixF(observation);
	printf("\nThe initial guess\n");
	printMatrixF(dG);

	printf("\n\n");
	#endif*/


	// clone emission queue
	for (i = 0; i < NB_CLUSTER; i++){
			//printf("[Cluster %d PE %d] cloning queue : %d \n",__k1_get_cluster_id(),__k1_get_cpu_id(), i);
			mppa_async_segment_clone(&paramQueue[i], i, NULL, 0, NULL);
			//printf("[Cluster %d PE %d] SUCCESS on clone queue : %d \n",__k1_get_cluster_id(),__k1_get_cpu_id(), i);
	}


	// create RDMA segment for each matrix
		// create RDMA seg. for the visibility mat.
	
	mppa_async_segment_create(&dataMatrixSegment, DATA_MATRIX_SEG_ID, data->data,
							nb_polarization*nb_antenna*(nb_antenna-1)*nb_time_freq*(sizeof(Complex_f)),0,0, NULL);
	mppa_async_segment_create(&jacobianMatrixSegment, JAC_MATRIX_SEG_ID, J0->data,
							nb_antenna*(nb_antenna-1)*nb_time_freq*nb_polarization*nb_polarization*sizeof(Complex_f) , 0,0, NULL);
	mppa_async_segment_create(&residualMatrixSegment, RES_MATRIX_SEG_ID, dG->data,
							 nb_antenna*nb_polarization*sizeof(Complex_f), 0,0, NULL);

	mppa_async_segment_create(&flagSegment, FLAG_SEG_ID, &flag, (NB_CLUSTER+1)*sizeof(long long), 0,0, NULL);
	// main loop

	InitThread initValue;
	initValue.nb_antenna = nb_antenna;
	initValue.nb_time_freq = nb_time_freq;
	initValue.nb_polarization = nb_polarization;

	// send the init values to fetch all the useful blocks
	for (i = 0; i < NB_CLUSTER ; i++){
		mppa_async_enqueue(&paramQueue[i], &initValue, sizeof(InitThread), 0, NULL);
	}

	IterationThread iterThread;



	for(i=0; i< nbIter; i++)
	{
		// Update the Jacobian 
		giveJfF(dG, sky_coherency, lbl, J0, nb_time_freq, nb_antenna); //TODO update it in the clusters
		
		// full memory barrier
		__builtin_k1_wpurge();
		__builtin_k1_fence();

		unsigned int data_sent = 0;
		#if TIME_DISPLAY == 1
		uint64_t start_iteration = __k1_read_dsu_timestamp();
		#endif
		while(__builtin_k1_ldu(&flag[NB_CLUSTER]) < nb_antenna)
		{
			for (j = 0; j < NB_CLUSTER; j++){

				if(flag[j]>0)continue;
				if(flag[j]==0)
				{
					flag[j]=1;
					iterThread.currentAntenna = data_sent;

					// send data
					mppa_async_enqueue(&paramQueue[j],&iterThread, sizeof(IterationThread),0, NULL);
					
					data_sent++;
				}
				if (data_sent >= nb_antenna) break; //no more jobs to send
			}

			//synchronization step 
			while(__builtin_k1_ldu(&flag[NB_CLUSTER]) != data_sent);

		}
		#if TIME_DISPLAY == 1
		uint64_t stop_iteration = __k1_read_dsu_timestamp();
		time_per_iteration += stop_iteration - start_iteration;
		#endif

		// clean stopFlag
		__builtin_k1_sdu(&flag[NB_CLUSTER], 0);

		// Invalidate cache memory -> force reading in MEM
		__builtin_k1_dinval();

		matNormF(dG); //TODO Parallelize on IOs
		//printf("[Cluster %d PE %d] end of iteration : %d \n",__k1_get_cluster_id(),__k1_get_cpu_id(), i);
	}
	#if TIME_DISPLAY == 1
	uint64_t program_stop = __k1_read_dsu_timestamp();
	time_spent = program_stop - program_start;

	printf("iteration time : %f ms\n", (float)((time_per_iteration)/(CHIP_FREQ*nbIter)));
	printf("execution time : %f s\n", (float)((time_spent)/(CHIP_FREQ*1000)));
	#endif
/*
	printf("NB_CLUSTER %d, NB_PE_CLUSTER %d, NB_PE %d\n",NB_CLUSTER, NB_PE_CLUSTER, NB_PE);
*/
	for (j = 0; j < NB_CLUSTER; j++)
	{
		iterThread.currentAntenna = -1;
		//printf("sending end flag to pe %d\n",j);
		mppa_async_enqueue(&paramQueue[j],&iterThread, sizeof(IterationThread),0, NULL);
	}

	#if TEST_MODE == 1

	subMatrixF(observation,dG);

	printf("\n\n if the guessed (and then computed) matrix - the original gain tend to be 0, the convergence is ok \n");
	printMatrixF(observation);

	
	printf("\nThe result should be \n\n");
	printf(" 0.000000 +0.000000*j\n 0.001944 +0.030689*j\n-0.005808 +0.061909*j\n-0.006643 -0.125064*j\n");	
	printf("-0.169734 -0.211088*j\n-0.158894 -0.192926*j\n-0.145033 -0.172928*j\n-0.162924 -0.075779*j\n");
 	printf(" 0.034903 -0.082082*j\n 0.043399 -0.104007*j\n 0.046430 -0.125429*j\n-0.003213 +0.016504*j\n");
	#endif
	
	// destroy RDMA
	mppa_async_segment_destroy(&dataMatrixSegment);
	mppa_async_segment_destroy(&jacobianMatrixSegment);
	mppa_async_segment_destroy(&residualMatrixSegment);

	mppa_async_segment_destroy(&flagSegment);


	freeMatrixF(observation);
	freeMatrixF(data);
	freeMatrixF(dG);
	freeMatrixF(J0);
	freeMatrixI(lbl);
}

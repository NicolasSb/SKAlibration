#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <utask.h>
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



#ifndef NB_PE_CLUSTER
#define NB_PE_CLUSTER 1
#endif



#ifdef __k1__

int running = 1;
long long job [NB_PE_CLUSTER-1];
int nb_running_thread = 1;

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
*	\brief update the current number of threads
*/
void setNbThread(int nb)
{
	int min = nb > NB_PE_CLUSTER ? NB_PE_CLUSTER : nb;
	nb_running_thread = min;
}





static void updateJobFlag(int index, long long new_job)
{
	__builtin_k1_sdu(&job[index], new_job);
}

static long long getCurrentJob(int index)
{
	return __builtin_k1_ldu(&job[index]);
}


/**
* 	\brief update all the PE flags
*/
static void setJobs(long long new_job, int nb_thread)
{
	int i = 0;
	int min = nb_thread > NB_PE_CLUSTER ? NB_PE_CLUSTER : nb_thread;
	for(i=0; i< min-1; i++)
		__builtin_k1_sdu(&job[i], new_job);
}


static int syncronizePE(int nbThread)
{
	int i;
	long long current_job = 0;

	for (i = 0; i< nbThread-1; i++)
	{
		current_job += __builtin_k1_ldu(&job[i]);
	}

	if(current_job != 0)
		return 0;
	return 1;
}	


/**
*	\brief update the data structure that will be sent to the PE
*
*/
static void updatePEdata(Pe_Data *pe_data, Matrix_f* mat0, Matrix_f * mat1, Matrix_f *mat2, float f, int nb_thread)
{
	if(pe_data)
	{
		if(mat0)
		{
			pe_data->mat[0] = mat0;
		}
		if(mat1)
		{
			pe_data->mat[1] = mat1;
		}
		if(mat2)
		{
			pe_data->mat[2] = mat2;
		}
		pe_data->scaleArg = f;
		pe_data->running_threads = nb_thread;

		// full memory barrier
		__builtin_k1_wpurge();
		__builtin_k1_fence();
	}
}


static void masterGiveJobs(int job_id, Pe_Data * pe_data)
{
	__builtin_k1_dinval();
	switch(job_id)
	{
			case HESSIAN:
			{
				setJobs((long long) HESSIAN, NB_PE_CLUSTER);
				//start main PE computing
				giveHessianPar(pe_data);
				break;
			}
			case INVERT:
			{
				//setJobs((long long) INVERT, 1);	
				matInvF(pe_data);
				break;
			}
			case TMP_STEP:
			{
				setJobs((long long) TMP_STEP, 4);
				tmpStep(pe_data); /*it only is the product of J0 transposed and conjugate with the visibility vector*/
				break;
			}
			case DOT_PRODUCT:
			{
				setJobs((long long) DOT_PRODUCT, 4);
				dotProductF(pe_data);
				break;
			}
			case SCALE_ADD:
			{
				setJobs((long long) SCALE_ADD, 4);
				scaleMatAddT(pe_data);
				break;
			}
			default:
				printf("unsupported Job function\n");
				break;
	}
	while(!syncronizePE(nb_running_thread));
}

/**
 * @brief run_test the calibration algorithms that computes for an antenna
 * @param arg a structure containing matrices and indexes
 */
void update_step(void * arg)
{
	/* Init */
	Cluster_data * cluster_data = (Cluster_data *)arg;
	
	Complex_f *dG = cluster_data->dG, *J0 = cluster_data->jac;

	float f = 0.5f;

	Pe_Data *pe_data = cluster_data->p;

	Matrix_f *jac = (Matrix_f *)malloc(sizeof(Matrix_f));
	Matrix_f *residual = (Matrix_f *)malloc(sizeof(Matrix_f));

	pe_data->id = cluster_data->current_antenna;

	residual->data = dG;
	residual->nbRow = cluster_data->nb_polarization;
	residual->nbColumn = 1;

	jac->data = J0;
	jac->nbRow = cluster_data->nb_polarization*cluster_data->nb_time_freq*(cluster_data->nb_antenna-1);
	jac->nbColumn = cluster_data->nb_polarization;

	pe_data->nb_antenna = cluster_data->nb_antenna;
	pe_data->nb_pola = cluster_data->nb_polarization;
	pe_data->nb_time_freq = cluster_data->nb_time_freq;

	/* Treatement */

	setJobs((long long) DO_NOTHING, NB_PE_CLUSTER);

	setNbThread(NB_PE_CLUSTER);
	updatePEdata(pe_data, jac, cluster_data->H,  NULL, 0.f, nb_running_thread);
	masterGiveJobs(HESSIAN, pe_data);

	setNbThread(1);
	updatePEdata(pe_data, NULL, NULL, cluster_data->Hinv, 0.f, nb_running_thread);
	masterGiveJobs(INVERT, pe_data);

	setNbThread(4);
	updatePEdata(pe_data, NULL, cluster_data->data,  cluster_data->tmpMat, 0.f, nb_running_thread);	
	masterGiveJobs(TMP_STEP, pe_data);

	setNbThread(4);
	updatePEdata(pe_data, cluster_data->Hinv, cluster_data->dg0,  NULL, 1/(f+1), nb_running_thread);
	masterGiveJobs(DOT_PRODUCT, pe_data);

	setNbThread(4);
	updatePEdata(pe_data, residual, NULL,  NULL, f/(f+1), nb_running_thread);
	masterGiveJobs(SCALE_ADD, pe_data);

	//printMatrixF( cluster_data->H);

	//free allocated memory
	free(jac);
	free(residual);
}



void * pe_main(void * args)
{
	Pe_Data *pe_data = (Pe_Data *) args;
	while (running)
	{
		int cpu_id = __k1_get_cpu_id()-1;
		long long current_job = getCurrentJob(cpu_id);
		switch(current_job)
		{
			case DO_NOTHING:
			{
				break;
			}
			case HESSIAN:
			{
				__builtin_k1_dinval();
				giveHessianPar(pe_data);
				updateJobFlag(cpu_id, (long long) DO_NOTHING);
				break;
			}
			case INVERT:
			{
				//__builtin_k1_dinval();
				updateJobFlag(cpu_id, (long long) DO_NOTHING);
				break;
			}
			case TMP_STEP:
			{
				__builtin_k1_dinval();
				tmpStep(pe_data);
				updateJobFlag(cpu_id, (long long) DO_NOTHING);
				break;
			}
			case DOT_PRODUCT:
			{
				__builtin_k1_dinval();
				dotProductF(pe_data);
				updateJobFlag(cpu_id, (long long) DO_NOTHING);
				break;
			}
			case SCALE_ADD:
			{
				__builtin_k1_dinval();
				scaleMatAddT(pe_data);
				updateJobFlag(cpu_id, (long long) DO_NOTHING);
				break;
			}
			default:
				printf("unsupported Job function\n");
				break;
		}
	}
	return 0;
}


/**
     * @brief computeHessianF computes the Hessian of a function
     * @param L the half Jacobian of a matrix
     * @param Lb an empty matrix (allocated) to store the conjugate transpose of the half jacobian
     * @param H an empty matrix(initialized) to store the Hessian
     *
     *
     * As we know that we've got the half jacobian and that the hessian is the product of it and it's conjugate transpose
     * we can avoid a huge matrix product and simply multiply 4*4 blocks for each antenna
     *
     */
void giveHessianPar(Pe_Data *pe_data)
{
	if(pe_data)
    {
    	// /!\ do not modify the indexes
    	Matrix_f *jac = pe_data->mat[0];
    	Matrix_f *H = pe_data->mat[1];
        int l, i;
        int id = __k1_get_cpu_id();
        int nb = pe_data->nb_pola*pe_data->nb_pola;
        int jacHeight = (pe_data->nb_pola*pe_data->nb_time_freq*(pe_data->nb_antenna-1));
        int threads = pe_data->running_threads;
        if (nb%threads==0)
        	nb /= threads;
    	else
        	nb = nb/threads +1;

        int offset = id*nb;
        Complex_f tmp, tmpa, tmpb;


        for (i=offset; i<offset+nb; i++)
       	{
       		if(i>= (pe_data->nb_pola*pe_data->nb_pola))
       		{
       			break;
       		}
       		tmp.re=0; tmp.im =0;
			for(l=0; l<jacHeight; l++)
			{
				tmpa = matrixGetF(jac, i/pe_data->nb_pola, l);
				tmpb = matrixGetF(jac, i%pe_data->nb_pola, l);
				tmpa = multiplyF(tmpb, conjugateF(tmpa));
				tmp = addF(tmp, tmpa);
			}
			matrixSetF(H, 0, i, tmp);
		}
		__builtin_k1_wpurge();
    	__builtin_k1_fence();
    }
}




/**
     * @brief matInvF to invert any square matrix
     * @param a
     * @param ainv
     */
void matInvF(Pe_Data *pe_data)
{
    /*if(pe_data)
    {
    	// /!\ do not modify the indexes
    	Matrix_f *a = pe_data->mat[1];
    	Matrix_f *ainv = pe_data->mat[2];
        Complex_f f, tmp;
        int res = 0;
        int i, j;
		while(i++ < __k1_get_cluster_id()*10000000);
		printf("Hessian \n");
		printMatrixF(a);
        createIdentityF(ainv);
        for(i =0; i<(int)a->nbColumn;i++)
        {
            tmp = matrixGetF(a,i,i);
            if( !equals(tmp.re, 0.0f, 1e-9) || !equals(tmp.re, 0.0f, 1e-9))
            {
                f = invF(tmp);
                scaleLineF(ainv, i, f);
            }
            for(j=0; j<(int)a->nbRow; j++)
            {
                if(i!=j)
                {
                    f = matrixGetF(a,j,i);
                    if(!equals(f.re, 0.0f, 1e-9) || !equals(f.re, 0.0f, 1e-9))
                    {
                        res = subXLinesF(ainv,j,i,f);
                        if (res==1)
                        {
                            createIdentityF(ainv);
                        }
                    }
                }
            }
        }
		while(i++ < __k1_get_cluster_id()*10000000);
		printf("Hessian inv \n");
		printMatrixF(ainv);
    }*/
   		Matrix_f *a = pe_data->mat[1];
    	Matrix_f *ainv = pe_data->mat[2];
        Complex_f f, tmp;
        int i;
        createIdentityF(ainv);
        for(i =0; i<(int)a->nbColumn;i++)
        {
            tmp = matrixGetF(a,i,i);
            if( !equals(tmp.re, 0.0f, 1e-9) || !equals(tmp.re, 0.0f, 1e-9))
            {
                f = invF(tmp);
                matrixSetF(ainv, i,i,f);
            }
        }
}


/**
     * @brief tmpStep used for the update step, this functions is used to avoid doing another full matrix multiplication
     * @param J0 jacobian matrix
     * @param ytmp visibilities
     * @param tmpMat to store the result
     */
void tmpStep(Pe_Data *pe_data)
{
	// /!\ do not modify the indexes
    Matrix_f *jac = pe_data->mat[0];
    Matrix_f *ytmp = pe_data->mat[1];
    Matrix_f *tmpMat = pe_data->mat[2];
    int i, j;
    int nb_antenna = pe_data->nb_antenna;
    int nb_time_freq = pe_data->nb_time_freq;
    int id = __k1_get_cpu_id();
    int nb_elements = nb_time_freq*NB_POLA*(nb_antenna-1);
	int threads = pe_data->running_threads; 
    Complex_f tmp, tmpa;

    int nb = pe_data->nb_pola;
    if (nb%threads==0)
    	nb /= threads;
    else
    	nb = nb/threads +1;

    int offset = id*nb; 
    for(i=offset; i< offset + nb ; i++)
    {
    	if(i >= pe_data->nb_pola)
	    		break;
    	tmp.re = tmp.im = 0;
	    for(j=0; j<nb_elements; j++)
	    {
		    tmpa = conjugateF(matrixGetF(jac,i,j));
		    tmpa = multiplyF(matrixGetF(ytmp,j,0),tmpa); // because the polarisations are consecutive
		    tmp = addF(tmpa, tmp);
	    }
	    //printf("setting id : %d  with value : %f %fi at adress : %p\n", i, tmp.re, tmp.im, tmpMat);
	    tmpMat->data[i] = tmp;

	    //matrixSetF(tmpMat, 0,i,tmp);
	}
	__builtin_k1_wpurge();
   	__builtin_k1_fence();
}


/**
     * @brief dotProductF compute C = A*B
     * @param A complex matrix
     * @param B complex matrix
     * @param C complex matrix to store the result
     */
void dotProductF(Pe_Data *pe_data)
{
    if(pe_data)
    {   
    	// /!\ do not modify the indexes
    	Matrix_f *A = pe_data->mat[0];	
    	Matrix_f *B = pe_data->mat[2];
    	Matrix_f *C = pe_data->mat[1];
    	float scale = pe_data->scaleArg;
        int i, j;
        int threads = pe_data->running_threads; 
        int id = __k1_get_cpu_id();
        Complex_f tmp, tmpa;

        int nb = pe_data->nb_pola;
   	 	if (nb%threads==0)
    		nb /= threads;
    	else
    		nb = nb/threads +1;

    	int offset = id*nb; 
    	for(i=offset; i<nb+offset; i++)
    	{
    		if(i>=pe_data->nb_pola)
	        		break;
    		tmp.re=tmp.im =0;
	        for(j=0; j<pe_data->nb_pola; j++)
	        {
		        tmpa = matrixGetF(A, i, j);
		        tmpa = multiplyF(tmpa,matrixGetF(B, j, 0));
		        tmp = addF(tmp, tmpa);
	        }
	        tmp = scaleF(tmp, scale);
	        //printf("[c : %d, pe : %d] setting index %d at value %f %fi \n", __k1_get_cluster_id(), id, i, tmp.re, tmp.im);
	        matrixSetF(C, 0, i, tmp);
	    }
	    __builtin_k1_wpurge();
    	__builtin_k1_fence();
	}
}


/**
     * @brief scaleMatrixF multiply each data by a constant (complex) and add the matrix BT
     * @param a the matrix to scale
     * @param c the complex
     * @todo Parallel programming
     */
void scaleMatAddT(Pe_Data *pe_data)
{
    if(pe_data)
    {
    	Matrix_f *a = pe_data->mat[0];
    	Matrix_f *b = pe_data->mat[1];
    	float f = pe_data->scaleArg;
        int i;
        Complex_f tmp, tmpadd;
        int threads = pe_data->running_threads; 
        int id = __k1_get_cpu_id();

        int nb = pe_data->nb_pola;
   	 	if (nb%threads==0)
    		nb /= threads;
    	else
    		nb = nb/threads +1;

    	int offset = id*nb; 
        for (i=offset; i<nb+offset; i++)
        {
        	if(i>=NB_POLA)
        		break;
            tmp = matrixGetF(a, i, 0);
            tmpadd = matrixGetF(b, i, 0);
            matrixSetF(a, i, 0, addF(scaleF(tmp, f), tmpadd));
        }
        __builtin_k1_wpurge();
    	__builtin_k1_fence();
    }
}

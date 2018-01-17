/*
 * Copyright (C) 2015-2017 Kalray SA. All rights reserved.
 */

#include "mOS_common_types_c.h"

#include "mOS_vcore_u.h"
#include "mOS_segment_manager_u.h"
#include "stdlib.h"
#include "stdio.h"
#include "vbsp.h"
#include "utask.h"
#include <HAL/hal/core/cpu.h>
#include <HAL/hal/core/mp.h>
#include <math.h>
#include <stdlib.h>
#include <mppa_power.h>
#include <mppa_async.h>
#include <mppa_remote.h>
#include <vbsp.h>
#include <string.h>
#include <assert.h>

#include "calibration.h"

/* Main executed on PE0 */
int main(int argc __attribute__((unused)), const char *argv[] __attribute__((unused)))
{
	/**************************************************************************************
	*
	*						  CLUSTER INIT SECTION DO NOT MODIFY
	*
	**************************************************************************************/
	mppa_rpc_client_init();
	mppa_async_init();
	mppa_remote_client_init();
	mppa_rpc_barrier_all(); /* synchronize all booted compute cluster */
	/* runtime library initialisation please do not do anthing before this */

	/**************************************************************************************
	*
	*							 END OF CLUSTER INIT SECTION
	*
	**************************************************************************************/


	//printf("[Cluster %d PE %d] Cluster Start\n",__k1_get_cluster_id(),__k1_get_cpu_id());




	cluster_main(NULL);


	/**************************************************************************************
	*
	*						  CLUSTER CLOSE SECTION DO NOT MODIFY
	*
	**************************************************************************************/


	mppa_rpc_barrier_all(); /* synchronize all booted compute cluster */
	//printf("[Cluster %d] Goodbye\n", cid);
	mppa_rpc_barrier_all(); /* synchronize all booted compute cluster */
	mppa_async_final();
	//while(1); // for now....
	return 0;
}

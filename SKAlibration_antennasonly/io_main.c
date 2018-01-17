/**
 * Copyright (C) 2015-2016 Kalray SA.
 *
 * All rights reserved.
 */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <mppa_power.h>
#include <mppa_routing.h>
#include <mppa_remote.h>
#include <mppa_rpc.h>
#include <mppa_async.h>
#include <assert.h>
#include <utask.h>
#include <HAL/hal/board/boot_args.h>


#include "calibration.h"
#include <math.h>


#define CHIP_FREQ ((float)__bsp_frequency/1000) // frequence en Hz du MPPA

static void*
mppa_rpc_task(void *arg __attribute__ ((unused)) )
{
	while(1) mppa_rpc_server();
	return NULL;
}

int main() {
	int i;
	mppadesc_t pcie_fd = 0;
	if (__k1_spawn_type() == __MPPA_PCI_SPAWN) {
		pcie_fd = pcie_open();
		pcie_queue_init(pcie_fd);
		//pcie_register_console(pcie_fd, stdin, stdout);
	}
	mppa_rpc_server_init(1/*rm where to run the server*/, 0/* ddr indirection*/, NB_CLUSTER);
	mppa_async_server_init();
	mppa_remote_server_init(pcie_fd, NB_CLUSTER);
	mppa_remote_server_enable_scall();
	
	for(i=0;i<NB_CLUSTER;i++)
	{
		if (mppa_power_base_spawn(i, "cluster_bin", NULL, NULL, MPPA_POWER_SHUFFLING_ENABLED) == -1){
			printf("# [IODDR0] Fail to Spawn cluster %d\n", i);
		}
	}
	pthread_t thread_rpc;
	//pthread_create(&thread_rpc, NULL, (void*) mppa_rpc_server_start, NULL);
	pthread_create(&thread_rpc, NULL, (void*) mppa_rpc_task, NULL);

	//printf("IO%d Hello\n", __k1_get_cpu_id());
	
	calibration(0, NULL);

	if(__k1_spawn_type() == __MPPA_PCI_SPAWN){
		int status = 0;
		int remote_status;
		pcie_queue_barrier(pcie_fd, status, &remote_status);
		pcie_queue_barrier(pcie_fd, status, &remote_status);
		mppa_remote_server_disable_scall();
		mppa_remote_server_exit();
	}
	/* Send an exit message on pcie interface */
	if (__k1_spawn_type() == __MPPA_PCI_SPAWN) {
		//int status = 0;
		//int remote_status;
		//pcie_unregister_console(pcie_fd);
		pcie_queue_exit(pcie_fd, 0, NULL);
	}
	#if 0
	uint64_t start = __k1_read_dsu_timestamp(); 
	uint64_t end = __k1_read_dsu_timestamp(); 
	#define CHIP_FREQ ((float)__bsp_frequency/1000) // frequence en Hz du MPPA
	float ms = (float)(end-start)/CHIP_FREQ;
	#endif
	return 0;
}


/*
 * Copyright (C) 2015-2017 Kalray SA. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <pcie.h>
#include <unistd.h>

int main(int argc, char **argv)
{
	mppadesc_t fd;
	int mppa_ret;

	if(!(fd = pcie_open_device(0))) 
		exit(-1);

	if (argc >= 3) {
		if ( pcie_load_io_exec_mb(fd, argv[1], argv[2], PCIE_LOAD_FULL ) ) {
			printf ("Boot of MPPA failed\n");
			exit(1);
		}
	}
	printf("# [HOST] pcie queue init\n");
	pcie_queue_init(fd);
	/* pcie_queue init needs to be called to enable pcie communication via queues */ 
	printf("# [HOST] init queue ok\n");	
	//pcie_register_console(fd, stdin, stdout);
	printf("# [HOST] pcie_register_console ok\n");
	int status;
	int local_status = 0;
	printf("# [HOST] waits\n");	
	pcie_queue_barrier(fd, local_status, &status);
	#define SLEEP_SEC (0)
	printf("# [HOST] sleep %d\n", SLEEP_SEC);	
	sleep(SLEEP_SEC);
	pcie_queue_barrier(fd, local_status, &status);
	pcie_queue_exit(fd, 0, &mppa_ret);
	if(mppa_ret != 0)
		return -1;
	printf("# [HOST] MPPA exited with status %d\n", status);
	pcie_close(fd);
	printf("# [HOST] Goodbye\n");
	return status;
}

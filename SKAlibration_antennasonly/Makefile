#
# Copyright (C) 2015-2016 Kalray SA. All rights reserved.
#
ARCH ?= k1b
nb_cluster ?= 1
nb_pe_cluster ?= 1
time_display ?= 0
test_mode ?= 0
cluster_system ?= bare

arch := $(ARCH)
board := developer

system-name := bare
COMPILE_OPTI := -g
#EXTRA_FLAGS := -ffunction-sections -fdata-sections -DMPPA_TRACE_ENABLE
INCLUDE_FLAGS := -I .

cflags := -Wall -Wextra -Winit-self -Wswitch-default -Wfloat-equal -Wundef -Wshadow -Wuninitialized -Werror
# Cluster rules
cluster-bin := cluster_bin
cluster-system := $(cluster_system)
cluster_bin-srcs := cluster.c calibration.c matrix.c complex.c parallel_matrix.c
cluster-cflags := -g -DNB_CLUSTER=$(nb_cluster) -DNB_PE_CLUSTER=$(nb_pe_cluster) -DTEST_MODE=$(test_mode) -DTIME_DISPLAY=$(time_display) ${COMPILE_OPTI} -mhypervisor $(INCLUDE_FLAGS) $(EXTRA_FLAGS) -Werror


cluster-lflags := -g -mhypervisor -lm -Wl,--defsym=USER_STACK_SIZE=0x2000 -Wl,--defsym=KSTACK_SIZE=0x1000

ifeq ($(cluster_system), bare)
cluster-lflags += -lvbsp -lutask -lmppa_remote -lmppa_async -lmppa_request_engine -lmppapower -lmppanoc -lmpparouting -Wl,--defsym=_LIBNOC_DISABLE_FIFO_FULL_CHECK=0 -Wl,--gc-sections
endif
ifeq ($(cluster_system), nodeos)
cluster-lflags += -pthread -lmppa_remote -lmppa_async -lmppa_request_engine -lmppapower -lmppanoc -lmpparouting -Wl,--defsym=_LIBNOC_DISABLE_FIFO_FULL_CHECK=0 -Wl,--gc-sections
endif

io-bin := io_bin
io_bin-srcs := io_main.c matrix.c complex.c calibration.c parallel_matrix.c
io_bin-cflags := $(INCLUDE_FLAGS) -DNB_CLUSTER=$(nb_cluster) -DNB_PE_CLUSTER=$(nb_pe_cluster) -DTEST_MODE=$(test_mode) -DTIME_DISPLAY=$(time_display) -g ${COMPILE_OPTI} -mhypervisor $(EXTRA_FLAGS)-Werror


#io_bin-lflags := -lvbsp -lmppa_remote -lmppa_async -lmppa_request_engine -lpcie_queue -lutask  -lmppapower -lmppanoc -lmpparouting -mhypervisor -Wl,--defsym=_LIBNOC_DISABLE_FIFO_FULL_CHECK=0
io_bin-lflags := -lutask -lvbsp  -lvbsp -lpcie_queue -lutask -lmppa_remote -lmppa_async -lmppa_request_engine  -lmppapower -lmppanoc -lmpparouting -mhypervisor -Wl,--defsym=_LIBNOC_DISABLE_FIFO_FULL_CHECK=0 -Wl,--defsym=USER_STACK_SIZE=0x4000 -Wl,--defsym=KSTACK_SIZE=0x1000
io_bin-lflags += -Wl,--defsym=_UTASK_SMEM_HEAP_SIZE=0x5000 -lm

mppa-bin := multibin_bin
multibin_bin-objs = io_bin cluster_bin


host-bin := host_bin
host_bin-srcs := host_main.c
host-cflags := ${COMPILE_OPTI}
#host-lflags := -lpthread -lm -lrt -lmppa_remote -lpcie
host-lflags := -lpthread -lrt -lmppa_remote -lpcie
host-bin    := host_bin

include $(K1_TOOLCHAIN_DIR)/share/make/Makefile.kalray

run_sim: all
	$(K1_TOOLCHAIN_DIR)/bin/k1-mppa --march=bostan -- ${O}/bin/multibin_bin.mpk

run_pcie_sim: all
	$(K1_TOOLCHAIN_DIR)/bin/k1-pciesim-runner ./${O}/bin/host_bin ./${O}/bin/multibin_bin.mpk io_bin

ifeq (${board}, emb01)
# When running on emb01, the following variable must be defined.
# k1-remote-runner will also need a password. We STRONGLY recommend
# to do not use --passwd option for that.
# Instead, please consider providing the proper user ssh public key
# from developer to the EMB01 (you may use ssh-copyid standard tool).
# It will allow connecting without providing a password.
EMB01_ADDRESS := $(shell echo $$EMB01_ADDRESS)
EMB01_USER := $(shell echo $$EMB01_USER)

run_jtag: all
	@[ "" != "${EMB01_ADDRESS}" ]  || ( echo "EMB01_ADDRESS not defined"  ; false)
	@[ "" != "${EMB01_USER}" ]     || ( echo "EMB01_USER not defined"     ; false)
	cd $(BIN_DIR);k1-remote-runner --remote ${EMB01_ADDRESS} --user ${EMB01_USER} --host_exec host_bin --multibinary multibin_bin.mpk -- ./host_bin multibin_bin.mpk io_bin

run_hw: all
	@[ "" != "${EMB01_ADDRESS}" ]  || ( echo "EMB01_ADDRESS not defined"  ; false)
	@[ "" != "${EMB01_USER}" ]     || ( echo "EMB01_USER not defined"     ; false)
	cd $(BIN_DIR) 
else

run_jtag: all
	$(K1_TOOLCHAIN_DIR)/bin/k1-jtag-runner --multibinary=./${O}/bin/multibin_bin.mpk --exec-multibin=IODDR0:io_bin

run_jtag_debug: all
	$(K1_TOOLCHAIN_DIR)/bin/k1-jtag-runner --debug-pcie-application --multibinary=./${O}/bin/multibin_bin.mpk --exec-multibin=IODDR0:io_bin -D

run_hw: all
	./${O}/bin/host_bin ./${O}/bin/multibin_bin.mpk io_bin

run_hw_debug: all
	MPPA_LOAD_BYPASS=1 ./${O}/bin/host_bin ./${O}/bin/multibin_bin.mpk io_bin

run_hw_trace: all
	$(K1_TOOLCHAIN_DIR)/bin/k1-trace-util --acquire --clusters=ioddr0,node0,node1,node2,node3 -- ./${O}/bin/host_bin ./${O}/bin/multibin_bin.mpk io_bin

display_traces:
	$(K1_TOOLCHAIN_DIR)/bin/k1-stv --multibinary ${O}/bin/multibin_bin.mpk --hwtrace cluster_bin:node*.DSU.bin

endif

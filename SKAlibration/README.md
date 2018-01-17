# SKAlibration

## Calibration algorithms obtimized on a MPPA platform
These algorithms have been optimized for an up to date MPPA platform.
They are calibrating 64 antennas over 4 polarizations and up to 150 time frequency blocks studied while observing a source at the center of the image.

The algorithms are for now calibrating only test values. It will be able to calibrate real data.

## MeerKat
The aim of this algorithm is to help the SKA project (especially MeerKAT located in South Africa) and to study the portability and the strenght of the MPPA platform.

## Compiling
To build and run the project simply type :

"*make run_hw -B -j8 nb_cluster=16 nb_pe_cluster=16*"

the "*make run_hw -B -j8*"  command will build and run the program having forced it to recompile (B option) faster (j8 option)

it can take two parameters : 

_*nb_cluster*_ an integer from 1 to 16 corresponding to the number of clusters to start (1 by default)

_*nb_pe_cluster*_ an integer from 1 to 16 corresponding to the number of  processing elements per clusters to start (1 by default)

_*test_mode*_ a boolean to use a test vector and compare its result with the theoritical one

_*time_display*_ a boolean to show or not the execution times of the program
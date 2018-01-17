# SKAlibration
The SKAlibration project is aimed to give several solutions of antenna calibration depending on parameters like the number of antenna, the number of polarisations etc. These algorithms are optimized for a KALRAY MPPA platform and should be low power and fast ! 

## Getting Started

###### If you use Linux : 
You have to create an empty repository and to initialize it as a git repository before cloning the whole project.
run the commands 
>*git init* -
and then
>*git pull https://github.com/NicolasSb/SKAlibration.git*

###### If you use windows : 
I'm sorry, I can't do anything for you..

### Prerequisites

To compile and run this project, you do not have to install anything. You just need an up-to-date MPPA.
To have a better comprehension of the algorithms, I will provide some matlab files.
To read and test them, you'll need the Matlab software or the Octave software that you can dowload from the
command line : 
>*apt-get install octave*

### Installing

To build and run the program you should connect to the MPPA and run the following command in the repository containig
the code you would like to test.

> *make run_hw -B -j8 nb_cluster=16 nb_pe_cluster=16  test_mode=0 time_display=0*


the "*make run_hw -B -j8*"  command will build and run the program having forced it to recompile (B option) faster (j8 option)

it can take two parameters : 

_*nb_cluster*_ an integer from 1 to 16 corresponding to the number of clusters to start (1 by default)

_*nb_pe_cluster*_ an integer from 1 to 16 corresponding to the number of  processing elements per clusters to start (1 by default)

_*test_mode*_ a boolean to use a test vector and compare its result with the theoritical one

_*time_display*_ a boolean to show or not the execution times of the program

## Repositories

#### SKAlibration 
This rep. contains the first calibration algorithm. It can be used to calibrate several polarised antennas.

#### SKAlibration_antennaonly
This rep contains an algorithm to calibrate unpolarized data. 

## The code 

As the MPPA has a lot of cores and a very special architecture, it leads us to multiply the files.

The entry point of every project is the host_main.c file.
It will call the io_main.c file that will initialize the RPC server and then, start the clusters.
It will then call functions in the calibration.c file.
The functions calibration() and TestAllF() are running on the IOS and the cluster_main() function is the main function of each cluster.  
The cluster are then working in the "parallel_matrix" file and using the functions of a library (matrix and complex files)


## Authors

* **Nicolas Sourbier** 


## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc


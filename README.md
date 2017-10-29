# GiBUU

This is the repository of the Giessen Boltzmann-Uehling-Uhlenbeck project (GiBUU),
a hadronic transport model for the simulation of nuclear reactions and heavy-ion collisions.
For documentation see the website https://gibuu.hepforge.org.


### Prerequisites

In order to build GiBUU from source, the following tools are required:
* a Fortran compiler (gfortran, ifort, sunf95, ...)
* GNU make
* perl or makedepf90
* libbzip2

For further details and version requirements see https://gibuu.hepforge.org/trac/wiki/tools.


### Compiling

To compile the GiBUU code, it is sufficient to type

    make

After the build process finishes, you should find an executable called "GiBUU.x" in the directory "testRun/".
For further details and additional compilation options, see https://gibuu.hepforge.org/trac/wiki/compiling.


### Running

To run the program, just type

    cd testRun
    ./GiBUU.x < "Name of your job card"

Sample job cards are provided in the directory "testRun/jobCards/".
See also https://gibuu.hepforge.org/trac/wiki/running.

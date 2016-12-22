# navier-stokes
Navier Stokes Equation Solver

In order to compile, you need to edit the configure script used to setup the
environment. Follow these steps:
1. Set the NS_ROOT variable to the directory where the Navier-Stokes project is.
2. Make sure you update the environment variables found under
'Environment paths ...' section to reflect the situation on your system, meaning
the paths to needed libraries are correctly set.
3. While in the project directory, run in a bash shell the following command:
"source ./configure" -- this will setup the environment and you're ready to go !
4. To effectively compile, run in the same terminal the following command:
"ns_build --clean && ns_build --build" -- this will make a clean build.

Note: Steps 1 and 2 are performed only once, when you first setup the
environment !

In order to run the binary resulted from compilation, run the following command
in the same terminal where you compiled: "./run".
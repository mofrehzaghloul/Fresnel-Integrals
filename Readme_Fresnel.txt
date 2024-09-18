Fortran90 Package for Fresnel Sine and Cosine Integrals
This package provides a multiple precision (single, double, and quad) implementation for calculating the Fresnel sine and cosine integrals for both complex and real arguments.

Installation
-------------
Save all the accompanying files in a single directory with a name that is compatible with your Fortran compiler.

Compilation
---------------
1- Using Makefile
Use the provided makefile to compile the package:
make -f makefile

2- Manual Compilation
Alternatively, you can compile the package manually using the terminal. For the gfortran compiler, use:
> gfortran -O3 set_rk.f90 FS_FC_Z_MOD_RK.f90 FS_FC_Z_DRIVER_RK.f90 -o FS_FC_Z_DRIVER_RK

For the Intel Fortran compiler (ifort), replace gfortran with ifort.

Running the Driver Code
---------------------------
Run the driver code by typing the following in the terminal:
>FS_FC_Z_DRIVER_RK
The code will execute automatically and display the results on the screen.

Selecting Precision
----------------------
To select or change the precision, set the integer rk in the set_rk.f90 file to the desired value 

File Descriptions
=============
Makefile
makefile: A script for compiling the package.

Source Files (.f90)
---------------
set_rk.f90                            : An auxiliary module to select the precision by setting the value of the integer rk.
parameters_FC.f90              : Contains numerical constants and parameters used in calculating the Fresnel Cosine Integral.
parameters_FS.f90               : Contains numerical constants and parameters used in calculating the Fresnel Sine Integral.
FS_FC_Z_MOD_RK.f90     : A Fortran90 module that provides generic interfaces for the Fresnel Sine (FresnelS_z(z) & FresnelS_x(x)) 
                                                and Fresnel Cosine (FresnelC_z(z) & FresnelC_x(x)) functions for complex (z = x + iy) or real (z = x) arguments. 
                                               The precision is determined by the integer rk in the subsidiary module set_rk.
FS_FC_Z_DRIVER_RK.f90 : An example of a Fortran driver code or main program.

Text Files
-----------
Readme_Fresnel.txt             : The present file describing the contents of the package.
FresnelC_ref_values2.txt     : Externally generated data using Maple for accuracy checking of the Fresnel Cosine with real arguments.
FresnelS_ref_values2.txt     : Externally generated data using Maple for accuracy checking of the Fresnel Sine with real arguments.
FresnelC_308_zref.txt          : Externally generated data using Maple for accuracy checking of the Fresnel Cosine with complex arguments.
FresnelS_308_zref.txt           : Externally generated data using Maple for accuracy checking of the Fresnel Sine with complex arguments.
Disclaimer_and_License.txt  : A file containing a disclaimer and the license agreement for using the package.



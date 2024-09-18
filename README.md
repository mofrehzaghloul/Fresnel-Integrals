\documentclass{article}
\usepackage{amsmath}
\usepackage{hyperref}
\usepackage{geometry}
\geometry{margin=1in}

\title{Fortran90 Package for Fresnel Sine and Cosine Integrals}
\author{}
\date{}

\begin{document}

\maketitle

\section*{Introduction}
This package provides a multiple precision (single, double, and quad) implementation for calculating the Fresnel sine and cosine integrals for both complex and real arguments.

\section*{Installation}
Save all the accompanying files in a single directory with a name that is compatible with your Fortran compiler.

\section*{Compilation}
\subsection*{1. Using Makefile}
Use the provided makefile to compile the package:
\begin{verbatim}
make -f makefile
\end{verbatim}

\subsection*{2. Manual Compilation}
Alternatively, you can compile the package manually using the terminal. For the gfortran compiler, use:
\begin{verbatim}
gfortran -O3 set_rk.f90 FS_FC_Z_MOD_RK.f90 FS_FC_Z_DRIVER_RK.f90 -o FS_FC_Z_DRIVER_RK
\end{verbatim}
For the Intel Fortran compiler (ifort), replace \texttt{gfortran} with \texttt{ifort}.

\section*{Running the Driver Code}
Run the driver code by typing the following in the terminal:
\begin{verbatim}
FS_FC_Z_DRIVER_RK
\end{verbatim}
The code will execute automatically and display the results on the screen.

\section*{Selecting Precision}
To select or change the precision, set the integer \texttt{rk} in the \texttt{set\_rk.f90} file to the desired value.

\section*{File Descriptions}
\subsection*{Makefile}
\begin{itemize}
    \item \texttt{makefile}: A script for compiling the package.
\end{itemize}

\subsection*{Source Files (.f90)}
\begin{itemize}
    \item \texttt{set\_rk.f90}: An auxiliary module to select the precision by setting the value of the integer \texttt{rk}.
    \item \texttt{parameters\_FC.f90}: Contains numerical constants and parameters used in calculating the Fresnel Cosine Integral.
    \item \texttt{parameters\_FS.f90}: Contains numerical constants and parameters used in calculating the Fresnel Sine Integral.
    \item \texttt{FS\_FC\_Z\_MOD\_RK.f90}: A Fortran90 module that provides generic interfaces for the Fresnel Sine (\texttt{FresnelS\_z(z)} \& \texttt{FresnelS\_x(x)}) and Fresnel Cosine (\texttt{FresnelC\_z(z)} \& \texttt{FresnelC\_x(x)}) functions for complex (\texttt{z = x + iy}) or real (\texttt{z = x}) arguments. The precision is determined by the integer \texttt{rk} in the subsidiary module \texttt{set\_rk}.
    \item \texttt{FS\_FC\_Z\_DRIVER\_RK.f90}: An example of a Fortran driver code or main program.
\end{itemize}

\subsection*{Text Files}
\begin{itemize}
    \item \texttt{Readme\_Fresnel.txt}: The present file describing the contents of the package.
    \item \texttt{FresnelC\_ref\_values2.txt}: Externally generated data using Maple for accuracy checking of the Fresnel Cosine with real arguments.
    \item \texttt{FresnelS\_ref\_values2.txt}: Externally generated data using Maple for accuracy checking of the Fresnel Sine with real arguments.
    \item \texttt{FresnelC\_308\_zref.txt}: Externally generated data using Maple for accuracy checking of the Fresnel Cosine with complex arguments.
    \item \texttt{FresnelS\_308\_zref.txt}: Externally generated data using Maple for accuracy checking of the Fresnel Sine with complex arguments.
    \item \texttt{Disclaimer\_and\_License.txt}: A file containing a disclaimer and the license agreement for using the package.
\end{itemize}

\end{document}



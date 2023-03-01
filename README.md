# Software for the Computation of Analytic Sensitivities of Eigenpairs in Isogeometric Analysis

This repository contains the Software for the Computation of Analytic Sensitivities of Eigenpairs in Isogeometric Analysis as shown in the Paper "On the computation of analytic sensitivities of eigenpairs in isogeometric analysis" by Anna Ziegler, Melina Merkel, Peter Gangl and Sebastian Schöps.

The repository contains: 
-The code to compute the higher order derivatives of the system matrices of the Maxwell Eigenvalue Problem. 
-The geometry descriptions of the two considered pillbox cavities presented in the benchmark example. 
-The geometry descriptions of the TESLA cavity and the pillbox cavity for the numerical example on the shape morphing.

## Requirements
The code was developed using Matlab R2022a and GeoPDEs in Version 3.2.0 [1, 2] and the NURBS package in Version 1.4.1. [3]
The GeoPDEs package needs to be installed seperately and added to the working directory.
The Symbolic Math Toolbox [4] is used to compute the derivatives of the terms C(t) and A(t) for the assembly of the derivatives of the system matrices. If the Toolbox is not installed, the derivatives are read as function handles from a file. 

[1] R. Vázquez, A new design for the implementation of isogeometric analysis in Octave and Matlab: GeoPDEs 3.0, Computers & Mathematics with Applications, Volume 72, Issue 3, 2016, Pages 523-554.

[2] C. de Falco, A. Reali, and R. Vazquez. GeoPDEs: A research tool for isogeometric analysis of PDEs. Advances in Engineering Software, 42(12):1020-1034, 2011.

[3] M. Spink, D. Claxton, C. de Falco, R. Vazquez, The NURBS toolbox, http://octave.sourceforge.net/nurbs/index.html.

[4] The Mathworks, Inc., Natick, Massachusetts, Symbolic Math Toolbox (R2022a) (2022).

## Run me
Running the file "runme.me" computes the derivatives of the examples shown in the paper and displays the geometries.
The code to also compute the higher derivatives of the eigenpair [5] is not made available in this repository. It belongs to the chair of Theoretische Elektrotechnik of the Technische Universität Berlin.

[5] P. Jorkowski, Zur numerischen Berechnung von parametrischen und nichtlinearen Eigenwertproblemen in der elektromagnetischen Feldsimulation, Dissertation, Technische Universität Berlin (2020).

## Acknowledgment
The authors thank Philipp Jorkowski and Rolf Schuhmann for their support and the many fruitful discussions, as well as for supplying us with their implementation of the eigenvalue derivative computation. Furthermore, the authors thank Annalisa Buffa and Rafael Vázquez for the fruitful discussions.
This work is supported by the Graduate School CE within the Centre for Computational Engineering at Technische Universität Darmstadt, the Federal Ministry of Education and Research (BMBF) and the state of Hesse as part of the NHR Program and the SFB TRR361/F90 CREATOR (grant number 492661287) funded by the German and Austrian Research Foundations DFG and FWF. Moreover, support by the FWF funded project P32911 and DFG funded project SCHO 1562/6-1 is acknowledged.

## License
This repository is licensed under the GNU General Public License v3.0.



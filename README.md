# Abinitio DGA
                                                                                                                                                                       
AbinitioDGA is an implementation of a DGA approach for electronic structure calculations of materials. It allows the inclusion of non-local correlations beyond DMFT.

## Getting Started

### Prerequisites

* HDF5 library (>= 1.8)
* LAPACK library
* BLAS library
* Fortran90 Compiler (with MPI support)

### Installation

There are 3 custom library include files provided which work on the according clusters. These can be found in **_make\_configs/make\_config\_\*_** and are detected by **_make\_configs/make\_config\_auto.sh_**. You can add your custom cluster include file to the **_make\_configs_** folder and extend the detection script or simply add your include file in the root folder of the repository under the name **_make\_config_**.

### More Information

For more detailed instructions please visit **documentation/README.pdf** and **documentation/configspec**.

### Post-Processing

The default output of ADGA is a partially compressed HDF5 file. For that reason we strongly recommend using Python for post-processing purposes. Our recommended Python libraries are:

* h5py
* numpy
* scipy
* matplotlib

These packages are however not necessary for the purely computational part of ADGA. If one does not want to deal with Python it is still possible to output all data in form of text-files (see **documentation/configspec**). The instructions in **documentation/python\_intro.pdf** and the scripts found in **documentation/scripts/** provide a first starting point for extracting data and making plots.

## Test Case SrVO3

The test data and results within ADGA for our usual testbed SrVO3 (strongly correlated metal) can be found in **srvo3-testdata/** and **documentation/README.pdf** respectively.

## Authors
Anna Galler, Patrik Thunström, Josef Kaufmann, Matthias Pickem, Jan M. Tomczak and Karsten Held.

## License
This project is licensed under the GNU General Public License v3 which can be found in **LICENSE**.

## Publications
[Ab initio dynamical vertex approximation - PRB](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.115107)

[Ab initio dynamical vertex approximation - arXiv](https://arxiv.org/abs/1610.02998)

[Towards ab initio Calculations with the Dynamical Vertex Approximation - arXiv](https://arxiv.org/abs/1709.02663)

[The AbinitioDGA Project v1.0: Non-local correlations beyond and susceptibilities within dynamical mean-field-theory - arXiv](https://arxiv.org/abs/1710.06651)

## Acknowledgments

This project has been supported by European Research Council under the European Union's Seventh Framework Program (FP/2007-2013) through Grant agreement No. 306447.

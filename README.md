#Molecular Docking Tool

This tool follows closely [AutoDock 4](http://autodock.scripps.edu) implementation. The goal is to develop heterogenous parallel program to accelerate molecular docking.

![Binding Mode](https://raw.github.com/ekaakurniawan/FpgaNeuroDock/master/Images/Molecule/ProteinSS_hsg1_ind.png)

##Python-OpenCL Implementation

Tested on:
* Python 2.7.3
* NumPy 1.7.1
* PyOpenCL 2013.1

To run:
* Go to **PyNeuroDock** directory
* Execute **python NeuroDock.py**

Benchmark:
* On Unix shell (Bash), execute: ./Benchmark.sh
* Python: Sequential processing run on 2.3GHz Intel Core i7
* Python-OpenCL GPU: Parallel processing run on NVIDIA GeForce GT 650M 1GB (384 Cores)
* Python-OpenCL CPU: Parallel processing run on 2.3GHz Intel Core i7 (8 Cores)

![Pyton-OpenCL Benchmark](https://raw.github.com/ekaakurniawan/FpgaNeuroDock/master/Images/Benchmark/Python-OpenCL_500Gens.png)

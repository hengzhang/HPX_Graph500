HPX_Graph500
============

Graph500 benchmark rewritten in C++ using HPX.

Introduction
============
[HPX](https://github.com/STEllAR-GROUP/hpx) is a general purpose C++ runtime system for parallel and distributed applications of any scale. The goal of HPX is to create a high quality, freely available, open source implementation of the ParalleX model for conventional systems

Graph500 is the benchmark used to model important factors of many modern parallel analytical workloads. The current Graph500 benchmarks are implemented using OpenMP and MPI. 

Porting Graph500 to [HPX](https://github.com/STEllAR-GROUP/hpx), which is well suited for the fine-grain and irregular workloads of graph applications.

The project replaces the inherent barrier synchronization with asynchronous communications of HPX.

The motivation is to produce a new benchmark for the HPC community as well as an addition to the HPX benchmark suite and compare HPX’s performance with reference implementations.


Build Instructions
==================

Before starting to build HPX_GRAPH500, please install [HPX](https://github.com/STEllAR-GROUP/hpx) C++ runtime system.

## Linux
### Install the HPX C++ runtime system
You should follow the instruction on the [HPX](https://github.com/STEllAR-GROUP/hpx) github website.

### Install HPX_GRAPH500
1. Clone the HPX_GRAPH500 master git repository:


		git clone git://github.com/hengzhang/HPX_GRAPH500.git

2. Invoke CMake from your build directory. 
	
		cmake -DHPX_ROOT=/where_hpx_be_installed \
			  ...	

3. Invoke GNU make. `-jN` flag to your make invocation, where `N` is the number of cores on your machine plus one:
		
		make -j5

4. To complete the build and install HPX_GRAPH500:
		
		make install
		ldconfig 	

5. Run hpx_graph500 on the bin directory.




-----------------------

Source code in github about Breath First Search algorithm ([BFS on my github](https://github.com/hengzhang/hpx_bfs_test)),  implemented using HPX, is referenced to this repository.
# 15-418 Final Project
# Rachel Yuan and May Li

# Background

The graph coloring problem involves finding a color for each node such that none of its neighbors have that same color. Practical applications of graph coloring include scheduling, register allocation for compiler optimization, and more. In a graph with N nodes, the graph coloring problem can be trivially solved by assigning each node their own unique color, resulting in N colors being used. The minimum graph coloring problem involves finding an assignment of colors to nodes that results in the least number of colors used. 

In our approach to parallelize graph coloring, we will implement a parallel algorithm using OpenMP and OpenMPI. See more about our proposal, milestone report, and final report below for results.

# Implementations

## Sequential

The sequential implementation of graph coloring is located in src/seq-coloring.cpp. In order to compile this, run `make` in the project directory, and run `./color-release -seq -f [file-name]`.

## OpenMP

There are three OpenMP implementations, which are located in the files src/openmp-coloring.cpp, src/jpopenmp-coloring.cpp, and src/half-jpopenmp-coloring.cpp.

The first of these implementations uses a shared color map, and each thread colors a node and writes the color to this map. Due to race conditions, coloring conflicts may arise. To resolve these conflicts, we assign a newly generated unique number to all the coloring conflicts. This results in high speedup and non-optimal coloring.

The second of these implementations implements the Jones-Plassman algorithm, coloring an independent set of vertices in paralell for each iteration until all nodes are colored. This colors the nodes optimally but achieves poor speedup.

The third of these implementations implements a combination of the two, in which it takes the approach of the first one and resolves conflicts with the Jones-Plassman algorithm. The implementation achieves high speedup and near-optimal coloring.

To run these, instead of the `-seq` flag, run with either the `-openmp`, `-jpop`, or `-half` flags respectively.

## OpenMPI

There are 4 OpenMP implemenentations, and to compile them, run `make CONFIGURATION=MPI` or `make CONFIGURATION=MPI2` to run the synchronous MPI coloring.

One implementation is synch-openmpi-coloring.cpp, which attempts to color the graph by requesting colors from its neighbors asynchronously and coloring as long as a node has the colors of its neighbors. This implementation achieves optimal coloring but poor speedup.

## Reports

Project Proposal: https://docs.google.com/document/d/1tic5K9rLctJeTJyd1os6lVys0_OXbuj4-iLJWPYas-Y/edit?usp=sharing

Milestone Report: https://docs.google.com/document/d/12bovf26HVuqaxgXckGCHJ0vLp09WyjbSFthsIoTKr-k/edit?usp=sharing

Final Report: https://docs.google.com/document/d/1pUw0mvtp_b0GKbmuk9YUBuub1qia7bKdhqIsGOnGYxc/edit?usp=sharing




# 15-418 Final Project
# Rachel Yuan and May Li

https://docs.google.com/document/d/1tic5K9rLctJeTJyd1os6lVys0_OXbuj4-iLJWPYas-Y/edit?usp=sharing

## Progress

The sequential implementation is done, which follows a naive greedy coloring algorithm. This is implemented in src/seq-coloring.cpp.
The parallel implementation using OpenMP is done as well and is implemented src/openmp-coloring.cpp.
The parallel implementation using OpenMPI is currently in development.
Test cases are under the tests folder with a couple in the src/main.cpp file in order to reduce storage space.

The sequential implementation currently colors much more optimally than the parallel implementation with OpenMP.
The OpenMP implementation is able to achieve speedup however, but due to the dependencies to color the graph, we've had to make compromises in how the graph is colored.

Our milestone report is here: https://docs.google.com/document/d/12bovf26HVuqaxgXckGCHJ0vLp09WyjbSFthsIoTKr-k/edit?usp=sharing


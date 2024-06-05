#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <utility>
#include <omp.h>
#include<bits/stdc++.h>
#include "common.h"
#include <random>
#include <math.h>
#include <mpi.h>
#include "parallel.h"
#include <cassert>
using namespace std;

#define MPI_SAFE_CALL( call ) do {                               \
    int err = call;                                              \
    if (err != MPI_SUCCESS) {                                    \
        fprintf(stderr, "MPI error %d in file '%s' at line %i",  \
               err, __FILE__, __LINE__);                         \
        exit(1);                                                 \
    } } while(0)



void usage(int argc, char** argv)
{
    if(argc < 6) {
        cerr << "please have enough command line arguments";
        exit(EXIT_FAILURE);
    } 
}



int main(int argc, char* argv[]){
    int num_procs = 0, rank = 0;
    MPI_Init(&argc, &argv);
    MPI_SAFE_CALL(MPI_Comm_size(MPI_COMM_WORLD, &num_procs));
    MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

    if(rank == 0) {
        cout << "num of processes: " << num_procs << endl;
    }
    int xSize = 2, ySize = 2, zSize = 2, x = 0, y = 0, z = 0;

    usage(argc, argv);

    if(argc > 1){
        xSize = atoi(argv[1]);
        ySize = atoi(argv[2]);
        zSize = atoi(argv[3]);
        x = atoi(argv[4]);
        y = atoi(argv[5]);
        z = atoi(argv[6]);
    }

    
    int totalVertices = xSize * ySize * zSize;
    int src = x * ySize * zSize + y * zSize + z;
    assert(src < totalVertices);
    
    vector<int> serial(totalVertices, 0);
    vector<int> parallel(totalVertices, 0);
    vector<int> MPIparallel(totalVertices, 0);
    vector<int> spath(totalVertices, -1);
    vector<int> ppath(totalVertices, -1);
    vector<int> MPIpath(totalVertices, -1);
   
    
    // Create a graph with the appropriate number of vertices
    Graph graph(totalVertices);
    random_device rd;
    mt19937 gen(rd());

    // Define the distribution for integers in the range [1, 10]
    uniform_real_distribution<float> distribution(1, 10);

    // Generate and print a random integer from the distribution
    int randomn = 0;
    if (rank == 0) {
        for (int x = 0; x < xSize; ++x) {
            for (int y = 0; y < ySize; ++y) {
                for (int z = 0; z < zSize; ++z) {
                    int from = x * ySize * zSize + y * zSize + z;
                    // Connect to the right neighbor
                    if (x < xSize - 1) {
                        int to = (x + 1) * ySize * zSize + y * zSize + z;
                        randomn = (int)distribution(gen);
                        graph.addEdge(from, to, randomn);

                    }
                    // Connect to the bottom neighbor
                    if (y < ySize - 1) {
                        int to = x * ySize * zSize + (y + 1) * zSize + z;
                        randomn = (int)distribution(gen);
                        graph.addEdge(from, to, randomn);

                    }
                    // Connect to the forward neighbor
                    if (z < zSize - 1) {
                        int to = x * ySize * zSize + y * zSize + (z + 1);
                        randomn = (int)distribution(gen);
                        graph.addEdge(from, to, randomn);

                    }
                }
            }
        }
    }
    // Connect each vertex to its neighbors
    //if (rank == 0) graph.printGraph();
    
    uint64_t start_t;
    uint64_t end_t;
    InitTSC();
    ///serial 
    if (rank == 0) {
        start_t = ReadTSC();
        BellmanFord(graph, src, serial, spath);
        end_t = ReadTSC();
        cout << "Time for serial: " << ElapsedTime(end_t - start_t) << " seconds" << endl;
        //printArr(serial, totalVertices, xSize);
    }
    //printpath(spath, xSize);

    ///openMP
    if (rank == 0) {
        cout << "parallel start" << endl;
        double startt = MPI_Wtime();
        BF_openMP(graph, src, parallel, ppath);
        double endt = MPI_Wtime();
        cout << "Time for openMP: " << endt - startt << " seconds" << endl;
        //printArr(parallel, totalVertices, xSize);
    }
    //printpath(ppath, xSize);
    //parprintpath(ppath, xSize);

    ///openMPI
    double startt = MPI_Wtime();
    BF_MPI(graph, x, y, z, xSize, MPIparallel, MPIpath);
    double endt = MPI_Wtime();
    if (rank == 0) cout << "Time for MPI: " << endt - startt << " seconds" <<endl;
    //if (rank ==0) printArr(MPIparallel, totalVertices, xSize);
    if(rank==0) verify(serial, MPIparallel, totalVertices);
    MPI_Finalize();
}

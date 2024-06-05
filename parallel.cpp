#include <mpi.h>
#include <vector>
#include "parallel.h"
#include "common.h"
#include <climits>
#include <iostream>
using namespace std;
#define MPI_SAFE_CALL( call ) do {                               \
    int err = call;                                              \
    if (err != MPI_SUCCESS) {                                    \
        fprintf(stderr, "MPI error %d in file '%s' at line %i",  \
               err, __FILE__, __LINE__);                         \
        exit(1);                                                 \
    } } while(0)

void verify(const vector<int>& serial, const vector<int>& parallel, int V) {
    int err = 0;
    for(int i=0; i<V; i++) {
        if(serial[i] != parallel[i]) err++; 
    }
    if (err != 0) cout << "num of error: " << err << endl;
    else cout << "pass" << endl;
}
void printArr(const vector<int>& dist, int n, int size) {
    cout << " Vertex Distance from Source "<< endl;
    for (int i = 0; i < n; ++i) {
        int X_axis = i/(size*size);
        int Y_axis = (i%(size*size))/size;
        int Z_axis = (i%(size*size))%size;
        printf("Grid[%d][%d][%d] \t\t %d\n", X_axis, Y_axis, Z_axis, dist[i]);
        
    }
}
void printpath(const vector<int>& path, int size) {
    vector<string> finalpath(5);
    int iter = 0, X_axis = 0, Y_axis = 0, Z_axis = 0;
    for (int i = 1; i < 20; i=i*2) {
        int current = i;
        X_axis = i/(size*size);
        Y_axis = (i%(size*size))/size;
        Z_axis = (i%(size*size))%size;
        printf("starting point Grid[%d][%d][%d]:\n", X_axis, Y_axis, Z_axis);
        while(path[current] != -1) {
            X_axis = path[current]/(size*size);
            Y_axis = (path[current]%(size*size))/size;
            Z_axis = (path[current]%(size*size))%size;
            finalpath[iter] += "-";
            finalpath[iter] += "[";
            finalpath[iter] += to_string(X_axis);
            finalpath[iter] += "]";
            finalpath[iter] += "[";
            finalpath[iter] += to_string(Y_axis);
            finalpath[iter] += "]";
            finalpath[iter] += "[";
            finalpath[iter] += to_string(Z_axis);
            finalpath[iter] += "]";
            current = path[current];
        }
        cout << finalpath[iter]<<endl;
        cout << endl;
        iter++;
    } 
}
void parprintpath(const vector<int>& path, int size) {
    vector<string> finalpath(2);
    int first = 0, second = 0, iter = 0, X_axis = 0, Y_axis = 0, Z_axis = 0;
    first = 12*size*size + 12*size + 12;
    second = 37*size*size + 37*size + 37;
    while(path[first] != -1) {
            X_axis = path[first]/(size*size);
            Y_axis = (path[first]%(size*size))/size;
            Z_axis = (path[first]%(size*size))%size;
            finalpath[iter] += "-";
            finalpath[iter] += "[";
            finalpath[iter] += to_string(X_axis);
            finalpath[iter] += "]";
            finalpath[iter] += "[";
            finalpath[iter] += to_string(Y_axis);
            finalpath[iter] += "]";
            finalpath[iter] += "[";
            finalpath[iter] += to_string(Z_axis);
            finalpath[iter] += "]";
            first = path[first];
    }
    cout << "path from [12][12][12] to source: " << endl;
    cout << finalpath[iter]<<endl;
    iter++;
    while(path[second] != -1) {
            X_axis = path[second]/(size*size);
            Y_axis = (path[second]%(size*size))/size;
            Z_axis = (path[second]%(size*size))%size;
            finalpath[iter] += "-";
            finalpath[iter] += "[";
            finalpath[iter] += to_string(X_axis);
            finalpath[iter] += "]";
            finalpath[iter] += "[";
            finalpath[iter] += to_string(Y_axis);
            finalpath[iter] += "]";
            finalpath[iter] += "[";
            finalpath[iter] += to_string(Z_axis);
            finalpath[iter] += "]";
            second = path[second];
    }
    cout << "path from [37][37][37] to source: " << endl;
    cout << finalpath[iter]<<endl;
    
}
    
void BellmanFord(Graph graph, int src, vector<int>& serial, vector<int>& spath) {
    int V = graph.V;
    vector<vector<pair<int, int>>> adjList = graph.adjList;
    vector<int> dist(V, INT_MAX);
    dist[src] = 0;
    for(int k=1; k<V; k++) {
        for(int i=0; i<V; i++) {
            for(auto neighbor:adjList[i]) {
                if(dist[i] != INT_MAX && dist[i] + neighbor.second < dist[neighbor.first]) {
                    dist[neighbor.first] = dist[i] + neighbor.second;
                    spath[neighbor.first] = i;   
                }
            }   
        }
    }
    
    /*for(int i=0; i<V; i++) {
        for(auto neighbor:adjList[i]) {
            if(dist[i] != INT_MAX && dist[i] + neighbor.second < dist[neighbor.first]) {
                cout << "there is negative cycle" << endl;
                return;
            }
        }   
    }*/
    serial = dist;
    
}

void BF_openMP(Graph graph, int src, vector<int>& parallel, vector<int>& ppath) {
    int V = graph.V;
    vector<vector<pair<int, int>>> adjList = graph.adjList;
    vector<int> dist(V, INT_MAX);
    dist[src] = 0;
    #pragma omp parallel for schedule(static)
    for(int k=1; k<V; k++) {
        for(int i=0; i<V; i++) {
            for(auto neighbor:adjList[i]) {
                if(dist[i] != INT_MAX && dist[i] + neighbor.second < dist[neighbor.first]) {
                    dist[neighbor.first] = dist[i] + neighbor.second;
                    ppath[neighbor.first] = i;   
                }
            }   
        }
        
    }
    
    /*for(int i=0; i<V; i++) {
        for(auto neighbor:adjList[i]) {
            if(dist[i] != INT_MAX && dist[i] + neighbor.second < dist[neighbor.first]) {
                cout << "there is negative cycle" << endl;
                return;
            }
        }   
    }*/

    parallel = dist;
    
    //printpath(path);
}

void BF_MPforMPI (const vector<vector<pair<int, int>>>& adjList, vector<int>& dist, vector<int>& ppath, int baserow, int gsize, int& updated, int rank, int num_procs) {
    int offset = baserow * gsize * gsize;
    int localstart = 0;
    int localend = dist.size();
    int cap = localend + offset;
    
    if (rank == 0) {
        localstart = 1*gsize*gsize;
        cap = cap - localstart;
        #pragma omp parallel for schedule(static)
        for(int i=localstart; i<localend; i++) {
            for(auto neighbor:adjList[i-localstart]) {
                if(dist[i] == INT_MAX || (neighbor.first >= cap)) {
                    continue;
                } 
                else if(dist[i] + neighbor.second < dist[(neighbor.first+localstart)]) {
                    dist[(neighbor.first+localstart)] = dist[i] + neighbor.second;
                    //ppath[(neighbor.first-offset)] = i;    
                    updated = 1;
                }
            }   
        }
    }
    else if (rank == num_procs-1) {
        localend -= gsize*gsize;
        cap -= gsize*gsize;
        #pragma omp parallel for schedule(static)
        for(int i=localstart; i<localend; i++) {
            for(auto neighbor:adjList[i+offset]) {
                if(dist[i] == INT_MAX || (neighbor.first >= cap) || neighbor.first < offset) {
                    continue;
                } 
                else if(dist[i] + neighbor.second < dist[(neighbor.first-offset)]) {
                    dist[(neighbor.first-offset)] = dist[i] + neighbor.second;
                    //ppath[(neighbor.first-offset)] = i;   
                    updated = 1;
                }
            }   
        }
    }
    else {
        #pragma omp parallel for schedule(static)
        for(int i=localstart; i<localend; i++) {
            for(auto neighbor:adjList[i+offset]) {
                if(dist[i] == INT_MAX || (neighbor.first >= cap) || neighbor.first < offset) {
                    continue;
                } 
                else if(dist[i] + neighbor.second < dist[(neighbor.first-offset)]) {
                    dist[(neighbor.first-offset)] = dist[i] + neighbor.second;
                    //ppath[(neighbor.first-offset)] = i;   
                    updated = 1;
                }
            }   
        }
    }
    
    //parallel = dist;
    //printpath(path);
}

void BF_MPI (Graph graph, int x, int y, int z, int gsize, vector<int>& parallel, vector<int>& ppath) {
    int num_procs = 0, rank = 0;
    MPI_SAFE_CALL(MPI_Comm_size(MPI_COMM_WORLD, &num_procs));
    MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
    ///MPI_type for pair<int, int>
    int num_in_struct = 2;
    int* ele_block_len = new int[num_in_struct];
    ele_block_len[0] = 1;
    ele_block_len[1] = 1;
    MPI_Aint displa[2] = {0, sizeof(int)};
    MPI_Datatype* pairtype = new MPI_Datatype[2];
    pairtype[0] = MPI_INT;
    pairtype[1] = MPI_INT;
    MPI_Datatype my_pair;
    MPI_Type_create_struct(num_in_struct, ele_block_len, displa, pairtype, &my_pair);
    MPI_Type_commit(&my_pair);//remember free
    

    ///MPI_cart 1d
    int dims[1] = {num_procs};
    int periods[1] = {0};
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 1, &cart_comm);
    int left_proc, right_proc;
    MPI_Cart_shift(cart_comm, 0, 1, &left_proc, &right_proc);

    ///local initialize
    int V = graph.V, offset = 0, displace = 0;
    vector<vector<pair<int, int>>> adjList = graph.adjList;
    vector<vector<pair<int, int>>> adj2d;
    adj2d.resize(V);
    vector<pair<int, int>> adj1dList;
    vector<int> localarrsize(num_procs, 0);
    vector<int> offsets(num_procs, 0);
    int remain = gsize%num_procs;
    int localsize = gsize/num_procs;
    if (rank == 0) {
        flattenvec(adjList, adj1dList);
    }
    int size1d = adj1dList.size();
    int srcrank = 0;
    MPI_SAFE_CALL(MPI_Bcast(&size1d, 1, MPI_INT, 0, MPI_COMM_WORLD));
    if (rank != 0) {
        adj1dList.resize(size1d);
    }
    MPI_SAFE_CALL(MPI_Bcast(&adj1dList[0], size1d, my_pair, 0, MPI_COMM_WORLD));
    unflattenvec(adj2d, adj1dList);
    vector<int> displacement(num_procs, 0);
    vector<int> recvbuf(num_procs, 0);
    for (int i=0; i<num_procs; i++) {
        localarrsize[i] = localsize;
        if (remain > 0) {
            localarrsize[i]++;
            remain--;
        }
        offsets[i] = offset; // which row starts
        displacement[i] = displace;
        offset += localarrsize[i];
        displace += localarrsize[i]*gsize*gsize;
        recvbuf[i] = localarrsize[i]*gsize*gsize;
    }
    for(int i=1; i<num_procs; i++) {
        if( x >= offsets[i-1] && x < offsets[i]) srcrank = i-1;
        else if (x >= offsets[i] && i == num_procs-1) srcrank = i;
    }
    //cout << "offset: "<<offsets[rank] << endl;
    int flocalsize = (localarrsize[rank]+2)*gsize*gsize;
    int sendsize = localarrsize[rank]*gsize*gsize;
    vector<int> sendarray(sendsize, 0);
    vector<int> localarray(flocalsize, 0); //ghostcell: + 2, for local dist                                                                
    for (int i = 0; i < localarrsize[rank]+2; ++i) {                                                                        
        for (int j = 0; j < gsize; ++j) {
            for (int k = 0; k < gsize; ++k) {
                localarray[i * gsize * gsize + j * gsize + k] = INT_MAX;
            }
        }
    }
    
    if (rank == srcrank) {
        x -= offsets[rank];
        localarray[(x+1)* gsize * gsize + y*gsize + z] = 0; /// src == 0
    }

    ///MPI_subarray
    int gsizes[2] = {gsize, gsize};
    int starts[2] = {0, 0};
    int subsizes[2] = {gsize, gsize};
    MPI_Datatype subarray; 
    MPI_Type_create_subarray(2, gsizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
    MPI_Type_commit(&subarray);//remember free

    int updated;
    int maxiter = V/num_procs;
    for (int i=1; i<maxiter; i++) {
        updated = 0;
        ///MPI_Sendrecv
        MPI_Sendrecv(&localarray[localarrsize[rank]*gsize*gsize], 1, subarray, right_proc, 0,
                     &localarray[0], gsize*gsize, MPI_INT, left_proc, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&localarray[1*gsize*gsize], 1, subarray, left_proc, 0,
                 &localarray[(localarrsize[rank]+1)*gsize*gsize], gsize*gsize, MPI_INT, right_proc, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
        if (rank == 0) {
            int baserow = offsets[rank];
            BF_MPforMPI(adj2d, localarray, ppath, baserow, gsize, updated, rank, num_procs);
        }
        else {
            int baserow = offsets[rank]-1;
            BF_MPforMPI(adj2d, localarray, ppath, baserow, gsize, updated, rank, num_procs);
        }
        int globalUpdated;
        MPI_Allreduce(&updated, &globalUpdated, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        if (!globalUpdated) break;
    }
    
    int end = localarrsize[rank];
    preparedata(localarray, sendarray, rank, end, gsize);
    int sendcount = sendarray.size();
    
    //MPI_SAFE_CALL(MPI_Gather(sendarray.data(), sendcount, MPI_INT, parallel.data(), sendcount, MPI_INT, 0, MPI_COMM_WORLD));
    MPI_SAFE_CALL(MPI_Gatherv(sendarray.data(), sendcount, MPI_INT, parallel.data(), &recvbuf[0], &displacement[0], MPI_INT, 0, MPI_COMM_WORLD));
    
    MPI_Type_free(&subarray);
    MPI_Type_free(&my_pair);
}

void flattenvec(const vector<vector<pair<int, int>>>& vec2d, vector<pair<int, int>>& vec1d) {
    for(auto row:vec2d) {
        vec1d.insert(vec1d.end(), row.begin(), row.end());
        vec1d.push_back({-1,-1});
    }
}

void unflattenvec(vector<vector<pair<int, int>>>& vec2d, vector<pair<int, int>>& vec1d) {
    int i = 0;
    for(auto index:vec1d) {
        if(index.first == -1) {
            i++;
            continue;
        }
        vec2d[i].push_back(index);
    }
}

void preparedata(const vector<int>& localarray, vector<int>& sendarray, int rank, int flocalsize, int gsize) {
    int start = gsize*gsize;
    int end = (flocalsize+1)*gsize*gsize;
    for(int i=start; i<end; i++) {
        sendarray[i-start] = localarray[i]; 
    }
}
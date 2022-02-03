//Parallel and Distributed Systems
//
//Spaias Georgios
//AEM: 8910

#include "cuda_runtime.h"
#include "device_launch_parameters.h"


#include "stdio.h"
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>
#include "Utils.cuh"


void sequential_Ising(int8_t*** G, int8_t*** G_trans, int n, int k);

void calculate_grid(int* numOfBlocks, int* numOfThreads, int pointsPerThread, int lattice_size);

__host__ void Ising_comp(int8_t** G, int8_t** G_trans, int* size, int* steps, int* numofBlocks, int* numofThreads, int* pointsPerThread, int tile_size);

__host__ void swap_dptr(void** x, void** y);

__host__ __device__ int8_t sgn(int8_t x);

__global__ void compute_State(int8_t* a, int8_t* b, int* size, int* dev_pointsPerThread);


int main(int argc, char* argv[]) {

    //n:Lattice Dimension & k:Number of Iterations
    int n =80, k =1, cnvrg = 0;
    int* dev_n;
    int* dev_cnvrg;

    int numOfBlocks = 1;
    int numOfThreads = 1;
    int pointsPerThread = 5;

    dev_n = (int*)malloc(sizeof(int));
    dev_cnvrg = (int*)malloc(sizeof(int));

    int8_t** lattice;
    int8_t* dev_lattice = (int8_t*)malloc(n * n * sizeof(int8_t));
    int8_t* dev_lattice_trans = (int8_t*)malloc(n * n * sizeof(int8_t));
    int8_t** lat_trans;

    cudaError_t cudaStatus;
    cudaEvent_t start, stop;
    float milliseconds = 0;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    FILE* f_ptr;
    if (EVALUATION_MODE)
    {

        f_ptr = fopen("eval.bin", "rb");
        if (f_ptr == NULL)
        {
            fprintf(stderr, "main: Failed to open eval.bin .\n");
            exit(1);
        }
        fread(&n, sizeof(int), 1, f_ptr);
        fread(&k, sizeof(int), 1, f_ptr);
        printf("Generated Lattice with n=%d and k=%d\n", n, k);
        fclose(f_ptr);
    }

    //Allocate memory for the host lattice
    lattice = lattice_init(n);
    lat_trans = lattice_init(n);

    //Initialize lattice with a random starting state
    start_state(lattice, n);
    printf("GPU ISING COMPUTATION V2\n\n");
    //printLattice(lattice, &n);

    cudaStatus = latticeInitCuda(&lattice, &dev_lattice, &dev_lattice_trans, &n, &cnvrg, &dev_n, &dev_cnvrg);
    if (cudaStatus != cudaSuccess)
    {
        fprintf(stderr, "latticeInitCuda Failed:\t%s\n", cudaGetErrorString(cudaStatus));
        return 1;
    }
    cudaDeviceSynchronize();

    calculate_grid(&numOfBlocks, &numOfThreads, pointsPerThread, n);
    int tile_size = (((numOfThreads*pointsPerThread)/n+1)*n+2*n) * sizeof(int8_t);
    cudaEventRecord(start);

    Ising_comp(&dev_lattice, &dev_lattice_trans, dev_n, &k, &numOfBlocks, &numOfThreads, &pointsPerThread, tile_size);

    cudaDeviceSynchronize();

    cudaEventRecord(stop);

    cudaEventElapsedTime(&milliseconds, start, stop);
    printf("1~GPU Ising Computation: Done\n");
    //printf("lat Address %p G Address %p |lat_trans Address %p  G_trans Address %p\n", &dev_lattice, dev_lattice, &dev_lattice_trans, dev_lattice_trans);


    sequential_Ising(&lattice, &lat_trans, n, k);
    printf("2~Sequential Ising Computation: Done\n");


    for (size_t i = 0; i < n; i++)
    {
        cudaStatus = cudaMemcpy(*(lattice + i), (dev_lattice + n * i), n * sizeof(int8_t), cudaMemcpyDeviceToHost);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed1!");
            return cudaStatus;
        }
    }
    //printLattice(lattice, &n);
    int flag = 1;
    for (int i = 0; i < n; i++)
    {
         for (int j = 0; j < n; j++)
         {
             if (lattice[i][j] != lat_trans[i][j])
             {
                 //printf("\nInconsistent outcome between parallel and Sequential computation on element: %d/%d!~~~~~~\n\n ", i * n + j, n * n);
                 flag = 0;
             }
         }
    }
    if (flag) printf("3~GPU and Sequential Ising Computation are consistent\n");
    
    //printLattice(lat_trans, &n);

    cudaEventElapsedTime(&milliseconds, start, stop);
    printf("\nIterations(k)\tDim(n)\tMaxThreadNumPerBlock\tNumofBlocks\tnumOfThreads\tpointsPerThread\tElapsed time (ms)\n");
    printf("%d\t%d\t1024\t%d\t%d\t%d\t%f\n", k, n, numOfBlocks, numOfThreads, pointsPerThread, milliseconds);

    cudaFree(dev_lattice);
    cudaFree(dev_lattice_trans);
    cudaFree(dev_n);
    cudaFree(dev_cnvrg);

    free(lattice);
    free(lat_trans);


    return 0;
}







__host__ void Ising_comp(int8_t** lat, int8_t** lat_trans, int* size, int* steps, int* numofBlocks, int* numofThreads, int* pointsPerThread ,int tile_size)
{
    cudaError_t cudaStatus;
    int k = *steps;
    
    int* dev_pointsperThread = 0;
    cudaStatus = cudaMalloc(&dev_pointsperThread, sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "CUDA MALLOC Failed:\t%s\n", cudaGetErrorString(cudaStatus));
    }
    cudaStatus = cudaMemcpy((void*)dev_pointsperThread, (void*)pointsPerThread, sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "MemCopy Failed:\t%s\n", cudaGetErrorString(cudaStatus));
    }
    cudaDeviceSynchronize();
    for (int i = 0; i < k; i++)
    {
        //cudaDeviceSynchronize();
        compute_State <<< *numofBlocks, *numofThreads , tile_size >>> (*lat, *lat_trans, size, dev_pointsperThread);
        fprintf(stderr, "kernel Launch Error:\t%s\n", cudaGetErrorString(cudaGetLastError()));
        cudaDeviceSynchronize();
        printf("%s\n", cudaGetErrorString(cudaThreadSynchronize()));
        swap_dptr((void**)lat, (void**)lat_trans);
        
        fprintf(stderr, "Sync Error:\t%s\n", cudaGetErrorString(cudaGetLastError()));

    }
    
    cudaFree(dev_pointsperThread);
    fprintf(stderr, "GPU Ising Computation finished with:\t%s\n", cudaGetErrorString(cudaGetLastError()));
}


__global__ void compute_State(int8_t* G, int8_t* G_trans, int* size, int* pointsPerThread)
{   
    extern __shared__ int8_t shrd_tile[];
    
    int sum,i, j;
    int ppt = *pointsPerThread;
    int n = *size;
    int tile_size = ((((int) blockDim.x )* ppt) / n + 1) * n +2*n;
    int k = tile_size / n;
    int threadIndex = (blockIdx.x * blockDim.x + threadIdx.x) * ppt;
    int index_start = ((threadIndex/n)+1)* n + threadIndex % n;
    int index_end = index_start + blockDim.x * ppt;
    printf("%d Part 0 with TileSize: %d |blockDim %d |blockIdx %d | ppt %d\n", threadIndex, tile_size, blockDim.x, blockIdx.x, ppt);
    for (int k = 0; k < ppt; k++)
    {
        shrd_tile[index_start + k] = G[threadIndex + k];
    }
    printf("Part 1 \n");
    if (threadIdx.x==0)
    {
       
        for (int k = 0; k < (tile_size- index_end); k++)
        {  
            //printf("\n%d\n", (blockDim.x * ppt + i) % (n * n));
            shrd_tile[index_end+k] =  G[(blockDim.x * ppt  + k) % (n * n)];
        }
    }
    //printf("Part 2 \n");
    if (threadIdx.x == 0)
    {   
        for (int k = 0; k < index_start ; k++)
        {
            shrd_tile[k] = G[(threadIndex - index_start+ n*n+ k) % (n * n)];
        }
    }
    //printf("Part 3 \n");
    __syncthreads();

    //if (threadIdx.x == 0)
    //{
    //    printf("TILE_SIZE: %d , index_start: %d , index_end: %d\n", tile_size, index_start, index_end);
    //    for (int i = 0; i <tile_size; i++)
    //    {
    //        if (i % n == 0) { printf("\n"); }
    //        printf("%d\t", shrd_tile[i]);
    //    }
    //}
    sum = 0;
    for (int index =0 ; index < ppt; index++)
    {
        i = (index_start + index) / n;
        j = (index_start + index) % n;

        if (threadIndex + index < (n*n))
        {
            sum = shrd_tile[i * n + j];
            //neighbor 1 -->(i,j+1)
            sum += shrd_tile[(i * n + (j + 1) % n)];
            //neighbor 2 -->(i,j-1)
            sum += shrd_tile[(i * n + (j - 1 + n) % n) ];
            //neighbor 3 -->(i-1,j)
       
            sum += shrd_tile[(((i - 1 + n) % n) * n + j) ];
            //neighbor 4 -->(i+1,j)
            
            sum += shrd_tile[(((i + 1) % n) * n + j) ];

            G_trans[(threadIndex+index)] = sgn(sum);
        }
    }
    
}


__forceinline__ __host__ __device__ int8_t sgn(int8_t x)
{
    x = (int8_t)(x > 0) - (int8_t)(x < 0);
    return x;
}

__host__ void swap_dptr(void** x, void** y)
{
    void* temp = *x;
    *x = *y;
    *y = temp;
}




void calculate_grid(int* numOfBlocks, int* numOfThreads, int pointsPerThread,  int lattice_size) {
    cudaDeviceProp cudaProperties;
    cudaGetDeviceProperties(&cudaProperties, 0);
    int maxThreadsPerBlock = cudaProperties.maxThreadsPerBlock;
    int lat_points = lattice_size * lattice_size;
    int pointsPerThreadSquared = pointsPerThread ;
    int block_num = 1;
    int thread_num = 1;

    block_num = lat_points / (pointsPerThreadSquared * maxThreadsPerBlock) + ((lat_points % (pointsPerThreadSquared * maxThreadsPerBlock)) > 0);
    thread_num = lat_points / (pointsPerThreadSquared * block_num) + ((lat_points % (pointsPerThreadSquared * block_num)) > 0);
    
    int tile_size = (((thread_num * pointsPerThread) / lattice_size + 1) * lattice_size + 2 * lattice_size) * sizeof(int8_t);
    printf("max thread num: %d ||block_num: % d || thread_num : % d|| tile size : %d\n", maxThreadsPerBlock, block_num, thread_num, tile_size);

    *numOfBlocks = block_num;
    *numOfThreads = thread_num;

}


void sequential_Ising(int8_t*** lat, int8_t*** lat_trans, int n, int k) {

    int8_t** G = *lat;
    int8_t** G_trans = *lat_trans;

    int cnvrg = 0;
    int n_i = 0, n_j = 0;
    for (int l = 0; l < k; l++)
    {
        cnvrg = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                G_trans[i][j] = 0;
                n_i = (i + 1) % n;              //neighbor 1 -->(i+1,j)
                n_j = j;
                G_trans[i][j] += G[n_i][n_j];
                n_i = i;                        //neighbor 2 -->(i,j+1)
                n_j = (j + 1) % n;
                G_trans[i][j] += G[n_i][n_j];
                n_i = (i - 1 + n) % n;          //neighbor 3 -->(i-1,j)
                n_j = j;
                G_trans[i][j] += G[n_i][n_j];
                n_i = i;                        //neighbor 4 -->(i,j-1)
                n_j = (j - 1 + n) % n;
                G_trans[i][j] += G[n_i][n_j];
                G_trans[i][j] += G[i][j];
                G_trans[i][j] = sgn(G_trans[i][j]);
                cnvrg += abs(sgn(G_trans[i][j] - G[i][j]));
            }
        }
        if (cnvrg == 0)
        {
            printf("\nSequential Ising model has converged on step %ld/%d \n", l, k);
            break;
        }
        swap_dptr((void**)&G, (void**)&G_trans);

    }
    swap_dptr((void**)&G, (void**)&G_trans);
    *lat = G;
    *lat_trans = G_trans;
}
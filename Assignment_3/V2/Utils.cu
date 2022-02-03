#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "stdio.h"
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>
#include "Utils.cuh"



//Allocate memory for the int ** pointer
int8_t** lattice_init(int n) {

    int8_t** ptr;
    ptr = (int8_t**)malloc(n * sizeof(int8_t*));
    for (size_t i = 0; i < n; i++)
    {
        *(ptr + i) = (int8_t*)calloc(n, sizeof(int8_t));
    }
    return ptr;
}

//Initialize lattice with a random starting state
void start_state(int8_t** lat, int n) {

    if (EVALUATION_MODE)
    {
        read_bin(lat, n);
    }
    else
    {
        srand((unsigned int)time(NULL));
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                lat[i][j] = 1 - 2 * (rand() % 2); //Random Moments of spin
                //lat[i][j] = 1 - 2 * ((i % 2 + j) % 2);
            }
        }
    }
    return;
}

//Reading initial state of known lattice for EVALUATION purposes. 
void read_bin(int8_t** lat, int n) {

    FILE* f_ptr = fopen("eval.bin", "rb");
    if (f_ptr == NULL)
    {
        fprintf(stderr, "read_bin: Failed to open eval.bin .\n");
        exit(1);
    }
    fseek(f_ptr, 2 * sizeof(int), SEEK_SET);


    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            fread(&(lat[i][j]), sizeof(int8_t), 1, f_ptr);
            if (lat[i][j] != 1 && lat[i][j] != -1)
            {
                printf("\nProblem afentiko!\n\n ");
            }
        }
    }
    fclose(f_ptr);
}



void resultCheck(int8_t** lat) {
    int n = 10;
    int8_t val;
    FILE* f_p = fopen("eval.bin", "rb");
    if (f_p == NULL)
    {
        fprintf(stderr, "resultCheck: Failed to open eval.bin .\n");
        exit(1);
    }
    fread(&n, sizeof(int), 1, f_p);
    fseek(f_p, 2 * sizeof(int) + n * n * sizeof(int8_t), SEEK_SET);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            fread(&val, sizeof(int8_t), 1, f_p);
            if (val != lat[i][j])
            {
                printf("\n~~~~~~~Incorrect Ising model %d/%d!~~~~~~\n\n ", i * n + j, n * n);

                fclose(f_p);
                return;
            }
        }
    }
    printf("\n~~~~~~~Correct Evaluation of Ising model!~~~~~~~\n\n ");
    fclose(f_p);
}



cudaError_t latticeInitCuda(int8_t*** lattice, int8_t** dev_lattice, int8_t** dev_lattice_trans, int* size, int* cnvrg, int** dev_n, int** dev_cnvrg) {

    cudaError_t cudaStatus;
    int n = *size;
    int8_t** lat = *lattice;

    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        return cudaStatus;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)dev_lattice, n * n * sizeof(int8_t));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "%s\n", cudaGetErrorString(cudaStatus));
        return cudaStatus;
    }
    for (size_t i = 0; i < n; i++)
    {
        cudaStatus = cudaMemcpy((*dev_lattice + (n * i)), *(lat + i), n * sizeof(int8_t), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed1!");
            return cudaStatus;
        }
    }

    cudaStatus = cudaMalloc((void**)dev_lattice_trans, n * n * sizeof(int8_t));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "%s\n", cudaGetErrorString(cudaStatus));
        return cudaStatus;
    }
    for (size_t i = 0; i < n; i++)
    {
        cudaStatus = cudaMemcpy((*dev_lattice_trans + n * i), *(lat + i), n * sizeof(int8_t), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpy failed1!");
            return cudaStatus;
        }
    }
    cudaStatus = cudaMalloc((void**)dev_cnvrg, sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        return cudaStatus;
    }
    cudaStatus = cudaMemcpy(*dev_cnvrg, cnvrg, sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        return cudaStatus;
    }
    cudaStatus = cudaMalloc((void**)dev_n, sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        return cudaStatus;
    }
    cudaStatus = cudaMemcpy(*dev_n, size, sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        return cudaStatus;
    }

    return cudaStatus;
}



void printLattice(int8_t** lat, int* size) {
    int n = *size;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++)
        {
            printf("%d\t", lat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
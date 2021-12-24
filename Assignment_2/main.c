#include <stdlib.h>
#include <stdio.h>
#include <mpi/mpi.h>
#include "medianMPI.h"


int quickSelect(double * Arr,int* Indexes,int left ,int right ,int k);

double  calc_distance(double * A , double * B ,int dim);

void distributeByMedian(double * distance, double *rec_dist,int* indexes);

void median_QuickSelect(double * array,double *rec_dist);

int main(int argc, char* argv []){
   if (argc != 2) {
       printf("USAGE: ./bin/main <matrix.bin>  ");
       exit(1);
   }
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    int64_t meta_data[2]={-1};
    MPI_Status * status;
   
    int pivot=-1,rc ,root_rank=0, dest,tag=1;

    //open the input file to read the .bin matrix
    FILE * f_ptr;
    f_ptr=fopen(argv[1],"rb");

    rc=fread(meta_data,sizeof(int64_t),2,f_ptr);
    if (world_rank==0) printf(".bin points| Dimensions: %ld | Number %ld\n",meta_data[0],meta_data[1]);
    
    
    double **points , *dist_from_pivot , *pivot_point , *median;
    median=(double*)malloc(sizeof(double));
    pivot_point=(double*)malloc(meta_data[0]*sizeof(double));
    dist_from_pivot= (double*)malloc(4*sizeof(double));
    points= (double**)malloc(4*sizeof(double*));
    for (size_t i = 0; i < 4; i++)
    {
        points[i]=(double *)malloc(meta_data[0]*sizeof(double));
    }
    root_rank=0;
    if (world_rank==root_rank)
    {   
        pivot=rand()%(world_size*4);
        
        // dest=pivot/world_size;
        // printf("pivot: %d  Destination %d\n",pivot,dest);
        // pivot=pivot%4;
        // rc=MPI_Send(&pivot,1,MPI_INT,dest,tag,MPI_COMM_WORLD);

        

    }
    rc=MPI_Bcast(&pivot,1,MPI_INT,0,MPI_COMM_WORLD);
   
    printf("Rank: %d out of %d processors and i recieved the pivot %d\n",world_rank, world_size,pivot);
    
    
    
    unsigned long counter, counter2 ;
    
    for (size_t i = 0; i < 4; i++)
    {
     counter=0;
     counter2=0 ;
    fseek(f_ptr, 2 * sizeof(int64_t)+(world_rank+4*i)*sizeof(double), SEEK_SET);
    for (size_t j = 0; j <meta_data[0]; j++)
    {  
       rc=fseek(f_ptr, meta_data[0] * sizeof(double), SEEK_CUR); 
       rc=fread(&points[i][j],sizeof(double),1,f_ptr);
       if(points[i][j]!=0)  counter++;
       counter2++;
       
    }
   // printf("\n\n %ld %ld rank %d out of %d processors\n",counter,counter2,world_rank, world_size);
    }
    fclose(f_ptr);
    int source;
    if(pivot<(world_rank+1)*4 && pivot>=world_rank*4 )
    {
        int pivot_index= pivot %4 ;
        //printf("Rank: %d out of %d I HAVE the pivot %d in index %d \n",world_rank, world_size,pivot,pivot_index);
        pivot_point=points[pivot_index];
        source=world_rank;
    }
    rc=MPI_Bcast(pivot_point,meta_data[0],MPI_DOUBLE,source,MPI_COMM_WORLD);
    
    for (size_t i = 0; i < 4; i++)
    {
        dist_from_pivot[i]=calc_distance(points[i],pivot_point,meta_data[0]);
        //printf("Rank: %d distance %f\n",world_rank,dist_from_pivot[i]);
    }
    double *rec_dist;
    rec_dist=(double*)malloc(world_size*4*sizeof(double));
    int * indexes;
    indexes=(int*)malloc(world_size*4*sizeof(int));
    for (size_t i = 0; i < world_size*4; i++) indexes[i]=i;
    
    
    root_rank=0;
    if(MPI_Barrier(MPI_COMM_WORLD)==MPI_SUCCESS) printf("Rank %d Waiting\n",world_rank);
    // printf("Rank %d Done Waiting\n",world_rank);
    
    


    int median_index=0;
    distributeByMedian(dist_from_pivot,rec_dist,indexes);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    

    
    MPI_Barrier(MPI_COMM_WORLD);
    free(rec_dist);
    free(median);
    free(points);
    free(pivot_point);
    free(dist_from_pivot);

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}


double  calc_distance(double * A , double * B ,int dim)
{
    double distance=0;
    for (size_t i = 0; i < dim; i++)
    {
        distance= distance + (A[i]-B[i])*(A[i]-B[i]) ;
    }
    return distance;
    
}

void distributeByMedian(double * distance,  double * rec_dist,int* indexes)
{   
    int size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    int source,median_index;
    
    MPI_Barrier(MPI_COMM_WORLD);
    source=0;
    MPI_Gather(distance,4,MPI_DOUBLE,rec_dist,4,MPI_DOUBLE,source,MPI_COMM_WORLD);
    
    int left=0;
    int right=4*size-1;
    int k=2*size;
    if(rank==0)
    {   
        //median_index=quickSelect(rec_dist,indexes,left,right,k);
        for (size_t i = 0; i < size*4; i++)
        {
        printf("Rank: %d distance %f\n",rank,rec_dist[i]);
        }
        
    }
    median_index=quickSelect(rec_dist,indexes,left,right,k);
    MPI_Barrier(MPI_COMM_WORLD);
    return;

}



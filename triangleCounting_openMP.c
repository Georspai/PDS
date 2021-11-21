//Parallel and Distributed Systems
//
//Triangle Counting using OpenMP for Graph Adjacency Matrix.
//Spaias Georgios
//AEM: 8910

// OpenMP header 
#include <omp.h> 

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <sys/time.h>


//CSC format for Sparce Matrices
typedef struct csc_format 
{
    uint32_t * values ;
    uint32_t * row_index ;
    uint32_t * col_ptr ;
    
}csc_format;

//COO format for Sparce Matrices
typedef struct coo_format 
{
    uint32_t * values ;
    uint32_t * rows ;
    uint32_t * cols ;
    
}coo_format;


//Initializes csc_format pointer 
csc_format * csc_init(uint32_t n , uint32_t nnz);

//Computes the number of triangles in adjacency Matrix
void countTriangles(uint32_t *col_ptr,uint32_t *row_index , uint32_t n, uint32_t nnz,int numThreads);

//Gets the row and col for each non-zero element from .txt
coo_format * read_mat (FILE* file_ptr, uint32_t * size ,uint32_t *nnz);

//gets the size (n*m) and number of non zero elements from .txt
void get_matrix_info (FILE* file_ptr, uint32_t * size ,uint32_t * nnz);

//Copied from pds-codebase
//Transforms COO formatted sparse matrix to CSC
void coo2csc(
  uint32_t       * const row,       /*!< CSC row start indices */
  uint32_t       * const col,       /*!< CSC column indices */
  uint32_t const * const row_coo,   /*!< COO row indices */
  uint32_t const * const col_coo,   /*!< COO column indices */
  uint32_t const         nnz,       /*!< Number of nonzero elements */
  uint32_t const         n,         /*!< Number of rows/columns */
  uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
);

struct timeval tic();

double toc(struct timeval begin);



int main(int argc, char* argv [])
{   
  if (argc != 3) 
    {
        printf("USAGE: ./bin/main <matrix.txt> <number of threads> ");
        exit(1);
    }
    int numofthreads= atoi(argv[2]);


    uint32_t n,nnz;
    csc_format * Adj_mat ;
    coo_format * coo_matrix ;

    //open the input file to read the Adjacency matrix
    FILE * f_ptr=fopen(argv[1],"r");

    //get the Size and number of non zero elements from input matrix
    get_matrix_info ( f_ptr, &n , &nnz);
    printf("\nn=   %d\nm=   %d \nnnz= %d",n,n,nnz);

    //read adjacency matrix in COO format
    coo_matrix=read_mat (f_ptr,&n,&nnz);
    if (coo_matrix == NULL) 
    {
        fprintf(stderr, "main: COO Matrix Read failed.\n");
        exit(1);
    }

    fclose(f_ptr);

    printf("\nRead COO formatted matrix from txt...");


    Adj_mat = csc_init(n ,nnz);
    if (Adj_mat == NULL) 
    {
        fprintf(stderr, "main: Adjacency Matrix Init failed.\n");
        exit(1);
    }
    
    coo2csc(Adj_mat->row_index,Adj_mat->col_ptr,coo_matrix->rows,coo_matrix->cols,nnz,n,0);
    printf("\nCOO->CSC...");

    printf("\nTriangle Counting...");
    countTriangles(Adj_mat->col_ptr,Adj_mat->row_index, n, nnz,numofthreads);
    

    
    free(Adj_mat);
    free(coo_matrix);

    return 0;
}




void countTriangles(uint32_t *col_ptr,uint32_t *row_index , uint32_t n, uint32_t nnz,int numThreads)
{   
    double elapsed_time;
    uint32_t counter= 0 ;
    uint32_t *num_difTriangles;
    num_difTriangles=( uint32_t*)malloc(numThreads*sizeof(uint32_t));
    for(size_t i=0;i<numThreads;i++) num_difTriangles[i]=0;

    
    // number of threads given by argv[2]
    omp_set_num_threads(numThreads);

    struct timeval t_start=tic();
    /*~~~~~~~~~~~~~~~~~~~~~~~~Start of openMP parallel Code~~~~~~~~~~~~~~~~~~~~~~~~~*/ 
    #pragma omp parallel 
    {
    int numthrds ,tid;
    int nonZeroflag= 0;
    uint32_t val_res=0;
    tid=omp_get_thread_num();
    numthrds=omp_get_num_threads();
    
    
    for (size_t i = tid; i < n ; i=i+numthrds)
    {   val_res=0;
        for (size_t j = col_ptr[i]; j < col_ptr[i+1]; j++)
        {   
            for (size_t k = col_ptr[row_index[j]]; k < col_ptr[row_index[j]+1]; k++)
            {   
                for (size_t l = col_ptr[i]; l < col_ptr[i+1]; l++)
                {
                    if (row_index[l]==row_index[k])
                    {
                        val_res++;
                        nonZeroflag=1;
                        break;
                    }
                }        
            }
            if (nonZeroflag)  nonZeroflag=0 ;      
        } 
    num_difTriangles[tid]+=val_res;             
    }

    } 
    /*~~~~~~~~~~~~~~~~~~~~~~~End of openMP parallel Code~~~~~~~~~~~~~~~~~~~~~~~~~~*/   
    elapsed_time=toc(t_start);
    
    for (size_t i = 0; i < numThreads; i++)
    {
        counter=counter+num_difTriangles[i];
    }
    counter=counter/6;
    
    printf("\n\nNumber of different triangles: %d",counter);
   
    free(num_difTriangles);
    return;
}




//Memory allocation for CSC Matrix struct
csc_format * csc_init(uint32_t n , uint32_t nnz)
{
    csc_format * ptr ;
    ptr = (csc_format *) malloc(sizeof(csc_format));

    ptr->values = (uint32_t *) malloc(nnz * sizeof(uint32_t));
    ptr->row_index = (uint32_t *) malloc(nnz * sizeof(uint32_t));
    ptr->col_ptr = (uint32_t *) malloc((n+1) * sizeof(uint32_t)); 
    return ptr;
}




void coo2csc(
  uint32_t       * const row,       /*!< CSC row start indices */
  uint32_t       * const col,       /*!< CSC column indices */
  uint32_t const * const row_coo,   /*!< COO row indices */
  uint32_t const * const col_coo,   /*!< COO column indices */
  uint32_t const         nnz,       /*!< Number of nonzero elements */
  uint32_t const         n,         /*!< Number of rows/columns */
  uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
) {

  // ----- cannot assume that input is already 0!
  for (uint32_t l = 0; l < n+1; l++) col[l] = 0;


  // ----- find the correct column sizes
  for (uint32_t l = 0; l < nnz; l++)
    col[col_coo[l] - isOneBased]++;

  // ----- cumulative sum
  for (uint32_t i = 0, cumsum = 0; i < n; i++) {
    uint32_t temp = col[i];
    col[i] = cumsum;
    cumsum += temp;
  }
  col[n] = nnz;
  // ----- copy the row indices to the correct place
  for (uint32_t l = 0; l < nnz; l++) {
    uint32_t col_l;
    col_l = col_coo[l] - isOneBased;

    uint32_t dst = col[col_l];
    row[dst] = row_coo[l] - isOneBased;

    col[col_l]++;
  }
  // ----- revert the column pointers
  for (uint32_t i = 0, last = 0; i < n; i++) {
    uint32_t temp = col[i];
    col[i] = last;
    last = temp;
  }

}



coo_format * read_mat (FILE* file_ptr, uint32_t * size ,uint32_t *nnz)
{
    
    float val=0;
    uint32_t n= *size;
    uint32_t non_zero=*nnz;
    
    coo_format * ptr;
    ptr= (coo_format *)malloc(sizeof(coo_format));
    
    ptr->cols= (uint32_t* )malloc(non_zero*sizeof(uint32_t));
    ptr->rows= (uint32_t* )malloc(non_zero*sizeof(uint32_t));
    

    //fscanf (file_ptr,"%d %d %d\n",n,&m,nz);
    //printf("\nn=%d m=%d nnz=%d",*n,m,*nz);
    for (int i=0;i<non_zero/2;i++)
    {  
        if(!fscanf(file_ptr,"%d %d %g\n",&ptr->rows[i],&ptr->cols[i],&val)){
          perror("Reading matrix from .txt");
        }
        ptr->cols[non_zero/2+i]=ptr->rows[i];
        ptr->rows[non_zero/2+i]=ptr->cols[i];  
    }
   
    return ptr;
          
}

void get_matrix_info (FILE* file_ptr, uint32_t * size ,uint32_t * nnz)
{
    uint32_t *nz,*n,m;
    n=size;
    nz=nnz;

    if(!fscanf (file_ptr,"%d %d %d\n",n,&m,nz)){
      perror("Reading matrix from .txt");
    }
    *nz=(*nz)*2;
    
    return;
          
}

struct timeval tic()
{ 
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return tv;
}

double toc(struct timeval begin)
{
  struct timeval end;
  gettimeofday(&end, NULL);
  double stime=((double)end.tv_sec-(double)begin.tv_sec)*1000 +
                ((double)end.tv_usec-(double)begin.tv_usec)/1000 ;
  stime= stime/1000;
  printf("\nElapsed time :%g s",stime);
  return stime;
}
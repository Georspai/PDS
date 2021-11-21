//Parallel and Distributed Systems
//
//Triangle Counting using Pthreads for Graph Adjacency Matrix.
//Spaias Georgios
//AEM: 8910

//Pthread header
#include <pthread.h>

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

//Struct to pass data to threads 
typedef struct pthread_data
{
    csc_format * adj_mat ;
    uint32_t     n ;
    uint32_t     nnz ;
    uint32_t    num_difTriangles;
    uint32_t    thr_col_index;

    pthread_mutex_t * index_lock;
    
}pthread_data;


//Initializes csc_format pointer 
csc_format * csc_init(uint32_t n , uint32_t nnz);

//Computes the number of triangles in adjacency Matrix
void * countTriangles(void * args);

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
    

    uint32_t n,nnz;
    csc_format * Adjacency_mat ;
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

    Adjacency_mat = csc_init(n ,nnz);
    if (Adjacency_mat == NULL) 
    {
        fprintf(stderr, "main: Adjacency Matrix Init failed.\n");
        exit(1);
    }
   
    
    coo2csc(Adjacency_mat->row_index,Adjacency_mat->col_ptr,coo_matrix->rows,coo_matrix->cols,nnz,n,0);
    printf("\nCOO->CSC...");

    //number of Threads given as user input
    int thread_num= atoi(argv[2]);
    pthread_t tid[thread_num];
    

    //Thread Input data
    pthread_data  pth_data;
    pth_data.n      =n;
    pth_data.nnz    =nnz;
    pth_data.adj_mat=Adjacency_mat;
    pth_data.thr_col_index=0;
    pth_data.num_difTriangles=0;
   
    pth_data.index_lock= (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(pth_data.index_lock,NULL);
    
 
    printf("\nTriangle Counting... \n");
    for (int i = 0; i < thread_num; i++) {
    if (pthread_create(&tid[i], NULL, countTriangles,(void *) &pth_data)) {
      printf("Error creating threads %d\n", i);
    }
    }
    for (int i = 0; i < thread_num; i++) {
    if (pthread_join(tid[i],NULL)) {
      printf("Error joining thread %d\n", i);
    }
    }

    pth_data.num_difTriangles=pth_data.num_difTriangles/6;
    printf("\n\n Number of different triangles: %d",pth_data.num_difTriangles);


    pthread_mutex_destroy(pth_data.index_lock);

    //Deallocating Memory
    free(Adjacency_mat);
    free(coo_matrix);

    return 0;

}




void * countTriangles(void * args )
{   
    pthread_data * data = (pthread_data *)args;

    csc_format * adj           =  data->adj_mat;
    uint32_t n                 =  data->n;
    uint32_t nnz               =  data->nnz;
    uint32_t *col_index        = &data->thr_col_index;
    uint32_t *num_difTriangles = &data->num_difTriangles;

    int nonZeroflag= 0;
    double elapsed_time;
    uint32_t val_res=0;
    uint32_t c=0;

   
    
    struct timeval t_start=tic();
    uint32_t i = 0;
    while ( i<n && (*col_index)<n )
    {   val_res=0;

        pthread_mutex_lock(data->index_lock);
        i= (*col_index)++ ;
        pthread_mutex_unlock(data->index_lock);

        for (size_t j = adj->col_ptr[i]; j < adj->col_ptr[i+1]; j++)
        {       
            for (size_t k = adj->col_ptr[adj->row_index[j]]; k < adj->col_ptr[adj->row_index[j]+1]; k++)
            {   
                for (size_t l = adj->col_ptr[i]; l < adj->col_ptr[i+1]; l++)
                {
                    if (adj->row_index[l]==adj->row_index[k])
                    {
                        val_res++;
                        nonZeroflag=1;
                        break;
                    }
                }          
            }
            if (nonZeroflag)  nonZeroflag=0 ;
        } 

     c+=val_res;        
    }
    pthread_mutex_lock(data->index_lock);
    *num_difTriangles+=c;
    pthread_mutex_unlock(data->index_lock);
    elapsed_time=toc(t_start);   
    
    return (NULL);
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
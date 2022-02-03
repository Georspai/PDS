/*Parallel and Distributed Systems
//Ising Model V0
/*Spaias Georgios
AEM: 8910

*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <sys/time.h>
#include <math.h>

int8_t ** lattice_init( int n);

void print_lattice(int8_t ** lat , int n);

void compute_State(int8_t *** a,int8_t *** b, int n );

void Ising_comp(int8_t *** G, int n ,int k);

void start_state(int8_t ** lat , int n);

int8_t sgn(int8_t x);



void swap_dptr(int8_t *** x , int8_t *** y);

int main(int argc, char* argv []){

    int k=7,n=4 ;
    int8_t ** lat ;

    /*Allocate memory for the int ** pointer*/
    lat=lattice_init(n);

    //Initialize lattice with a random starting state
    start_state(lat, n);

    //Print Lattice
    print_lattice( lat ,  n);

    Ising_comp(&lat,n,k);

    print_lattice( lat ,  n);
    printf("\n Lattice pointer Address: %p \n",lat);


    free(lat);
    return 0;
}


//Allocate memory for the int ** pointer
int8_t ** lattice_init(int n){

    int8_t ** ptr;
    ptr= (int8_t **)malloc(n*sizeof(int8_t *));
    for (size_t i = 0; i < n; i++)
    {
         *(ptr+i)=(int8_t *)calloc(n, sizeof(int8_t));
    }
    return ptr ;
}

//Print Lattice
void print_lattice(int8_t ** lat , int n)
{
    for (size_t i = 0; i < n; i++)
    {
       for (size_t j = 0; j < n; j++)
       {
           printf("%d\t",lat[i][j]);
       }
       printf("\n");
    }


}

//Initialize lattice with a random starting state
void start_state(int8_t ** lat , int n)
{   srand(time(NULL));
    for (size_t i = 0; i < n; i++)
    {
       for (size_t j = 0; j < n; j++)
       {
           //lat[i][j]=1-2*(rand()%2);
           lat[i][j]=1-2*((i%2+j)%2);
       }

    }

 return ;
}


//Compute spin state of lattice after k steps
void compute_State(int8_t *** a,int8_t *** b, int n )
{   int8_t ** G=*a;
    int8_t ** G_trans=*b;

    int n_i=0, n_j=0;


    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            G_trans[i][j]=0;
            //neighbor 1 -->(i+1,j)
            n_i=(i+1)%n;
            n_j= j ;
            G_trans[i][j]+=G[n_i][n_j];

            //neighbor 2 -->(i,j+1)
            n_i=i;
            n_j= (j+1)%n ;
            G_trans[i][j]+=G[n_i][n_j];

            //neighbor 3 -->(i-1,j)
            n_i=(i-1+n)%n;
            n_j= j ;
            G_trans[i][j]+=G[n_i][n_j];

            //neighbor 4 -->(i,j-1)
            n_i=i;
            n_j= (j-1+n)%n;
            G_trans[i][j]+=G[n_i][n_j];

            //+self -->(i,j)
            G_trans[i][j]+=G[i][j];

            //
            G_trans[i][j]=sgn(G_trans[i][j]);

        }

    }

}

void Ising_comp(int8_t *** G, int n ,int k)
{
    int8_t *** G_trans =(int8_t ***)malloc(sizeof(int8_t **));
    *G_trans=lattice_init(n);

    for (size_t i = 0; i < k; i++)
    {   printf("\n\n~~~ Step %ld ~~~\n",i);
        compute_State(G,G_trans, n );
        print_lattice(*G,n);
        swap_dptr(G,G_trans);
        printf("\nMem address G %p G_trans %p \n", *G, *G_trans);
    }

    free(G_trans);

}

int8_t sgn(int8_t x)
{
    int8_t one = (int8_t)(x>0)-(int8_t)(x<0);
    return one;
}

void swap_dptr(int8_t *** x , int8_t *** y)
{

    int8_t ** temp ;
    temp=*x ;
    *x=*y;
    *y=temp ;

}

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "medianMPI.h"




void swapElements(double* i, double* j)
{   //printf("Swaaap\n");
    double temp;
    temp = *i;
    *i = *j;
    *j = temp;
}
void swapIndex(int* i, int* j)
{   //printf("SwaaapPAPA\n");
    int temp;
    temp = *i;
    *i = *j;
    *j = temp;
}



int quickSelect(double * Arr,int * Indexes ,int left ,int right ,int k){
    
    int pivotIndex=-1;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    printf("rank%d in QS\n",rank);
    while (rank==0)
    {
        if(left==right)   return left;
        pivotIndex=pivot(Arr,Indexes,left,right);
        pivotIndex=partition(Arr,Indexes ,left,right,pivotIndex,k);

        if(pivotIndex==k){
           return k;
        }
        else if (k<pivotIndex){
            right = pivotIndex -1;
        }
        else{
            left = pivotIndex +1;
        }    
    }
    
    printf("rank%d DONE WITH QS\n",rank);
    return pivotIndex;
}




int partition(double * Arr,int * Indexes, int left , int right , int pivotIndex, int k){
    
    double pivotValue = Arr[pivotIndex];
    swapElements(Arr + pivotIndex,Arr+right);
    swapIndex(Indexes+pivotIndex,Indexes+right);

    int storeIndexEq=left,storeIndex=left;

    // Move all elements smaller than the pivot to the left of the pivot
    for (int i = left; i < right; i++)
    {
        if (Arr[i]<pivotValue)
        {
            swapElements(Arr+storeIndex,Arr+i);
            swapIndex(Indexes+storeIndex,Indexes+i);
            storeIndex++;
        }   
    }
    // Move all elements equal to the pivot right after the smaller elements
    storeIndexEq=storeIndex;
    for (int i = storeIndex; i < right; i++)
    {
        if (Arr[i]==pivotValue)
        {
           swapElements(Arr+storeIndexEq,Arr+i);
           swapIndex(Indexes+storeIndex,Indexes+i);
           storeIndexEq++;
        }    
    }
    swapElements(Arr+right,Arr+storeIndexEq);
    swapIndex(Indexes+right,Indexes+storeIndexEq);

    // Return location of pivot considering the desired location n
    if (k<storeIndex){
        return storeIndex;
    }
    if (k<=storeIndexEq){
        return k;
    }   
}




int pivot( double * Arr,int* Indexes,int left , int right)
{   
    // for 5 or less elements just get median
    if (right-left<5)
    {
        return insertionSort(Arr,Indexes,left,right);
    }
    
    // otherwise move the medians of five-element subgroups to the first n/5 positions
    for (int i = left; i < right; )
    {
        // get the median position of the i'th five-element subgroup
        int subRight = i + 4;
        if (subRight>right)    subRight=right;  
        int median5 =insertionSort(Arr,Indexes,i,subRight);
        swapElements(Arr+median5,Arr+left+(int)floor((double)((i-left)/5))); 
        swapIndex(Indexes+median5,Indexes+left+(int)floor((double)((i-left)/5)));
        i=i+5;
    }
    // compute the median of the n/5 medians-of-five
    int mid = (right - left) / 10 + left + 1;
    return quickSelect(Arr, Indexes, left , left+(int)floor((double)((right-left)/5)),mid);   
}



int insertionSort(double * Arr,int* Indexes, int left , int right )
{   
    int i = left + 1;
    while(i <= right){
        int j = i;
        while(j > left && Arr[j-1] > Arr[j]){
            swapElements(Arr+j-1 , Arr+j);
            swapIndex(Indexes+j-1,Indexes+j);
            j--;
        }
        i++;
    }
    return (int)floor((left + right) / 2);
}
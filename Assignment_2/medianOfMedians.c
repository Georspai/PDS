#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "median.h"

int Indexes[8];

int main(){

    double A[8]={2.0,4.0,13.0,9.0,3.0,6.0,8.0,94.0};
    for (size_t i = 0; i < 8; i++)
    {
        Indexes[i]=i;
    }
    
    int left=0;
    int right=7;
    int median=(int)((left+right)/2);
    int median_final=quickSelect(A,left,right,median);
    printf("Final median: Value %f index %d\n",A[median_final],median_final);
    for (int i = 0; i < right+1; i++)
    {
        printf(" Value %f index %d -----Indexes %d\n",A[i],i,Indexes[i]);
    }
    
    return 0;
}



void swapElements(double* i, double* j)
{   printf("Swap\n");
    double temp;
    temp = *i;
    *i = *j;
    *j = temp;
}
void swapIndex(int* i, int* j)
{   printf("Swap Indexes\n");
    int temp;
    temp = *i;
    *i = *j;
    *j = temp;
}

// function select(list, left, right, n)
//     loop
//         if left = right then
//             return left
//         pivotIndex := pivot(list, left, right)
//         pivotIndex := partition(list, left, right, pivotIndex, n)
//         if n = pivotIndex then
//             return n
//         else if n < pivotIndex then
//             right := pivotIndex - 1
//         else
//             left := pivotIndex + 1

int quickSelect(double * Arr,int left ,int right ,int k){
    printf("quickSelect\n");
    int pivotIndex;

    while (right>=left)
    {
        if(left==right)   return left;
        pivotIndex=pivot(Arr,left,right);
        pivotIndex=partition(Arr, left,right,pivotIndex,k);

        if(pivotIndex==k){
           return k;
        }
        else if (k<pivotIndex){
            right= pivotIndex -1;
        }
        else{
            left= pivotIndex +1;
        }    
    }
    return pivotIndex;
}


// function partition(list, left, right, pivotIndex, n)
//     pivotValue := list[pivotIndex]
//     swap list[pivotIndex] and list[right]  // Move pivot to end
//     storeIndex := left
//     // Move all elements smaller than the pivot to the left of the pivot
//     for i from left to right − 1 do
//         if list[i] < pivotValue then
//             swap list[storeIndex] and list[i]
//             increment storeIndex
//     // Move all elements equal to the pivot right after
//     // the smaller elements
//     storeIndexEq = storeIndex
//     for i from storeIndex to right − 1 do
//         if list[i] = pivotValue then
//             swap list[storeIndexEq] and list[i]
//             increment storeIndexEq
//     swap list[right] and list[storeIndexEq]  // Move pivot to its final place
//     // Return location of pivot considering the desired location n
//     if n < storeIndex then
//         return storeIndex  // n is in the group of smaller elements
//     if n ≤ storeIndexEq then
//         return n  // n is in the group equal to pivot


int partition(double * Arr, int left , int right , int pivotIndex, int k){
    printf("Partition\n");
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


// function pivot(list, left, right)
//     // for 5 or less elements just get median
//     if right − left < 5 then
//         return partition5(list, left, right)
//     // otherwise move the medians of five-element subgroups to the first n/5 positions
//     for i from left to right in steps of 5
//         // get the median position of the i'th five-element subgroup
//         subRight := i + 4
//         if subRight > right then
//             subRight := right
//         median5 := partition5(list, i, subRight)
//         swap list[median5] and list[left + floor((i − left)/5)]

//     // compute the median of the n/5 medians-of-five
//     mid := (right − left) / 10 + left + 1
//     return select(list, left, left + floor((right − left) / 5), mid)

int pivot( double * Arr,int left , int right)
{   printf("Pivot\n");
    // for 5 or less elements just get median
    if (right-left<5)
    {
        return insertionSort(Arr,left,right);
    }
    
    // otherwise move the medians of five-element subgroups to the first n/5 positions
    for (int i = left; i < right; )
    {
        // get the median position of the i'th five-element subgroup
        int subRight = i + 4;
        if (subRight>right)    subRight=right;  
        int median5 =insertionSort(Arr,i,subRight);
        swapElements(Arr+median5,Arr+left+(int)floor((double)((i-left)/5))); 
        swapIndex(Indexes+median5,Indexes+left+(int)floor((double)((i-left)/5)));
        i=i+5;
    }
    // compute the median of the n/5 medians-of-five
    int mid = (right - left) / 10 + left + 1;
    return quickSelect(Arr, left , left+(int)floor((double)((right-left)/5)),mid);   
}


// function partition5( list, left, right)
//     i := left + 1
//     while i ≤ right
//         j := i
//         while j > left and list[j−1] > list[j] do
//             swap list[j−1] and list[j]
//             j := j − 1
//         i :=  i + 1
            
//     return floor((left + right) / 2)

int insertionSort(double * Arr, int left , int right )
{   printf("InsertionSort\n");
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

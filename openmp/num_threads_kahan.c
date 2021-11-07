#include <stdio.h>
#include <omp.h>


int main (void){

    int num_threads = omp_get_max_threads();

    printf("Número de threads: %d\n", num_threads);


    return 0;
}
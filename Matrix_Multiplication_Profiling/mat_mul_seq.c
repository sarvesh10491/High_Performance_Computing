#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAT_SIZE 1024

static double matA[MAT_SIZE][MAT_SIZE];
static double matB[MAT_SIZE][MAT_SIZE];
static double matC[MAT_SIZE][MAT_SIZE];


void init_matrix(){
    for(int i=0; i<MAT_SIZE; i++)
        for(int j=0; j<MAT_SIZE; j++){
            matA[i][j] = rand() % 10;
            matB[i][j] = rand() % 10;
            matC[i][j] = 0.0;
        }
}


void mult_matrix(){
    for(int i=0; i<MAT_SIZE; i++)
        for(int j=0; j<MAT_SIZE; j++)
            for(int k=0; k<MAT_SIZE; k++)
                matC[i][j] += matA[i][k] * matB[k][j];
}


void print_matrix() {
    printf("Matrix A\n");
    for(int i=0; i<MAT_SIZE; i++){
        for(int j=0; j<MAT_SIZE; j++)
            printf(" %7.2f", matA[i][j]);
        printf("\n");
    }
    printf("\n\n");
    
    printf("Matrix B\n");
    for(int i=0; i<MAT_SIZE; i++){
        for(int j=0; j<MAT_SIZE; j++)
            printf(" %7.2f", matB[i][j]);
        printf("\n");
    }
    printf("\n\n");
    
    printf("Matrix C\n");
    for(int i=0; i<MAT_SIZE; i++){
        for(int j=0; j<MAT_SIZE; j++)
            printf(" %7.2f", matC[i][j]);
        printf("\n");
    }
    printf("\n\n");
}


static inline void timespec_diff(struct timespec *a, struct timespec *b, struct timespec *result){
    result->tv_sec  = b->tv_sec - a->tv_sec;
    result->tv_nsec = b->tv_nsec - a->tv_nsec;

    if(result->tv_nsec < 0){
        --result->tv_sec;
        result->tv_nsec += 1000000000L;
    }
}


int main(){
    time_t t;
    srand((unsigned) time(&t));

    struct timespec initStart, initEnd, initDiff;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &initStart);
    init_matrix();
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &initEnd);
    timespec_diff(&initStart, &initEnd, &initDiff);
    // print_matrix();
    printf("Matrix initialize time : ");
    printf("%lld.%.9ld sec\n\n", (long long)initDiff.tv_sec, initDiff.tv_nsec);

    struct timespec multStart, multEnd, multDiff;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &multStart);
    mult_matrix();
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &multEnd);
    timespec_diff(&multStart, &multEnd, &multDiff);
    print_matrix();
    printf("Matrix multiplication time : ");
    printf("%lld.%.9ld sec\n\n", (long long)multDiff.tv_sec, multDiff.tv_nsec);

    return 0;
}
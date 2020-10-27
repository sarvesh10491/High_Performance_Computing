#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>

#define MAT_SIZE    1024
#define NUM_THREADS 2

static double matA[MAT_SIZE][MAT_SIZE];
static double matB[MAT_SIZE][MAT_SIZE];
static double matC[MAT_SIZE][MAT_SIZE];

struct t_data
{
    int t_id;
    int t_num;
};

pthread_cond_t start;           // Condition variable for start sync
pthread_mutex_t start_mut;      // Mutex used in conjunction with condition variables

int chunk, row_start, row_end;
chunk = MAT_SIZE / NUM_THREADS;


void init_matrix(){
    for(int i=0; i<MAT_SIZE; i++)
        for(int j=0; j<MAT_SIZE; j++){
            matA[i][j] = rand() % 10;
            matB[i][j] = rand() % 10;
            matC[i][j] = 0.0;
        }
}


void *mult_matrix(void *argptr)
{
    // Thread function computations
    struct t_data *tptr;
    tptr = (struct t_data *)argptr;

    row_start = (tptr->t_num) * chunk;
    row_end = (tptr->t_num + 1) * chunk;

    printf("Periodic thread invoked : %d\n", tptr->t_num);
    printf("row_start : %d || row_end : %d\n", row_start, row_end);

    // Each thread wait here till conditional variable broadcasts
    pthread_mutex_lock(&start_mut);
    int cret = pthread_cond_wait(&start, &start_mut);       // All threads waiting on start condition variable to start simultaneously
    if(cret!=0)
        printf("Periodic condition variable error %d\n",cret);
    pthread_mutex_unlock(&start_mut);
    
    printf("Calculating chunk : %d\n", tptr->t_num);
    for(int i=row_start; i<row_end; i++)
        for(int j=0; j<MAT_SIZE; j++)
            for(int k=0; k<MAT_SIZE; k++)
                matC[i][j] += matA[i][k] * matB[k][j];
    
    printf("Periodic thread exiting : %d\n", tptr->t_num);   

    pthread_exit(NULL);
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


    int tret=-1;
    pthread_t tid[NUM_THREADS];
    struct t_data *dptr[NUM_THREADS];

    pthread_attr_t fifo_attr;
    pthread_attr_init(&fifo_attr);
    pthread_attr_setschedpolicy(&fifo_attr,SCHED_FIFO);

    int cret = pthread_cond_init(&start, NULL);           // initialising condition variables
    if(cret!=0)
        printf("Start condn init error %d\n",cret);
    pthread_mutex_init(&start_mut, NULL);                 // Initialising condition variable mutexes(thread mutex are different)


    for(int i=0; i<NUM_THREADS; i++)
    {
        dptr[i] = (struct t_data *)malloc(sizeof(struct t_data));
        dptr[i]->t_id = tid[i];
        dptr[i]->t_num = i;

        tret = pthread_create(&tid[i], &fifo_attr, mult_matrix, (void *)dptr[i]);
        if(tret)
        {
            printf("ERROR. Return code from pthread_create() is %d\n", tret);
            exit(-1);
        }
    }

    usleep(1000);       // allowing all threads to be created and waitinng on condition variable before broadcasting

    struct timespec multStart, multEnd, multDiff;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &multStart);
    pthread_cond_broadcast(&start);                    // all threads start running
    for(int i=0; i<NUM_THREADS; i++){
        pthread_join(tid[i], NULL);
        free(dptr[i]);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &multEnd);
    timespec_diff(&multStart, &multEnd, &multDiff);
    print_matrix();
    printf("Matrix multiplication time : ");
    printf("%lld.%.9ld sec\n\n", (long long)multDiff.tv_sec, multDiff.tv_nsec);

    pthread_cond_destroy(&start);
    printf("Terminated successfully.\n");

    return 0;
}
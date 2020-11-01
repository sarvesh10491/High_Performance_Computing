ifeq ($(mode),0)
CC = gcc
		TARGET = mat_mul_seq
endif

ifeq ($(mode),1)
CC = gcc
LIBS = -pthread -lm -lrt
		TARGET = mat_mul_pthread
endif

ifeq ($(mode),2)
CC = gcc
LIBS = -fopenmp -lm -lrt
		TARGET = mat_mul_openmp
endif

ifeq ($(mode),3)
CC = mpicc
LIBS = -lm -lrt
		TARGET = mat_mul_mpi
endif

CFLAGS  = -Wall

all: $(TARGET)

$(TARGET): $(TARGET).c $(DEPS)
		$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c $(LIBS)

clean:
		$(RM) mat_mul_seq mat_mul_pthread mat_mul_openmp mat_mul_mpi
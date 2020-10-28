ifeq ($(mode),0)
		TARGET = mat_mul_seq
endif

ifeq ($(mode),1)
LIBS = -pthread -lm -lrt
		TARGET = mat_mul_pthread
endif

ifeq ($(mode),2)
LIBS = -fopenmp -lm -lrt
		TARGET = mat_mul_openmp
endif

# the compiler: gcc for C program
CC = gcc
CFLAGS  = -Wall
# DEPS = queue_ops.h

all: $(TARGET)

$(TARGET): $(TARGET).c $(DEPS)
		$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c $(LIBS)

clean:
		$(RM) $(TARGET)
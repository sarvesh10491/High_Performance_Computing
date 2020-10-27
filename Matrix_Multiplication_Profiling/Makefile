# the compiler: gcc for C program
CC = gcc
# LIBS = -pthread -lm -lrt
CFLAGS  = -Wall
# DEPS = queue_ops.h
TARGET = mat_mul_seq

all: $(TARGET)

$(TARGET): $(TARGET).c $(DEPS)
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c $(LIBS)

clean:
	$(RM) $(TARGET)
ifeq ($(mode),0)
		TARGET = mat_mul_seq
else
		TARGET = mat_mul_pthread
endif

# the compiler: gcc for C program
CC = gcc
LIBS = -pthread -lm -lrt
CFLAGS  = -Wall
# DEPS = queue_ops.h

all: $(TARGET)

$(TARGET): $(TARGET).c $(DEPS)
		$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c $(LIBS)

clean:
		$(RM) $(TARGET)
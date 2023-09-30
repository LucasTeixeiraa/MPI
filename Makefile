CC   = gcc
MPICC=mpicc
CFLAGS=-g -O
LDFLAGS=-lm
OBJS=*.o
BINS= serial_stencil 

# Habilite essa opcao caso tenha compilacao com MPI
USE_MPI = 1

ifeq ($(USE_MPI),1)
	CC=mpicc
	CFLAGS+=-DHAVE_MPI
endif

all: $(BINS)

%.o: %.c 
	$(CC) $(CFLAGS) $< -c -o $@

# serial_stencil: serial_stencil.o save_array.o appctx.o appctx.h Makefile
# 	$(CC) $(CFLAGS) -o $@ serial_stencil.o save_array.o appctx.o $(LDFLAGS)

mpi_stencil: mpi_stencil.o save_array.o appctx.o appctx.h Makefile
	$(CC) $(CFLAGS) -o $@ serial_stencil.o save_array.o appctx.o $(LDFLAGS)


clean:
	rm -f $(BINS) $(OBJS)
	rm -f output*.bmp *.txt *.bin *.vti 

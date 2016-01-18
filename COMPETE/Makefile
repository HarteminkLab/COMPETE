CC=gcc 
CFLAGS=-O3 -funroll-loops -I./libconfig/libconfig-1.1_inst/include
LFLAGS=-lm -L./libconfig/libconfig-1.1_inst/lib -lconfig -lpthread

all: compete

compete: bc.o compete.o
	$(CC) $(CFLAGS) -o compete bc.o compete.o $(LFLAGS)

clean:
	rm -f compete compete.o bc.o

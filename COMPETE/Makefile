CC=gcc 
CFLAGS=-O3 -funroll-loops -I./libconfig/libconfig-1.1_inst/include
LFLAGS=-lm -L./libconfig/libconfig-1.1_inst/lib -lconfig -lpthread

all: model2dot compete

compete: bc.o compete.o
	$(CC) $(CFLAGS) -o compete bc.o compete.o $(LFLAGS)

model2dot: bc.o model2dot.o
	$(CC) $(CFLAGS) -o model2dot bc.o model2dot.o $(LFLAGS)

clean:
	rm -f compete compete.o model2dot model2dot.o bc.o

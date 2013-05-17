CC=gcc
CFLAGS=-Wall -I/usr/include -I.
LFLAGS=-lm -L/usr/lib -lgsl -lgslcblas -L.
modulos= correlation.o 


%.out:%.o $(modulos)
	$(CC) $^ $(LFLAGS) -o $@

clean:
	rm -rf *.o
	rm -rf *.out
	rm -rf *.dat	


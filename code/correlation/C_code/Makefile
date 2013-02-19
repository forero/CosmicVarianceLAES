CC=gcc
CFLAGS=-Wall -I/usr/include -I.
LFLAGS=-lm -L/usr/lib -lgsl -lgslcblas -L.
modulos= mass_distance_selector.o randomsample.o correlation.o mass_function.o LAE_catalog_generator.o correlation_statistic.o


%.out:%.o $(modulos)
	$(CC) $^ $(LFLAGS) -o $@

clean:
	rm -rf *.o
	rm -rf *.out
	rm -rf *.dat	x


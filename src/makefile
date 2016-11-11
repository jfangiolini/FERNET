CC = gcc
CLIBS = -lgsl -lgslcblas -lm -largtable2 -lconfig -ltiff
CFLAGS = -Wall -std=gnu99 -pedantic

fernet: fernet.o point.o multi.o line.o parseconfig.o raster.o stack.o spim.o parseargs.o photons.o orbit.o fernet.h
	$(CC) $(CFLAGS) -o fernet fernet.o multi.o point.o line.o raster.o stack.o spim.o parseconfig.o parseargs.o photons.o orbit.o $(CLIBS)

point.o: point.c fernet.h
	$(CC) $(CFLAGS) -c point.c

multi.o: multi.c fernet.h
	$(CC) $(CFLAGS) -c multi.c

line.o: line.c fernet.h
	$(CC) $(CFLAGS) -c line.c

raster.o: raster.c fernet.h
	$(CC) $(CFLAGS) -c raster.c

stack.o: stack.c fernet.h
	$(CC) $(CFLAGS) -c stack.c

spim.o: spim.c fernet.h
	$(CC) $(CFLAGS) -c spim.c

fernet.o: fernet.c fernet.h
	$(CC) $(CFLAGS) -c fernet.c

parseconfig.o: parseconfig.c fernet.h
	$(CC) $(CFLAGS) -c parseconfig.c

parseargs.o: parseargs.c fernet.h
	$(CC) $(CFLAGS) -c parseargs.c

photons.o: photons.c fernet.h
	$(CC) $(CFLAGS) -c photons.c

orbit.o: orbit.c fernet.h
	$(CC) $(CFLAGS) -c orbit.c

clean:
	-@rm -rf *.o fernet 2>/dev/null || true

install:
	sudo cp fernet /usr/local/bin

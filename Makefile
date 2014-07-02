CC=gcc
CFLAGS=-Imatrix -Imaterial-data/choi-okos -Imaterial-data/pasta -I. -ggdb -Wall
LDFLAGS=-lm
VPATH=matrix material-data material-data/pasta 

all: visco elastic 

# Make stuff from other projects using their makefile
matrix.a:
	$(MAKE) -C matrix matrix.a
	cp matrix/matrix.a .

material-data.a: matrix.a
	cp matrix/matrix.a material-data
	$(MAKE) -C material-data material-data.a
	cp material-data/material-data.a .

crank.o: drying.h
visco.o: drying.h
stress.o: drying.h
elastic.o: drying.h

visco: visco.o stress.o crank.o matrix.a material-data.a
	$(CC) -o $@ $(CFLAGS) $^ $(LDFLAGS)

elastic: elastic.o stress.o crank.o matrix.a material-data.a
	$(CC) -o $@ $(CFLAGS) $^ $(LDFLAGS)

doc: Doxyfile
	doxygen Doxyfile

clean:
	rm -rf *.o *.a visco elastic
	$(MAKE) -C material-data clean
	$(MAKE) -C matrix clean


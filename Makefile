CC=gcc
CFLAGS=-Imatrix -Imaterial-data/choi-okos -Imaterial-data/pasta -I. -ggdb -Wall
LDFLAGS=-lm
VPATH=matrix material-data material-data/pasta 

all: dry

# Make stuff from other projects using their makefile
matrix.a:
	$(MAKE) -C matrix matrix.a
	cp matrix/matrix.a .

material-data.a: matrix.a
	cp matrix/matrix.a material-data
	$(MAKE) -C material-data material-data.a
	cp material-data/material-data.a .

crank.o: drying.h
main.o: drying.h
stress.o: drying.h

dry: main.o stress.o crank.o matrix.a material-data.a
	$(CC) -o dry $(CFLAGS) $^ $(LDFLAGS)

doc: Doxyfile
	doxygen Doxyfile

clean:
	rm -rf *.o *.a dry
	$(MAKE) -C material-data clean
	$(MAKE) -C matrix clean


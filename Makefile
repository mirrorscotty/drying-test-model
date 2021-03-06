CC=gcc
CFLAGS=-Imatrix -Imaterial-data -I. -ggdb -Wall
LDFLAGS=-lm
VPATH=matrix material-data material-data/pasta viscoelastic

all: visco-flux visco-profile

# Make stuff from other projects using their makefile
matrix.a:
	$(MAKE) -C matrix matrix.a
	cp matrix/matrix.a .

material-data.a: matrix.a
	cp matrix/matrix.a material-data
	$(MAKE) -C material-data material-data.a
	cp material-data/material-data.a .

strain.o: drying.h
stress.o: drying.h
flux.o:
moisture.o: drying.h
profile.o:

crank.o: drying.h
force.o: drying.h
elastic.o: drying.h

visco-flux: flux.o strain.o stress.o force.o moisture.o crank.o matrix.a material-data.a
	$(CC) -o $@ $(CFLAGS) $^ $(LDFLAGS)

visco-profile: profile.o strain.o force.o moisture.o crank.o matrix.a material-data.a
	$(CC) -o $@ $(CFLAGS) $^ $(LDFLAGS)

elastic: elastic.o force.o crank.o matrix.a material-data.a
	$(CC) -o $@ $(CFLAGS) $^ $(LDFLAGS)

doc: Doxyfile
	doxygen Doxyfile

clean:
	rm -rf *.o *.a visco-flux visco-profile elastic
	$(MAKE) -C material-data clean
	$(MAKE) -C matrix clean


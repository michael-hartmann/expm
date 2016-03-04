PATH_LIBHADES=../../libhades/src/
SOURCES = ${PATH_LIBHADES}/libhades.c  ${PATH_LIBHADES}/expm.c  ${PATH_LIBHADES}/odeint.c  ${PATH_LIBHADES}/optimize.c ${PATH_LIBHADES}/parse_npy_dict.c
CFLAGS += -std=c99 -malign-double -Wall -Wmissing-prototypes -Wstrict-prototypes -I${PATH_LIBHADES}

clean:
	rm -f *.o
	
all:
	$(CC) $(OPT) $(CFLAGS) $(LDFLAGS) -c ${SOURCES} expm_wrapper.c
	ld -r *.o -o packed.o
	mv packed.o expm.o
	f95 -O3 -Wall test.f expm.o -lm -llapack -lblas

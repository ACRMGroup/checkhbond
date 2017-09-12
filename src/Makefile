#COPTS    = -I$(HOME)/include -L$(HOME)/lib -g -Wall -pedantic -ansi
COPTS     = -I$(HOME)/include -L$(HOME)/lib -O3 -Wall -pedantic -ansi
NOWARN    = -Wno-unused-but-set-variable
CHBCOMMON = residues.o orientate.o
HMCOMMON  = orientate.o cavallo_userfunc.o
BINDIR    = ../bin
CC	  = gcc
LIBS      = -lbiop -lgen -lm -lxml2

EXE = hydrogen_matrices \
     hydrogen_matrices_Ndonor \
     hydrogen_matrices_Oacceptor \
     hydrogen_matrices_SCMC \
     checkhbond \
     checkhbond_Ndonor \
     checkhbond_Oacceptor

all : $(EXE)


clean :
	\rm -f $(EXE) $(CHBCOMMON) $(HMCOMMON) \
	hydrogen_matrices.o \
	hydrogen_matrices_Ndonor.o \
	hydrogen_matrices_Oacceptor.o \
	hydrogen_matrices_SCMC.o \
	checkhbond.o \
	checkhbond_Ndonor.o \
	checkhbond_Oacceptor.o

hydrogen_matrices :  hydrogen_matrices.o $(HMCOMMON)
	$(CC) $(COPTS) -o $@  hydrogen_matrices.o $(HMCOMMON) $(LIBS)

hydrogen_matrices_Ndonor : hydrogen_matrices_Ndonor.o $(HMCOMMON)
	$(CC) $(COPTS) -o $@ hydrogen_matrices_Ndonor.o $(HMCOMMON) $(LIBS)

hydrogen_matrices_Oacceptor : hydrogen_matrices_Oacceptor.o $(HMCOMMON)
	$(CC) $(COPTS) -o $@ hydrogen_matrices_Oacceptor.o $(HMCOMMON) $(LIBS)

hydrogen_matrices_SCMC : hydrogen_matrices_SCMC.o $(HMCOMMON)
	$(CC) $(COPTS) -o $@ hydrogen_matrices_SCMC.o $(HMCOMMON) $(LIBS)

checkhbond : checkhbond.o $(CHBCOMMON)
	$(CC) $(COPTS) -o $@ checkhbond.o $(CHBCOMMON) $(LIBS)

checkhbond_Ndonor : checkhbond_Ndonor.o $(CHBCOMMON)
	$(CC) $(COPTS) -o $@ checkhbond_Ndonor.o $(CHBCOMMON) $(LIBS)

checkhbond_Oacceptor : checkhbond_Oacceptor.o $(CHBCOMMON)
	$(CC) $(COPTS) -o $@ checkhbond_Oacceptor.o $(CHBCOMMON) $(LIBS)


hydrogen_matrices.o : hydrogen_matrices.c orientate.h hbondmat2.h cavallo_userfunc.h
	$(CC) -c $(COPTS) $(NOWARN) -o $@ hydrogen_matrices.c

hydrogen_matrices_Ndonor.o : hydrogen_matrices.c orientate.h hbondmat2.h cavallo_userfunc.h
	$(CC) -D MCDONOR -c $(COPTS) -o $@ hydrogen_matrices.c

hydrogen_matrices_Oacceptor.o : hydrogen_matrices.c orientate.h hbondmat2.h cavallo_userfunc.h
	$(CC) -D MCACCEPTOR -c $(COPTS) $(NOWARN) -o $@ hydrogen_matrices.c

hydrogen_matrices_SCMC.o : hydrogen_matrices.c orientate.h hbondmat2.h cavallo_userfunc.h
	$(CC) -D SCMC -c $(COPTS) $(NOWARN) -o $@ hydrogen_matrices.c

checkhbond.o : checkhbond.c hbondmat2.h
	$(CC) -c $(COPTS) -o $@ checkhbond.c 

checkhbond_Ndonor.o : checkhbond.c hbondmat2.h
	$(CC) -D MCDONOR -c $(COPTS) -o $@ checkhbond.c 

checkhbond_Oacceptor.o : checkhbond.c hbondmat2.h
	$(CC) -D MCACCEPTOR -c $(COPTS) -o $@ checkhbond.c 

.c.o :
	$(CC) -c $(COPTS) -o $@ $<

install :
	mkdir -p $(BINDIR)
	cp hydrogen_matrices $(BINDIR)
	cp hydrogen_matrices_Ndonor $(BINDIR)
	cp hydrogen_matrices_Oacceptor $(BINDIR)
	cp hydrogen_matrices_SCMC $(BINDIR)
	cp checkhbond $(BINDIR)
	cp checkhbond_Ndonor $(BINDIR)
	cp checkhbond_Oacceptor $(BINDIR)

#COPTS    = -I$(HOME)/include -L$(HOME)/lib -g -Wall -pedantic -ansi
COPTS    = -I$(HOME)/include -L$(HOME)/lib -O3 -Wall -pedantic -ansi
CHBCOMMON = residues.o orientate.o
HMCOMMON = orientate.o cavallo_userfunc.o
BINDIR   = ../bin
CC	 = gcc
LIBS     = -lbiop -lgen -lm -lxml2

all: $(BINDIR)/hydrogen_matrices \
     $(BINDIR)/hydrogen_matrices_Ndonor \
     $(BINDIR)/hydrogen_matrices_Oacceptor \
     $(BINDIR)/hydrogen_matrices_SCMC \
     $(BINDIR)/checkhbond \
     $(BINDIR)/checkhbond_Ndonor \
     $(BINDIR)/checkhbond_Oacceptor

clean :
	\rm -f $(CHBCOMMON) $(HMCOMMON) \
	hydrogen_matrices.o \
	hydrogen_matrices_Ndonor.o \
	hydrogen_matrices_Oacceptor.o \
	hydrogen_matrices_SCMC.o \
	checkhbond.o \
	checkhbond_Ndonor.o \
	checkhbond_Oacceptor.o

$(BINDIR)/hydrogen_matrices :  hydrogen_matrices.o $(HMCOMMON)
	$(CC) $(COPTS) -o $@  hydrogen_matrices.o $(HMCOMMON) $(LIBS)

$(BINDIR)/hydrogen_matrices_Ndonor : hydrogen_matrices_Ndonor.o $(HMCOMMON)
	$(CC) $(COPTS) -o $@ hydrogen_matrices_Ndonor.o $(HMCOMMON) $(LIBS)

$(BINDIR)/hydrogen_matrices_Oacceptor : hydrogen_matrices_Oacceptor.o $(HMCOMMON)
	$(CC) $(COPTS) -o $@ hydrogen_matrices_Oacceptor.o $(HMCOMMON) $(LIBS)

$(BINDIR)/hydrogen_matrices_SCMC : hydrogen_matrices_SCMC.o $(HMCOMMON)
	$(CC) $(COPTS) -o $@ hydrogen_matrices_SCMC.o $(HMCOMMON) $(LIBS)

$(BINDIR)/checkhbond : checkhbond.o $(CHBCOMMON)
	$(CC) $(COPTS) -o $@ checkhbond.o $(CHBCOMMON) $(LIBS)

$(BINDIR)/checkhbond_Ndonor : checkhbond_Ndonor.o $(CHBCOMMON)
	$(CC) $(COPTS) -o $@ checkhbond_Ndonor.o $(CHBCOMMON) $(LIBS)

$(BINDIR)/checkhbond_Oacceptor : checkhbond_Oacceptor.o $(CHBCOMMON)
	$(CC) $(COPTS) -o $@ checkhbond_Oacceptor.o $(CHBCOMMON) $(LIBS)


hydrogen_matrices.o : hydrogen_matrices.c orientate.h hbondmat2.h cavallo_userfunc.h
	$(CC) -c $(COPTS) -o $@ hydrogen_matrices.c

hydrogen_matrices_Ndonor.o : hydrogen_matrices.c orientate.h hbondmat2.h cavallo_userfunc.h
	$(CC) -D MCDONOR -c $(COPTS) -o $@ hydrogen_matrices.c

hydrogen_matrices_Oacceptor.o : hydrogen_matrices.c orientate.h hbondmat2.h cavallo_userfunc.h
	$(CC) -D MCACCEPTOR -c $(COPTS) -o $@ hydrogen_matrices.c

hydrogen_matrices_SCMC.o : hydrogen_matrices.c orientate.h hbondmat2.h cavallo_userfunc.h
	$(CC) -D SCMC -c $(COPTS) -o $@ hydrogen_matrices.c

checkhbond.o : checkhbond.c hbondmat2.h
	$(CC) -c $(COPTS) -o $@ checkhbond.c 

checkhbond_Ndonor.o : checkhbond.c hbondmat2.h
	$(CC) -D MCDONOR -c $(COPTS) -o $@ checkhbond.c 

checkhbond_Oacceptor.o : checkhbond.c hbondmat2.h
	$(CC) -D MCACCEPTOR -c $(COPTS) -o $@ checkhbond.c 

.c.o :
	$(CC) -c $(COPTS) -o $@ $<


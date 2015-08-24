COPTS    = -I/home/bsm/martin/include -L/home/bsm/martin/lib -g
CHBFILES = checkhbond.o residues.o orientate.o
HMCOMMON = orientate.o cavallo_userfunc.o
BINDIR   = ../bin

all: $(BINDIR)/hydrogen_matrices \
     $(BINDIR)/hydrogen_matrices_Ndonor \
     $(BINDIR)/hydrogen_matrices_Oacceptor \
     $(BINDIR)/checkhbond

clean :
	\rm -f $(CHBFILES) $(HMFILES) \
        hydrogen_matrices.o \
        hydrogen_matrices_Ndonor.o \
        hydrogen_matrices_Oacceptor.o

$(BINDIR)/hydrogen_matrices :  hydrogen_matrices.o $(HMCOMMON)
	cc $(COPTS) -o $@  hydrogen_matrices.o $(HMCOMMON) -lbiop -lgen -lm

$(BINDIR)/hydrogen_matrices_Ndonor : hydrogen_matrices_Ndonor.o $(HMCOMMON)
	cc $(COPTS) -o $@ hydrogen_matrices_Ndonor.o $(HMCOMMON) -lbiop -lgen -lm

$(BINDIR)/hydrogen_matrices_Oacceptor : hydrogen_matrices_Oacceptor.o $(HMCOMMON)
	cc $(COPTS) -o $@ hydrogen_matrices_Oacceptor.o $(HMCOMMON) -lbiop -lgen -lm

$(BINDIR)/checkhbond : $(CHBFILES)
	cc $(COPTS) -o $@ $(CHBFILES) -lbiop -lgen -lm

hydrogen_matrices.o : hydrogen_matrices.c orientate.h hbondmat2.h cavallo_userfunc.h
	cc -c $(COPTS) -o $@ hydrogen_matrices.c

hydrogen_matrices_Ndonor.o : hydrogen_matrices_Ndonor.c orientate.h hbondmat2.h cavallo_userfunc.h
	cc -c $(COPTS) -o $@ hydrogen_matrices_Ndonor.c

hydrogen_matrices_Oacceptor.o : hydrogen_matrices_Oacceptor.c orientate.h hbondmat2.h cavallo_userfunc.h
	cc -c $(COPTS) -o $@ hydrogen_matrices_Oacceptor.c

checkhbond.o : checkhbond.c hbondmat2.h
	cc -c $(COPTS) -o $@ checkhbond.c 

.c.o :
	cc -c $(COPTS) -o $@ $<


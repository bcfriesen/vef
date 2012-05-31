FC = gfortran
FCFLAGS = -c -g -O0
LD = $(FC)
LDFLAGS = -g -O0
LAPACK = -L/home/friesen/lapack-3.4.1/build/lib -llapack -lblas

SOURCES = \
	main.f90
OBJECTS = $(SOURCES:.f90=.o)
EXECUTABLE = run

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(LD) $(LDFLAGS) $(OBJECTS) -o $@ $(LAPACK)

%.o: %.f90
	$(FC) $(FCFLAGS) $< -o $@

clean:
	rm -rf $(EXECUTABLE) *.o *.dat *.mod

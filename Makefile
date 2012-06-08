FC = nagfor
FCFLAGS = -C=all -g -O0 -r8 -c
LAPACK = -llapack -lblas
LDFLAGS = -C=all -g -O0 -r8

#FC = ifort
#FCFLAGS = -warn all -g -O0 -r8 -c
#LAPACK = -mkl
#LDFLAGS = -warn all -g -O0 -r8

LD = $(FC)

SOURCES = \
	interfaces.f90 \
	global.f90 \
	twerp.f90 \
	locate.f90 \
	stop_exit.f90 \
	gauleg.f90 \
	planck_fn.f90 \
	formal_soln.f90 \
	calc_moments.f90 \
	calc_feautrier_vars.f90 \
	calc_source_fn.f90 \
	calc_vefs.f90 \
	solve_scatt_prob.f90 \
	main.f90
OBJECTS = $(SOURCES:.f90=.o)
EXECUTABLE = run

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(LD) $(LDFLAGS) $(OBJECTS) -o $@ $(LAPACK)

%.o: %.f90
	$(FC) $(FCFLAGS) $< -o $@

clean:
	rm -rf $(EXECUTABLE) *.o *.dat *.mod *__genmod* fort.*

FC = nagfor
FCFLAGS = -C=all -g -g90 -O0 -c
LAPACK = -llapack -lblas
LDFLAGS = -C=all -g -g90 -O0

#FC = ifort
#FCFLAGS = -warn all -g -O0 -c
#LAPACK = -mkl
#LDFLAGS = -warn all -g -O0

LD = $(FC)

SOURCES = \
	precision_mod.f90 \
	interfaces.f90 \
	global.f90 \
	write_vefs.f90 \
	write_source_fn.f90 \
	write_moments.f90 \
	twerp.f90 \
	locate.f90 \
	stop_exit.f90 \
	gauleg.f90 \
	planck_fn.f90 \
	calc_moments.f90 \
	calc_source_fn.f90 \
	calc_vefs.f90 \
	solve_scatt_prob.f90 \
	solve_rte.f90 \
	main.f90
OBJECTS = $(SOURCES:.f90=.o)
EXECUTABLE = run

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(LD) $(LDFLAGS) $(OBJECTS) -o $@ $(LAPACK)

%.o: %.f90
	$(FC) $(FCFLAGS) $< -o $@

clean:
	rm -rf $(EXECUTABLE) *.o *.dat *.mod *__genmod* fort.* *.g90

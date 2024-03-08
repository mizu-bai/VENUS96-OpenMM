FC = gfortran
FFLAGS = -std=legacy -fno-automatic -O2

.PHONY: all clean

all: venus96.x venus96_openmm.x

venus96.x: venus96.f
	$(FC) $(FFLAGS) venus96.f -o venus96.x

venus96_openmm.x: venus96_openmm.f openmm_driver.f90
	$(FC) $(FFLAGS) venus96_openmm.f openmm_driver.f90 -o venus96_openmm.x

clean:
	rm -rf *.o *.mod *.x

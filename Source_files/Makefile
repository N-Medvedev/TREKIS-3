# !***************************************************************
# ! This file is part of TREKIS-3
# !***************************************************************
# Simple compilation:

default:
	@echo "Compiling TREKIS-3 code"

# For compilation with gfortran, the following options are required:
# -ffree-line-length-0 -fdec-format-defaults -cpp -fdefault-real-8 -fdefault-double-8
# -D OMP_inside

	ifort -qopenmp -O5 -fpp -cpp -ipo -real-size 64 -standard-semantics -assume nofpe_summary Universal_MC_for_SHI_MAIN.f90 -o TREKIS.x
	
	@echo "Executable: TREKIS.x"

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.yaml
	rm -f *.optrpt
	rm -f TREKIS.x

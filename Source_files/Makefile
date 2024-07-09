# !***************************************************************
# ! This file is part of TREKIS-3
# !***************************************************************
# Simple compilation:

default:
	@echo "Compiling TREKIS-3 code"

# For compilation with gfortran, the following options are required:
# -ffree-line-length-0 -fdec-format-defaults -cpp -fdefault-real-8 -fdefault-double-8
# -D OMP_inside

# mpiifort -O5 -fpp -cpp -ipo -real-size 64 -D MPI_USED -standard-semantics -assume nofpe_summary Universal_MC_for_SHI_MAIN.f90 -o TREKIS_MPI.x

# mpifort -O5 -fpp -cpp -ipo -real-size 64 -D MPI_USED -standard-semantics -assume nofpe_summary Universal_MC_for_SHI_MAIN.f90 -o TREKIS_MPI.x

ifeq ($(mpi),y)
	mpifort -O5 -fpp -cpp -ipo -real-size 64 -D MPI_USED -standard-semantics -assume nofpe_summary -ffree-line-length-0 -fdec-format-defaults -cpp -fdefault-real-8 -fdefault-double-8 Universal_MC_for_SHI_MAIN.f90 -o TREKIS_MPI.x
	@echo "----------------------------------------------"
	@echo "With MPI: Executable: TREKIS_MPI.x"
	@echo "Run the code as mpiexec -np [n] ./TREKIS_MPI.x"
	@echo "----------------------------------------------"
endif

ifeq ($(OMP),y) # OMP=no means no omp
	ifort -O5 -fpp -cpp -ipo -real-size 64 -standard-semantics -assume nofpe_summary Universal_MC_for_SHI_MAIN.f90 -o TREKIS.x
	@echo "----------------------------------------------"
	@echo "With OpenMP: Executable: TREKIS.x"
	@echo "Run the code as ./TREKIS.x"
	@echo "----------------------------------------------"
endif

ifeq ($(non),y) # OMP=no means no omp
	ifort -qopenmp -O5 -fpp -cpp -ipo -real-size 64 -standard-semantics -assume nofpe_summary Universal_MC_for_SHI_MAIN.f90 -o TREKIS.x
	@echo "----------------------------------------------"
	@echo "No parallelization: Executable: TREKIS.x"
	@echo "Run the code as ./TREKIS.x"
	@echo "----------------------------------------------"
endif

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.yaml
	rm -f *.optrpt
	rm -f *.x

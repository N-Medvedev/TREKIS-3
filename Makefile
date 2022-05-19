# !***************************************************************
# ! This file is part of TREKIS-3
# !***************************************************************
# Simple compilation:

default:
	@echo "Compiling TREKIS-3 code"
	
	ifort -qopenmp -O5 -fpp -ipo -real-size 64 -standard-semantics Universal_MC_for_SHI_MAIN.F90 -o TREKIS.x
	
	@echo "Executable: TREKIS.x"

clean:
	rm -f *.o
	rm -f *.mod
	rm -f TREKIS.x
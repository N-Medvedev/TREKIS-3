# !***************************************************************
# ! This file is part of TREKIS-3
# !***************************************************************

all:
	ifort -qopenmp -O5 -fpp -vec-report0 -ipo -real-siz 64 Universal_MC_for_SHI_MAIN.F90 -o TREKIS.x

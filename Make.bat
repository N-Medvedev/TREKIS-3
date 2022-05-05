:: !***************************************************************
:: ! This file is part of TREKIS-3
:: !***************************************************************

ifort.exe /F9999999999 /fpp /Qopenmp /O3 /Qvec-report0 /Qipo /real-size:64 Universal_MC_for_SHI_MAIN.F90 -o TREKIS.exe


del *.obj *.mod

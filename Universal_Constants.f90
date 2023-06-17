!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains all the global variables:universal parameters

MODULE Universal_Constants
 implicit none

real(8) :: g_Pi, g_e, g_re, g_me, g_lc, g_cvel, g_Mp, g_h, g_kb, g_e0, g_mu0, g_Ry, g_a0, g_v0, g_Na, g_alpha, g_G   ! universal constants 
real(8) :: g_me_eV, g_me_MeV, g_r0

complex :: g_CI

!UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
! Universal constants:
parameter (g_Pi = 3.1415926535897932384626433832795d0)            ! Pi
parameter (g_e = 1.602176487d-19)         ! Electron charge       [Coulomb]
parameter (g_me = 9.1093821545d-31)       ! Electron mass         [kg]
parameter (g_re = 2.8179403267d-5)        ! Classical electron radius [A]
parameter (g_cvel = 299792458.0d0)        ! Light velosity        [m/sec]
parameter (g_lc = 2.4263102389d-2)        ! Compton wavelength [A]
parameter (g_Mp =1836.1526724780d0*g_me)  ! Proton mass           [kg]
parameter (g_h = 1.05457162853d-34)       ! Plank constant        [J*sec]
parameter (g_kb = 11604.51928260d0)       ! Boltzmann constant    [K/eV]
parameter (g_Na = 6.0221413d+23)          ! Avogadro's number     [mol^-1]
parameter (g_e0 = 8.854187817620d-12)     ! Electrical constant   [F/m]
parameter (g_mu0 = 1.2566370614359d-6)    ! Magnetic constant     [H*A^-2]
parameter (g_Ry = 13.606d0)               ! Rydberg constant      [eV]
parameter (g_a0 = 0.5291772085936d0)      ! Bohr radius           [A]
parameter (g_v0 = dsqrt(2.0d0*g_Ry*g_e/g_me))    ! Bohr velosity         [m/s]
parameter (g_alpha = g_e*g_e/(g_h*g_cvel*4.0d0*g_Pi*g_e0))    ! Fine structure constant
parameter (g_G = 6.67384d-11)             ! Newtonian constant of gravitation [m3/kg/s^2]
parameter (g_r0		= g_alpha*g_alpha*g_a0)	! classical electron radius [A]

parameter (g_me_MeV	= 0.51099906d0)		! Electron mass	[MeV]
parameter (g_me_eV	= g_me_MeV*1.0d6)	! Electron mass in [eV]; or mc^2

parameter (g_CI = dcmplx(0.0d0,1.0d0) )	! complex unity
!UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU

END MODULE Universal_Constants

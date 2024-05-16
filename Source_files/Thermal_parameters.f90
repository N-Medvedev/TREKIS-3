!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains subroutines to calculate thermal parameters:
! - electron-ion (electron-phonon) coupling,
! - electron heat capacity,
! - electron thermal conductivity,
! as functions of the electronic temperature
!***************************************************************
! NOT READY YET
!--------------

module Thermal_parameters

use Universal_Constants
use Objects
use Reading_files_and_parameters, only : print_time_step, Find_in_array_monoton
use Variables, only: dashline, starline
use Cross_sections, only : Tot_EMFP, CS_from_MFP, Equilibrium_charge_Target


implicit none

!==============================================



private  ! hides items not listed on public statement

public :: Get_thermal_parameters

character(100), parameter :: m_thermal = 'OUTPUT_thermal'


contains


subroutine Get_thermal_parameters(Output_path, CDF_Phonon, Target_atoms, Matter, NumPar, Mat_DOS, File_names)
   character(*), intent(in) :: Output_path
   type(CDF), intent(in) :: CDF_Phonon      ! phononic part of the CDF
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects
   type(Solid), intent(in) :: Matter     ! all material parameters
   type(Flag), intent(in) :: NumPar         ! numerical parameters
   type(Density_of_states), intent(in) :: Mat_DOS       ! Density of states
   type(All_names), intent(inout) :: File_names    ! file names to use later for gnuplot printing
   !--------------------------
   type(diff_CS_single), dimension(:), allocatable :: ds_E
   real(8), dimension(:), allocatable :: Te_grid
   real(8) :: Eemax, E_min, dE, E_cur, Zt, Zeff, Sigma, dEdx, E_rate, G_coupling, Dens, hw, d_hw, V, DOS, Ffactor
   integer :: FN, i, Nsiz, Nsiz_ds, j, i_T, k, N_siz_DOS
   character(300) :: File_name, File_short
   logical :: file_exist, file_opened

   ! Check if the user wants to calculate the thermal parameters:
   if (.not. NumPar%get_thermal) then
        if (NumPar%verbose) call print_time_step('Thermal parameters calculation is skipped:', msec=.true.)
        return ! exit, if they don't
   endif
   if (NumPar%verbose) call print_time_step('Thermal parameters calculation started:', msec=.true.)

   ! make filename for output:
   call create_filename(Output_path, NumPar, Matter, File_name, File_short) ! below
   ! Save it for gnuplotting, if required:
   if (allocated(File_names%F)) then
      File_names%F(30) = trim(adjustl(File_short)) ! save for later use
   endif

   ! open the corresponding file:
   FN = 8756
   open(FN, file=trim(adjustl(File_name)))
   write(FN, '(a)') '# Te   Coupling'
   write(FN, '(a)') '# [K]  [W/(K*m^3)]'

   ! Get the thermal parameters:
   ! Create the grid of electron temperatures:
   call get_Te_grid(Te_grid)    ! below
   Nsiz = size(Te_grid)

   ! Define maximal relevant electron energy:
   E_min = 0.0d0    ! starting point [eV]
   Eemax = 10.0d0 * (Te_grid(Nsiz)*g_kb_EV) ! [eV]
   dE = 1.0d0   ! assuming step of 1 eV
   Nsiz_ds = ceiling(Eemax/dE)
   allocate(ds_E(Nsiz_ds))

   ! 1) Electron-ion (electron-phonon) coupling:
   ! 1.1) Precalculate electronic differential cross section of elastic scattering:
   do i = 1, Nsiz_ds
      allocate(ds_E(i)%dsdhw(10000))
      allocate(ds_E(i)%hw(10000))
      E_cur = E_min + dE*dble(i)    ! [eV]
      ds_E(i)%E = E_cur ! total energy [eV]

      ! Target mean atomic charge:
      Zt = SUM(target_atoms(:)%Zat*dble(target_atoms(:)%Pers))/dble(SUM(target_atoms(:)%Pers)) ! mean atomic number of target atoms
      select case(numpar%CDF_elast_Zeff)
      case (0) ! Barkas-like charge
         Zeff = 1.0d0 + Equilibrium_charge_Target(E_cur, g_me, Zt, (Zt-1.0e0), 0, 1.0d0) ! Equilibrium charge, see below
      case (1) ! one, as used in old CDF expression : Fourier of the unscreened Coulomb with charge Z=1
         Zeff = 1.0d0    ! electron charge
      case default ! full charge will be screened with the electronic CDF(q,w) inside the cross section
         Zeff = 1.0d0    ! so no additional multiplication needed here
      end select

      ! Get differential mean free path [A]:
      call Tot_EMFP(E_cur, Target_atoms, CDF_Phonon, Matter, Sigma, dEdx, NumPar, Mat_DOS, 'electron', Zeff, &
                    dsdhw=ds_E(i)%dsdhw, d_hw = ds_E(i)%hw)   ! module "Cross_sections"

      ! Convert it into CS:
      do j = 1, size(ds_E(i)%dsdhw)
         ds_E(i)%dsdhw(j) = CS_from_MFP(ds_E(i)%dsdhw(j), Matter%At_Dens) * 1.0d-20  ! [m^2] module "Cross_sections"
      enddo ! j = 1, size(ds_E(i)%dsdhw)
   enddo ! i = 1, Nsiz_ds

   if (NumPar%verbose) call print_time_step('Thermal parameters CS precalculated:', msec=.true.)

   ! atomic density:
   Dens = (Matter%At_Dens*1.0d6)    ! [1/m^3]

   N_siz_DOS = size(Mat_DOS%E)  ! size of the DOS array

   ! 1.2) Calculate the coupling parameter:
   do i_T = 1, Nsiz     ! for all electronic temperature points
      E_rate = 0.0d0    ! to start with

      !$omp parallel &
      !$omp private (i, dE, DOS, V, j, d_hw, k)
      !$omp do schedule(dynamic) reduction(+: E_rate )
      do i = 1, Nsiz_ds ! integral by energy

         if (i > 1) then
            dE = ds_E(i)%E - ds_E(i-1)%E
         else
            dE = ds_E(i)%E
         endif

         ! Get material DOS:
         if (ds_E(i)%E > Mat_DOS%E(N_siz_DOS)) then   ! free-electron DOS
            DOS = sqrt(ds_E(i)%E) / sqrt(Mat_DOS%E(N_siz_DOS))*maxval(Mat_DOS%DOS(:)) ! normalization of choice (maybe changed later)
         else   ! material DOS:
            call Find_in_array_monoton(Mat_DOS%E, ds_E(i)%E, k) ! module "Reading_files_and_parameters"
            DOS = Mat_DOS%DOS(k)    ! [electrons/atom/eV]
         endif

         ! Electron velosity:
         V = SQRT(2.0d0*ds_E(i)%E*g_e/g_me) ! [m/s]

         do j = 1, size(ds_E(i)%dsdhw)  ! integral by transferred energy
            if (j > 1) then
               d_hw = ds_E(i)%hw(j) - ds_E(i)%hw(j-1)
            else
               d_hw = ds_E(i)%hw(j)
            endif

            ! Combination of distribution functions:
            Ffactor = 1.0d0

            ! Energy exchange rate:
            E_rate = E_rate + Dens * ds_E(i)%dsdhw(j) * ds_E(i)%hw(j) * dE * V * DOS * Ffactor
         enddo ! j = 1, size(ds_E(i)%dsdhw)
      enddo ! i = 1, Nsiz_ds
      !$omp end do
      !$omp end parallel

      ! Convert into coupling parameter:
      G_coupling = (E_rate * g_e) * Dens / ( Te_grid(i_T) - Matter%temp )

      ! printout the calculated coupling parameter:
      write(FN, '(e20.7, e20.7)') Te_grid(i_T), G_coupling
   enddo ! i_T = 1, Nsiz

   if (NumPar%verbose) then
      call print_time_step('Thermal parameters calculations finished :', msec=.true.)
   endif
   write(*,'(a,a)') 'Thermal parameters are stored in the file ', trim(adjustl(File_name))

   ! clean up:
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)
end subroutine Get_thermal_parameters



subroutine get_Te_grid(Te_grid)
    real(8), dimension(:), allocatable, intent(inout) :: Te_grid
    !-----------------
    integer :: i, Nsiz
    real(8) :: Te, dTe, Tmin, Tmax

    Tmin = 300.0d0      ! [K]
    Tmax = 50000.0d0    ! [K]
    dTe = 1000.0d0      ! [K] grid step
    Tmin = max(Tmin,dTe)
    Nsiz = ceiling( (Tmax - Tmin)/dTe )

    allocate(Te_grid(Nsiz))

    Te_grid(1) = Tmin
    do i = 2, Nsiz
        Te_grid(i) = Te_grid(i-1) + dTe
    enddo
end subroutine get_Te_grid



subroutine create_filename(Output_path, NumPar, Matter, File_name, File_short)
    character(*), intent(in) :: Output_path
    type(Flag), intent(in) :: NumPar         ! numerical parameters
    type(Solid), intent(in) :: Matter     ! all material parameters
    character(*), intent(out) :: File_name, File_short
    !--------------------------
    character(100) :: temp_char, temp_char1, KCS

    ! Create a file to print thermal parameters into:
    KCS = ''
    select case (NumPar%kind_of_CDF_ph)
    case (0)    ! Ritchie-Howie
        KCS = '_CDF_'
    case (1)    ! single-pole CDF
        KCS = '_spCDF_'
    endselect

    ! Info about temperature
    write(temp_char, '(f7.2, a)') Matter%temp, '_K'

    ! Add info about effective charge:
    select case(NumPar%CDF_elast_Zeff)
    case default
        temp_char1 = 'Zeff_'//trim(adjustl(temp_char))
    case (1)
        temp_char1 = 'Z=1_'//trim(adjustl(temp_char))
    case (2)
        temp_char1 = 'Z_CDFe_FF_'//trim(adjustl(temp_char))
    case (3)
        temp_char1 = 'Z_CDFe_'//trim(adjustl(temp_char))
    endselect

    File_short = trim(adjustl(m_thermal))//trim(adjustl(KCS))//trim(adjustl(temp_char1))//'.dat'
    File_name = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(File_short))
end subroutine create_filename


end module Thermal_parameters

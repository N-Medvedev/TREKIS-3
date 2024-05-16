!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains subroutines to create Gnuplot scripts
! To be useful, it needs installed Gnuplot prior to execution of TREKIS
! Gnuplot is a freeware program distributed here:
! http://www.gnuplot.info/

module Gnuplotting_subs
use Objects, only : Flag, All_names, Atom, Ion, Solid

#ifndef __GFORTRAN__
   USE IFLPORT, only : system, chdir   ! library, allowing to operate with directories in intel fortran
#endif

use Dealing_with_EADL, only : Count_lines_in_file, Count_columns_in_file


implicit none
PRIVATE  ! hides items not listed on public statement

! this is a function that returns the order of a passed number:
interface find_order_of_number
   module procedure find_order_of_number_real
   module procedure find_order_of_number_int
end interface find_order_of_number

public :: Gnuplot_ion, Gnuplot_electron_hole, Gnuplot_transients


!----------------------------------------------
! Reminder:
! It is using the following order of file names
! type(All_names), intent(out) :: File_names:
! File_names%F(1) = Gnuplot path
! File_names%F(2) = electrons IMFP
! File_names%F(3) = holes IMFP
! File_names%F(4) = electrons EMFP
! File_names%F(5) = holes EMFP
! File_names%F(6) = SHI dEdx; MFP; Range
! File_names%F(7) = Photon IMFP
! File_names%F(8) = Photon EMFP
! File_names%F(9) = DOS
! File_names%F(10) = Output directory name for transients
! File_names%F(11) = Total_numbers
! File_names%F(12) = Total_energies
! File_names%F(13) = Electronic spectrum
! File_names%F(14) = Emitted electrons spectrum
! File_names%F(15) = Valence holes spectrum
! File_names%F(16) = Electron velosity angular sitribution
! File_names%F(17) = Valence hole velosity angular sitribution
! File_names%F(18) = Electron raidal sitribution
! File_names%F(19) = Electron energy raidal sitribution
! File_names%F(20) = Radial_electron_temperature
! File_names%F(21) = Radial_photon_density
! File_names%F(22) = Radial_photon_energy
! File_names%F(23) = Radial_Lattice_energy
! File_names%F(24) = Radial_Track_energy
! File_names%F(25) = Radial_Lattice_temperature
! File_names%F(26) = VB Holes density vs R
! File_names%F(27) = Radial_holes_pot_energy
! File_names%F(28) = Radial_holes_kin_energy
! File_names%F(29) = Radial_holes_temperature
! File_names%F(30) = Thermal_parameters

!----------------------------------------------

contains


!------------------------------
! Dynamical quantities:
subroutine Gnuplot_transients(Tim, NumPar, Matter, Target_atoms, File_names)
   real(8), intent(in) :: Tim ! total simulation time [fs]
   type(Flag), intent(in) :: NumPar
   type(Solid), intent(in) :: Matter   ! all material parameters
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   type(All_names), intent(in) :: File_names   ! all file names for printing out stuff
   !-----------------------------
   character(300) :: output_path, Filename, In_file
   character(10) :: call_slash, sh_cmd
   integer :: FN, ext_ind, leng
   logical :: file_opened

   ! Set the output directory path:
   output_path = trim(adjustl(File_names%F(10)))//trim(adjustl(NumPar%path_sep))

   ! Find the extension of the gnuplot scripts:
   call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"

   !--------------
   ! 1) Print total numbers:
   In_file = trim(adjustl(File_names%F(11)))
   leng = LEN(trim(adjustl(In_file)))
   Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-4)))//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_total_numbers(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar)  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it


   !--------------
   ! 2) Print total energies:
   In_file = trim(adjustl(File_names%F(12)))
   leng = LEN(trim(adjustl(In_file)))
   Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-4)))//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_total_NRG(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), trim(adjustl(File_names%F(11))), NumPar) ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it


   !--------------
   ! 3) Print electron spectrum:
   In_file = trim(adjustl(File_names%F(13)))
   leng = LEN(trim(adjustl(In_file)))
   Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-4)))//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_spectrum(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar)  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it

   ! 3.1) Print emitted electron spectrum:
   if (Matter%work_function .GT. 0.0d0) then
      In_file = trim(adjustl(File_names%F(14)))
      leng = LEN(trim(adjustl(In_file)))
      Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-4)))//trim(adjustl(sh_cmd))
      open(newunit=FN, FILE = trim(adjustl(Filename)))
      call gnuplot_spectrum(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar)  ! below
      inquire(unit=FN,opened=file_opened)    ! check if this file is opened
      if (file_opened) close(FN)             ! and if it is, close it
   endif

   ! 3.2) Print valence holes spectrum:
   In_file = trim(adjustl(File_names%F(15)))
   leng = LEN(trim(adjustl(In_file)))
   Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-4)))//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_spectrum(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar, .true.)  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it


   !--------------
   ! 4) Print electron velosity angular distribution:
   In_file = trim(adjustl(File_names%F(16)))
   leng = LEN(trim(adjustl(In_file)))
   Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-4)))//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_theta_distribution(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar)  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)

   ! 4.1) Print VB hole velosity angular distribution:
   In_file = trim(adjustl(File_names%F(17)))
   leng = LEN(trim(adjustl(In_file)))
   Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-4)))//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_theta_distribution(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar)  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)

   !--------------
   ! 5) Print electron density radial distribution:
   In_file = trim(adjustl(File_names%F(18)))
   leng = LEN(trim(adjustl(In_file)))
   Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-13)))//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_raidal_distribution(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar, 1.0d16, &
            'Electron density (1/cm^3)')  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)

   ! 5.1) Print VB holes density radial distribution:
   In_file = trim(adjustl(File_names%F(26)))
   leng = LEN(trim(adjustl(In_file)))
   Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-13)))//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_raidal_distribution(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar, 1.0d16, &
            'Valence holes density (1/cm^3)')  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)

   ! 5.2) Print photon density radial distribution:
   if (NumPar%include_photons) then
      In_file = trim(adjustl(File_names%F(21)))
      leng = LEN(trim(adjustl(In_file)))
      Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-13)))//trim(adjustl(sh_cmd))
      open(newunit=FN, FILE = trim(adjustl(Filename)))
      call gnuplot_raidal_distribution(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar, 1.0d16, &
            'Photons density (1/cm^3)')  ! below
      inquire(unit=FN,opened=file_opened)    ! check if this file is opened
      if (file_opened) close(FN)
   endif

   !--------------
   ! 6) Print electron energy density radial distribution:
   In_file = trim(adjustl(File_names%F(19)))
   leng = LEN(trim(adjustl(In_file)))
   Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-13)))//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_raidal_distribution(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar, 1.0d-6, &
            'Electrons energy density (eV/A^3)')  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)

   ! 6.1) Print VB hole potential energy density radial distribution:
   In_file = trim(adjustl(File_names%F(27)))
   leng = LEN(trim(adjustl(In_file)))
   Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-13)))//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_raidal_distribution(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar, 1.0d-6, &
            'Holes potential energy density (eV/A^3)')  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)
   ! 6.2) Print VB hole kinetic energy density radial distribution:
   In_file = trim(adjustl(File_names%F(28)))
   leng = LEN(trim(adjustl(In_file)))
   Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-13)))//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_raidal_distribution(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar, 1.0d-6, &
            'Holes kinetic energy density (eV/A^3)')  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)

   ! 6.3) Print atomic energy density radial distribution:
   In_file = trim(adjustl(File_names%F(23)))
   leng = LEN(trim(adjustl(In_file)))
   Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-13)))//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_raidal_distribution(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar, 1.0d-6, &
            'Atoms kinetic energy density (eV/A^3)')  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)
   ! 6.4) Print total atomic energy density radial distribution (scattering + nonthermal):
   In_file = trim(adjustl(File_names%F(24)))
   leng = LEN(trim(adjustl(In_file)))
   Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-13)))//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_raidal_distribution(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar, 1.0d-6, &
            'Track energy density (eV/A^3)')  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)

   ! 6.5) Print photon energy density radial distribution:
   if (NumPar%include_photons) then
      In_file = trim(adjustl(File_names%F(22)))
      leng = LEN(trim(adjustl(In_file)))
      Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-13)))//trim(adjustl(sh_cmd))
      open(newunit=FN, FILE = trim(adjustl(Filename)))
      call gnuplot_raidal_distribution(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar, 1.0d-6, &
            'Track energy density (eV/A^3)')  ! below
      inquire(unit=FN,opened=file_opened)    ! check if this file is opened
      if (file_opened) close(FN)
   endif

   !--------------
   ! 7) Print electron temperature radial distribution:
   In_file = trim(adjustl(File_names%F(20)))
   leng = LEN(trim(adjustl(In_file)))
   Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-4)))//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_raidal_distribution(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar, 1.0d1, &
            'Electron kinetic temperature (K)')  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)

   ! 7.1) Print VB holes temperature radial distribution:
   In_file = trim(adjustl(File_names%F(29)))
   leng = LEN(trim(adjustl(In_file)))
   Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-4)))//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_raidal_distribution(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar, 1.0d1, &
            'Valence holes kinetic temperature (K)')  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)

   ! 7.2) Print atomic temperature radial distribution:
   In_file = trim(adjustl(File_names%F(25)))
   leng = LEN(trim(adjustl(In_file)))
   Filename = trim(adjustl(Output_path))//trim(adjustl(In_file(1:leng-4)))//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_raidal_distribution(FN, Tim, Target_atoms, Filename, trim(adjustl(In_file)), NumPar, 1.0d1, &
            'Atomic kinetic temperature (K)')  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)

   !----------------
   ! Collect all gnuplot scripts together into one, and execute it:
   call collect_gnuplots(NumPar, trim(adjustl(Output_path)))   ! below
end subroutine Gnuplot_transients




subroutine gnuplot_raidal_distribution(FN, Tim, Target_atoms, Filename, file_data, NumPar, y_min, value_name)
   integer, intent(in) :: FN  ! file with gnuplot script
   real(8), intent(in) :: Tim, y_min ! total simulation time [fs], minimal value of y-axis
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   character(*), intent(in) :: Filename, file_data
   type(Flag), intent(in) :: NumPar
   character(*), intent(in) :: value_name
   !-----------------
   character(10) :: plot_extension, path_sep
   character(20) :: time_step
   integer :: ext_ind, File_num
   integer :: i, j, Nat, shl, col_count, VB_count, leng, N_col
   character(200) :: datafile, ymin, ymax, xmin, xmax
   character(3) :: col
   real(8) :: L_min, L_max, T_min, T_max, dt, x_tics
   character(8) :: temp, time_order
   logical :: x_log
   !-----------------

   ! number of columns (time instants) in the file to plot:
   N_col = size(NumPar%time_grid)

   plot_extension = trim(adjustl(NumPar%plot_extension))
   path_sep = trim(adjustl(NumPar%path_sep))

   ! Get index of file extension:
   call get_extension_index(plot_extension, ext_ind)   ! below

   !L_max = 10.0d0   ! maximal
   L_min = y_min      ! minimal
   !write(ymax,'(i10)') ceiling(L_max)
   write(ymin,'(es12.2)') L_min

   ! Prepare gnuplot script header:
   T_min = 1.0d0
   x_tics = 10.0d0  ! for log scale, assume base 10
   T_max = Tim
   write(xmin,'(f12.5)') T_min
   if (Tim < 5.0) then
      write(xmax,'(i10)') 1000
   elseif (Tim < 50.0) then
      write(xmax,'(i10)') 10000
   else
      write(xmax,'(i10)') 100000
   endif
   x_log = .true.

   ! File with the data:
   datafile = trim(adjustl(file_data))
   leng = LEN(trim(adjustl(datafile)))

   call write_gnuplot_script_header_new(FN, ext_ind, 3.0e0, x_tics, 'Distribution', 'Radius (A)', trim(adjustl(value_name)), &
         trim(adjustl(datafile(1:leng-3)))//trim(adjustl(plot_extension)), path_sep, 0, &
         set_x_log=x_log, set_y_log=.true., fontsize=14) ! below

   ! Prepare the plotting line:
   if (path_sep .EQ. '\') then	! if it is Windows
      ! All time instants:
      do i = 1, N_col    ! all time-steps
         col_count = 1 + i
         write(col,'(i3)') col_count
         write(time_step,'(f12.2)')  NumPar%time_grid(i)

         if (i == 1) then  ! Start:
            write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
               trim(adjustl(ymin))//':] "'// &
               trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "'//&
               trim(adjustl(time_step))//' fs" ,\'
         elseif (i == N_col) then    ! last
            write(FN, '(a)') ' "'//trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
               trim(adjustl(time_step))//' fs" '
         else
            write(FN, '(a)') ' "'//trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
               trim(adjustl(time_step))//' fs" ,\'
         endif
      enddo ! i

   else ! It is linux
      ! All time instants:
      do i = 1, N_col    ! all time-steps
         col_count = 1 + i
         write(col,'(i3)') col_count
         write(time_step,'(f12.2)')  NumPar%time_grid(i)

         if (i == 1) then  ! Start:
            write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
               trim(adjustl(ymin))//':] \"'// &
               trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"'//&
               trim(adjustl(time_step))//' fs\" ,\'
         elseif (i == N_col) then    ! last
            write(FN, '(a)') ' \"'//trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
               trim(adjustl(time_step))//' fs\" '
         else
            write(FN, '(a)') ' \"'//trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
               trim(adjustl(time_step))//' fs\" ,\'
         endif
      enddo ! i
   endif

   ! Prepare the ending:
   call write_gnuplot_script_ending_new(FN, Filename, path_sep)  ! below

end subroutine gnuplot_raidal_distribution





subroutine gnuplot_theta_distribution(FN, Tim, Target_atoms, Filename, file_data, NumPar, holes)
   integer, intent(in) :: FN  ! file with gnuplot script
   real(8), intent(in) :: Tim ! total simulation time [fs]
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   character(*), intent(in) :: Filename, file_data
   type(Flag), intent(in) :: NumPar
   logical, optional :: holes
   !-----------------
   character(10) :: plot_extension, path_sep
   character(20) :: time_step
   integer :: ext_ind, File_num
   integer :: i, j, Nat, shl, col_count, VB_count, leng, N_col
   character(200) :: datafile, ymin, ymax, xmin, xmax
   character(3) :: col
   real(8) :: L_min, L_max, T_min, T_max, dt, x_tics
   character(8) :: temp, time_order
   logical :: x_log, holes_spectrum
   !-----------------

   if (present(holes)) then
      holes_spectrum = holes
   else ! default, electron"
      holes_spectrum = .false.
   endif

   ! number of columns (time instants) in the file to plot:
   N_col = size(NumPar%time_grid)

   plot_extension = trim(adjustl(NumPar%plot_extension))
   path_sep = trim(adjustl(NumPar%path_sep))

   ! Get index of file extension:
   call get_extension_index(plot_extension, ext_ind)   ! below

   !L_max = 10.0d0   ! maximal
   L_min = 1.0d-3      ! minimal
   !write(ymax,'(i10)') ceiling(L_max)
   write(ymin,'(es12.2)') L_min

   ! Prepare gnuplot script header:
   T_min = 0.0
   x_tics = 20.0d0  ! for log scale, assume base 10
   T_max = 180.0d0  ! degrees
   write(xmin,'(f12.5)') T_min
   write(xmax,'(i10)') ceiling(T_max)
   x_log = .false.

   ! File with the data:
   datafile = trim(adjustl(file_data))
   leng = LEN(trim(adjustl(datafile)))

   call write_gnuplot_script_header_new(FN, ext_ind, 3.0e0, x_tics, 'Theta', 'Angle (deg)', 'Spectrum (arb.units)', &
         trim(adjustl(datafile(1:leng-3)))//trim(adjustl(plot_extension)), path_sep, 0, &
         set_x_log=x_log, set_y_log=.true., fontsize=14) ! below

   ! Prepare the plotting line:
   if (path_sep .EQ. '\') then	! if it is Windows
      ! All time instants:
      do i = 1, N_col    ! all time-steps
         col_count = 1 + i
         write(col,'(i3)') col_count
         write(time_step,'(f12.2)')  NumPar%time_grid(i)

         if (i == 1) then  ! Start:
            write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
               trim(adjustl(ymin))//':] "'// &
               trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "'//&
               trim(adjustl(time_step))//' fs" ,\'
         elseif (i == N_col) then    ! last
            write(FN, '(a)') ' "'//trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
               trim(adjustl(time_step))//' fs" '
         else
            write(FN, '(a)') ' "'//trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
               trim(adjustl(time_step))//' fs" ,\'
         endif
      enddo ! i

   else ! It is linux
      ! All time instants:
      do i = 1, N_col    ! all time-steps
         col_count = 1 + i
         write(col,'(i3)') col_count
         write(time_step,'(f12.2)')  NumPar%time_grid(i)

         if (i == 1) then  ! Start:
            write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
               trim(adjustl(ymin))//':] \"'// &
               trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"'//&
               trim(adjustl(time_step))//' fs\" ,\'
         elseif (i == N_col) then    ! last
            write(FN, '(a)') ' \"'//trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
               trim(adjustl(time_step))//' fs\" '
         else
            write(FN, '(a)') ' \"'//trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
               trim(adjustl(time_step))//' fs\" ,\'
         endif
      enddo ! i
   endif

   ! Prepare the ending:
   call write_gnuplot_script_ending_new(FN, Filename, path_sep)  ! below

end subroutine gnuplot_theta_distribution





subroutine gnuplot_spectrum(FN, Tim, Target_atoms, Filename, file_data, NumPar, holes)
   integer, intent(in) :: FN  ! file with gnuplot script
   real(8), intent(in) :: Tim ! total simulation time [fs]
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   character(*), intent(in) :: Filename, file_data
   type(Flag), intent(in) :: NumPar
   logical, optional :: holes
   !-----------------
   character(10) :: plot_extension, path_sep
   character(20) :: time_step
   integer :: ext_ind, File_num
   integer :: i, j, Nat, shl, col_count, VB_count, leng, N_col
   character(200) :: datafile, ymin, ymax, xmin, xmax
   character(3) :: col
   real(8) :: L_min, L_max, T_min, T_max, dt, x_tics
   character(8) :: temp, time_order
   logical :: x_log, holes_spectrum
   !-----------------

   if (present(holes)) then
      holes_spectrum = holes
   else ! default, electron"
      holes_spectrum = .false.
   endif

   ! number of columns (time instants) in the file to plot:
   N_col = size(NumPar%time_grid)

   plot_extension = trim(adjustl(NumPar%plot_extension))
   path_sep = trim(adjustl(NumPar%path_sep))

   ! Get index of file extension:
   call get_extension_index(plot_extension, ext_ind)   ! below

   !L_max = 10.0d0   ! maximal
   L_min = 1.0d-6      ! minimal
   !write(ymax,'(i10)') ceiling(L_max)
   write(ymin,'(es12.2)') L_min

   ! Prepare gnuplot script header:
   if (holes_spectrum) then   ! hole uses linear scale:
      T_min = 0.0
      x_tics = 2.0d0  ! for log scale, assume base 10
      T_max = 20.0
      write(xmin,'(f12.5)') T_min
      write(xmax,'(i10)') ceiling(T_max)
      x_log = .false.
   else ! electron uses log scale:
      T_min = 1.0d0
      x_tics = 10.0d0  ! for log scale, assume base 10
      T_max = Tim
      write(xmin,'(f12.5)') T_min
      write(xmax,'(i10)') 100000
      x_log = .true.
   endif

   ! File with the data:
   datafile = trim(adjustl(file_data))
   leng = LEN(trim(adjustl(datafile)))

   call write_gnuplot_script_header_new(FN, ext_ind, 3.0e0, x_tics, 'Spectrum', 'Energy (eV)', 'Spectrum (arb.units)', &
         trim(adjustl(datafile(1:leng-3)))//trim(adjustl(plot_extension)), path_sep, 0, &
         set_x_log=x_log, set_y_log=.true., fontsize=14) ! below

   ! Prepare the plotting line:
   if (path_sep .EQ. '\') then	! if it is Windows
      ! All time instants:
      do i = 1, N_col    ! all time-steps
         col_count = 1 + i
         write(col,'(i3)') col_count
         write(time_step,'(f12.2)')  NumPar%time_grid(i)

         if (i == 1) then  ! Start:
            write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
               trim(adjustl(ymin))//':] "'// &
               trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "'//&
               trim(adjustl(time_step))//' fs" ,\'
         elseif (i == N_col) then    ! last
            write(FN, '(a)') ' "'//trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
               trim(adjustl(time_step))//' fs" '
         else
            write(FN, '(a)') ' "'//trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
               trim(adjustl(time_step))//' fs" ,\'
         endif
      enddo ! i

   else ! It is linux
      ! All time instants:
      do i = 1, N_col    ! all time-steps
         col_count = 1 + i
         write(col,'(i3)') col_count
         write(time_step,'(f12.2)')  NumPar%time_grid(i)

         if (i == 1) then  ! Start:
            write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
               trim(adjustl(ymin))//':] \"'// &
               trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"'//&
               trim(adjustl(time_step))//' fs\" ,\'
         elseif (i == N_col) then    ! last
            write(FN, '(a)') ' \"'//trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
               trim(adjustl(time_step))//' fs\" '
         else
            write(FN, '(a)') ' \"'//trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
               trim(adjustl(time_step))//' fs\" ,\'
         endif
      enddo ! i
   endif

   ! Prepare the ending:
   call write_gnuplot_script_ending_new(FN, Filename, path_sep)  ! below

end subroutine gnuplot_spectrum




subroutine gnuplot_total_NRG(FN, Tim, Target_atoms, Filename, file_NRG, file_Numbers, NumPar)
   integer, intent(in) :: FN  ! file with gnuplot script
   real(8), intent(in) :: Tim ! total simulation time [fs]
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   character(*), intent(in) :: Filename, file_NRG, file_Numbers
   type(Flag), intent(in) :: NumPar
   !-----------------
   character(10) plot_extension, path_sep
   integer :: ext_ind   ! file extension index
   integer :: i, j, Nat, shl, col_count, VB_count, leng
   character(200) :: datafile, ymin, ymax, xmin, xmax
   character(3) :: col
   real(8) :: L_min, L_max, T_min, T_max, dt, x_tics
   character(8) :: temp, time_order
   logical :: x_log
   !-----------------

   plot_extension = trim(adjustl(NumPar%plot_extension))
   path_sep = trim(adjustl(NumPar%path_sep))

   ! Get index of file extension:
   call get_extension_index(plot_extension, ext_ind)   ! below

   !L_max = 10.0d0   ! maximal
   L_min = 0.0d0      ! minimal
   !write(ymax,'(i10)') ceiling(L_max)
   write(ymin,'(i10)') floor(L_min)

   ! Prepare gnuplot script header:
   if (NumPar%dt_flag <= 0) then ! linear time scale used:
      T_min = 0.0d0
      x_log = .false.
      ! Define the time ticks:
      call order_of_time((T_max - T_min), time_order, temp, x_tics)  ! module "Little_subroutines"
   else ! log-scale
      T_min = 0.01d0
      x_log = .true.
      x_tics = 10.0d0  ! for log scale, assume base 10
   endif
   T_max = Tim
   write(xmin,'(f12.5)') T_min
   write(xmax,'(i10)') ceiling(T_max)

   ! File with the data:
   datafile = trim(adjustl(file_NRG))
   leng = LEN(trim(adjustl(datafile)))

   call write_gnuplot_script_header_new(FN, ext_ind, 3.0e0, x_tics, 'Energies', 'Time (fs)', 'Energy (eV)', &
         trim(adjustl(datafile(1:leng-3)))//trim(adjustl(plot_extension)), path_sep, 2, &
         set_x_log=x_log, set_y_log=.false., fontsize=14) ! below

   Nat = size(Target_atoms)
   col_count = 1  ! to start with
   ! Prepare the plotting line:
   if (path_sep .EQ. '\') then	! if it is Windows

      ! Total:
      write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
                  trim(adjustl(ymin))//':] "'// &
                  trim(adjustl(file_Numbers)) // '"u 1:4 w l lw LW title "Total" ,\'
      ! Electrons:
      write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:2 w l lw LW title "Electrons" ,\'
      ! Atoms:
      write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:3 w l lw LW title "Atoms" ,\'
      ! Photons:
      write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:4 w l lw LW title "Photons" ,\'
      ! Fields
      ! write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:5 w l lw LW title "Fields" ,\' ! SKIP IT< BOT READY
      ! Valence:
      VB_count = 5 + size(Target_atoms(1)%Ip)
      write(col,'(i3)') VB_count ! add valence band
      write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "Valence holes" ,\'
      ! Core holes:
      col_count = 5
      do i = 1, Nat   ! all atoms
         shl = size(Target_atoms(i)%Ip)
         do j = 1, shl
            col_count = col_count + 1  ! column number to print
            write(col,'(i3)') col_count

            if ((i == 1) .and. (j == shl)) then    ! VB
               !VB_count = col_count ! save column number for VB to plot before last line
            elseif ((i == Nat) .and. (j == shl)) then    ! last one
               write(FN, '(a)') ' "'//trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" '

               !write(col,'(i3)') VB_count ! add valence band
               !write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "Valence holes" ,\'
            else
               write(FN, '(a)') ' "'//trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" ,\'
            endif
         enddo ! j
      enddo ! i

   else ! It is linux
      ! Total:
      write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
                  trim(adjustl(ymin))//':] \"'// &
                  trim(adjustl(file_Numbers))//'\"u 1:4 w l lw \"$LW\" title \"Total\" ,\'
      ! Electrons:
      write(FN, '(a)') '\"'//trim(adjustl(datafile))//'\"u 1:2 w l lw \"$LW\" title \"Electrons\" ,\'
      ! Atoms:
      write(FN, '(a)') '\"'//trim(adjustl(datafile))//'\"u 1:3 w l lw \"$LW\" title \"Atoms\" ,\'
      ! Photons:
      write(FN, '(a)') '\"'//trim(adjustl(datafile))//'\"u 1:4 w l lw \"$LW\" title \"Photons\" ,\'
      ! Fields:
      !write(FN, '(a)') '\"'//trim(adjustl(datafile))//'\"u 1:5 w l lw \"$LW\" title \"Fields\" ,\'   ! SKIP IT, NOT READT
      ! Valence:
      VB_count = 5 + size(Target_atoms(1)%Ip)
      write(col,'(i3)') VB_count ! add valence band
      write(FN, '(a)') '\"'//trim(adjustl(datafile))//'\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Valence holes\" ,\'
      ! Core holes:
      col_count = 5
      do i = 1, Nat   ! all atoms
         shl = size(Target_atoms(i)%Ip)
         do j = 1, shl
            col_count = col_count + 1  ! column number to print
            write(col,'(i3)') col_count

            if ((i == 1) .and. (j == shl)) then    ! VB
               !VB_count = col_count ! save column number for VB to plot before last line
            elseif ((i == Nat) .and. (j == shl)) then    ! last one
               write(FN, '(a)') '\"'//trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" '

            else
               write(FN, '(a)') '\"'//trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'
            endif
         enddo ! j
      enddo ! i
   endif

   ! Prepare the ending:
   call write_gnuplot_script_ending_new(FN, Filename, path_sep)  ! below

end subroutine gnuplot_total_NRG



subroutine gnuplot_total_numbers(FN, Tim, Target_atoms, Filename, file_data, NumPar)
   integer, intent(in) :: FN  ! file with gnuplot script
   real(8), intent(in) :: Tim ! total simulation time [fs]
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   character(*), intent(in) :: Filename, file_data
   type(Flag), intent(in) :: NumPar
   !-----------------
   character(10) plot_extension, path_sep
   integer :: ext_ind   ! file extension index
   integer :: i, j, Nat, shl, col_count, VB_count, leng
   character(200) :: datafile, ymin, ymax, xmin, xmax
   character(3) :: col
   real(8) :: L_min, L_max, T_min, T_max, dt, x_tics
   character(8) :: temp, time_order
   logical :: x_log
   !-----------------

   plot_extension = trim(adjustl(NumPar%plot_extension))
   path_sep = trim(adjustl(NumPar%path_sep))

   ! Get index of file extension:
   call get_extension_index(plot_extension, ext_ind)   ! below

   !L_max = 10.0d0   ! maximal
   L_min = 0.0d0      ! minimal
   !write(ymax,'(i10)') ceiling(L_max)
   write(ymin,'(i10)') floor(L_min)

   ! Prepare gnuplot script header:
   if (NumPar%dt_flag <= 0) then ! linear time scale used:
      T_min = 0.0d0
      x_log = .false.
      ! Define the time ticks:
      call order_of_time((T_max - T_min), time_order, temp, x_tics)  ! module "Little_subroutines"
   else ! log-scale
      T_min = 0.01d0
      x_log = .true.
      x_tics = 10.0d0  ! for log scale, assume base 10
   endif
   T_max = Tim
   write(xmin,'(f12.5)') T_min
   write(xmax,'(i10)') ceiling(T_max)

   ! File with the data:
   datafile = trim(adjustl(file_data))
   leng = LEN(trim(adjustl(datafile)))

   call write_gnuplot_script_header_new(FN, ext_ind, 3.0e0, x_tics, 'Numbers', 'Time (fs)', 'Number (arb.units)', &
         trim(adjustl(datafile(1:leng-3)))//trim(adjustl(plot_extension)), path_sep, 1, &
         set_x_log=x_log, set_y_log=.false., fontsize=14) ! below

   Nat = size(Target_atoms)
   col_count = 1  ! to start with
   ! Prepare the plotting line:
   if (path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
         trim(adjustl(ymin))//':] "'// &
         trim(adjustl(datafile)) // '"u 1:2 w l lw LW title "Excited e-" ,\'

      write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:3 w l lw LW title "Emitted e-" ,\'
      write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:6 w l lw LW title "Photons" '

   else ! It is linux
      write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
         trim(adjustl(ymin))//':] \"'// &
         trim(adjustl(datafile)) // '\"u 1:2 w l lw \"$LW\" title \"Excited e-\" ,\'

      write(FN, '(a)') '\"'//trim(adjustl(datafile))//'\"u 1:3 w l lw \"$LW\" title \"Emitted e-\" ,\'
      write(FN, '(a)') '\"'//trim(adjustl(datafile))//'\"u 1:6 w l lw \"$LW\" title \"Photons\" '
   endif

   ! Prepare the ending:
   call write_gnuplot_script_ending_new(FN, Filename, path_sep)  ! below

end subroutine gnuplot_total_numbers




!------------------------------
! Electron, hole, photon MFPs; DOS:
subroutine Gnuplot_electron_hole(NumPar, Target_atoms, File_names, Output_path)   ! From modlue "Gnuplotting_subs"
   type(Flag), intent(in) :: NumPar
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   type(All_names), intent(in) :: File_names   ! all file names for printing out stuff
   character(*), intent(in) :: Output_path
   !----------------
   integer :: FN, ext_ind
   logical :: file_opened
   character(200) :: Filename, File_short
   character(10) :: call_slash, sh_cmd

   ! Find the extension of the gnuplot scripts:
   call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"


   !----------------
   ! 1) Plot electron MFPs:
   Filename = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//'Gnuplot_electron_MFP'//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_electron_MFP(FN, Target_atoms, Filename, File_names%F(2), File_names%F(4), &
         NumPar%plot_extension, trim(adjustl(NumPar%path_sep)))  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it

   !----------------
   ! 2) Plot hole MFPs:
   Filename = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//'Gnuplot_hole_MFP'//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_hole_MFP(FN, Target_atoms, Filename, File_names%F(3), File_names%F(5), &
         NumPar%plot_extension, trim(adjustl(NumPar%path_sep)))  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it

   !----------------
   ! 3) Plot photon MFPs:
   if (NumPar%include_photons) then
      Filename = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//'Gnuplot_photon_MFP'//trim(adjustl(sh_cmd))
      open(newunit=FN, FILE = trim(adjustl(Filename)))
      call gnuplot_photon_MFP(FN, Target_atoms, Filename, File_names%F(7), NumPar%plot_extension, trim(adjustl(NumPar%path_sep)))  ! below
      inquire(unit=FN,opened=file_opened)    ! check if this file is opened
      if (file_opened) close(FN)             ! and if it is, close it
   endif

   !----------------
   ! 4) Plot DOS things:
   Filename = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//'Gnuplot_DOS'//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_DOS(FN, Target_atoms, Filename, File_names%F(9), NumPar%plot_extension, trim(adjustl(NumPar%path_sep)))  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it

   !----------------
   ! Collect all gnuplot scripts together into one, and execute it:
   call collect_gnuplots(NumPar, trim(adjustl(Output_path)))   ! below

end subroutine Gnuplot_electron_hole



subroutine gnuplot_DOS(FN, Target_atoms, Filename, file_DOS, plot_extension, path_sep)
   integer, intent(in) :: FN  ! file with gnuplot script
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   character(*), intent(in) :: Filename, file_DOS
   character(*), intent(in) :: plot_extension   ! file extension
   character(1), intent(in) :: path_sep
   !-----------------
   integer :: ext_ind   ! file extension index
   integer :: i, j, Nat, shl, col_count, VB_count
   character(100) :: datafile, ymin, ymax, xmin, xmax
   character(3) :: col
   real(8) :: L_min, L_max, E_min, E_max

   ! Get index of file extension:
   call get_extension_index(plot_extension, ext_ind)   ! below

   L_max = 10.0d0   ! maximal DOS
   L_min = 0.0d0      ! minimal DOS
   write(ymax,'(i10)') ceiling(L_max)
   write(ymin,'(i10)') floor(L_min)
   E_min = 0.0d0
   E_max = 20.0d0
   write(xmin,'(f12.5)') E_min
   write(xmax,'(i10)') ceiling(E_max)
   ! File with the data:
   datafile = trim(adjustl(file_DOS))


   ! Prepare gnuplot script header:
   call write_gnuplot_script_header_new(FN, ext_ind, 3.0e0, 10.0e0, 'DOS', 'Energy (eV)', 'DOS or mass (arb.units)', &
      'DOS.'//trim(adjustl(plot_extension)), path_sep, 0, set_x_log=.false., set_y_log=.false., fontsize=14) ! below

   Nat = size(Target_atoms)
   col_count = 1  ! to start with
   ! Prepare the plotting line:
   if (path_sep .EQ. '\') then	! if it is Windows
      !write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
         !trim(adjustl(ymin))//':'//trim(adjustl(ymax))//'] "'// &
      write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':]['// &
         trim(adjustl(ymin))//':] "'// &
         trim(adjustl(datafile)) // '"u 1:($3*10) w l lw LW title "DOS" ,\'

      write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:5 w l lw LW title "Effective mass [me]" '

   else ! It is linux
      !write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
         !trim(adjustl(ymin))//':'//trim(adjustl(ymax))//'] \"'// &
      write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':]['// &
         trim(adjustl(ymin))//':] \"'// &
         trim(adjustl(datafile)) // '\"u 1:(\$3*10) w l lw \"$LW\" title \"DOS\" ,\'

      write(FN, '(a)') '\"'//trim(adjustl(datafile))//'\"u 1:5 w l lw \"$LW\" title \"Effective mass [me]\" '
   endif

   ! Prepare the ending:
   call write_gnuplot_script_ending_new(FN, Filename, path_sep)  ! below

end subroutine gnuplot_DOS




subroutine gnuplot_electron_MFP(FN, Target_atoms, Filename, file_IMFP, file_EMFP, plot_extension, path_sep)
   integer, intent(in) :: FN  ! file with gnuplot script
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   character(*), intent(in) :: Filename, file_IMFP, file_EMFP  ! file with data to plot
   character(*), intent(in) :: plot_extension   ! file extension
   character(1), intent(in) :: path_sep
   !-----------------
   integer :: ext_ind   ! file extension index
   integer :: i, j, Nat, shl, col_count, VB_count
   character(100) :: datafile_IMFP, datafile_EMFP, ymin, ymax, xmin, xmax
   character(3) :: col
   real(8) :: L_min, L_max, E_min, E_max

   ! Get index of file extension:
   call get_extension_index(plot_extension, ext_ind)   ! below

   L_max = 10000.0d0   ! maximal MFP
   L_min = 1.0d0      ! minimal MFP
   write(ymax,'(i10)') ceiling(L_max)
   write(ymin,'(i10)') floor(L_min)
   E_min = 0.01d0
   E_max = 1.0d5
   write(xmin,'(f12.5)') E_min
   write(xmax,'(i10)') ceiling(E_max)
   ! File with the data:
   datafile_IMFP = trim(adjustl(file_IMFP))
   datafile_EMFP = trim(adjustl(file_EMFP))

   ! Prepare gnuplot script header:
   call write_gnuplot_script_header_new(FN, ext_ind, 3.0e0, 10.0e0, 'Electron MFP', 'Electron energy (eV)', 'Mean free path (A)', &
      'Electron_MFPs.'//trim(adjustl(plot_extension)), path_sep, 1, set_x_log=.true., set_y_log=.true., fontsize=14) ! below

   Nat = size(Target_atoms)
   col_count = 1  ! to start with
   ! Prepare the plotting line:
   if (path_sep .EQ. '\') then	! if it is Windows
      ! IMFPs:
      do i = 1, Nat   ! all atoms
         shl = size(Target_atoms(i)%Ip)
         do j = 1, shl
            col_count = col_count + 1  ! column number to print
            write(col,'(i3)') col_count

            if ((i == 1) .and. (j == 1)) then    ! first shell
               write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
                  trim(adjustl(ymin))//':'//trim(adjustl(ymax))//'] "'// &
                  trim(adjustl(datafile_IMFP)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                  trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" ,\'
            elseif ((i == 1) .and. (j == shl)) then    ! VB
               VB_count = col_count ! save column number for VB to plot before last line
            else
               write(FN, '(a)') ' "'//trim(adjustl(datafile_IMFP)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" ,\'
            endif
            ! Add valcence and total:
            if ((i == Nat) .and. (j == shl)) then    ! last one
               write(col,'(i3)') VB_count ! add valence band
               write(FN, '(a)') ' "'//trim(adjustl(datafile_IMFP))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "Valence" ,\'

               write(col,'(i3)') col_count+1 ! last column is the total MFP
               write(FN, '(a)') ' "'//trim(adjustl(datafile_IMFP))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "Total IMFP" ,\'
            endif
         enddo ! j
      enddo ! i

      ! EMFPs:
      write(FN, '(a)') ' "'//trim(adjustl(datafile_EMFP))//'"u 1:2 w l lw LW title "Total EMFP" '

   else ! It is linux
      do i = 1, Nat   ! all atoms
         shl = size(Target_atoms(i)%Ip)
         do j = 1, shl
            col_count = col_count + 1  ! column number to print
            write(col,'(i3)') col_count

            if ((i == 1) .and. (j == 1)) then    ! first shell
               write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
                  trim(adjustl(ymin))//':'//trim(adjustl(ymax))//'] \"'// &
                  trim(adjustl(datafile_IMFP)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                  trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'
            elseif ((i == 1) .and. (j == shl)) then    ! VB
               VB_count = col_count ! save column number for VB to plot before last line
            else
               write(FN, '(a)') '\"'//trim(adjustl(datafile_IMFP)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'
            endif
            ! Add valence and total:
            if ((i == Nat) .and. (j == shl)) then    ! last one
               write(col,'(i3)') VB_count ! add valence band
               write(FN, '(a)') '\"'//trim(adjustl(datafile_IMFP))//'\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Valence\" ,\'

               write(col,'(i3)') col_count+1 ! last column is the total MFP
               write(FN, '(a)') '\"'//trim(adjustl(datafile_IMFP))//'\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Total IMFP\" ,\'
            endif
         enddo ! j
      enddo ! i

      ! EMFPs:
      write(FN, '(a)') '\"'//trim(adjustl(datafile_EMFP))//'\"u 1:2 w l lw \"$LW\" title \"Total EMFP\" '
   endif

   ! Prepare the ending:
   call write_gnuplot_script_ending_new(FN, Filename, path_sep)  ! below
end subroutine gnuplot_electron_MFP



subroutine gnuplot_hole_MFP(FN, Target_atoms, Filename, file_IMFP, file_EMFP, plot_extension, path_sep)
   integer, intent(in) :: FN  ! file with gnuplot script
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   character(*), intent(in) :: Filename, file_IMFP, file_EMFP  ! file with data to plot
   character(*), intent(in) :: plot_extension   ! file extension
   character(1), intent(in) :: path_sep
   !-----------------
   integer :: ext_ind   ! file extension index
   integer :: i, j, Nat, shl, col_count, VB_count
   character(100) :: datafile_IMFP, datafile_EMFP, ymin, ymax, xmin, xmax
   character(3) :: col
   real(8) :: L_min, L_max, E_min, E_max

   ! Get index of file extension:
   call get_extension_index(plot_extension, ext_ind)   ! below

   L_max = 10000.0d0   ! maximal MFP
   L_min = 1.0d0      ! minimal MFP
   write(ymax,'(i10)') ceiling(L_max)
   write(ymin,'(i10)') floor(L_min)
   E_min = 0.0d0
   E_max = 30.0d0
   write(xmin,'(f12.5)') E_min
   write(xmax,'(i10)') ceiling(E_max)
   ! File with the data:
   datafile_IMFP = trim(adjustl(file_IMFP))
   datafile_EMFP = trim(adjustl(file_EMFP))

   ! Prepare gnuplot script header:
   call write_gnuplot_script_header_new(FN, ext_ind, 3.0e0, 10.0e0, 'Valence hole MFP', 'Valence hole energy (eV)', 'Mean free path (A)', &
      'Hole_MFPs.'//trim(adjustl(plot_extension)), path_sep, 0, set_x_log=.false., set_y_log=.true., fontsize=14) ! below

   Nat = size(Target_atoms)
   col_count = 1  ! to start with
   ! Prepare the plotting line:
   if (path_sep .EQ. '\') then	! if it is Windows
      ! IMFPs:
      do i = 1, Nat   ! all atoms
         shl = size(Target_atoms(i)%Ip)
         do j = 1, shl
            col_count = col_count + 1  ! column number to print
            write(col,'(i3)') col_count

            if ((i == 1) .and. (j == 1)) then    ! first shell
               write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
                  trim(adjustl(ymin))//':'//trim(adjustl(ymax))//'] "'// &
                  trim(adjustl(datafile_IMFP)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                  trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" ,\'
            elseif ((i == 1) .and. (j == shl)) then    ! VB
               VB_count = col_count ! save column number for VB to plot before last line
            else
               write(FN, '(a)') ' "'//trim(adjustl(datafile_IMFP)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" ,\'
            endif
            ! Add valence and total:
            if ((i == Nat) .and. (j == shl)) then    ! last one
               write(col,'(i3)') VB_count ! add valence band
               write(FN, '(a)') ' "'//trim(adjustl(datafile_IMFP))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "Valence" ,\'

               write(col,'(i3)') col_count+1 ! last column is the total MFP
               write(FN, '(a)') ' "'//trim(adjustl(datafile_IMFP))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "Total IMFP" ,\'
            endif
         enddo ! j
      enddo ! i

      ! EMFPs:
      write(FN, '(a)') ' "'//trim(adjustl(datafile_EMFP))//'"u 1:2 w l lw LW title "Total EMFP" '

   else ! It is linux
      do i = 1, Nat   ! all atoms
         shl = size(Target_atoms(i)%Ip)
         do j = 1, shl
            col_count = col_count + 1  ! column number to print
            write(col,'(i3)') col_count

            if ((i == 1) .and. (j == 1)) then    ! first shell
               write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
                  trim(adjustl(ymin))//':'//trim(adjustl(ymax))//'] \"'// &
                  trim(adjustl(datafile_IMFP)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                  trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'
            elseif ((i == 1) .and. (j == shl)) then    ! VB
               VB_count = col_count ! save column number for VB to plot before last line
            else
               write(FN, '(a)') '\"'//trim(adjustl(datafile_IMFP)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'
            endif
            ! Add valcence and total:
            if ((i == Nat) .and. (j == shl)) then    ! last one
               write(col,'(i3)') VB_count ! add valence band
               write(FN, '(a)') '\"'//trim(adjustl(datafile_IMFP))//'\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Valence\" ,\'

               write(col,'(i3)') col_count+1 ! last column is the total MFP
               write(FN, '(a)') '\"'//trim(adjustl(datafile_IMFP))//'\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Total IMFP\" ,\'
            endif
         enddo ! j
      enddo ! i

      ! EMFPs:
      write(FN, '(a)') '\"'//trim(adjustl(datafile_EMFP))//'\"u 1:2 w l lw \"$LW\" title \"Total EMFP\" '
   endif

   ! Prepare the ending:
   call write_gnuplot_script_ending_new(FN, Filename, path_sep)  ! below
end subroutine gnuplot_hole_MFP




subroutine gnuplot_photon_MFP(FN, Target_atoms, Filename, file_IMFP, plot_extension, path_sep)
   integer, intent(in) :: FN  ! file with gnuplot script
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   character(*), intent(in) :: Filename, file_IMFP ! file with data to plot
   character(*), intent(in) :: plot_extension   ! file extension
   character(1), intent(in) :: path_sep
   !-----------------
   integer :: ext_ind   ! file extension index
   integer :: i, j, Nat, shl, col_count, VB_count
   character(100) :: datafile_IMFP, datafile_EMFP, ymin, ymax, xmin, xmax
   character(3) :: col
   real(8) :: L_min, L_max, E_min, E_max

   ! Get index of file extension:
   call get_extension_index(plot_extension, ext_ind)   ! below

   L_max = 1.0d6   ! maximal MFP
   L_min = 10.0d0      ! minimal MFP
   write(ymax,'(i10)') ceiling(L_max)
   write(ymin,'(i10)') floor(L_min)
   E_min = 1.0d0
   E_max = 1.0d5
   write(xmin,'(f12.5)') E_min
   write(xmax,'(i10)') ceiling(E_max)
   ! File with the data:
   datafile_IMFP = trim(adjustl(file_IMFP))

   ! Prepare gnuplot script header:
   call write_gnuplot_script_header_new(FN, ext_ind, 3.0e0, 10.0e0, 'Photon MFP', 'Photon energy (eV)', 'Mean free path (A)', &
      'Photon_MFPs.'//trim(adjustl(plot_extension)), path_sep, 1, set_x_log=.true., set_y_log=.true., fontsize=14) ! below

   Nat = size(Target_atoms)
   col_count = 1  ! to start with
   ! Prepare the plotting line:
   if (path_sep .EQ. '\') then	! if it is Windows
      ! IMFPs:
      do i = 1, Nat   ! all atoms
         shl = size(Target_atoms(i)%Ip)
         do j = 1, shl
            col_count = col_count + 1  ! column number to print
            write(col,'(i3)') col_count

            if ((i == 1) .and. (j == 1)) then    ! first shell
               write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
                  trim(adjustl(ymin))//':'//trim(adjustl(ymax))//'] "'// &
                  trim(adjustl(datafile_IMFP)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                  trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" ,\'
            elseif ((i == 1) .and. (j == shl)) then    ! VB
               VB_count = col_count ! save column number for VB to plot before last line
            else
               write(FN, '(a)') ' "'//trim(adjustl(datafile_IMFP)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" ,\'
            endif
            ! Add valence and total:
            if ((i == Nat) .and. (j == shl)) then    ! last one
               write(col,'(i3)') VB_count ! add valence band
               write(FN, '(a)') ' "'//trim(adjustl(datafile_IMFP))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "Valence" ,\'

               write(col,'(i3)') col_count+1 ! last column is the total MFP
               write(FN, '(a)') ' "'//trim(adjustl(datafile_IMFP))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "Total" '
            endif
         enddo ! j
      enddo ! i

   else ! It is linux
      do i = 1, Nat   ! all atoms
         shl = size(Target_atoms(i)%Ip)
         do j = 1, shl
            col_count = col_count + 1  ! column number to print
            write(col,'(i3)') col_count

            if ((i == 1) .and. (j == 1)) then    ! first shell
               write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
                  trim(adjustl(ymin))//':'//trim(adjustl(ymax))//'] \"'// &
                  trim(adjustl(datafile_IMFP)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                  trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'
            elseif ((i == 1) .and. (j == shl)) then    ! VB
               VB_count = col_count ! save column number for VB to plot before last line
            else
               write(FN, '(a)') '\"'//trim(adjustl(datafile_IMFP)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'
            endif
            ! Add valence band and total:
            if ((i == Nat) .and. (j == shl)) then    ! last one
               write(col,'(i3)') VB_count ! add valence band
               write(FN, '(a)') '\"'//trim(adjustl(datafile_IMFP))//'\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Valence\" ,\'

               write(col,'(i3)') col_count+1 ! last column is the total MFP
               write(FN, '(a)') '\"'//trim(adjustl(datafile_IMFP))//'\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Total\"'
            endif
         enddo ! j
      enddo ! i

   endif

   ! Prepare the ending:
   call write_gnuplot_script_ending_new(FN, Filename, path_sep)  ! below
end subroutine gnuplot_photon_MFP


!------------------------
! SHI plotting:

subroutine Gnuplot_ion(NumPar, SHI, Target_atoms, File_names, Output_path_SHI)   ! From modlue "Gnuplotting_subs"
   type(Flag), intent(in) :: NumPar
   type(Ion), intent(in) :: SHI   ! declare SHI as an object with atributes "Ion"
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   type(All_names), intent(in) :: File_names   ! all file names for printing out stuff
   character(*), intent(in) :: Output_path_SHI
   !----------------
   integer :: FN, ext_ind
   character(200) :: Filename
   logical :: file_opened
   character(10) :: call_slash, sh_cmd

   ! Find the extension of the gnuplot scripts:
   call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"


   !----------------
   ! 1) Plot ion MFPs:
   Filename = trim(adjustl(Output_path_SHI))//trim(adjustl(NumPar%path_sep))//'Gnuplot_'//trim(adjustl(File_names%F(6)))// &
               '_MFP'//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_SHI_MFP(FN, SHI, Target_atoms, Filename, File_names%F(6), &
         NumPar%plot_extension, trim(adjustl(NumPar%path_sep)))  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it

   !----------------
   ! 2) Plot ion dEdx:
   Filename = trim(adjustl(Output_path_SHI))//trim(adjustl(NumPar%path_sep))//'Gnuplot_'//trim(adjustl(File_names%F(6)))// &
               '_Se'//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_SHI_dEdx(FN, SHI, Target_atoms, Filename, File_names%F(6), &
         NumPar%plot_extension, trim(adjustl(NumPar%path_sep)))  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it

   !----------------
   ! 3) Plot ion range:
   Filename = trim(adjustl(Output_path_SHI))//trim(adjustl(NumPar%path_sep))//'Gnuplot_'//trim(adjustl(File_names%F(6)))// &
               '_Range'//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_SHI_Range(FN, SHI, Target_atoms, Filename, File_names%F(6), &
         NumPar%plot_extension, trim(adjustl(NumPar%path_sep)))  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it


   !----------------
   ! 4) Plot ion effective charge:
   Filename = trim(adjustl(Output_path_SHI))//trim(adjustl(NumPar%path_sep))//'Gnuplot_'//trim(adjustl(File_names%F(6)))// &
               '_Zeff'//trim(adjustl(sh_cmd))
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_SHI_Zeff(FN, SHI, Target_atoms, Filename, File_names%F(6), &
         NumPar%plot_extension, trim(adjustl(NumPar%path_sep)))  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it

   ! 5) Collect all gnuplot scripts together into one, and execute it:
   call collect_gnuplots(NumPar, trim(adjustl(Output_path_SHI)))   ! below

end subroutine Gnuplot_ion



subroutine gnuplot_SHI_Zeff(FN, SHI, Target_atoms, Filename, file_ion_MFP, plot_extension, path_sep)
   integer, intent(in) :: FN  ! file with gnuplot script
   type(Ion), intent(in) :: SHI   ! declare SHI as an object with atributes "Ion"
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   character(*), intent(in) :: Filename, file_ion_MFP  ! file with data to plot
   character(*), intent(in) :: plot_extension   ! file extension
   character(1), intent(in) :: path_sep
   !-----------------
   integer :: ext_ind   ! file extension index
   integer :: i, j, Nat, shl, col_count, VB_count
   character(100) :: datafile, ymax, xmin, xmax
   character(3) :: col
   real(8) :: Se_max, E_min, E_max

   ! Get index of file extension:
   call get_extension_index(plot_extension, ext_ind)   ! below

   Se_max = dble(SHI%Zat)
   write(ymax,'(i10)') ceiling(Se_max)
   E_min = 1.0d0*SHI%Zat/100.0d0
   E_max = 1.0d5*(SHI%Zat/100.0d0)**1.5
   write(xmin,'(f12.5)') E_min
   write(xmax,'(i10)') ceiling(E_max)

   ! File with the data:
   datafile = trim(adjustl(file_ion_MFP))//'_effective_charges.dat'

   ! Prepare gnuplot script header:
   call write_gnuplot_script_header_new(FN, ext_ind, 3.0e0, 10.0e0, 'SHI Se', 'Ion energy (MeV)', 'Effective charge (Z)', &
       trim(adjustl(file_ion_MFP))//'_Zeff.'//trim(adjustl(plot_extension)), path_sep, 1, &
       set_x_log=.true., set_y_log=.false., fontsize=14) ! below

   Nat = size(Target_atoms)
   col_count = 1  ! to start with
   ! Prepare the plotting line:
   if (path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//'][0.0:'//trim(adjustl(ymax))//'] "'// &
         trim(adjustl(datafile))//'"u 1:2 w l lw LW title "Zeff"'
   else ! It is linux
      write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//'][0.0:'//trim(adjustl(ymax))//'] \"'// &
         trim(adjustl(datafile))//'\"u 1:2 w l lw \"$LW\" title \"Zeff\"'
   endif

   ! Prepare the ending:
   call write_gnuplot_script_ending_new(FN, Filename, path_sep)  ! below
end subroutine gnuplot_SHI_Zeff




subroutine gnuplot_SHI_Range(FN, SHI, Target_atoms, Filename, file_ion_MFP, plot_extension, path_sep)
   integer, intent(in) :: FN  ! file with gnuplot script
   type(Ion), intent(in) :: SHI   ! declare SHI as an object with atributes "Ion"
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   character(*), intent(in) :: Filename, file_ion_MFP  ! file with data to plot
   character(*), intent(in) :: plot_extension   ! file extension
   character(1), intent(in) :: path_sep
   !-----------------
   integer :: ext_ind   ! file extension index
   integer :: i, j, Nat, shl, col_count, VB_count
   character(100) :: datafile, ymax, xmax, ymin, xmin
   character(3) :: col
   real(8) :: Se_max, R_max, R_min

   ! Get index of file extension:
   call get_extension_index(plot_extension, ext_ind)   ! below

   Se_max = 10000.0d0*(SHI%Zat/100.0d0)**0.5   ! maximal Se estimate, to plot up to
   R_min = 10.d6*SHI%Zat/100.0d0
   R_max = max(1.0d9*(SHI%Zat/100.0d0)**1.5,1e7)
   write(ymax,'(i10)') ceiling(Se_max)
   write(xmin,'(f16.3)') R_min
   write(xmax,'(f16.3)') R_max

   ! File with the data:
   datafile = trim(adjustl(file_ion_MFP))//'_Range.dat'

   ! Prepare gnuplot script header:
   call write_gnuplot_script_header_new(FN, ext_ind, 3.0e0, 10.0e0, 'SHI Range', 'Residual range (A)', 'Stopping power, Se (eV/A)', &
       trim(adjustl(file_ion_MFP))//'_Range.'//trim(adjustl(plot_extension)), path_sep, 0, &
       set_x_log=.true., set_y_log=.false., fontsize=14) ! below

   Nat = size(Target_atoms)
   col_count = 1  ! to start with
   ! Prepare the plotting line:
   if (path_sep .EQ. '\') then	! if it is Windows
      !write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//'][0.0:] "'// trim(adjustl(datafile)) // &
      write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':][0.0:] "'// trim(adjustl(datafile)) // &
         '"u 3:2 w l lw LW title "Range"'
   else ! It is linux
      !write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//'][0.0:] \"'// trim(adjustl(datafile)) // &
      write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':][0.0:] \"'// trim(adjustl(datafile)) // &
         '\"u 3:2 w l lw \"$LW\" title \"Range\"'
   endif

   ! Prepare the ending:
   call write_gnuplot_script_ending_new(FN, Filename, path_sep)  ! below
end subroutine gnuplot_SHI_Range




subroutine gnuplot_SHI_dEdx(FN, SHI, Target_atoms, Filename, file_ion_MFP, plot_extension, path_sep)
   integer, intent(in) :: FN  ! file with gnuplot script
   type(Ion), intent(in) :: SHI   ! declare SHI as an object with atributes "Ion"
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   character(*), intent(in) :: Filename, file_ion_MFP  ! file with data to plot
   character(*), intent(in) :: plot_extension   ! file extension
   character(1), intent(in) :: path_sep
   !-----------------
   integer :: ext_ind   ! file extension index
   integer :: i, j, Nat, shl, col_count, VB_count
   character(100) :: datafile, ymax, xmin, xmax
   character(3) :: col
   real(8) :: Se_max, E_min, E_max

   ! Get index of file extension:
   call get_extension_index(plot_extension, ext_ind)   ! below

   Se_max = 10000.0d0*(SHI%Zat/100.0d0)**0.5   ! maximal Se estimate, to plot up to
   write(ymax,'(i10)') ceiling(Se_max)
   E_min = 1.0d0*SHI%Zat/100.0d0
   E_max = 1.0d5*(SHI%Zat/100.0d0)**1.5
   write(xmin,'(f12.5)') E_min
   write(xmax,'(i10)') ceiling(E_max)

   ! File with the data:
   datafile = trim(adjustl(file_ion_MFP))//'_dEdx.dat'

   ! Prepare gnuplot script header:
   call write_gnuplot_script_header_new(FN, ext_ind, 3.0e0, 10.0e0, 'SHI Se', 'Ion energy (MeV)', 'Stopping power, Se (eV/A)', &
       trim(adjustl(file_ion_MFP))//'_Se.'//trim(adjustl(plot_extension)), path_sep, 0, &
       set_x_log=.true., set_y_log=.false., fontsize=14) ! below

   Nat = size(Target_atoms)
   col_count = 1  ! to start with
   ! Prepare the plotting line:
   if (path_sep .EQ. '\') then	! if it is Windows
      do i = 1, Nat   ! all atoms
         shl = size(Target_atoms(i)%Ip)
         do j = 1, shl
            col_count = col_count + 1  ! column number to print
            write(col,'(i3)') col_count

            if ((i == 1) .and. (j == 1)) then    ! first shell
               write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//'][0.0:] "'// &
                  trim(adjustl(datafile))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                  trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" ,\'
            elseif ((i == 1) .and. (j == shl)) then    ! VB
               VB_count = col_count ! save column number for VB to plot before last line
            else
               write(FN, '(a)') ' "'//trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" ,\'
            endif
            ! Add valence and total:
            if ((i == Nat) .and. (j == shl)) then    ! last one
               write(col,'(i3)') VB_count ! add valence band
               write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "Valence" ,\'

               write(col,'(i3)') col_count+1 ! last column is the total MFP
               write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "Total"'
            endif
         enddo ! j
      enddo ! i
   else ! It is linux
      do i = 1, Nat   ! all atoms
         shl = size(Target_atoms(i)%Ip)
         do j = 1, shl
            col_count = col_count + 1  ! column number to print
            write(col,'(i3)') col_count

            if ((i == 1) .and. (j == 1)) then    ! first shell
               write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//'][0.0:] \"'// &
                  trim(adjustl(datafile))//'\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                  trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'
            elseif ((i == 1) .and. (j == shl)) then    ! VB
               VB_count = col_count ! save column number for VB to plot before last line
            else
               write(FN, '(a)') '\"'//trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'
            endif
            ! Add valence and total:
            if ((i == Nat) .and. (j == shl)) then    ! last one
               write(col,'(i3)') VB_count ! add valence band
               write(FN, '(a)') '\"'//trim(adjustl(datafile))//'\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Valence\" ,\'

               write(col,'(i3)') col_count+1 ! last column is the total MFP
               write(FN, '(a)') '\"'//trim(adjustl(datafile))//'\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Total\"'
            endif
         enddo ! j
      enddo ! i
   endif

   ! Prepare the ending:
   call write_gnuplot_script_ending_new(FN, Filename, path_sep)  ! below
end subroutine gnuplot_SHI_dEdx




subroutine gnuplot_SHI_MFP(FN, SHI, Target_atoms, Filename, file_ion_MFP, plot_extension, path_sep)
   integer, intent(in) :: FN  ! file with gnuplot script
   type(Ion), intent(in) :: SHI   ! declare SHI as an object with atributes "Ion"
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   character(*), intent(in) :: Filename, file_ion_MFP  ! file with data to plot
   character(*), intent(in) :: plot_extension   ! file extension
   character(1), intent(in) :: path_sep
   !-----------------
   integer :: ext_ind   ! file extension index
   integer :: i, j, Nat, shl, col_count, VB_count
   character(100) :: datafile, ymin, ymax, xmin, xmax
   character(3) :: col
   real(8) :: L_min, L_max, E_min, E_max

   ! Get index of file extension:
   call get_extension_index(plot_extension, ext_ind)   ! below

   L_max = 100.0d0/SHI%Zat*100.0d0   ! maximal Se estimate, to plot up to
   L_min = 0.01d0/SHI%Zat*100.0d0   ! maximal Se estimate, to plot up to
   write(ymax,'(i10)') ceiling(L_max)
   write(ymin,'(f12.5)') L_min
   E_min = 1.0d0*SHI%Zat/100.0d0
   E_max = 1.0d6*(SHI%Zat/100.0d0)**2
   write(xmin,'(f12.5)') E_min
   write(xmax,'(i10)') ceiling(E_max)
   ! File with the data:
   datafile = trim(adjustl(file_ion_MFP))//'_IMFP.dat'

   ! Prepare gnuplot script header:
   call write_gnuplot_script_header_new(FN, ext_ind, 3.0e0, 10.0e0, 'SHI MFP', 'Ion energy (MeV)', 'Ion mean free path (A)', &
       trim(adjustl(file_ion_MFP))//'_MFP.'//trim(adjustl(plot_extension)), path_sep, 1, &
       set_x_log=.true., set_y_log=.true., fontsize=14) ! below

   Nat = size(Target_atoms)
   col_count = 1  ! to start with
   ! Prepare the plotting line:
   if (path_sep .EQ. '\') then	! if it is Windows
      do i = 1, Nat   ! all atoms
         shl = size(Target_atoms(i)%Ip)
         do j = 1, shl
            col_count = col_count + 1  ! column number to print
            write(col,'(i3)') col_count

            if ((i == 1) .and. (j == 1)) then    ! first shell
               write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
                  trim(adjustl(ymin))//':'//trim(adjustl(ymax))//'] "'// &
                  trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                  trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" ,\'
            elseif ((i == 1) .and. (j == shl)) then    ! VB
               VB_count = col_count ! save column number for VB to plot before last line
            else
               write(FN, '(a)') ' "'//trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" ,\'
            endif
            ! Add valence and total:
            if ((i == Nat) .and. (j == shl)) then    ! last one
               write(col,'(i3)') VB_count ! add valence band
               write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "Valence" ,\'

               write(col,'(i3)') col_count+1 ! last column is the total MFP
               write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "Total"'
            endif
         enddo ! j
      enddo ! i
   else ! It is linux
      do i = 1, Nat   ! all atoms
         shl = size(Target_atoms(i)%Ip)
         do j = 1, shl
            col_count = col_count + 1  ! column number to print
            write(col,'(i3)') col_count

            if ((i == 1) .and. (j == 1)) then    ! first shell
               write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
                  trim(adjustl(ymin))//':'//trim(adjustl(ymax))//'] \"'// &
                  trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                  trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'
            elseif ((i == 1) .and. (j == shl)) then    ! VB
               VB_count = col_count ! save column number for VB to plot before last line
            else
               write(FN, '(a)') '\"'//trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'
            endif
            ! Add valence and total:
            if ((i == Nat) .and. (j == shl)) then    ! last one
               write(col,'(i3)') VB_count ! add valence band
               write(FN, '(a)') '\"'//trim(adjustl(datafile))//'\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Valence\" ,\'

               write(col,'(i3)') col_count+1 ! last column is the total MFP
               write(FN, '(a)') '\"'//trim(adjustl(datafile))//'\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Total\"'
            endif
         enddo ! j
      enddo ! i
   endif

   ! Prepare the ending:
   call write_gnuplot_script_ending_new(FN, Filename, path_sep)  ! below
end subroutine gnuplot_SHI_MFP



!=========================================
! General gnuplot subroutines:

subroutine write_gnuplot_script_header_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, path_sep, setkey, set_x_log, set_y_log, fontsize)
   integer, intent(in) :: FN, ind
   real(8), intent(in) :: LW, x_tics
   character(1), intent(in) :: path_sep ! path separator defines which system it is
   character(*), intent(in) :: labl, xlabl, ylabl, Out_file
   integer, intent(in), optional :: setkey, fontsize
   logical, intent(in), optional :: set_x_log, set_y_log
   !---------------------------
   integer :: font_size, set_key
   logical :: x_log, y_log

   if (present(fontsize)) then   ! user-set font size
      font_size = fontsize
   else  ! default font size
      font_size = 14
   endif

   if (present(setkey)) then
      set_key = setkey
   else  ! default
      set_key = 1
   endif

   if (present(set_x_log)) then
      x_log = set_x_log
   else  ! default
      x_log = .false.
   endif

   if (present(set_y_log)) then
      y_log = set_y_log
   else  ! default
      y_log = .false.
   endif

   if (path_sep .EQ. '\') then	! if it is Windows
      call write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, setkey, x_log, y_log, font_size)
   else ! it is linux
      call write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, setkey, x_log, y_log, font_size)
   endif
end subroutine write_gnuplot_script_header_new


subroutine write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, setkey, set_x_log, set_y_log, font_size)
   integer, intent(in) :: FN, ind
   real(8), intent(in) :: LW, x_tics
   integer, intent(in), optional :: setkey, font_size
   logical, intent(in), optional :: set_x_log, set_y_log
   character(*), intent(in) :: labl, xlabl, ylabl, Out_file
   character(20) :: temp, temp2

   if (present(font_size)) then
      write(temp2,'(i0)') font_size
   else
      write(temp2,'(i0)') 14
   endif

   select case (ind)
   case(1:) ! any file format
      write(FN, '(a)') '#!/bin/bash'
      write(FN, '(a)') ''
      write(FN, '(a)') 'NAME='//trim(adjustl(Out_file))
   end select
   write(FN, '(a,f3.1)') 'LW=', LW
   write(FN, '(a)') 'LABL="'//trim(adjustl(labl))//'"'
   write(temp, '(f12.2)') x_tics
   write(FN, '(a)') 'TICSIZ='//trim(adjustl(temp))
   write(FN, '(a)') 'echo " '
   select case (ind)
      case (1)  ! eps
         write(FN, '(a)') 'set terminal postscript enhanced \"Helvetica\" 16 color '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (2)  ! jpeg
         write(FN, '(a)') 'set terminal jpeg font \"arial,'//trim(adjustl(temp2))//'\" '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (3)  ! gif
         write(FN, '(a)') 'set terminal gif font \"arial,'//trim(adjustl(temp2))//'\" '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (4)  ! png
         !write(FN, '(a)') 'set terminal png font \"arial,14\" '
         write(FN, '(a)') 'set terminal pngcairo dashed font \"arial,'//trim(adjustl(temp2))//'\" '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (5)  ! pdf
         write(FN, '(a)') 'set terminal pdf color font \"arial,'//trim(adjustl(temp2))//'\" '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (6)  ! animated gif
         write(FN, '(a)') 'set terminal gif animate delay 10 font \"arial,'//trim(adjustl(temp2))//'\" '
         write(FN, '(a)') 'set output \"$NAME\"'
      case (0)
         write(FN, '(a)') 'set terminal x11 persist'
         write(FN, '(a)') 'unset label'
   endselect
!    write(FN, '(a)') 'set xlabel \"'//trim(adjustl(xlabl))//' \"        font \"Helvetica,20\" '
!    write(FN, '(a)') 'set ylabel \"'//trim(adjustl(ylabl))//' \"      font \"Helvetica,20\" '
   write(FN, '(a)') 'set xlabel \"'//trim(adjustl(xlabl))//'\" font \"arial,18\" '
   write(FN, '(a)') 'set ylabel \"'//trim(adjustl(ylabl))//'\" font \"arial,18\" '

   !write(FN, '(a)') 'set label \"$LABL\" at 150,-8 font \"Helvetica,22\" '
   if (present(setkey)) then
      select case(setkey)
      case (1)
         write(FN, '(a)') 'set key right bottom '
      case (2)
         write(FN, '(a)') 'set key left top '
      case (3)
         write(FN, '(a)') 'set key left bottom '
      case (4)
         write(FN, '(a)') 'unset key '
      case default
         write(FN, '(a)') 'set key right top '
      endselect
   else
      write(FN, '(a)') 'set key right top '
   endif

   if (present(set_x_log)) then
      if (set_x_log) then
         write(FN, '(a)') "set logscale x"
         write(FN, '(a)') 'set format x \"10^{\%L}\"'
      endif
   endif

   if (present(set_y_log)) then
      if (set_y_log) then
         write(FN, '(a)') "set logscale y"
         write(FN, '(a)') 'set format y \"10^{\%L}\"'
      endif
   endif

   write(FN, '(a)') 'set xtics \"$TICSIZ\" '
end subroutine write_gnuplot_script_header_linux_new



subroutine write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, setkey, set_x_log, set_y_log, font_size)
   integer, intent(in) :: FN, ind
   real(8), intent(in) :: LW, x_tics
   integer, intent(in), optional :: setkey, font_size
   logical, intent(in), optional :: set_x_log, set_y_log
   character(*), intent(in) :: labl, xlabl, ylabl, Out_file
   character(20) :: temp, temp2

   if (present(font_size)) then
      write(temp2,'(i0)') font_size
   else
      write(temp2,'(i0)') 14
   endif

   select case (ind)
   case(1:)
      write(FN, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
   end select
   write(FN, '(a,f3.1)') 'LW=', LW

    select case (ind)
      case (1)  ! eps
         write(FN, '(a)') 'set terminal postscript enhanced "Helvetica,'//trim(adjustl(temp2))//'" color '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (2)  ! gpeg
         write(FN, '(a)') 'set terminal jpeg large font "arial,'//trim(adjustl(temp2))//'" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (3)  ! gif
         write(FN, '(a)') 'set terminal gif large font "arial,'//trim(adjustl(temp2))//'" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (4)  ! png
         !write(FN, '(a)') 'set terminal png font "arial,14" '
         write(FN, '(a)') 'set terminal pngcairo dashed font "arial,'//trim(adjustl(temp2))//'" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (5)  ! pdf
         write(FN, '(a)') 'set terminal pdf color font "arial,'//trim(adjustl(temp2))//'" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (6)  ! animated gif
         write(FN, '(a)') 'set terminal gif animate delay 10 font "arial,'//trim(adjustl(temp2))//'" '
         write(FN, '(a)') 'set output "'//trim(adjustl(Out_file))//'"'
      case (0)
         write(FN, '(a)') 'set terminal x11 persist'
         write(FN, '(a)') 'unset label'
   endselect
   write(FN, '(a)') 'set xlabel "'//trim(adjustl(xlabl))//'" font "arial,18"'
   write(FN, '(a)') 'set ylabel "'//trim(adjustl(ylabl))//'" font "arial,18"'

   !write(FN, '(a)') 'set label \"$LABL\" at 150,-8 font \"Helvetica,22\" '
   if (present(setkey)) then
      select case(setkey)
      case (1)
         write(FN, '(a)') 'set key right bottom '
      case (2)
         write(FN, '(a)') 'set key left top '
      case (3)
         write(FN, '(a)') 'set key left bottom '
      case (4)
         write(FN, '(a)') 'unset key '
      case default
         write(FN, '(a)') 'set key right top '
      endselect
   else
      write(FN, '(a)') 'set key right top '
   endif

   if (present(set_x_log)) then
      if (set_x_log) then
         write(FN, '(a)') "set logscale x"
         write(FN, '(a)') 'set format x "10^{%L}"'
      endif
   endif

   if (present(set_y_log)) then
      if (set_y_log) then
         write(FN, '(a)') "set logscale y"
         write(FN, '(a)') 'set format y "10^{%L}"'
      endif
   endif

   !write(FN, '(a,f6.2)') 'set xtics ', x_tics
   write(temp, '(f12.2)') x_tics
   write(FN, '(a,a)') 'set xtics ', trim(adjustl(temp))
end subroutine write_gnuplot_script_header_windows_new


subroutine write_gnuplot_script_ending_new(FN, File_name, path_sep)
   integer, intent(in) :: FN
   character(*), intent(in) :: File_name
   character(1), intent(in) :: path_sep ! path separator defines which system it is
   integer :: iret

   if (path_sep .EQ. '\') then	! if it is Windows
      ! no need to add anything here
   else ! it is linux
      write(FN, '(a)') 'reset'
      write(FN, '(a)') '" | gnuplot '
      !call system('chmod +x '//trim(adjustl(File_name))) ! make the output-script executable
      iret = system('chmod +x '//trim(adjustl(File_name))) ! make the output-script executable
   endif
end subroutine write_gnuplot_script_ending_new



subroutine get_extension_index(text_ext, ind)
   character(*), intent(in) :: text_ext   ! extension
   integer, intent(out) :: ind   ! internal index
   !--------------------
   select case (trim(adjustl(text_ext)))
   case ('eps', 'EPS')  ! eps
      ind = 1
   case ('JPG', 'jpg', 'JPEG', 'jpeg')  ! jpeg
      ind = 2
   case ('GIF', 'gif')  ! gif
      ind = 3
   case ('PNG', 'png')  ! png
      ind = 4
   case ('PDF', 'pdf')  ! pdf
      ind = 5
   case ('animated_gif')  ! animated gif
      ind = 6
   case default ! default
      print*, 'Using default gnuplot format for plots: jpeg'
      ind = 2 ! exclude
   endselect
end subroutine get_extension_index



subroutine order_of_time(tim, text, gnu_text, x_tics)
   real(8), intent(in) :: tim ! time to find its order
   character(*), intent(out) :: text ! fs, ps, ns, mks, ms, s
   character(*), intent(out), optional :: gnu_text ! culomn to set in gnuplot
   real(8), intent(out), optional :: x_tics ! tics for gnuplot
   integer :: time_ord
   time_ord = find_order_of_number(tim) ! below
   if (present(x_tics)) then
      x_tics = 10.0d0**(time_ord) ! set tics for gnuplot
      if (tim/dble(x_tics) > 0.5) then
         x_tics = 10.0d0**(time_ord-1) ! set tics for gnuplot
      else if (tim/dble(x_tics) > 0.2) then
         x_tics = 0.5d0*10.0d0**(time_ord-1) ! set tics for gnuplot
      else
         x_tics = 10.0d0**(time_ord-2) ! set tics for gnuplot
      endif
   endif

   if (time_ord > 1e15) then ! s
      text = '(s)'
      if (present(gnu_text)) gnu_text = '($1/1e15)'
   else if (time_ord > 1e12) then ! ms
      text = '(ms)'
      if (present(gnu_text)) gnu_text = '($1/1e12)'
   else if (time_ord > 1e9) then ! mks
      text = '(mks)'
      if (present(gnu_text)) gnu_text = '($1/1e9)'
   else if (time_ord > 1e6) then ! ns
      text = '(ns)'
      if (present(gnu_text)) gnu_text = '($1/1e6)'
   else if (time_ord > 1e3) then ! ps
      text = '(ps)'
      if (present(gnu_text)) gnu_text = '($1/1e3)'
   else ! fs
      text = '(fs)'
      if (present(gnu_text)) gnu_text = '($1)'
   endif
end subroutine order_of_time



pure function find_order_of_number_real(num)
   integer find_order_of_number_real
   real(8), intent(in) :: num
   character(64) :: temp
   write(temp,'(i8)') CEILING(num) ! make it a string
   find_order_of_number_real = LEN(TRIM(adjustl(temp))) ! find how many characters in this string
end function find_order_of_number_real

pure function find_order_of_number_int(num)
   integer find_order_of_number_int
   integer, intent(in) :: num
   character(64) :: temp
   write(temp,'(i8)') num ! make it a string
   find_order_of_number_int = LEN(TRIM(adjustl(temp))) ! find how many characters in this string
end function find_order_of_number_int


!===================================================
! Gnuplotting all the scripts:
subroutine collect_gnuplots(numpar, out_path)
   type(Flag), intent(in) :: numpar        ! all numerical parameters
   character(*), intent(in) :: out_path    ! folder with the cmd-files
   !------------------------
   character(200) :: File_name, command, Gnuplot_all_file
   integer :: FN, N_f, i, n_slash
   integer :: open_status, iret, idir, leng
   character(200), dimension(:), allocatable :: All_files
   character(300) :: output_path
   character(5) ::  call_slash, sh_cmd

   ! In which directory to collect all gnu scripts:
   output_path = out_path

   ! Create a temporary file:
   Gnuplot_all_file = 'Gnuplot_all'

   ! Find the extension of the gnuplot scripts:
   call cmd_vs_sh(numpar%path_sep, call_slash, sh_cmd)  ! module "Gnuplotting"
   ! Include the extension of the script:
   Gnuplot_all_file = 'Gnuplot_all'//trim(adjustl(sh_cmd))

   ! Include the path to the directory:
   File_name = trim(adjustl(output_path))//numpar%path_sep//trim(adjustl(Gnuplot_all_file))

   ! Save the names of all gnuplot scripts into this file:
   if (numpar%path_sep == '\') then	! if it is Windows
      command = 'dir '//trim(adjustl(output_path))//'\*'//trim(adjustl(sh_cmd))//' /b >'//trim(adjustl(File_name))
   else ! linux:
      command = "ls -t "//trim(adjustl(output_path))//" | grep '"//trim(adjustl(sh_cmd))//"' >"//trim(adjustl(File_name))
   endif

   iret = system(trim(adjustl(command)))   ! execute the command to save file names in the temp file
   !call system(trim(adjustl(command))) ! execute the command to save file names in the temp file

   ! Open the files with gnuplot script names:
   open(NEWUNIT=FN, file=trim(adjustl(File_name)), iostat=open_status)
   if ( open_status /= 0 ) then
      print *, 'Could not open ',trim(adjustl(File_name)),' for gnuplotting.', ' Unit = ', FN
   endif

   ! Find out how many there are:
   call Count_lines_in_file(FN, N_f) ! below

   ! Allocate array with them:
   allocate(All_files(N_f)) ! array with all relevant file names

   ! Read file names:
   do i = 1,N_f
      read(FN,*) All_files(i)
   enddo

   ! Rewind file to overwrite including the calls:
   rewind(FN)
   ! Make the script executable:
   if (numpar%path_sep == '\') then	! if it is Windows
      write(FN,'(a)') '@echo off'
   else
      write(FN,'(a)') '#!/bin/bash'
   endif
   do i = 1,N_f
      if (trim(adjustl(All_files(i))) /= trim(adjustl(Gnuplot_all_file))) then ! exclude the file itself
         if (numpar%path_sep == '\') then	! if it is Windows
            write(FN,'(a)') trim(adjustl(call_slash))//' '//trim(adjustl(All_files(i)))
         else
            leng = LEN(trim(adjustl(All_files(i))))
            if (trim(adjustl(All_files(i)(leng-2:))) == '.sh') then  ! to exclude other files possible containing "sh"
               write(FN,'(a)') trim(adjustl(call_slash))//trim(adjustl(All_files(i)))
            endif
         endif
      endif
   enddo
   close (FN)
   if (numpar%path_sep /= '\') then	! if it is Linux
      iret = system('chmod +x '//trim(adjustl(File_name))) ! make the output-script executable
   endif
   !pause 'Execute all'

   !--------------
   ! Execute all the gnuplot scripts:
   idir = chdir(trim(adjustl(output_path))) ! go into the directory with output files
   !call chdir(trim(adjustl(output_path))) ! go into the directory with output files

   if (numpar%path_sep == '\') then	! if it is Windows
      iret = system( '@echo off' )   ! create the folder
      iret = system(trim(adjustl(call_slash))//' '//trim(adjustl(Gnuplot_all_file)))   ! create the folder
      !call system( '@echo off' )
      !call system(trim(adjustl(call_slash))//' '//trim(adjustl(Gnuplot_all_file)))
   else ! linux:
      iret = system( '#!/bin/bash' )
      iret = system(trim(adjustl(call_slash))//trim(adjustl(Gnuplot_all_file)))
      !call system( '#!/bin/bash' )
      !call system(trim(adjustl(call_slash))//trim(adjustl(Gnuplot_all_file)))
   endif

   ! Count how many times the system has to go out of the directory to get back into the original directory:
   ! Defined by the number of slashes in the path given:
   n_slash = count( (/ (trim(adjustl(output_path(i:i))), i=1,len_trim(output_path)) /) == trim(adjustl(numpar%path_sep)) )
   do i = 1, n_slash+1  ! go up in the directory tree as many times as needed
      idir = chdir("../")    ! exit the directory with output files
      !call chdir("../")    ! exit the directory with output files
   enddo
end subroutine collect_gnuplots



pure subroutine cmd_vs_sh(path_sep, call_slash, sh_cmd)
   character(*), intent(in) :: path_sep
   character(*), intent(out) :: call_slash, sh_cmd
   if (path_sep .EQ. '\') then	! if it is Windows
      call_slash = 'call '
      sh_cmd = '.cmd'
   else ! It is linux
      call_slash = './'
      sh_cmd = '.sh'
   endif
end subroutine cmd_vs_sh




!------------------------------------
! Obsolete subroutines:

subroutine Gnuplot_execute_old(Folder, File_name, n)
   character(*), intent(in) :: Folder, File_name    ! folder, where the script is; script name
   integer, intent(in), optional :: n
   integer i, idir
   !call chdir(trim(adjustl(Folder))) ! go to the directory with the created Gnuplot file
   idir = chdir(trim(adjustl(Folder))) ! go to the directory with the created Gnuplot file
   !call system(trim(adjustl(File_name)))  ! execute there the Gnuplot script
   idir = system(trim(adjustl(File_name)))  ! execute there the Gnuplot script
   if (present(n)) then
      do i = 1, n
         !call chdir("../")    ! go back into the parent directory
         idir = chdir("../")    ! go back into the parent directory
      enddo
   else
      !call chdir("../")    ! go back into the parent directory
      idir = chdir("../")    ! go back into the parent directory
   endif
end subroutine Gnuplot_execute_old


subroutine Gnuplot_ion_old(path_sep, Gnuplot_path, Output_path, Filename, y_column)   ! From modlue "Gnuplotting_subs"
   character(*), intent(in) :: path_sep, Gnuplot_path, Output_path, Filename
   integer, intent(in) :: y_column
   integer FN
   character(200) :: File_name
   logical :: do_short
   do_short = .true.
   File_name = 'Gnuplot_Ion_losses.cmd'
   open(newunit=FN, FILE = trim(adjustl(Output_path))//trim(adjustl(path_sep))//trim(adjustl(File_name)))
   call Gnuplot_header(FN, Gnuplot_path, 'Ion_dEdx', 1, 0, 'Energy, (MeV)', 'dE/dx, (eV/A)', 0)
   call Gnuplot_plotting(FN=FN, File_name=trim(adjustl(Filename))//'_dEdx.dat', x0=1.0d0, x1=1d5, y0=1d0, y1=5d3, x_column=1, x_mult=1d-6, y_column=y_column, lw=3, title='Ion energy losses', go_on=.false.)

   call Gnuplot_header(FN, Gnuplot_path, 'Ion_MFP', 1, 1, 'Energy, (MeV)', 'Mean free path, (A)', 0, do_short)
   call Gnuplot_plotting(FN=FN, File_name=trim(adjustl(Filename))//'_IMFP.dat', x0=1.0d0, x1=1d5, y0=1d-2, y1=1d1, x_column=1, x_mult=1d-6, y_column=y_column, lw=3, title='Ion mean free path', go_on=.false.)

   call Gnuplot_header(FN, Gnuplot_path, 'Ion_range', 1, 0, 'Range, (microns)', 'dE/dx, (eV/A)', 0, do_short)
   call Gnuplot_plotting(FN=FN, File_name=trim(adjustl(Filename))//'_Range.dat', x0=1.0d0, x1=1d5, y0=1d0, y1=5d3, x_column=3, x_mult=1d-4, y_column=2, lw=3, title='Ion range', go_on=.false.)
   close(FN)
   call Gnuplot_execute_old(Output_path, File_name,2) ! execute the created gnuplot script
end subroutine Gnuplot_ion_old


subroutine Gnuplot_electrons_MFP_old(path_sep, Gnuplot_path, Output_path, Filename, y_column)
   character(*), intent(in) :: path_sep, Gnuplot_path, Output_path, Filename
   integer, intent(in) :: y_column
   integer FN
   character(200) :: File_name

   File_name = 'Gnuplot_electron_IMFPs.cmd'
   open(newunit=FN, FILE = trim(adjustl(Output_path))//trim(adjustl(path_sep))//trim(adjustl(File_name)))

   call Gnuplot_header(FN, Gnuplot_path, 'Electron_IMFP_total', 1, 1, 'Energy, (eV)', 'Mean free path, (A)', 1)
   call Gnuplot_plotting(FN=FN, File_name=trim(adjustl(Filename)), x0=1.0d0, x1=10d4, y0=1d0, y1=1d3, x_column=1, y_column=y_column, lw=3, title='Total electron MFP', go_on=.false.)
   close(FN)

   call Gnuplot_execute_old(Output_path, File_name) ! execute the created gnuplot script
end subroutine Gnuplot_electrons_MFP_old



subroutine Gnuplot_header(FN, Gnuplot_path, Fig_name, x_scale, y_scale, x_label, y_label, set_key, G_short)
   integer, intent(in) :: FN    ! file number, which to write into
   character(*), intent(in) :: Gnuplot_path ! e.g. C:\Program Files (x86)\gnuplot\bin
   character(*), intent(in) :: Fig_name     ! name of the figure to be created
   integer, intent(in) :: x_scale, y_scale  ! 0=linear, 1=log
   character(*), intent(in) :: x_label  ! what to write on X axis
   character(*), intent(in) :: y_label  ! what to write on Y axis
   integer, intent(in) :: set_key ! where to show the legend: 0=right-top; 1=right-bottom; 2=left-bottom; 3=left-top; 4=don't
   logical, intent(in), optional :: G_short   ! short or full?
   if (present(G_short)) then
      if (.not. G_short) then
         write(FN, '(a,a,a)') '@echo off & set PATH="', trim(adjustl(Gnuplot_path)), '"'
         write(FN, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
      endif
   else
      write(FN, '(a,a,a)') '@echo off & set PATH="', trim(adjustl(Gnuplot_path)), '"'
      write(FN, '(a,a,a)') '@echo off & call gnuplot.exe -e "echo=', "'#';", 'set macros" "%~f0" & goto :eof'
   endif
   write(FN, '(a)') 'set terminal gif enhanced 16'
   write(FN, '(a,a,a)') "set output '", trim(adjustl(Fig_name)) ,".gif'"
   select case (x_scale)
   case (1)
      write(FN, '(a)') "set logscale x"
   case default
      write(FN, '(a)') "unset logscale x"
   endselect
   select case (y_scale)
   case (1)
      write(FN, '(a)') "set logscale y"
   case default
      write(FN, '(a)') "unset logscale y"
   endselect
   write(FN, '(a,a,a)') 'set xlabel "', trim(adjustl(x_label)), '"'
   write(FN, '(a,a,a)') 'set ylabel "', trim(adjustl(y_label)), '"'
   select case (set_key)
   case (1)
      write(FN, '(a)') 'set key right bottom'
   case (2)
      write(FN, '(a)') 'set key left bottom'
   case (3)
      write(FN, '(a)') 'set key left top'
   case (4)
      write(FN, '(a)') 'unset key'
   case default
      write(FN, '(a)') 'set key right top'
   end select
end subroutine Gnuplot_header


subroutine Gnuplot_plotting(FN, File_name, x0, x1, y0, y1, x_column, x_mult, y_column, y_mult, lw, title, go_on)
   integer, intent(in) :: FN    ! file-number, to which file to write it
   character(*), intent(in) :: File_name    ! from which file to take data to plot
   real(8), intent(in), optional :: x0, x1, y0, y1  ! initial and final x and y points
   integer, intent(in) :: x_column, y_column    ! which columns from the file to use as X and Y 
   real(8), intent(in), optional :: x_mult, y_mult
   integer, intent(in), optional :: lw    ! which line width to use
   character(*), intent(in) :: title    ! which title to give to this curve
   logical, intent(in) :: go_on ! is it the last curve for this plot? true/false

   if (present(x0)) then
      write(FN, '(a,es)', advance='no') 'p[', x0
   else
      write(FN, '(a)', advance='no') 'p[:'
   endif
   if (present(x1)) then
      write(FN, '(a,es,a)', advance='no') ':', x1, ']['
   else
      write(FN, '(a)', advance='no') ']['
   endif
   if (present(y0)) then
      write(FN, '(es,a)', advance='no') y0, ':'
   else
      write(FN, '(a)', advance='no') ':'
   endif
   if (present(y1)) then
      write(FN, '(es,a)', advance='no') y1, '] '
   else
      write(FN, '(a)', advance='no') '] '
   endif
   write(FN, '(a,a,a)', advance='no') '"', trim(adjustl(File_name)) ,'" u '
   
   if (present(x_mult)) then
      write(FN, '(a,i2,a,es,a)', advance='no') '($', x_column, '*', x_mult, ') :'
   else
      write(FN, '(i2,a)', advance='no') x_column, ':'
   endif
   if (present(y_mult)) then
      write(FN, '(a,i2,a,es,a)', advance='no') '($', y_column,'*', y_mult, ') w l '
   else
      write(FN, '(i2,a)', advance='no') y_column, ' w l '
   endif
   
   if (present(lw)) write(FN, '(a,i)', advance='no') 'lw ', lw
   write(FN, '(a,a,a)', advance='no') ' title "', trim(adjustl(title)), '"'
   if (go_on) then
      write(FN, '(a)', advance='no') ','
   else
      write(FN, '(a)') ''
   endif
end subroutine Gnuplot_plotting

end module Gnuplotting_subs

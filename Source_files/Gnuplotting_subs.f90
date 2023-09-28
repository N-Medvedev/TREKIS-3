!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains subroutines to create Gnuplot scripts
! To be useful, it needs installed Gnuplot prior to execution of TREKIS
! Gnuplot is a freeware program distributed here:
! http://www.gnuplot.info/

module Gnuplotting_subs
use Objects, only : Flag, All_names, Atom, Ion
!use ifport  ! library, allowing to operate with directories in intel fortran

implicit none
PRIVATE  ! hides items not listed on public statement

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

public :: Gnuplot_ion, Gnuplot_electrons


contains



subroutine Gnuplot_ion(NumPar, SHI, Target_atoms, File_names, Output_path_SHI)   ! From modlue "Gnuplotting_subs"
   type(Flag), intent(in) :: NumPar
   type(Ion), intent(in) :: SHI   ! declare SHI as an object with atributes "Ion"
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   type(All_names) :: File_names   ! all file names for printing out stuff
   character(*), intent(in) :: Output_path_SHI
   !----------------
   integer :: FN, ext_ind
   character(200) :: Filename
   logical :: file_opened

   !----------------
   ! 1) Plot ion MFPs:
   Filename = trim(adjustl(Output_path_SHI))//trim(adjustl(NumPar%path_sep))//'Gnuplot_'//trim(adjustl(SHI%Name))//'_MFP.cmd'
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_SHI_MFP(FN, SHI, Target_atoms, Filename, File_names%F(6), &
         NumPar%plot_extension, trim(adjustl(NumPar%path_sep)))  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it

   !----------------
   ! 2) Plot ion dEdx:
   Filename = trim(adjustl(Output_path_SHI))//trim(adjustl(NumPar%path_sep))//'Gnuplot_'//trim(adjustl(SHI%Name))//'_Se.cmd'
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_SHI_dEdx(FN, SHI, Target_atoms, Filename, File_names%F(6), &
         NumPar%plot_extension, trim(adjustl(NumPar%path_sep)))  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it

   !----------------
   ! 3) Plot ion range:
   Filename = trim(adjustl(Output_path_SHI))//trim(adjustl(NumPar%path_sep))//'Gnuplot_'//trim(adjustl(SHI%Name))//'_Range.cmd'
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_SHI_Range(FN, SHI, Target_atoms, Filename, File_names%F(6), &
         NumPar%plot_extension, trim(adjustl(NumPar%path_sep)))  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it


   !----------------
   ! 4) Plot ion effective charge:
   Filename = trim(adjustl(Output_path_SHI))//trim(adjustl(NumPar%path_sep))//'Gnuplot_'//trim(adjustl(SHI%Name))//'_Zeff.cmd'
   open(newunit=FN, FILE = trim(adjustl(Filename)))
   call gnuplot_SHI_Zeff(FN, SHI, Target_atoms, Filename, File_names%F(6), &
         NumPar%plot_extension, trim(adjustl(NumPar%path_sep)))  ! below
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it



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
   E_max = 1.0d6*(SHI%Zat/100.0d0)**1.5
   write(xmin,'(f12.5)') E_min
   write(xmax,'(i10)') ceiling(E_max)

   ! File with the data:
   datafile = trim(adjustl(file_ion_MFP))//'_effective_charges.dat'

   ! Prepare grnplot script header:
   call write_gnuplot_script_header_new(FN, ext_ind, 3.0, 10.0, 'SHI Se', 'Ion energy (MeV)', 'Effective charge (Z)', &
      trim(adjustl(SHI%Name))//'_Zeff.'//trim(adjustl(plot_extension)), path_sep, 0, set_x_log=.true., set_y_log=.false., fontsize=14) ! below

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

   Se_max = 20000.0d0*(SHI%Zat/100.0d0)**1.5   ! maximal Se estimate, to plot up to
   R_min = 10.d6*SHI%Zat/100.0d0
   R_max = max(1.0d10*(SHI%Zat/100.0d0)**2,1e7)
   write(ymax,'(i10)') ceiling(Se_max)
   write(xmin,'(f16.3)') R_min
   write(xmax,'(f16.3)') R_max

   ! File with the data:
   datafile = trim(adjustl(file_ion_MFP))//'_Range.dat'

   ! Prepare grnplot script header:
   call write_gnuplot_script_header_new(FN, ext_ind, 3.0, 10.0, 'SHI Range', 'Residual range (A)', 'Stopping power, Se (eV/A)', &
      trim(adjustl(SHI%Name))//'_Range.'//trim(adjustl(plot_extension)), path_sep, 0, set_x_log=.true., set_y_log=.false., fontsize=14) ! below

   Nat = size(Target_atoms)
   col_count = 1  ! to start with
   ! Prepare the plotting line:
   if (path_sep .EQ. '\') then	! if it is Windows
      write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//'][0.0:'//trim(adjustl(ymax))//'] "'// trim(adjustl(datafile)) // &
         '"u 3:2 w l lw LW title "Range"'
   else ! It is linux
      write(FN, '(a)') 'p [1e5:1e8][0.0:'//trim(adjustl(ymax))//'] \"'// trim(adjustl(datafile)) // &
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

   Se_max = 20000.0d0*(SHI%Zat/100.0d0)**1.5   ! maximal Se estimate, to plot up to
   write(ymax,'(i10)') ceiling(Se_max)
   E_min = 1.0d0*SHI%Zat/100.0d0
   E_max = 1.0d6*(SHI%Zat/100.0d0)**1.5
   write(xmin,'(f12.5)') E_min
   write(xmax,'(i10)') ceiling(E_max)

   ! File with the data:
   datafile = trim(adjustl(file_ion_MFP))//'_dEdx.dat'

   ! Prepare grnplot script header:
   call write_gnuplot_script_header_new(FN, ext_ind, 3.0, 10.0, 'SHI Se', 'Ion energy (MeV)', 'Stopping power, Se (eV/A)', &
      trim(adjustl(SHI%Name))//'_Se.'//trim(adjustl(plot_extension)), path_sep, 0, set_x_log=.true., set_y_log=.false., fontsize=14) ! below

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
               write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//'][0.0:'//trim(adjustl(ymax))//'] "'// &
                  trim(adjustl(datafile))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                  trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" ,\'
            elseif ((i == 1) .and. (j == shl)) then    ! VB
               VB_count = col_count ! save column number for VB to plot before last line
            elseif ((i == Nat) .and. (j == shl)) then    ! last one
               write(FN, '(a)') ' "'//trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" ,\'

               write(col,'(i3)') VB_count ! add valence band
               write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "Valence" ,\'

               write(col,'(i3)') col_count+1 ! last column is the total MFP
               write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "Total"'
            else
               write(FN, '(a)') ' "'//trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" ,\'
            endif
         enddo ! j
      enddo ! i
   else ! It is linux
      do i = 1, Nat   ! all atoms
         shl = size(Target_atoms(i)%Ip)
         do j = 1, shl
            col_count = col_count + 1  ! column number to print
            if ((i == 1) .and. (j == 1)) then    ! first shell
               write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//'][0.0:'//trim(adjustl(ymax))//'] \"'// &
                  trim(adjustl(datafile))//'\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                  trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'
            elseif ((i == 1) .and. (j == shl)) then    ! VB
               VB_count = col_count ! save column number for VB to plot before last line
            elseif ((i == Nat) .and. (j == shl)) then    ! last one
               write(FN, '(a)') ' "'//trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'

               write(col,'(i3)') VB_count ! add valence band
               write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Valence\" ,\'

               write(col,'(i3)') col_count+1 ! last column is the total MFP
               write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Total\"'
            else
               write(FN, '(a)') '\"'//trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'
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

   ! Prepare grnplot script header:
   call write_gnuplot_script_header_new(FN, ext_ind, 3.0, 10.0, 'SHI MFP', 'Ion energy (MeV)', 'Ion mean free path (A)', &
      trim(adjustl(SHI%Name))//'_MFP.'//trim(adjustl(plot_extension)), path_sep, 1, set_x_log=.true., set_y_log=.true., fontsize=14) ! below

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
            elseif ((i == Nat) .and. (j == shl)) then    ! last one
               write(FN, '(a)') ' "'//trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" ,\'

               write(col,'(i3)') VB_count ! add valence band
               write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "Valence" ,\'

               write(col,'(i3)') col_count+1 ! last column is the total MFP
               write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:'//trim(adjustl(col))//' w l lw LW title "Total"'
            else
               write(FN, '(a)') ' "'//trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw LW title "' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '" ,\'
            endif
         enddo ! j
      enddo ! i
   else ! It is linux
      do i = 1, Nat   ! all atoms
         shl = size(Target_atoms(i)%Ip)
         do j = 1, shl
            col_count = col_count + 1  ! column number to print
            if ((i == 1) .and. (j == 1)) then    ! first shell
               write(FN, '(a)') 'p ['//trim(adjustl(xmin))//':'//trim(adjustl(xmax))//']['// &
                  trim(adjustl(ymin))//':'//trim(adjustl(ymax))//'] \"'// &
                  trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                  trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'
            elseif ((i == 1) .and. (j == shl)) then    ! VB
               VB_count = col_count ! save column number for VB to plot before last line
            elseif ((i == Nat) .and. (j == shl)) then    ! last one
               write(FN, '(a)') ' "'//trim(adjustl(datafile)) // '"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'

               write(col,'(i3)') VB_count ! add valence band
               write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Valence\" ,\'

               write(col,'(i3)') col_count+1 ! last column is the total MFP
               write(FN, '(a)') ' "'//trim(adjustl(datafile))//'"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"Total\"'
            else
               write(FN, '(a)') '\"'//trim(adjustl(datafile)) // '\"u 1:'//trim(adjustl(col))//' w l lw \"$LW\" title \"' // &
                     trim(adjustl(Target_atoms(i)%Name))//' '//trim(adjustl(Target_atoms(i)%Shell_name(j))) // '\" ,\'
            endif
         enddo ! j
      enddo ! i
   endif

   ! Prepare the ending:
   call write_gnuplot_script_ending_new(FN, Filename, path_sep)  ! below
end subroutine gnuplot_SHI_MFP


subroutine Gnuplot_electrons()
end subroutine Gnuplot_electrons



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
      if (set_x_log) write(FN, '(a)') "set logscale x"
   endif

   if (present(set_y_log)) then
      if (set_y_log) write(FN, '(a)') "set logscale y"
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
      if (set_x_log) write(FN, '(a)') "set logscale x"
   endif

   if (present(set_y_log)) then
      if (set_y_log) write(FN, '(a)') "set logscale y"
   endif

   !write(FN, '(a,f6.2)') 'set xtics ', x_tics
   write(temp, '(f12.2)') x_tics
   write(FN, '(a,a)') 'set xtics ', trim(adjustl(temp))
end subroutine write_gnuplot_script_header_windows_new


subroutine write_gnuplot_script_ending_new(FN, File_name, path_sep)
   integer, intent(in) :: FN
   character(*), intent(in) :: File_name
   character(1), intent(in) :: path_sep ! path separator defines which system it is

   if (path_sep .EQ. '\') then	! if it is Windows
      ! no need to add anything here
   else ! it is linux
      write(FN, '(a)') 'reset'
      write(FN, '(a)') '" | gnuplot '
      call system('chmod +x '//trim(adjustl(File_name))) ! make the output-script executable
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




!------------------------------------
! Old subroutines:

subroutine Gnuplot_execute_old(Folder, File_name, n)
   character(*), intent(in) :: Folder, File_name    ! folder, where the script is; script name
   integer, intent(in), optional :: n
   integer i
   call chdir(trim(adjustl(Folder))) ! go to the directory with the created Gnuplot file
   call system(trim(adjustl(File_name)))  ! execute there the Gnuplot script
   if (present(n)) then
      do i = 1, n
         call chdir("../")    ! go back into the parent directory
      enddo
   else
      call chdir("../")    ! go back into the parent directory
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

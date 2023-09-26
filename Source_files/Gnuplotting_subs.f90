!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains subroutines to create Gnuplot scripts
! To be useful, it needs installed Gnuplot prior to execution of TREKIS
! Gnuplot is a freeware program distributed here:
! http://www.gnuplot.info/

module Gnuplotting_subs
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

public :: Gnuplot_ion, Gnuplot_electrons_MFP


contains


subroutine Gnuplot_ion(path_sep, Gnuplot_path, Output_path, Filename, y_column)   ! From modlue "Gnuplotting_subs"
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
   call Gnuplot_execute(Output_path, File_name,2) ! execute the created gnuplot script
end subroutine Gnuplot_ion


subroutine Gnuplot_electrons_MFP(path_sep, Gnuplot_path, Output_path, Filename, y_column)
   character(*), intent(in) :: path_sep, Gnuplot_path, Output_path, Filename
   integer, intent(in) :: y_column
   integer FN
   character(200) :: File_name
   
   File_name = 'Gnuplot_electron_IMFPs.cmd'
   open(newunit=FN, FILE = trim(adjustl(Output_path))//trim(adjustl(path_sep))//trim(adjustl(File_name)))

   call Gnuplot_header(FN, Gnuplot_path, 'Electron_IMFP_total', 1, 1, 'Energy, (eV)', 'Mean free path, (A)', 1)
   call Gnuplot_plotting(FN=FN, File_name=trim(adjustl(Filename)), x0=1.0d0, x1=10d4, y0=1d0, y1=1d3, x_column=1, y_column=y_column, lw=3, title='Total electron MFP', go_on=.false.)
   close(FN)

   call Gnuplot_execute(Output_path, File_name) ! execute the created gnuplot script
end subroutine Gnuplot_electrons_MFP


subroutine Gnuplot_execute(Folder, File_name, n)
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
end subroutine Gnuplot_execute



subroutine write_gnuplot_script_header_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, path_sep, setkey, fontsize)
   integer, intent(in) :: FN, ind
   real(8), intent(in) :: LW, x_tics
   character(1), intent(in) :: path_sep ! path separator defines which system it is
   character(*), intent(in) :: labl, xlabl, ylabl, Out_file
   integer, intent(in), optional :: setkey, fontsize
   integer :: font_size

   if (present(fontsize)) then   ! user-set font size
      font_size = fontsize
   else  ! default font size
      font_size = 14
   endif

   if (present(setkey)) then
      if (path_sep .EQ. '\') then	! if it is Windows
         call write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, setkey, font_size)
      else ! it is linux
         call write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, setkey, font_size)
      endif
   else
      if (path_sep .EQ. '\') then	! if it is Windows
         call write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, font_size)
      else ! it is linux
         call write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, font_size)
      endif
   endif
end subroutine write_gnuplot_script_header_new


subroutine write_gnuplot_script_header_linux_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, setkey, font_size)
   integer, intent(in) :: FN, ind
   real(8), intent(in) :: LW, x_tics
   integer, intent(in), optional :: setkey, font_size
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
   write(FN, '(a)') 'set xtics \"$TICSIZ\" '
end subroutine write_gnuplot_script_header_linux_new



subroutine write_gnuplot_script_header_windows_new(FN, ind, LW, x_tics, labl, xlabl, ylabl, Out_file, setkey, font_size)
   integer, intent(in) :: FN, ind
   real(8), intent(in) :: LW, x_tics
   integer, intent(in), optional :: setkey, font_size
   character(*), intent(in) :: labl, xlabl, ylabl, Out_file
   character(20) :: temp, temp2


   if (present(font_size)) then
      write(temp2,'(i0)') font_size
   else
      write(temp2,'(i0)') 14
   endif


   select case (ind)
   case(1:)	! eps
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





!------------------------------------
! Old subroutines:
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

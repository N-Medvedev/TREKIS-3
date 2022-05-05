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

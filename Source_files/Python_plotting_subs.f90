!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains subroutines to call Python scripts;
! Requires Python3 to be installed in you system, and the
! TREKIS-3 python package to exist - by default in the folder "TREKIS3_python_plots"

module Python_plotting_subs


use Objects, only : Flag


#ifndef __GFORTRAN__
   USE IFLPORT, only : system, chdir   ! library, allowing to operate with directories in intel fortran
#endif

implicit none
PRIVATE  ! hides items not listed on public statement


public :: Python_plots




contains

!------------------------------
! Plots with the python package TREKIS3_python_plots (alternative to gnuplot):
subroutine Python_plots(NumPar, folder, recursive)
   type(Flag), intent(in) :: NumPar
   character(*), intent(in) :: folder   ! folder with the output files to be plotted
   logical, intent(in) :: recursive     ! also plot all sub-folders?
   !---------------------
   character(1000) :: command
   character(300) :: script
   integer :: iret
   logical :: script_exists

   ! Launcher of the python package:
   if (len_trim(NumPar%python_script) > 0) then
      script = trim(adjustl(NumPar%python_script))
   else ! assume default folder with TREKIS3 python plotting package
      script = 'TREKIS3_python_plots'//trim(adjustl(NumPar%path_sep))//'run.py'
   endif
   inquire(file=trim(script), exist=script_exists)
   if (.not. script_exists) then
      print*, 'Python plotting skipped: cannot find ', trim(adjustl(script))
      print*, "Give its location with 'python_script <path>' in the input file"
      return
   endif

   ! Command line:
   command = trim(adjustl(NumPar%python_cmd))//' "'//trim(script)//'" "'//trim(adjustl(folder))//'"'
   command = trim(command)//' -f '//trim(adjustl(NumPar%python_formats))
   if (.not. recursive) command = trim(command)//' --no-recursive'

   call execute_command_line(trim(command), wait=.true., exitstat=iret)
   if (iret /= 0) then
      print*, 'Python plotting returned error code ', iret, ' for: ', trim(command)
   endif
end subroutine Python_plots





end module Python_plotting_subs

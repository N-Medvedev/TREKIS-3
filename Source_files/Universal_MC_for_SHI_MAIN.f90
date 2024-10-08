!---------------------------------------------------------------
! The code TREKIS-3: https://doi.org/10.5281/zenodo.8394462
! Time Resolved Electron Kinetics in solids after SHI Impact
! was developed by
! N. Medvedev, R.A. Rymzhanov, A.E. Volkov
! With contributions from D. Zainutdinov, F. Akhmetov, P. Babaev, S. Gorbunov
!
! The code simulates a Swift Heavy Ion impact in solids,
! where the properties of any solid target are provided as input files.
!
! The main references describing the model are:
! 1) N.A. Medvedev, R.A. Rymzhanov, A.E. Volkov, J. Phys. D. Appl. Phys. 48 (2015) 355303
! 2) R.A. Rymzhanov, N.A. Medvedev, A.E. Volkov, Nucl. Instrum. Methods B 388 (2016) 41
!
! Should you have any questions, address them to the author:
! nikita.medvedev@fzu.cz
!---------------------------------------------------------------
! TREKIS-3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Although we endeavour to ensure that the code TREKIS-3 and results delivered are correct,
! no warranty is given as to its accuracy. We assume no responsibility for possible errors or omissions.
! We shall not be liable for any damage arising from the use of this code or its parts
! or any results produced with it, or from any action or decision taken
! as a result of using this code or any related material.
!
! This code is distributed as is for non-commercial peaceful purposes only,
! such as research and education. It is explicitly *prohibited* to use the code,
! its parts, its results or any related material for military-related and other than peaceful purposes.
!
! By using this code or its materials, you agree with these terms and conditions.
!---------------------------------------------------------------
! (OPTIONAL) CONVENTIONS OF PROGRAMMING:
! 1) All global variables start with "g_", e.g. g_numpar, and all defined in the module "Variables"
! 2) All modular variable names are defined starting as "m_", e.g. "m_number"
! 3) All local variables used within subrounies should NOT start with "g_" or "m_"
! 4) Add a comment after each subroutine and function specifying in which module it can be found
! 5) Leave comments describing what EACH LINE of the code is doing
! 6) Each end(smth) statement should be commented to which block it belongs, e.g.: if (i<k) then ... endif ! (i<k)
!---------------------------------------------------------------



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Include all the separate file with modules to use in the main program:
#include 'Universal_Constants.f90'
#include 'Objects.f90'
#include 'MPI_subroutines.f90'
#include 'Variables.f90'
#include 'Dealing_with_EADL.f90'
#include 'Gnuplotting_subs.f90'
#include 'Reading_files_and_parameters.f90'
#include 'Cross_sections.f90'
#include 'Analytical_IMFPs.f90'
#include 'Monte_Carlo.f90'
#include 'Thermal_parameters.f90'
#include 'Sorting_output_data.f90'
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

PROGRAM Universal_MC_for_SHI
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
! Initiate only modules used here:
use Universal_Constants
use Objects
use Variables
use Gnuplotting_subs, only: Gnuplot_ion, Gnuplot_electron_hole, Gnuplot_transients
use Reading_files_and_parameters, only: Read_input_file, get_num_shells, Find_VB_numbers, print_time_step, &
                                    get_add_data, set_default_numpar
use Sorting_output_data, only: TREKIS_title, Radius_for_distributions, Allocate_out_arrays, Save_output, &
                            Deallocate_out_arrays, parse_time, print_parameters
use Cross_sections, only: SHI_TotIMFP, Equilibrium_charge_SHI, get_single_pole
use Analytical_IMFPs, only: Analytical_electron_dEdx, Analytical_ion_dEdx, printout_optical_CDF
use Monte_Carlo, only : do_Monte_Carlo
use Thermal_parameters, only : Get_thermal_parameters

use MPI_subroutines, only : initialize_MPI, initialize_random_seed, MPI_barrier_wrapper

#ifdef MPI_USED
use mpi
#endif
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
implicit none

!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
#ifdef _OPENMP
integer OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS ! OMP-related variables must be declared here
#endif
!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
Error_message%Err = .false. ! no errors at the beginning, turn it into 'true' if any occurs

!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
! MPI initialization:
call initialize_MPI(MPI_param, Error_message)   ! module "MPI_subroutines"
if (Error_message%Err) goto 2012  ! if we couldn't get MPI initialized, stop the program
! Print the program name:
if (MPI_param%process_rank == 0) then   ! only MPI master process does it
    call TREKIS_title(6)   ! module "Sorting_output_data"
    ! Get the starting date and time:
    call date_and_time(values=c1) ! standard FORTRAN time and date
    ctim=c1
    write(*, 1005) ctim(5), ctim(6), ctim(7), ctim(3), ctim(2), ctim(1)
endif

!------------------------------------------------------
! Synchronize MPI processes to make sure the title is printed before anything else
call MPI_barrier_wrapper(MPI_param)  ! module "MPI_subroutines"
!------------------------------------------------------

!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
call get_path_separator(Numpar%path_sep, Error_message, read_well, MPI_param)   ! module "Variables"

! One needs to make sure each process is using a different random seed:
call initialize_random_seed(MPI_param)   ! module "MPI_subroutines"


!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! Set default values of flags, proir to reading them from screen or file:
call set_default_numpar(Numpar) ! module "Reading_files_and_parameters"

! Get additional options provided by the user in the command line:
call get_add_data(Numpar, MPI_param) ! module "Reading_files_and_parameters"


!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! Reading input file:
call Read_input_file(Target_atoms, CDF_Phonon, Matter, Mat_DOS, SHI, Tim, dt, Output_path, Output_path_SHI, Material_name, &
        NMC, Num_th, Error_message, read_well, DSF_DEMFP, DSF_DEMFP_H, NumPar, File_names, aidCS, MPI_param)  ! module  'Reading_files_and_parameters'

if (.not. read_well) goto 2012  ! if we couldn't read the input files, there is nothing else to do, go to end
call get_num_shells(Target_atoms, Nshtot) ! from module 'Reading_files_and_parameters'


!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! Set OpenMP parallel threading parameters:
#ifdef _OPENMP
    call OMP_SET_DYNAMIC(0) ! standard openmp subroutine
    call OMP_SET_NUM_THREADS(Num_th)    ! start using threads with openmp: Num_th is the number of threads, defined in the input file
#else ! if you set not to use OpenMP in compiling: 'make OMP=no'
!   print*, 'No openmp to deal with...'
#endif

!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! If single-pole approximation is required, make it:
call get_single_pole(Target_atoms, NumPar, CDF_Phonon, Matter, Error_message, MPI_param)   ! module "Cross_sections"
if (Error_message%Err) goto 2012  ! if we couldn't get the cross section, cannot continue

! Optical CDF for reconstructed CDF from Ritchie-Howie fitted loss function:
if (MPI_param%process_rank == 0) then   ! only MPI master process does it
    call printout_optical_CDF(Output_path, Target_atoms, Matter, NumPar, Mat_DOS)  ! module "Analytical_IMFPs"
endif

!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! Print parameters on screen:
if (MPI_param%process_rank == 0) then   ! only MPI master process does it
    call print_parameters(6, SHI, Material_name, Target_atoms, Matter, NumPar, CDF_Phonon, &
                        Tim, dt, NMC, Num_th, .false., .true., MPI_param)   ! module "Sorting_output_data"
endif

if (SHI%Zat .LE. 0) goto  3012  ! if ion is to be skipped, skip ion:
!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Write down the analytical dEdx of SHI, and IMFP of an electron, valence hole, photon:
if (NumPar%verbose) call print_time_step('Starting SHI mean-free-paths calculations:', MPI_param, msec=.true.)

call Analytical_ion_dEdx(Output_path_SHI, Material_name, Target_atoms, SHI, SHI_MFP, Error_message, read_well, NumPar, Matter, Mat_DOS, File_names, MPI_param)  ! precalculate SHI mean free path and energy loss
if (.not. read_well) goto 2012  ! if we couldn't read the input files, there is nothing else to do, go to end
!if (allocated(File_names%F)) call Gnuplot_ion_old(NumPar%path_sep, File_names%F(1), Output_path_SHI, File_names%F(6), Nshtot+2)   ! From module "Gnuplotting_subs"


call Equilibrium_charge_SHI(SHI, Target_atoms)  ! get Barcas' equilibrium charge from module Cross_sections

! Plot SHI's parameters if requisted:
if (NumPar%do_gnuplot) then
    if (MPI_param%process_rank == 0) then   ! only MPI master process does it
        call Gnuplot_ion(NumPar, SHI, Target_atoms, File_names, Output_path_SHI)   ! module "Gnuplotting_subs"
    endif
endif
3012 continue ! if the ion skipped, go on from here:

! 1) Electron MFPs:
if (NumPar%verbose) call print_time_step('Starting electron mean-free-paths calculations:', MPI_param, msec=.true.)
kind_of_particle = 'Electron'
call Analytical_electron_dEdx(Output_path, Material_name, Target_atoms, CDF_Phonon, Matter, Total_el_MFPs, &
        Elastic_MFP, Error_message, read_well, DSF_DEMFP, Mat_DOS, NumPar, aidCS, kind_of_particle, File_names, MPI_param) ! from module Analytical_IMPS

! 2) Hole MFPs:
if (NumPar%verbose) call print_time_step('Starting VB hole mean-free-paths calculations:', MPI_param, msec=.true.)
kind_of_particle = 'Hole'
call Analytical_electron_dEdx(Output_path, Material_name, Target_atoms, CDF_Phonon, Matter, Total_Hole_MFPs, & 
          Elastic_Hole_MFP, Error_message, read_well, DSF_DEMFP_H, Mat_DOS, NumPar, aidCS, kind_of_particle, File_names, MPI_param) ! from module Analytical_IMPS

! 3) Photon MFPs:
if (NumPar%include_photons) then ! only if we include photons:
    if (NumPar%verbose) call print_time_step('Starting photon mean-free-paths calculations:', MPI_param, msec=.true.)
    kind_of_particle = 'Photon'
    call Analytical_electron_dEdx(Output_path, Material_name, Target_atoms, CDF_Phonon, Matter, Total_Photon_MFPs, &
            Elastic_MFP, Error_message, read_well, DSF_DEMFP, Mat_DOS, NumPar, aidCS, kind_of_particle, File_names, MPI_param) ! from module Analytical_IMPS
else
   allocate(Total_Photon_MFPs(0))
endif

! Plot electron, hole, and photon MFPs, as well as DOS, if requisted:
if (NumPar%do_gnuplot) then
    if (MPI_param%process_rank == 0) then   ! only MPI master process does it
        call Gnuplot_electron_hole(NumPar, Target_atoms, File_names, Output_path)   ! module "Gnuplotting_subs"
    endif
endif


!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! Thermal parameters (electron-phonon coupling, heat capacity, conductivity), if requested (unfinished):
call Get_thermal_parameters(Output_path, CDF_Phonon, Target_atoms, Matter, NumPar, Mat_DOS, File_names, MPI_param) ! module "Thermal_parameters"
!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

if ((NMC .LE. 0) .OR. (Tim .LE. 0.0d0)) then
    if (MPI_param%process_rank == 0) then   ! only MPI master process does it
        write(*,'(a)') trim(adjustl(dashline))
        write(*,'(a)') 'No Monte Carlo routine will be performed since'
        write(*,'(a,i6)') 'Number of MC iterations = ', NMC
        write(*,'(a, ES16.7)') 'Time of MC analysis = ', Tim
    endif
    goto 2012 ! skip MC at all
endif

! if we couldn't read the input files, there is nothing else to do, go to the end; or if skip ion:
if ((.not. read_well) .OR. (SHI%Zat .LE. 0)) goto 2012

! Prepare differential SHI MFP for the given energy:
if (NumPar%verbose) call print_time_step('Calculating SHI mean free paths for given energy:', MPI_param, msec=.true.)
do i = 1, size(Target_atoms)
    do j = 1, size(SHI_MFP(i)%ELMFP)
        call SHI_TotIMFP(SHI, Target_atoms, i, j, Sigma, dEdx, Matter, Mat_DOS, NumPar, diff_SHI_MFP)   ! module "Analytical_IMFPs"
        SHI_MFP(i)%ELMFP(j)%E(:) = SHI_MFP(1)%ELMFP(1)%E(:) ! [eV] energy
        Total_el_MFPs(i)%ELMFP(j)%E(:) = Total_el_MFPs(1)%ELMFP(1)%E(:) ! [eV] energy
        Total_Hole_MFPs(i)%ELMFP(j)%E(:) = Total_Hole_MFPs(1)%ELMFP(1)%E(:) ! [eV] energy
    enddo
enddo


!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
if (NumPar%verbose) call print_time_step('Preparing output directory and files:', MPI_param, msec=.true.)

! Find the number of atom and shell which correspond to the valence band:
call Find_VB_numbers(Target_atoms, Lowest_Ip_At, Lowest_Ip_Shl)    ! module Reading_files_and_parameters
! Fill the array of cylindrical radii:
call Radius_for_distributions(Target_atoms, Out_Distr, Out_V, Tim, dt, NumPar%dt_flag, Matter%Layer)      ! module Reading_files_and_parameters
!BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB


!------------------------------------------------------
! Synchronize MPI to make sure all the parameters are known to all processes
! (may not be necessary, just a precaution)
call MPI_barrier_wrapper(MPI_param)  ! module "MPI_subroutines"
!------------------------------------------------------


!MCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMC
! Allocate arrays for MC iterations:
call Allocate_out_arrays(target_atoms, Out_Distr, Mat_DOS, Out_tot_Ne, Out_tot_Nphot, Out_tot_E, Out_E_e, Out_E_phot, Out_E_at, Out_E_h, &
      Out_R, Out_V, Out_ne, Out_Ee, Out_nphot, Out_Ephot, Out_Ee_vs_E, Out_Eh_vs_E, Out_Elat, Out_nh, Out_Eh, Out_Ehkin, Out_Eat_dens, &
      Out_theta, Out_theta_h, Out_theta1, Out_field_all, Out_Ne_Em, Out_E_Em, Out_Ee_vs_E_Em, Out_E_field, Out_diff_coeff) ! Module 'Sorting_output_data.f90'

if (MPI_param%process_rank == 0) then   ! only MPI master process does it
    write(*,'(a)') ' '
    write(*,'(a)') trim(adjustl(dashline))
    if (NumPar%verbose) call print_time_step('Starting MC iterations:', MPI_param, msec=.true.)
endif

! The subroutine for MC is stored in the module Monte_Carlo_modelling:
call do_Monte_Carlo(NMC, SHI, SHI_MFP, diff_SHI_MFP, Target_atoms, Lowest_Ip_At, Lowest_Ip_Shl, CDF_Phonon, &
     Total_el_MFPs, Elastic_MFP, Total_Hole_MFPs, Elastic_Hole_MFP, Total_Photon_MFPs, Mat_DOS, Tim, dt, Matter, NumPar, aidCS, &
     Out_R, Out_V, Out_ne, Out_Ee, Out_nphot, Out_Ephot, Out_Ee_vs_E, Out_Eh_vs_E, Out_Elat, &
     Out_nh, Out_Eh, Out_Ehkin, Out_tot_Ne, Out_tot_Nphot, Out_tot_E, &
     Out_E_e, Out_E_phot, Out_E_at, Out_E_h, Out_Eat_dens, Out_theta, Out_theta_h, Out_theta1, Out_Ne_Em, Out_E_Em, Out_Ee_vs_E_Em, &
     Error_message, DSF_DEMFP, DSF_DEMFP_H, Out_field_all, Out_E_field, Out_diff_coeff, MPI_param) ! module "Monte_Carlo"

if (NumPar%verbose) call print_time_step('Preparing MC output data:', MPI_param, msec=.true.)

! Average the distributions over the statistics:
if (MPI_param%process_rank == 0) then   ! only MPI master process does it
    Out_tot_Ne = Out_tot_Ne/dble(NMC)
    Out_tot_Nphot = Out_tot_Nphot/dble(NMC)
    Out_tot_E = Out_tot_E/dble(NMC)
    Out_E_e = Out_E_e/dble(NMC)
    Out_E_phot = Out_E_phot/dble(NMC)
    Out_E_at = Out_E_at/dble(NMC)
    Out_E_h = Out_E_h/dble(NMC)
    Out_Eat_dens = Out_Eat_dens/dble(NMC)
    Out_Distr%ne = Out_ne/dble(NMC)
    Out_Distr%Ee = Out_Ee/dble(NMC)
    Out_nphot = Out_nphot/dble(NMC)
    Out_Ephot = Out_Ephot/dble(NMC)
    Out_Ee_vs_E = Out_Ee_vs_E/dble(NMC) ! electron spectrum
    Out_Eh_vs_E = Out_Eh_vs_E/dble(NMC) ! VB holes spectrum
    Out_Elat = Out_Elat/dble(NMC)
    Out_theta = Out_theta/dble(NMC)
    Out_theta_h = Out_theta_h/dble(NMC)
    Out_field_all = Out_field_all/dble(NMC)
    Out_E_field = Out_E_field/dble(NMC)
    Out_Ne_Em = Out_Ne_Em/dble(NMC)
    Out_E_Em = Out_E_Em/dble(NMC)
    Out_Ee_vs_E_Em = Out_Ee_vs_E_Em/dble(NMC)
    Out_diff_coeff = Out_diff_coeff/dble(NMC)
    do k = 1, size(Out_Distr%Atom(1)%Shl(1)%nh,1) ! for all the timesteps:
        do i = 1, size(Out_Distr%Atom)    ! for all atoms
            do j = 1, size(Out_Distr%Atom(i)%Shl)    ! for all shells of atom
                Out_Distr%Atom(i)%Shl(j)%nh(k,:) = Out_nh(k,:,i,j)/dble(NMC)
                Out_Distr%Atom(i)%Shl(j)%Eh(k,:) = Out_Eh(k,:,i,j)/dble(NMC)
                Out_Distr%Atom(i)%Shl(j)%Ehkin(k,:) = Out_Ehkin(k,:,i,j)/dble(NMC)
            enddo ! j
        enddo   ! i
    enddo   ! k
endif ! (MPI_param%process_rank == 0)
!OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
call Save_output(Output_path_SHI, File_names, ctim, NMC, Num_th, Tim, dt, Material_name, Matter, Target_atoms, Mat_DOS, CDF_Phonon, &
                SHI, Out_R, Out_tot_Ne, Out_tot_Nphot, Out_tot_E, Out_E_e, Out_E_phot, Out_nphot, Out_Ephot,&
                Out_Ee_vs_E, Out_Eh_vs_E, Out_E_at, Out_E_h, Out_Eat_dens, Out_Distr, &
                Out_Elat, Out_theta, Out_theta_h, Out_field_all, Out_Ne_Em, Out_E_Em, Out_Ee_vs_E_Em, &
                NumPar, Out_E_field, Out_diff_coeff, MPI_param) ! module 'Sorting_output_data.f90'

! clean up after using temporary arrays:
call Deallocate_out_arrays(Out_tot_Ne, Out_tot_Nphot, Out_tot_E, Out_E_e, Out_E_phot, Out_E_at, Out_E_h, Out_R, Out_V, &
    Out_ne, Out_Ee, Out_nphot, Out_Ephot, Out_Ee_vs_E, Out_Eh_vs_E, &
    Out_Elat, Out_nh, Out_Eh, Out_Ehkin, Out_Eat_dens, Out_theta, Out_theta_h, Out_field_all, Out_Ne_Em, Out_E_Em, &
    Out_Ee_vs_E_Em, Out_E_field)       ! module 'Sorting_output_data.f90'

! Gnuplot the data, if requested:
if (MPI_param%process_rank == 0) then   ! only MPI master process does it
    call Gnuplot_transients(Tim, NumPar, Matter, Target_atoms, File_names)   ! module "Gnuplotting_subs"
endif
!MCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of the program:
! Printing out the duration of the program, starting and ending time and date:
2012 continue

if (MPI_param%process_rank == 0) then   ! only MPI master process does it
    call date_and_time(values=c1)	    ! For calculation of the time of execution of the program
    as1=dble(24*60*60*(c1(3)-ctim(3))+3600*(c1(5)-ctim(5))+60*(c1(6)-ctim(6))+(c1(7)-ctim(7))+(c1(8)-ctim(8))*0.001) ! sec

    call parse_time(as1, Text_var) ! module "Sorting_output_data.f90"
    print*, '   '
    write(*,'(a, a)') 'Duration of execution of the program: ', trim(adjustl(Text_var))
    print*, '   '
    write(*, 1001) ctim(5),ctim(6),ctim(7), ctim(3), ctim(2), ctim(1)
    write(*, 1002) c1(5), c1(6), c1(7), c1(3),c1(2), c1(1)
    print*, '   '
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Closing the remaining opened files:
#ifdef MPI_USED
   ! https://www.open-mpi.org/doc/v3.0/man3/MPI_File_close.3.php
   !call MPI_FILE_CLOSE(Error_message%File_Num, MPI_param%ierror)    ! module "MPI"
   !if (MPI_param%ierror /= MPI_SUCCESS) then
   !   write(Error_message%Err_descript, '(A,I0,A)') '[MPI process #', MPI_param%process_rank, '] Failure in closing Error-log file'
   !   print*, trim(adjustl(Error_message%Err_descript)) ! print it also on the sreen
   !   call MPI_Abort(MPI_COMM_WORLD, -1, MPI_param%ierror)   ! module "MPI"
   !endif

   inquire(Error_message%File_Num, opened=file_opened)
   if (file_opened) then
      if (Error_message%Err) then ! if error occured (thus it's = "true") then save the error log file:
         close(Error_message%File_Num)
      else ! if there was no error, no need to keep the file, delete it:
         close(Error_message%File_Num, status='delete')
      endif
   endif

#else
   if (Error_message%Err) then ! if error occured (thus it's = "true") then save the error log file:
      close(Error_message%File_Num)
   else ! if there was no error, no need to keep the file, delete it:
      close(Error_message%File_Num, status='delete')
   endif
#endif


!-------------------------------------
#ifdef MPI_USED
	!write(*,'(a,i0,a)') '[MPI process #', MPI_param%process_rank, '] test 3'
	if (MPI_param%process_rank == 0) then   ! only MPI master process does it
        print*, 'Finilizing MPI'
	endif
	call MPI_FINALIZE(MPI_param%ierror)
	if (MPI_param%ierror /= 0) then
		write(*, *) 'Error finalizing MPI!'
	endif
#endif
!-------------------------------------

if (MPI_param%process_rank == 0) then   ! only MPI master process does it
    print*, 'TREKIS has completed its calculations...'
endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Formats defined for printing out on the screen:
1001 format ('Beginning: ', i2.2, ':', i2.2, ':', i2.2, '  ', i2.2, '/', i2.2, '/', i4.4)
1002 format ('The end :  ', i2.2, ':', i2.2, ':', i2.2, '  ', i2.2, '/', i2.2, '/', i4.4)
1005 format ('Start at : ', i2.2, ':', i2.2, ':', i2.2, '  ', i2.2, '/', i2.2, '/', i4.4)
1006 format ('Step at: ', i2.2, ':', i2.2, ':', i2.2, '  ', i2.2, '/', i2.2, '/', i4.4)
1007 format ('Step at: ', i2.2, ':', i2.2, ':', i2.2, ':', i3.3, '  ', i2.2, '/', i2.2, '/', i4.4)
1008 format (a, i4, a, i6, a, i2.2, ':', i2.2, ':', i2.2)

STOP

contains
! currently nothing

END PROGRAM Universal_MC_for_SHI

!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains subroutines for reading and interpreting input files:

module Reading_files_and_parameters
  use Universal_Constants   ! let it use universal constants
  use Objects   ! since it uses derived types, it must know about them from module 'Objects'
  use Dealing_with_EADL, only : Decompose_compound, check_atomic_parameters, Find_element_name, define_PQN, &
                                 Count_lines_in_file, m_atomic_folder, m_atomic_data_file, SkipCount_lines_in_file
  use Variables, only: dashline, starline
  
  use MPI_subroutines, only : MPI_barrier_wrapper, MPI_fileopen_wrapper, MPI_fileclose_wrapper, Save_error_details, MPI_error_wrapper


  ! Open_MP related modules from external libraries:
  !USE IFLPORT, only : system
#ifdef _OPENMP
   USE OMP_LIB, only : omp_get_max_threads
#endif

#ifdef MPI_USED
   use mpi !, only : MPI_FILE_OPEN, MPI_COMM_WORLD, MPI_MODE_RDWR, MPI_Abort, MPI_SUCCESS
#endif

  implicit none
private  ! hides items not listed on public statement


! this interface finds by itself which of the two subroutine to use depending on the dimensions of the array passed:
interface Find_in_array ! search cheaking one by one
    module procedure Find_in_1D_array
    module procedure Find_in_2D_array
end interface Find_in_array
! this interface finds by itself which of the two subroutine to use depending on the dimensions of the array passed:
interface Find_in_array_monoton ! search with bisection method
    module procedure Find_in_monotonous_1D_array
    module procedure Find_in_monotonous_2D_array
end interface Find_in_array_monoton

interface Linear_approx
    module procedure Linear_approx_2d
    module procedure Linear_approx_2x1d
end interface Linear_approx

interface Integrate_function ! integrate function
    module procedure Integrate_function_one
    module procedure Integrate_function_save
end interface Integrate_function

interface Trapeziod
    module procedure Trapeziod_one
    module procedure Trapeziod_save
end interface Trapeziod

public :: Find_in_array, Find_in_array_monoton, Linear_approx, get_file_stat, get_num_shells, print_time_step
public :: Read_input_file, Linear_approx_2x1d_DSF, Find_VB_numbers, read_file_here, read_SHI_MFP, get_add_data, m_INPUT_file, &
          Find_in_monoton_array_decreasing, set_default_numpar, &
          broadcast_SHI_MFP, broadcast_el_MFPs, broadcast_el_aidCS_electrons, broadcast_el_aidCS_holes, broadcast_Elastic_MFP


character(25), parameter :: m_form_factors_file = 'Atomic_form_factors.dat'
character(25), parameter :: m_INPUT_file = 'INPUT_PARAMETERS.txt'
character(10), parameter :: m_INPUT_CDF = 'INPUT_CDF'
character(10), parameter :: m_INPUT_DOS = 'INPUT_DOS'
character(10), parameter :: m_INPUT_DSF = 'INPUT_DSF'
!character(10), parameter :: m_INPUT_EADL = 'INPUT_EADL'

contains



subroutine get_file_stat(File_name, device_ID, Inode_number, File_mode, Number_of_links, O_uid, O_gid, where_located, &
                         File_size, Last_access_time, Last_modification_time, Last_status_change, blocks_allocated)
! See description here: https://software.intel.com/en-us/node/526830
   character(*), intent(in) :: File_name ! which file we are checking?
   integer, intent(out), optional :: device_ID ! Device the file resides on
   integer, intent(out), optional :: Inode_number ! File inode number
   integer, intent(out), optional :: File_mode ! Access mode of the file
   integer, intent(out), optional :: Number_of_links ! Number of hard links to the file
   integer, intent(out), optional :: O_uid ! User ID of owner
   integer, intent(out), optional :: O_gid ! Group ID of owner
   integer, intent(out), optional :: where_located ! Raw device the file resides on
   integer, intent(out), optional :: File_size ! Size of the file
   integer, intent(out), optional :: Last_access_time ! Time when the file was last accessed (*)
   integer, intent(out), optional :: Last_modification_time ! Time when the file was last modified(*)
   integer, intent(out), optional :: Last_status_change ! Time of last file status change (*)
   integer, intent(out), optional :: blocks_allocated ! Blocksize for file system I/O operations
   !(*) Times are in the same format returned by the TIME function (number of seconds since 00:00:00 Greenwich mean time, January 1, 1970).
   !=====================
   ! The preprocessor option defining compilation with Gfortran: https://gcc.gnu.org/onlinedocs/gfortran/Preprocessing-Options.html
#ifdef __GFORTRAN__
   ! for gfortran compiler:
   INTEGER :: info_array(13)
#else
   ! for intel fortran compiler:
   INTEGER :: info_array(12)
#endif

   ! Get the statistics on the file:
   call STAT(trim(adjustl(File_name)), info_array) ! intrinsec fortran subroutine

   if (present(device_ID)) device_ID = info_array(1)  ! Device the file resides on
   if (present(Inode_number)) Inode_number = info_array(2) ! File inode number
   if (present(File_mode)) File_mode = info_array(3) ! Access mode of the file
   if (present(Number_of_links)) Number_of_links = info_array(4) ! Number of hard links to the file
   if (present(O_uid)) O_uid = info_array(5) ! User ID of owner
   if (present(O_gid)) O_gid = info_array(6) ! Group ID of owner
   if (present(where_located)) where_located = info_array(7) ! Raw device the file resides on
   if (present(File_size)) File_size = info_array(8) ! Size of the file
   if (present(Last_access_time)) Last_access_time = info_array(9) ! Time when the file was last accessed (*)
   if (present(Last_modification_time)) Last_modification_time = info_array(10) ! Time when the file was last modified(*)
   if (present(Last_status_change)) Last_status_change = info_array(11) ! Time of last file status change (*)
   if (present(blocks_allocated)) blocks_allocated = info_array(12) ! Blocksize for file system I/O operations
end subroutine get_file_stat



subroutine set_default_numpar(Numpar)
   type(Flag), intent(inout) :: NumPar
   !------------------------
   ! Default to start with:
   NumPar%redo_IMFP = .false. ! don't recalculate inelastic MFPs, if possible
   NumPar%redo_EMFP = .false. ! don't recalculate elastic MFPs, if possible
   NumPar%redo_IMFP_SHI = .false. ! don't recalculate elastic MFPs, if possible
   NumPar%include_photons = .false. ! no photons by default (unless user includes them)
   NumPar%plasmon_Emax = .false. ! do not include plasmon integration limit in inelastic CDF
   NumPar%field_include = .false.   ! no fields (bc NOT READY!)
   NumPar%print_CDF = .false. ! don't print CDF file out
   NumPar%print_CDF_optical = .false.  ! don't print optical CDF
   NumPar%do_gnuplot = .true. ! gnuplot by default
   NumPar%plot_extension = 'jpeg' ! default jpeg-files
   NumPar%get_thermal = .false.   ! default: no thermal parameters
   NumPar%CS_method = 1    ! choice of the method of CS integration (integration grid): default - tabulated files
   NumPar%DOS_file = ''    ! no optional name, to use default
   NumPar%CDF_file = ''    ! no optional name, to use default
   NumPar%kind_of_EMFP = 1 ! kind of inelastic mean free path (0=atomic; 1=CDF, 2=DSF)
   NumPar%kind_of_CDF = 1   ! kind of CDF used for inelastic CS: 0=Ritchie; 1=Single=pole
   NumPar%kind_of_CDF_ph = 1  ! kind of CDF used for elastic CS: 0=Ritchie; 1=Single=pole
   NumPar%kind_of_DR = 0    ! target electron dispersion relation used in CDF calculations
   NumPar%dt_flag = 1      ! kind of time-grid (0=linear;1=log)
   NumPar%CDF_elast_Zeff = 0 ! kind of effective charge of target atoms (1=1, 0=Barkas-like Zeff)
   NumPar%out_dim = 1         ! dimensionality of the output plots: 0 = eV/A^3 (old); 1=eV/atom
   ! Printout for testing:
   NumPar%verbose = .false.
   NumPar%very_verbose = .false.
   ! Flags for marking parts of user-defined CDF:
   NumPar%VB_CDF_defined = .false.
   NumPar%phonon_CDF_defined = .false.
   ! Parameters for MD (unfinished)
   NumPar%MD_dt = 1.0d0
   NumPar%MD_grid = 0
   NumPar%Num_Z_points = 0
   NumPar%Zout_min = 0.0d0
   NumPar%Zout_max = 0.0d0
   NumPar%Zout_dz = 0.0d0
   NumPar%field_dt = 0.0d0    ! time-grid for fields update
end subroutine set_default_numpar



subroutine Read_input_file(Target_atoms, CDF_Phonon, Matter, Mat_DOS, SHI, Tim, dt, Output_path, Output_path_SHI, &
           Material_name, NMC, Num_th, Error_message, read_well, DSF_DEMFP, DSF_DEMFP_H, NumPar, File_names, aidCS, MPI_param)
   type(Atom), dimension(:), allocatable, intent(inout) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   type(CDF), intent(inout) :: CDF_Phonon   ! CDF parameters for phonon to be read from a file
   type(Solid), intent(inout) :: Matter   ! all material parameters
   type(Density_of_states), intent(inout) :: Mat_DOS  ! materail DOS
   type(Differential_MFP), dimension(:), allocatable, intent(inout) :: DSF_DEMFP, DSF_DEMFP_H
   class(Ion), intent (out) :: SHI  ! we'll read an information about SHI from input-file
   real(8), intent(out) :: Tim !  [fs] total duration of the analysis
   real(8), intent(out) :: dt  !  [fs] timestep
   integer, intent(out) :: NMC ! number of MC iterations
   integer, intent(out) :: Num_th   ! number of threads for parralel calculations with openmp
   character(100), intent(out) :: Output_path, Output_path_SHI, Material_name   ! path for the output files
   type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
   logical, intent(out) :: read_well    ! did we read the file without an error?
   type(All_names), intent(out) :: File_names
   type(Flag), intent(inout) :: NumPar
   type(All_diff_CS), intent(inout) :: aidCS    ! all integrated differential cross sections
   type(Used_MPI_parameters), intent(inout) :: MPI_param ! MPI parameters
   !---------------------
   real(8) M, SUM_DOS
   integer FN, Reason, i, temp, temp1, temp2(1), temper, FN2, err_ind, access_mode
   character(100)   Material_file   ! file with material parameters
   character(100)   Short_material_file ! short version of the 'material parameters file'
   character(100)   File_name   ! file name if needed to use
   character(100)   DOS_file    ! files with material DOS and CDF
   character(100)   DSF_file, DSF_file_h    ! file with DSF differential cross-sections
   character(100)   Error_descript  ! to write a description of an error, if any
   character(100)   Temp_char, temp_char1
   character(3) Name    ! for reading elements names
   character(30) Full_Name    ! for reading full elements names
   character(100) command, temp_ch, File_name_INPUT, text   ! to pass to cmd a command
   logical file_exist    ! to check where file to be open exists
   logical file_opened   ! to check if a file is still opened
   !----------------
   ! initialize variables:
   file_exist = .false. ! to start with
   file_opened = .false.
   Material_file = ''
   Short_material_file = ''
   File_name = ''
   DOS_file = ''
   DSF_file = ''
   DSF_file_h = ''
   Error_descript = ''
   Temp_char = ''
   temp_char1 = ''
   Name = ''
   Full_Name = ''
   command = ''
   temp_ch = ''
   File_name_INPUT = ''
   text = ''

   !----------------
   ! Reading the input file:
   ! In MPI version, it must be done by one process (master process),
   ! and then the data read must be broadcasted to all other processes.
   File_name_INPUT = trim(adjustl(m_INPUT_file))

   FN = 200
   if (MPI_param%process_rank == 0) then   ! only MPI master process will read the input file
      inquire(file=trim(adjustl(File_name_INPUT)),exist=file_exist)     ! check if input file excists
   else
      goto 3333   ! other processes should not try to read the file
   endif

   if (file_exist) then
      open(unit = FN, FILE = trim(adjustl(File_name_INPUT)), status = 'old', IOSTAT=err_ind, readonly)   ! if yes, open it and read
      !call MPI_fileopen_wrapper(MPI_param, trim(adjustl(File_name_INPUT)), FN, Error_message, readonly=.true.)  ! WRONG

      if (err_ind /= 0) then ! some error while opening file
         Error_descript = 'File '//trim(adjustl(File_name_INPUT))//' could not be opened!'    ! description of an error
#ifdef MPI_USED
         Error_descript = trim(adjustl(Error_descript))//' [MPI process #'//trim(adjustl(MPI_param%rank_ch))//']'
#endif
         !call Save_error_details(Error_message, 2, Error_descript, MPI_param) ! write it into the error-log file
         print*, trim(adjustl(Error_descript)) ! print it also on the sreen
         read_well = .false.   ! it didn't read well the input file...
         goto 2013 ! go to the end of the subroutine, there is nothing else we could do
      endif
   else ! if no, save error message about it:
      File_name = 'OUTPUT_Error_log.dat'
      Error_message%File_Num = 100	! file number with error-log
      ! Open file with error log:
      open(unit = Error_message%File_Num, FILE = trim(adjustl(File_name)), IOSTAT=err_ind, action="write")

      Error_descript = 'File '//trim(adjustl(File_name_INPUT))//' is not found!'    ! description of an error
#ifdef MPI_USED
      Error_descript = trim(adjustl(Error_descript))//' [MPI process #'//trim(adjustl(MPI_param%rank_ch))//']'
#endif
      ! If file opened ok, write into it:
      if (err_ind == 0) then
         call Save_error_details(Error_message, 1, Error_descript, MPI_param) ! write it into the error-log file
      else
         print*, 'File '//trim(adjustl(File_name))//' could not be opened [MPI process #'//trim(adjustl(MPI_param%rank_ch))//']'
      endif
      print*, trim(adjustl(Error_descript)) ! print it also on the screen
      read_well = .false.   ! it didn't read well the input file...
      goto 2013 ! go to the end of the subroutine, there is nothing else we could do
   endif
   !----------------------------------------------
   ! Start reading the input file:
   i = 0
   READ(FN,*,IOSTAT=Reason) Material_name
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013

   READ(FN,*,IOSTAT=Reason) SHI%Zat   ! atomic number
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013
   
   if (SHI%Zat .GT. 0) then ! when 0, skip ion at all
       call Find_element_name(ABS(SHI%Zat), Name, Full_Name, M) ! from module 'Dealing_with_EADL'
       SHI%Name = Name  ! name of the element
       SHI%Full_Name = Full_Name  ! full name of the element
       SHI%Mass = M  ! mass of the element in the proton-mass units
   endif
   
   READ(FN,*,IOSTAT=Reason)  SHI%E ! [MeV] SHI energy
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013
   call Particle_event(SHI, E=SHI%E*1.0d6, t0=0.0d0, tn=0.0d0, X=0.0d0, Y=0.0d0, Z=0.0d0) ! convert [MeV] into [eV] and set initial data
   
   READ(FN,*,IOSTAT=Reason)  M   ! [proton mass] define mass if it's nonstandard, put negative value if it's standard
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013
   if (M .GT. 0.0d0) SHI%Mass = M  ! overwrite the standard mass only if needed!!
   
   READ(FN,*,IOSTAT=Reason)  Tim   ! [fs] total time
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013
      
   READ(FN,*,IOSTAT=Reason)  dt, NumPar%dt_flag    ! [fs] timestep, and kind of grid
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013
   if ((dt .GT. Tim) .AND. (NumPar%dt_flag .LE. 0)) dt = Tim    ! timestep can't be larger than the total time
   if ((dt .LT. 2) .AND. (NumPar%dt_flag .GE. 1)) dt = 2    ! in logarithmic grid, it shouldn't be smaller than at least 2
   
   READ(FN,*,IOSTAT=Reason) Matter%cut_off      ! cut-off energy (electrons with lower energies are excluded from calculation)
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013
   
   READ(FN,*,IOSTAT=Reason) Matter%Layer    ! [A] thikness of layer
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013
   
   READ(FN,*,IOSTAT=Reason) Matter%temp    ! [K] temperature of the target
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013

   temper = int(Matter%temp)    ! material temperature
   write(temp_char1, '(i)') temper
   ! File where DSF differential EMFP for this material are storred:
   DSF_file = trim(adjustl(m_INPUT_DSF))//trim(adjustl(NumPar%path_sep))// &
              trim(adjustl(Material_name))//trim(adjustl(NumPar%path_sep))// &
              trim(adjustl(Material_name))//'_Electron_DSF_Differential_EMFPs_'//trim(adjustl(temp_char1))//'K.dat'
   ! Same for holes:
   DSF_file_h = trim(adjustl(m_INPUT_DSF))//trim(adjustl(NumPar%path_sep))// &
                trim(adjustl(Material_name))//trim(adjustl(NumPar%path_sep))// &
                trim(adjustl(Material_name))//'_Hole_DSF_Differential_EMFPs_'//trim(adjustl(temp_char1))//'K.dat'
   
   READ(FN,*,IOSTAT=Reason) SHI%Kind_Zeff, SHI%fixed_Zeff   ! 0=Barkas; 1=Bohr; 2=Nikolaev-Dmitriev; 3=Schiwietz-Grande, 4 - fixed value;
   call read_file(Reason, i, read_well) ! reports if everything read well
   if ((SHI%fixed_Zeff .LE. 0.0) .OR. (SHI%fixed_Zeff .GT. SHI%Zat)) SHI%fixed_Zeff = SHI%Zat
   if (.not. read_well) goto 2013
   
   READ(FN,*,IOSTAT=Reason) SHI%Kind_ion   ! 0=point-like charge; 1=Brandt-Kitagawa ion
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013
   
   ! kind of elastic cross-section: -1=no elastic, 0=Mott; 1=CDF phonons; 2=DSF
   ! and kind of effective target charge (used in CDF): 0=effective Barkas-like; 1=Z=1 (old); 2=Z^2/CDF_e_e
   READ(FN,*,IOSTAT=Reason) NumPar%kind_of_EMFP, NumPar%CDF_elast_Zeff
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013
   
   READ(FN,*,IOSTAT=Reason) NumPar%kind_of_DR, Matter%El_eff_mass
   call read_file(Reason, i, read_well)                                     ! reports if everything read well
   if ((NumPar%kind_of_DR .LE. 0) .OR. (NumPar%kind_of_DR .GT. 4)) NumPar%kind_of_DR = 1         ! (hk)^2/2m dispersion relation by default
   if (Matter%El_eff_mass .LT. 0) Matter%El_eff_mass = 1.0d0                ! default value is electron mass
   if (.not. read_well) goto 2013
   
   ! Use maximal plasmon energy as upper integration limit in cross-sections calculations
   READ(FN,*,IOSTAT=Reason) temp1       ! 1 - include, everything else - exclude
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013
   if (temp1 .NE. 1) then
      NumPar%plasmon_Emax = .false. ! do not include
   else
      NumPar%plasmon_Emax = .true.  ! include
   endif
   
   READ(FN,*,IOSTAT=Reason) Matter%hole_mass   ! parameter of effective mass of valence hole, [me]
   !if (Matter%hole_mass .LE. 0.0d0) Matter%hole_mass = 1.0d25 ! infinite mass, immobile holes
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013
   
   READ(FN,*,IOSTAT=Reason) temp    ! include photons or not?
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013
   if (temp .NE. 1) then
      NumPar%include_photons = .false. ! no photons
   else
      NumPar%include_photons = .true.  ! include radiative decays of holes and all that follows
   endif
   

   ! This option is not ready, so it is excluded from release:
!    READ(FN,*,IOSTAT=Reason) NumPar%field_include, NumPar%field_dt
!    if (NumPar%field_include .NE. 1.0d0) NumPar%field_include = 0.0d0 ! fields included only when = 1
!    if (NumPar%field_dt .LT. 0.0d0) NumPar%field_dt = 0.01
!    call read_file(Reason, i, read_well) ! reports if everything read well
!    if (.not. read_well) goto 2013
   

   READ(FN,*,IOSTAT=Reason) Matter%work_function, Matter%bar_length, Matter%bar_height       ! Work function to specify surface barrier form
   if (Matter%work_function .LE. 0.0d0) Matter%work_function = 0.0d0 ! Emission included when > 0
   if (Matter%bar_length .LE. 0.0d0) Matter%work_function = 0.0d0
   if (Matter%bar_height .LE. 0.0d0) Matter%work_function = 0.0d0
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013
   
   READ(FN,*,IOSTAT=Reason) NMC ! number of MC iterations
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013
   
   READ(FN,*,IOSTAT=Reason) Num_th  ! number of threads for parallel openmp calculations
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013

#ifdef _OPENMP
   if (Num_th < 1) then ! use default: maximum number of available threads
      Num_th = omp_get_max_threads() ! number of processors available by default
   endif
#else
   Num_th = 1  ! no OMP => no threads
#endif

   !------------------------------------------------------
   ! Read optional parameters:
   read_well = .true.   ! to start with
   RDID: do while (read_well)
      read(FN,'(a)',IOSTAT=Reason) text
      call read_file(Reason, i, read_well, no_end=.true.)   ! below
      if (.not. read_well) exit RDID  ! end of file, stop reading

      call interpret_additional_data_INPUT(text, NumPar, MPI_param) ! below
   enddo RDID
   ! If gnuplotting is required:
   if (NumPar%do_gnuplot) then
      allocate(File_names%F(100))
   endif

3333 continue ! here go the processes that are not reading the input file by themselves

   !------------------------------------------------------
   ! Synchronize MPI processes: make sure all processes are here to recieve the input data from the master process
   call MPI_barrier_wrapper(MPI_param)  ! module "MPI_subroutines"
   !------------------------------------------------------
   ! MPI master process shares info with all worker-processes:
   call broadcast_input_data(NumPar, Matter, Material_name, SHI, Tim, dt, NMC, Num_th, MPI_param)   ! below
   !------------------------------------------------------

   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      if (NumPar%verbose) print*, 'The input file is read: '//trim(adjustl(File_name_INPUT))
   endif

   !------------------------------------------------------
   ! Create an output folder:
   Output_path = 'OUTPUT_'//trim(adjustl(Material_name)) ! that should be a folder with output
   if (SHI%Zat .GT. 0) Output_path_SHI = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//'OUTPUT_'//trim(adjustl(SHI%Name))//'_in_'//trim(adjustl(Material_name)) ! that should be a folder with output

   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
#ifndef __GFORTRAN__
      ! for intel fortran compiler:
      inquire(DIRECTORY=trim(adjustl(Output_path)),exist=file_exist)    ! check if input file excists
#else
      ! for gfortran compiler:
      inquire(FILE==trim(adjustl(Output_path)),exist=file_exist)    ! check if input file excists
#endif

      if (file_exist) then
         write(*,'(a,a,a)') ' Folder ', trim(adjustl(Output_path)), ' already exists, save output files into it'
      else
         write(*,'(a,a,a)') ' Folder ', trim(adjustl(Output_path)), ' was created, save output files into it'
         command='mkdir '//trim(adjustl(Output_path)) ! to create a folder use this command
         CALL system(command)  ! create the folder
      endif
   endif

   !------------------------------------------------------
   ! Synchronize MPI processes after creating the output directory:
   call MPI_barrier_wrapper(MPI_param)  ! module "MPI_subroutines"
   !------------------------------------------------------

   Error_message%File_Num = 101 ! this is the file where error-log will be storred
   if (SHI%Zat .GT. 0) then
      File_name = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(SHI%Name))//'_in_'//trim(adjustl(Material_name))//'_error_log.txt'
   else
      File_name = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(Material_name))//'_error_log.txt'
   endif

   ! Opening file for (possibly) parallel i/o in it:
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      open(unit = Error_message%File_Num, FILE = trim(adjustl(File_name)), action="write")
   endif
   !call MPI_fileopen_wrapper(MPI_param, trim(adjustl(File_name)), Error_message%File_Num, Error_message)  ! module "MPI_subroutines"

   !------------------------------------------------------
   ! Read CDF-material parameters:
   if (LEN(trim(adjustl(NumPar%CDF_file))) > 0 ) then ! user provided the name:
      Material_file = trim(adjustl(m_INPUT_CDF))//trim(adjustl(NumPar%path_sep))//trim(adjustl(NumPar%CDF_file))
      Short_material_file = trim(adjustl(m_INPUT_CDF))//trim(adjustl(NumPar%path_sep))//trim(adjustl(NumPar%CDF_file))
   else ! assume default name:
      Material_file = trim(adjustl(m_INPUT_CDF))//trim(adjustl(NumPar%path_sep))//trim(adjustl(Material_name))//'.cdf'
      Short_material_file = trim(adjustl(m_INPUT_CDF))//trim(adjustl(NumPar%path_sep))//trim(adjustl(Material_name))//'.scdf'
   endif

   call reading_material_parameters(Material_file, Short_material_file, Target_atoms, NumPar, CDF_Phonon, Matter, &
                                    aidCS, Error_message, read_well, MPI_param) ! below
   !------------------------------------------------------
   ! MPI master process shares if the file read well or not:
#ifdef MPI_USED
   call mpi_bcast(read_well, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, 'reading_material_parameters:{read_well}') ! module "MPI_subroutines"
#endif
   !------------------------------------------------------
   if (.not. read_well) goto 2015
   !------------------------------------------------------
   ! MPI master process shares (CDF) material info with all worker-processes:
   call broadcast_material_parameters(Target_atoms, NumPar, CDF_Phonon, Matter, aidCS, MPI_param)   ! below
   !------------------------------------------------------
   
   do i = 1, size(Target_atoms) ! all atoms
       temp2(1) = maxval(Target_atoms(i)%KOCS_SHI(:)) ! check if there is at least one shell for which we use BEB instead of CDF
       if (temp2(1) .GE. 2) then ! BEB cross sections are used, then:
         if (MPI_param%process_rank == 0) then   ! only MPI master process does it
            print*, 'Cannot use Brandt-Kitagawa model with BEB cross section, using point charge instead'
         endif
         SHI%Kind_ion = 0
         exit
      endif
   enddo

   !-----------------
   ! Printout warning for too high (relativistic) or too low energies:
   if (SHI%Zat .GT. 0) then ! if there is an SHI:
       ! If SHI energy is too small:
       if (SHI%E .LE. (SHI%Mass*g_Mp + g_me)*(SHI%Mass*g_Mp + g_me)/(SHI%Mass*g_Mp*g_me)*Target_atoms(1)%Ip(size(Target_atoms(1)%Ip))/4.0d0) then
           if (MPI_param%process_rank == 0) then   ! only MPI master process does it
               write(*,'(a,e,a,e,a)') 'The SHI energy is ', SHI%E, ' [eV] which is smaller than the minimum allowed energy of ', (SHI%Mass*g_Mp + g_me)*(SHI%Mass*g_Mp + g_me)/(SHI%Mass*g_Mp*g_me)*Target_atoms(1)%Ip(size(Target_atoms(1)%Ip))/4.0d0, ' [eV]'
               print*, 'Calculations cannot be performed, sorry. Try to increase the ion energy.'
           endif
           read_well = .false.  ! under such conditions, calculations cannot be performed
           goto 2015
       endif
       ! if SHI energy is too high (relativistic):
       if (SHI%E .GE. 175.0d6/2.0d0*SHI%Mass) then
           if (MPI_param%process_rank == 0) then   ! only MPI master process does it
               write(*,'(a,f17.3,a,f17.3,a)') 'The SHI energy is ', SHI%E, ' [eV] which is higher than the maximum allowed energy of ', 175.0d6/2.0d0*SHI%Mass, ' [eV]'
               print*, 'Calculations cannot be performed, sorry. Try to decrease the ion energy.'
           endif
           read_well = .false.  ! under such conditions, calculations cannot be performed
           goto 2015
       else if (SHI%E .GE. 175.0d6/4.0d0*SHI%Mass) then
           if (MPI_param%process_rank == 0) then   ! only MPI master process does it
               write(*,'(a,f17.3,a,f17.3,a)') 'The SHI energy is ', SHI%E, ' [eV] which is higher than the energy of ', 17.50d6*SHI%Mass, ' [eV]'
               print*, 'Note that the calculations might not be very reliable, since NO relativistic effects are included!'
           endif
       else
         ! all good
       endif
   endif

   !------------------------------------------------------
   ! Read VB DOS file:
   if (LEN(trim(adjustl(NumPar%DOS_file))) > 0 ) then ! user provided the name:
      DOS_file = trim(adjustl(m_INPUT_DOS))//trim(adjustl(NumPar%path_sep))//trim(adjustl(NumPar%DOS_file))
   else ! assume default name:
      DOS_file = trim(adjustl(m_INPUT_DOS))//trim(adjustl(NumPar%path_sep))//trim(adjustl(Material_name))//'.dos'
   endif
   call reading_material_DOS(DOS_file, Mat_DOS, Matter, Target_atoms, NumPar, Error_message, read_well, MPI_param)    ! read material DOS
   !------------------------------------------------------
   ! MPI master process shares if the file read well or not:
#ifdef MPI_USED
   call mpi_bcast(read_well, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, 'reading_material_DOS:{read_well}') ! module "MPI_subroutines"
#endif
   !------------------------------------------------------
   if (.not. read_well) goto 2015
   !------------------------------------------------------
   ! MPI master process shares (CDF) material info with all worker-processes:
   call broadcast_material_DOS(Mat_DOS, MPI_param)   ! below
   !------------------------------------------------------


   ! Printout DOS parameters:
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      temp_char1 = 'OUTPUT_'//trim(adjustl(Material_name))//'_DOS_analysis.dat'
      File_names%F(9) = trim(adjustl(temp_char1))  ! save for plotting later
      open(25, file = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(temp_char1)))
      write(25, *) '#E[eV]   k[1/m]   DOS[a.u.]    Int_DOS[a.u.]   Eff_mass[me]'
      if (Target_atoms(1)%Ip(size(Target_atoms(1)%Ip)) .LT. 0.2d0) then      ! Metal: Use DOS from bottom of CB to calculate dispersion relation
        do i = 1, size(Mat_DOS%E)
            write(25,'(e,e,e,e,e)') Mat_DOS%E(i), Mat_DOS%k_inv(i), Mat_DOS%DOS_inv(i), Mat_DOS%int_DOS_inv(i), Mat_DOS%Eff_m_inv(i)
        enddo
      else         !! Insulator or semoconductor: Use DOS from top of CB to calculate dispersion relation
        do i = 1, size(Mat_DOS%k)
            write(25,'(e,e,e,e,e)') Mat_DOS%E(i), Mat_DOS%k(i), Mat_DOS%DOS(i), Mat_DOS%int_DOS(i), Mat_DOS%Eff_m(i)
        enddo
      endif
      close(25)
   endif ! (MPI_param%process_rank == 0)


   ! DSF files:
   if (NumPar%kind_of_EMFP .EQ. 2) then
       call reading_DSF_cross_sections(DSF_file, DSF_DEMFP, NumPar, Error_message, read_well, MPI_param) ! below
       !------------------------------------------------------
       ! MPI master process shares (CDF) material info with all worker-processes:
       call broadcast_DSF_cross_sections(DSF_DEMFP, read_well, NumPar, MPI_param, 'DSF_DEMFP')   ! below
       !------------------------------------------------------
       if (.not. read_well) goto 2015

       ! Find out when this file was last modified:
       if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         call get_file_stat(trim(adjustl(DSF_file)), Last_modification_time=NumPar%Last_mod_time_DSF) ! above
       endif
       !print*, 'DSF file last modified on:', NumPar%Last_mod_time_DSF

       call reading_DSF_cross_sections(DSF_file_h, DSF_DEMFP_H, NumPar, Error_message, read_well, MPI_param)   ! below
       !------------------------------------------------------
       ! MPI master process shares (CDF) material info with all worker-processes:
       call broadcast_DSF_cross_sections(DSF_DEMFP_H, read_well, NumPar, MPI_param, 'DSF_DEMFP_H')   ! below
       !------------------------------------------------------
       if (.not. read_well) goto 2015
   else
      allocate(DSF_DEMFP(0))
      allocate(DSF_DEMFP_H(0))
   endif


   if (NumPar%CDF_elast_Zeff == 2) then   ! elastic cross section requires form factors:
      ! Independent reads for each process, no need for MPI-specific reads
      if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         FN2 = 26
         open(FN2, file = trim(adjustl(m_atomic_folder))//trim(adjustl(NumPar%path_sep))//trim(adjustl(m_form_factors_file)))
         call read_form_factors(FN2, Matter%form_factor, Error_message, read_well, MPI_param)  ! below
         close(26)
      endif
      !------------------------------------------------------
      ! MPI master process shares (CDF) material info with all worker-processes:
      call broadcast_form_factors(Matter%form_factor, read_well, MPI_param)   ! below
      !------------------------------------------------------
      if (.not. read_well) goto 2015
   endif
   
!2013 if (.not. read_well) print*, 'Error in INPUT_PARAMETERS.txt file or input files. See log!!'
2013 if (.not. read_well) then
      Error_descript = 'Error in '//trim(adjustl(m_INPUT_file)) //' file or input files. See log!!'
#ifdef MPI_USED
      Error_descript = trim(adjustl(Error_descript))//' [MPI process #'//trim(adjustl(MPI_param%rank_ch))//']'
#endif
      if (MPI_param%process_rank == 0) print*, trim(adjustl(Error_descript))
   endif

2015 continue   ! if we must skip to the end for some reason
   if (MPI_param%process_rank == 0) then   ! only MPI master process does it
      inquire(unit=FN,opened=file_opened)    ! check if this file is opened
      if (file_opened) close(FN)             ! and if it is, close it
   endif
end subroutine Read_input_file




subroutine broadcast_input_data(NumPar, Matter, Material_name, SHI, Tim, dt, NMC, Num_th, MPI_param)
   type(Flag), intent(inout) :: NumPar ! numerical parameters
   type(Solid), intent(inout) :: Matter   ! material properties
   character(*), intent(inout) :: Material_name   ! path for the output files
   class(Ion), intent (inout) :: SHI  ! we'll read an information about SHI from input-file
   real(8), intent(inout) :: Tim !  [fs] total duration of the analysis
   real(8), intent(inout) :: dt  !  [fs] timestep
   integer, intent(inout) :: NMC ! number of MC iterations
   integer, intent(inout) :: Num_th   ! number of threads for parralel calculations with openmp
   type(Used_MPI_parameters), intent(inout) :: MPI_param ! MPI parameters
   !-------------------------------------------
   integer :: Nsiz
   character(100):: error_part

#ifdef MPI_USED

   error_part = 'ERROR in broadcast_input_data:'

   ! Mater thread shares info read from the input file with all the other processes:

   Nsiz = LEN(Material_name)
   call mpi_bcast(Material_name, Nsiz, MPI_CHARACTER, 0,MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Material_name}') ! module "MPI_subroutines"
   !print*, '[process #', MPI_param%process_rank, ']', trim(adjustl(Material_name)), LEN(trim(adjustl(Material_name)))

   call mpi_bcast(SHI%Zat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {SHI%Zat}') ! module "MPI_subroutines"

   Nsiz = LEN(SHI%Name)
   call mpi_bcast(SHI%Name, Nsiz, &
                  MPI_CHARACTER, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {SHI%Name}') ! module "MPI_subroutines"

   Nsiz = LEN(SHI%Full_Name)
   call mpi_bcast(SHI%Full_Name, Nsiz, &
                  MPI_CHARACTER, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {SHI%Full_Name}') ! module "MPI_subroutines"

   call mpi_bcast(SHI%Mass, 1, &
                  MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {SHI%Mass}') ! module "MPI_subroutines"

   call mpi_bcast(SHI%E, 1, &
                  MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {SHI%E}') ! module "MPI_subroutines"

   ! For processes that didn't do it yet:
   if (MPI_param%process_rank /= 0) then
      call Particle_event(SHI, E=SHI%E, t0=0.0d0, tn=0.0d0, X=0.0d0, Y=0.0d0, Z=0.0d0) ! set initial data
   endif

   call mpi_bcast(Tim, 1, &
                  MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Tim}') ! module "MPI_subroutines"

   call mpi_bcast(dt, 1, &
                  MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {dt}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%dt_flag, 1, &
                  MPI_INTEGER, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%dt_flag}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%cut_off, 1, &
                  MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Matter%cut_off}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%Layer, 1, &
                  MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Matter%Layer}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%temp, 1, &
                  MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Matter%temp}') ! module "MPI_subroutines"

   call mpi_bcast(SHI%Kind_Zeff, 1, &
                  MPI_INTEGER, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {SHI%Kind_Zeff}') ! module "MPI_subroutines"

   call mpi_bcast(SHI%fixed_Zeff, 1, &
                  MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {SHI%fixed_Zeff}') ! module "MPI_subroutines"

   call mpi_bcast(SHI%Kind_ion, 1, &
                  MPI_INTEGER, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {SHI%Kind_ion}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%kind_of_EMFP, 1, &
                  MPI_INTEGER, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%kind_of_EMFP}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%CDF_elast_Zeff, 1, &
                  MPI_INTEGER, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%CDF_elast_Zeff}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%kind_of_DR, 1, &
                  MPI_INTEGER, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%kind_of_DR}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%El_eff_mass, 1, &
                  MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Matter%El_eff_mass}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%plasmon_Emax, 1, &
                  MPI_LOGICAL, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%plasmon_Emax}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%hole_mass, 1, &
                  MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Matter%hole_mass}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%include_photons, 1, &
                  MPI_LOGICAL, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%include_photons}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%work_function, 1, &
                  MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Matter%work_function}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%bar_length, 1, &
                  MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Matter%bar_length}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%bar_height, 1, &
                  MPI_DOUBLE_PRECISION, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Matter%bar_height}') ! module "MPI_subroutines"

   call mpi_bcast(NMC, 1, &
                  MPI_INTEGER, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NMC}') ! module "MPI_subroutines"

   call mpi_bcast(Num_th, 1, &
                  MPI_INTEGER, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Num_th}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%CS_method, 1, &
                  MPI_INTEGER, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%CS_method}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%out_dim, 1, &
                  MPI_INTEGER, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%out_dim}') ! module "MPI_subroutines"

   Nsiz = LEN(NumPar%CDF_file)
   call mpi_bcast(NumPar%CDF_file, Nsiz, &
                  MPI_CHARACTER, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%CDF_file}') ! module "MPI_subroutines"

   Nsiz = LEN(NumPar%DOS_file)
   call mpi_bcast(NumPar%DOS_file, Nsiz, &
                  MPI_CHARACTER, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%DOS_file}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%redo_IMFP, 1, &
                  MPI_LOGICAL, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%redo_IMFP}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%redo_EMFP, 1, &
                  MPI_LOGICAL, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%redo_EMFP}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%redo_IMFP_SHI, 1, &
                  MPI_LOGICAL, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%redo_IMFP_SHI}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%get_thermal, 1, &
                  MPI_LOGICAL, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%get_thermal}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%print_CDF, 1, &
                  MPI_LOGICAL, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%print_CDF}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%print_CDF_optical, 1, &
                  MPI_LOGICAL, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%print_CDF_optical}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%verbose, 1, &
                  MPI_LOGICAL, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%verbose}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%very_verbose, 1, &
                  MPI_LOGICAL, 0, &
                  MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {NumPar%very_verbose}') ! module "MPI_subroutines"
#endif

end subroutine broadcast_input_data




subroutine read_form_factors(FN, form_factor, Err, read_well, MPI_param)
   integer, intent(in) :: FN  ! file number with the database
   real(8), dimension(:,:), allocatable, intent(inout) :: form_factor   ! coefficients used to construct form factors
   type(Error_handling), intent(inout) :: Err	! error log
   logical, intent(inout) :: read_well
   type(Used_MPI_parameters), intent(inout) :: MPI_param ! MPI parameters
   !----------------------------------
   real(8), dimension(:), allocatable :: read_vec	! to read the coeffs
   integer :: N_line, i, j, Reason, count_lines
   character(200) :: Error_descript, temp

   ! Count how many lines are in the file:
   call SkipCount_lines_in_file(FN, N_line, skip_lines=1)	! module "Dealing_with_files"
   ! Allocate array of copefficients:
   allocate(form_factor(N_line,5))

   allocate(read_vec(5))  ! temp
   count_lines = 1
   read(FN,*)  ! skip the first line with description

   do i = 1, N_line	! read the pair creation coeffs
      read(FN,*,IOSTAT=Reason) read_vec(:)   ! read into temporary array
      call read_file(Reason, count_lines, read_well)	! module "Dealing_with_files"

      if (read_well) then
         ! Save in the array:
         form_factor(i,:) = read_vec(:)   ! save into working array
      else
         write(Error_descript,'(a,i3)') 'In read_form_factors could not read line ', count_lines
#ifdef MPI_USED
         Error_descript = trim(adjustl(Error_descript))//' [MPI process #'//trim(adjustl(MPI_param%rank_ch))//']'
#endif
         call Save_error_details(Err, 2, Error_descript, MPI_param)	! module "Objects"
         return
      endif
   enddo
   rewind(FN)  ! for future use of the file

   ! Clean up at the end:
   deallocate(read_vec)
end subroutine read_form_factors


subroutine broadcast_form_factors(form_factor, read_well, MPI_param)   ! below
   real(8), dimension(:,:), allocatable, intent(inout) :: form_factor   ! coefficients used to construct form factors
   logical, intent(inout) :: read_well
   type(Used_MPI_parameters), intent(inout) :: MPI_param ! MPI parameters
   !---------------------------------------
   integer :: N1, N2
   character(100):: error_part

#ifdef MPI_USED
   error_part = 'ERROR in broadcast_form_factors:'
   call mpi_bcast(read_well, 1, MPI_LOGICAL, 0,MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {form_factor:read_well}') ! module "MPI_subroutines"

   ! Mater thread shares info read from the input file with all the other processes:
   if (MPI_param%process_rank == 0) then
      N1 = size(form_factor,1)
      N2 = size(form_factor,2)
   else  ! unknown sze yet
      N1 = 0
      N2 = 0
   endif
   call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {form_factor:N1}') ! module "MPI_subroutines"
   call mpi_bcast(N2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {form_factor:N2}') ! module "MPI_subroutines"

   if (MPI_param%process_rank /= 0) then
      if (.not.allocated(form_factor)) allocate(form_factor(N1,N2))
   endif

   call mpi_bcast(form_factor, N1*N2, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {form_factor}') ! module "MPI_subroutines"
#endif
end subroutine broadcast_form_factors




subroutine interpret_additional_data_INPUT(text_in, NumPar, MPI_param)
   character(*), intent(in) :: text_in ! text read from the file
   type(Flag), intent(inout) :: NumPar ! numerical parameters
   type(Used_MPI_parameters), intent(in) :: MPI_param ! MPI parameters
   !------------------------------------
   character(100) :: text, text2, ch_temp
   integer :: i, Reason, i_read
   logical :: read_well

   i = 0 ! to start with

   ! Read from the variable, in case there are more then one flag in one line:
   read(text_in, *, IOSTAT=Reason) text
   call read_file(Reason, i, read_well, do_silent=.true.) ! reports if everything read well
   if (.not. read_well) then  ! nothing more to do, skip this line
      return
   endif

   select case (text)
   case ('grid', 'GRID', 'Grid')
      ! Try to read which file extension to use:
      read(text_in, *, IOSTAT=Reason) text, i_read
      call read_file(Reason, i, read_well, do_silent=.true.) ! reports if everything read well
      if (.not. read_well) then  ! by default, use tabulated files
         !NumPar%CS_method = 1   ! was specified above in the default NumPar
         if (MPI_param%process_rank == 0) then   ! only MPI master process does it
            print*, 'Could not interprete grid index in line: ', trim(adjustl(text_in)), ', using default'
         endif
      else ! user provided grid index
         NumPar%CS_method = i_read
      endif

   case ('UNITS', 'Units', 'units')
      read(text_in, *, IOSTAT=Reason) text, i_read
      call read_file(Reason, i, read_well, do_silent=.true.) ! reports if everything read well
      if (.not. read_well) then  ! default
         if (MPI_param%process_rank == 0) then   ! only MPI master process does it
            print*, 'Could not interprete units index in line: ', trim(adjustl(text_in)), ', using default'
         endif
      else  ! use the provided name for DOS file:
         NumPar%out_dim = i_read
      endif

   case ('CDF', 'Cdf', 'cdf')
      read(text_in, *, IOSTAT=Reason) text, text2
      call read_file(Reason, i, read_well, do_silent=.true.) ! reports if everything read well
      if (.not. read_well) then  ! default
         if (MPI_param%process_rank == 0) then   ! only MPI master process does it
            print*, 'Could not interprete CDF-file name in line: ', trim(adjustl(text_in)), ', using default'
         endif
      else  ! use the provided name for DOS file:
         NumPar%CDF_file = trim(adjustl(text2))
      endif

   case ('DOS', 'Dos', 'dos')
      read(text_in, *, IOSTAT=Reason) text, text2
      call read_file(Reason, i, read_well, do_silent=.true.) ! reports if everything read well
      if (.not. read_well) then  ! default
         if (MPI_param%process_rank == 0) then   ! only MPI master process does it
            print*, 'Could not interprete DOS-file name in line: ', trim(adjustl(text_in)), ', using default'
         endif
      else  ! use the provided name for DOS file:
         NumPar%DOS_file = trim(adjustl(text2))
      endif

   case ('gnuplot', 'plot', 'gnu', 'GNUPLOT', 'PLOT', 'GNU')
      NumPar%do_gnuplot = .true. ! use gnuplot to create plots

      ! Try to read which file extension to use:
      read(text_in, *, IOSTAT=Reason) text, text2
      call read_file(Reason, i, read_well, do_silent=.true.) ! reports if everything read well
      if (.not. read_well) then  ! by default, use this extension
         NumPar%plot_extension = 'jpeg'
      else ! use the provided extension
         NumPar%plot_extension = trim(adjustl(text2))
      endif
      ! Check, if gnuplot was switched off:
      select case (trim(adjustl(NumPar%plot_extension)))
      case ('NO', 'No', 'no')
         NumPar%do_gnuplot = .false. ! no gnuplot to create plots
         if (MPI_param%process_rank == 0) then   ! only MPI master process does it
            print*, 'No gnuplot scripts will be created'
         endif
      case default
         if (MPI_param%process_rank == 0) then   ! only MPI master process does it
            print*, "Gnuplot scripts will be created with extension '."//trim(adjustl(NumPar%plot_extension))//"'"
         endif
      end select

   case ('redo_MFP', 'REDO_MFP', 'Redo_MFP', 'redo_mfp')
      NumPar%redo_IMFP = .true. ! recalculate inelastic MFPs, if possible
      NumPar%redo_EMFP = .true. ! recalculate elastic MFPs, if possible
      NumPar%redo_IMFP_SHI = .true. ! recalculate SHI inelastic MFPs, if possible
      if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         print*, 'Mean free paths will be recalculated'
      endif

   case ('redo_MFP_SHI', 'REDO_MFP_SHI', 'Redo_MFP_SHI', 'redo_mfp_shi', 'redo_IMFP_SHI', 'REDO_IMFP_SHI', 'Redo_IMFP_SHI', 'redo_imfp_shi')
      NumPar%redo_IMFP_SHI = .true. ! recalculate SHI inelastic MFPs, if possible
      if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         print*, 'Inelastic Mean free paths will be recalculate'
      endif

   case ('redo_IMFP', 'REDO_IMFP', 'Redo_IMFP', 'redo_imfp')
      NumPar%redo_IMFP = .true. ! recalculate inelastic MFPs, if possible
      if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         print*, 'Inelastic Mean free paths will be recalculated'
      endif

   case ('redo_EMFP', 'REDO_EMFP', 'Redo_EMFP', 'redo_emfp')
      NumPar%redo_EMFP = .true. ! recalculate elastic MFPs, if possible
      if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         print*, 'Elastic Mean free paths will be recalculated'
      endif

   case ('get_thermal', 'thermal', 'make_thermal', 'Get_thermal', 'Thermal', 'Make_thermal')
      NumPar%get_thermal = .true.
      if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         print*, 'Thermal parameters will be calculated'
      endif

   case ('print_CDF', 'Print_CDF', 'print_cdf', 'PRINT_CDF')
      NumPar%print_CDF = .true.
      if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         print*, 'File with CDF parameters will be printed out'
      endif

   case ('print_optical', 'Print_optical', 'Print_Optical', 'print_optical_cdf', 'PRINT_OPTICAL_CDF')
      NumPar%print_CDF_optical = .true.
      if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         print*, 'File with optical CDF will be printed out'
      endif

   case ('Verbose', 'verbose', 'VERBOSE')
      NumPar%verbose = .true.
      if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         print*, 'Verbose option is on, TREKIS will print a lot of extra stuff...'
      endif

   case ('Very_verbose', 'very_verbose', 'VERY_VERBOSE', 'Very_Verbose', 'Very-verbose', 'very-verbose', 'VERY-VERBOSE', 'Very-Verbose')
      NumPar%verbose = .true.
      NumPar%very_verbose = .true.
      if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         print*, 'Very-verbose option is on, TREKIS will print A LOT of extra stuff...'
      endif

   case ('INFO', 'Info', 'info')
      if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         print*, trim(adjustl(starline))
         print*, trim(adjustl(starline))
         print*, 'TREKIS-3 stands for: Time-Resolved Electron Kinteics in SHI-Irradiated Solids'
         print*, 'For all details and instruction, address the files !READ_ME_TREKIS_3.doc or !READ_ME_TREKIS_3.pdf'
         print*, 'DISCLAIMER'
         print*, 'Although we endeavour to ensure that the code TREKIS-3 and results delivered are correct, no warranty is given as to its accuracy. We assume no responsibility for possible errors or omissions. We shall not be liable for any damage arising from the use of this code or its parts or any results produced with it, or from any action or decision taken as a result of using this code or any related material.', &
         'This code is distributed as is for non-commercial peaceful purposes only, such as research and education. It is explicitly prohibited to use the code, its parts, its results or any related material for military-related and other than peaceful purposes. By using this code or its materials, you agree with these terms and conditions.'
         print*, 'HOW TO CITE'
         print*, 'The use of the code is at your own risk. Should you choose to use it, please cite the following works.'
         print*, 'The code itself: ', &
         'N. Medvedev, R. Rymzhanov, A.E. Volkov, (2023). TREKIS-3 [Computer software]. https://doi.org/10.5281/zenodo.8394462', &
         ' 1) N. Medvedev, R. A. Rymzhanov, A. E. Volkov, J. Phys. D. Appl. Phys. 48 (2015) 355303', &
         ' 2) R. A. Rymzhanov, N. Medvedev, A. E. Volkov, Nucl. Instrum. Methods B 388 (2016) 41', &
         'Should you use this code to create initial conditions for further molecular dynamics simulations of atomic response to the electronic excitation by a swift heavy ion (e.g. with LAMMPS), the following citation is required:', &
         ' 3) R. Rymzhanov, N. Medvedev, A. E. Volkov, J. Phys. D. Appl. Phys. 50 (2017) 475301', &
         'In a publication, we recommend that at least the following parameters should be mentioned for reproducibility of the results: material, its structure, density, speed of sound, the used CDF coefficients, which processes were included (active) in the simulation, ion type, its energy, the model for SHI charge, number of MC iterations.'
         print*, trim(adjustl(starline))
         print*, trim(adjustl(starline))
      endif
   endselect
end subroutine interpret_additional_data_INPUT




! Reads additional data from the command line passed along with the XTANT:
subroutine get_add_data(NumPar, MPI_param)
   type(Flag), intent(inout) :: NumPar ! numerical parameters
   type(Used_MPI_parameters), intent(in) :: MPI_param ! MPI parameters
   !---------------
   character(1000) :: string
   integer :: i_arg, count_args, N_arg

   ! Count how many arguments the user provided:
   N_arg = COMMAND_ARGUMENT_COUNT() ! Fortran intrinsic function

   count_args = 0 ! to start with

   ALLARG:do i_arg = 1, N_arg ! read all the arguments passed
      ! Read the argument provided:
      call GET_COMMAND_ARGUMENT(i_arg, string)  ! intrinsic

      ! Act on the command passed:
      call interpret_additional_data_INPUT(string, NumPar, MPI_param) ! above

   enddo ALLARG
end subroutine get_add_data



subroutine print_time_step(text, MPI_param, num, msec, print_to)
   CHARACTER(len=*) :: text   ! text to print out
   type(Used_MPI_parameters), intent(in) :: MPI_param ! MPI parameters
   real(8), intent(in), optional :: num   ! to print out this number
   logical, intent(in), optional :: msec  ! print msec or not?
   integer, intent(in), optional :: print_to ! file number to print to
   !--------------------------------
   character(len=100) :: var
   integer :: FN, c1(8) ! time stamps

   if (MPI_param%process_rank /= 0) then   ! only MPI master process does it
      return   ! other processes do not do anything here
   endif

   if (present(print_to)) then   ! file number
      FN = print_to
   else  ! screen
      FN = 6
   endif

   call date_and_time(values=c1) ! standard FORTRAN time and date
   if (present(num) .and. present(msec)) then
      write(var,'(f8.2)') num
      write(FN, 1002) trim(adjustl(text))//' '//trim(adjustl(var))//' at ', c1(5), c1(6), c1(7), c1(8), c1(3), c1(2), c1(1)
   elseif (present(msec)) then
      write(FN, 1002) trim(adjustl(text))//' at ', c1(5), c1(6), c1(7), c1(8), c1(3), c1(2), c1(1)
   elseif (present(num)) then
      write(var,'(f8.2)') num
      write(FN, 1001) trim(adjustl(text))//' '//trim(adjustl(var))//' at ', c1(5), c1(6), c1(7), c1(3), c1(2), c1(1)
   else
      write(FN, 1001) trim(adjustl(text))//' at ', c1(5), c1(6), c1(7), c1(3), c1(2), c1(1)
   endif
   ! Formats for printing:
   1001 format (a, i2.2, ':', i2.2, ':', i2.2, ' on ', i2.2, '/', i2.2, '/', i4.4)
   1002 format (a, i2.2, ':', i2.2, ':', i2.2, ':', i3.3, ' on ', i2.2, '/', i2.2, '/', i4.4)
end subroutine print_time_step



subroutine reading_material_parameters(Material_file, Short_material_file, Target_atoms, NumPar, CDF_Phonon, Matter, &
                                       aidCS, Error_message, read_well, MPI_param)
   character(100), intent(in) :: Material_file, Short_material_file  ! files with material parameters (full or short)
   type(Atom), dimension(:), allocatable, intent(inout) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   type(Flag), intent(inout) :: NumPar ! numerical parameters
   type(CDF), intent(inout) :: CDF_Phonon   ! CDF parameters for phonon to be read from a file
   type(Solid), intent(inout) :: Matter   ! all material parameters
   type(All_diff_CS), intent(inout) :: aidCS    ! all integrated differential cross sections
   type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
   logical, intent(inout) :: read_well  ! did we read the file well?
   type(Used_MPI_parameters), intent(inout) :: MPI_param ! MPI parameters
   !-------------------------
   real(8) :: M, temp_r(5)
   integer FN2, Reason, i, j, k, l, N, Shl, CDF_coef, Shl_num, comment_start, Nsiz
   character(200) Error_descript, temp, read_line
   character(3) Name
   character(11) Shell_name
   character(30) Full_Name
   logical file_opened, file_exist, file_exist2
   type(CDF) :: Ritchi_temp  ! Ritchi CDF
   
   read_well = .true.   ! no problem here

   if (MPI_param%process_rank /= 0) then
      return   ! nothing to do for non-master processes
   endif


   FN2 = 201

   inquire(file=trim(adjustl(Material_file)),exist=file_exist)    ! check if the full input file exists
   inquire(file=trim(adjustl(Short_material_file)),exist=file_exist2)    ! check if the short input file exists
   if (file_exist) then
      ! Independent read-only, no need for MPI-specific fileopen:
      open(unit = FN2, FILE = trim(adjustl(Material_file)), status = 'old', readonly)   ! yes, open full file and read
      !call MPI_fileopen_wrapper(MPI_param, trim(adjustl(Material_file)), FN2, Error_message, readonly=.true.)  ! module "MPI_subroutines"

      ! Save CDF-file name for output:
      NumPar%CDF_file = trim(adjustl(Material_file))
      ! Find out when this file was last modified:
      call get_file_stat(trim(adjustl(Material_file)), Last_modification_time=NumPar%Last_mod_time_CDF) ! above
      !print*, 'CDF file last modified on:', NumPar%Last_mod_time_CDF
   else if (file_exist2) then ! if at least short version exists, the rest can be used within atomic approximation
      ! Independent read-only, no need for MPI-specific fileopen:
      open(unit=FN2, FILE = trim(adjustl(Short_material_file)), status = 'old', readonly)   ! yes, open full file and read
      !call MPI_fileopen_wrapper(MPI_param, trim(adjustl(Short_material_file)), FN2, Error_message, readonly=.true.)  ! module "MPI_subroutines"

      ! Save CDF-file name for output:
      NumPar%CDF_file = trim(adjustl(Short_material_file))
      ! Find out when this file was last modified:
      call get_file_stat(trim(adjustl(Short_material_file)), Last_modification_time=NumPar%Last_mod_time_CDF) ! above
      call read_short_scdf(FN2, Target_atoms, NumPar, CDF_Phonon, Matter, Error_message, read_well, MPI_param) ! below
      goto 2014 ! short version is done, skip reading the long one below
   else ! if no, save error message about it:
      !if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         Error_descript = 'File '//trim(adjustl(Material_file))//' is not found!'    ! description of an error
#ifdef MPI_USED
         Error_descript = trim(adjustl(Error_descript))//' [MPI process #'//trim(adjustl(MPI_param%rank_ch))//']'
#endif
         print*, trim(adjustl(Error_descript)) ! print it also on the sreen
         call Save_error_details(Error_message, 2, Error_descript, MPI_param) ! write it into the error-log file
      !endif
      read_well = .false.   ! it didn't read well the input file...
      goto 2014 ! go to the end of the subroutine, there is nothing else we could do
   endif
   
   ! Defaults:
   NumPar%VB_CDF_defined = .false.
   !NumPar%phonon_CDF_defined = .false. ! unused

   i = 0
   !READ(FN2,*,IOSTAT=Reason) Matter%Target_name ! first line is the full material name
   Matter%Target_name = '' ! to start with
   read_line = '' ! to start with
   READ(FN2,'(a)',IOSTAT=Reason) read_line   ! read the full line
   comment_start = index(read_line, '!')  ! to exclude comment, if any
   ! Save the name into the variable:
   if (comment_start > 1) then   ! trim comment
      Matter%Target_name = read_line(1:comment_start-1)
   else  ! no comment
      Matter%Target_name = trim(adjustl(read_line))
   endif
   ! Also, remove TABs if any:
   comment_start = index(read_line, CHAR(9))  ! to exclude TABs, if any
   if (comment_start > 1) then   ! trim comment
      Matter%Target_name = read_line(1:comment_start-1)
   elseif (comment_start == 1) then   ! omit leading TAB
      Matter%Target_name = read_line(2:)
   else  ! no comment
      Matter%Target_name = trim(adjustl(read_line))
   endif
   Matter%Target_name = trim(adjustl(Matter%Target_name))

!    print*, LEN(trim(adjustl(Matter%Target_name))), trim(adjustl(Matter%Target_name))
!    print*, '::'//Matter%Target_name//'::'
!    pause 'Matter%Target_name'


   Matter%Chem = ''  ! to start with
   
   READ(FN2,*,IOSTAT=Reason) N   ! number of elements in this compound
   IF (Reason .GT. 0)  THEN ! if it's a formula:
      backspace(FN2) ! read this line again:
      temp = ''   ! to restart
      READ(FN2,*,IOSTAT=Reason) temp   ! chemical formula of this compound
      Matter%Chem = trim(adjustl(temp))   ! save chemical formula of the materials
      call read_file(Reason, i, read_well) ! reports if everything read well
      if (.not. read_well) goto 2014
      call Decompose_compound(temp, Numpar%path_sep, Target_atoms, read_well) ! from module "Dealing_with_EADL"
      if (.not.read_well) goto 2014
      N = size(Target_atoms) ! that's how many atoms we have

      if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         print*, 'The material formula ', trim(adjustl(temp)), ' was interpreted as:'
         do j = 1, N  ! read for each element it's basic data:
            write(*,'(a,a,a,i3,a,f5.1,a,f7.2)') 'Atom ', Target_atoms(j)%Name, ' number #', Target_atoms(j)%Zat, ', composition: ', Target_atoms(j)%Pers, ', mass', Target_atoms(j)%Mass
         enddo
      endif
   else ! if it is a number:
      call read_file(Reason, i, read_well) ! reports if everything read well
      if (.not. read_well) goto 2014
      if (.not. allocated(Target_atoms)) allocate(Target_atoms(N)) ! that's how many atom kinds we have
      do j = 1, N  ! read for each element it's basic data:
        READ(FN2,*,IOSTAT=Reason) Target_atoms(j)%Zat, Target_atoms(j)%Pers ! atmoic number and its persentage in the compound's molecule
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2014
        call Find_element_name(Target_atoms(j)%Zat, Name, Full_Name, M) ! from module 'Dealing_with_EADL'
        Target_atoms(j)%Name = Name  ! name of the element
        Target_atoms(j)%Full_Name = Full_Name  ! full name of the element
        Target_atoms(j)%Mass = M  ! mass of the element in the proton-mass units
        ! Construct chemical formula:
        if ( abs(Target_atoms(j)%Pers-1.0d0) > 1.0d-6 ) then
         write(temp,'(f10.2)') Target_atoms(j)%Pers
        else
         temp = ''
        endif
        Matter%Chem = trim(adjustl(Matter%Chem))//trim(adjustl(Target_atoms(j)%Name))//trim(adjustl(temp))
      enddo
   endif

   ! material density [g/cm^3], speed of sound in the material [m/s], Fermi energy [eV]; bandgap [eV]
   temp = ''   ! to restart
   READ(FN2,'(a)',IOSTAT=Reason) temp
   READ(temp,*,IOSTAT=Reason) Matter%Dens, Matter%Vsound, Matter%E_F, Matter%Egap
   if (Reason /= 0) then   ! try to read only 3 numbers
      READ(temp,*,IOSTAT=Reason) Matter%Dens, Matter%Vsound, Matter%E_F
      Matter%Egap = 0.0d0  ! assume metal
   endif
   if (Matter%Egap .LT. 1.0d-1) Matter%Egap = 1.0d-1 ! [eV], introduce at least a minimum "band gap"
   call read_file(Reason, i, read_well) ! reports if everything read well

   if (.not. read_well) goto 2014
   !Matter%At_Dens = 1.0d-3*Matter%Dens/(g_Mp*SUM(Target_atoms(:)%Mass)/size(Target_atoms)) ! [1/cm^3] atomic density
   Matter%At_Dens = 1.0d-3*Matter%Dens/(g_Mp*SUM(Target_atoms(:)%Mass*Target_atoms(:)%Pers)/SUM(Target_atoms(:)%Pers))
   Matter%v_f = sqrt(2.0d0*Matter%E_f/g_me) ! units as used in one-pole approximation (not SI!)
  

   !-----------
   ! Check if there is a flag, or a number:
   temp = ''   ! to restart
   READ(FN2,'(a)',IOSTAT=Reason) temp
   READ(temp,*,IOSTAT=Reason) Shl ! number of shells in the first atom kind:
   call read_file(Reason, i, read_well, do_silent=.true.) ! reports if everything read well

   SP_CDF:if (.not. read_well) then ! check if there is a CDF
      ! Defaults:
      NumPar%kind_of_CDF = 1  ! single-pole CDF
      NumPar%kind_of_CDF_ph = 1 ! use single-pole approximation of phonon CDF

      ! Check if there are partial CDF defined:
      ! markers "VALENCE" or "PHONON"
      backspace(FN2) ! read this line again:
      Reason = 0  ! restart
      do while (Reason == 0)
         READ(FN2,'(a)',IOSTAT=Reason) temp
         IF (Reason == 0)  THEN ! normal reading, try to interprete
            if (LEN(trim(adjustl(temp))) > 0) then ! there is something written
               select case( trim(adjustl(temp(1:3))) )
               !--------------
               case ('VAL', 'Val', 'val') ! valence CDF
                  NumPar%VB_CDF_defined = .true.
                  !print*, 'VB', NumPar%VB_CDF_defined

                  READ(FN2,*,IOSTAT=Reason) CDF_coef, Shl_num, temp_r(1), temp_r(2), temp_r(3)
                  if ((Reason == 0) .and. (CDF_coef > 0)) then ! there is CDF to be used in a cross section
                     allocate(Ritchi_temp%E0(CDF_coef))
                     allocate(Ritchi_temp%A(CDF_coef))
                     allocate(Ritchi_temp%Gamma(CDF_coef))
                     if (NumPar%kind_of_DR == 4) then
                        allocate(Ritchi_temp%alpha(CDF_coef))
                     endif
                     do l = 1, CDF_coef    ! read all the CDF parameters:
                        READ(FN2,*,IOSTAT=Reason) Ritchi_temp%E0(l), Ritchi_temp%A(l), Ritchi_temp%Gamma(l)
                        call read_file(Reason, i, read_well) ! reports if everything read well
                        if (.not. read_well) then
                           write(*,'(a)') ' Could not interprete VB CDF parameters in '//trim(adjustl(Material_file))//'. Using single-pole approximation.'
                           NumPar%VB_CDF_defined = .false.
                           exit
                        endif
                        ! For Delta-CDF:
                        if (NumPar%kind_of_DR == 4) then
                           call define_alpha(Ritchi_temp%A(l), Ritchi_temp%Gamma(l), Ritchi_temp%E0(l), Matter%Egap, Ritchi_temp%alpha(l))  ! below
                        endif
                     enddo
                  else   ! something went wrong, use single-pole CDF
                     write(*,'(a)') ' Could not interprete VB CDF parameters in '//trim(adjustl(Material_file))//'. Using single-pole approximation.'
                     NumPar%VB_CDF_defined = .false.
                  endif

               !--------------
               case ('PHO', 'Pho', 'pho') ! phonon CDF
                  NumPar%kind_of_CDF_ph = 0 ! user provided phonon CDF
                  !print*, 'Phonons:', NumPar%kind_of_CDF_ph

                  READ(FN2,*,IOSTAT=Reason) CDF_coef
                  if (.not. allocated(CDF_Phonon%A)) allocate(CDF_Phonon%A(CDF_coef))
                  if (.not. allocated(CDF_Phonon%E0)) allocate(CDF_Phonon%E0(CDF_coef))
                  if (.not. allocated(CDF_Phonon%Gamma)) allocate(CDF_Phonon%Gamma(CDF_coef))
                  if (NumPar%kind_of_DR .EQ. 4) then
                     if (.not. allocated(CDF_Phonon%alpha)) allocate(CDF_Phonon%alpha(CDF_coef))   ! that's how many CDF functions
                  endif
                  do l = 1, CDF_coef    ! read all the CDF parameters for phonon peak:
                     READ(FN2,*,IOSTAT=Reason) CDF_Phonon%E0(l), CDF_Phonon%A(l), CDF_Phonon%Gamma(l)
                     call read_file(Reason, i, read_well) ! reports if everything read well
                     if (.not. read_well) then
                        write(*,'(a)') ' Could not interprete phonon CDF parameters in '//trim(adjustl(Material_file))//'. Using single-pole approximation.'
                        deallocate(CDF_Phonon%A)
                        deallocate(CDF_Phonon%E0)
                        deallocate(CDF_Phonon%Gamma)
                        if (allocated(CDF_Phonon%alpha)) deallocate(CDF_Phonon%alpha)
                        NumPar%kind_of_CDF_ph = 1 ! use single-pole
                        exit
                     endif
                     ! For Delta-CDF:
                     if (NumPar%kind_of_DR .EQ. 4) then
                     call define_alpha(CDF_Phonon%A(l), CDF_Phonon%Gamma(l), CDF_Phonon%E0(l), 0.0d0, CDF_Phonon%alpha(l))  ! below
                     endif
                  enddo

               end select
            endif ! (LEN(trim(adjustl(temp))) > 0)
         ENDIF ! (Reason == 0)
      enddo ! while (Reason == 0)
      !pause 'TEST CDF'

      if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         if (.not.NumPar%VB_CDF_defined) then
            write(*,'(a)') ' No CDF parameters found in the file '//trim(adjustl(Material_file))//'. Using single-pole approximation.'
         else
            write(*,'(a)') ' Only VB CDF found in the file '//trim(adjustl(Material_file))//'. For the rest, using single-pole approximation.'
         endif
      endif

      ! Read the atomic parameters from EPICS-database:
      call check_atomic_parameters(NumPar, Target_atoms, Error_message=Error_message, read_well=read_well, MPI_param=MPI_param) ! from module 'Dealing_with_EADL'

      ! Construct the valence band form the valence atomic shells:
      call make_valence_band(Target_atoms, NumPar, Matter, Error_message, read_well) ! below

      ! Having all the atomic parameters, allocate CDFs:
      do j = 1, N  ! for each element
         Shl = size(Target_atoms(j)%Ip)

         if (.not. allocated(Target_atoms(j)%Ritchi)) then
            allocate(Target_atoms(j)%Ritchi(Shl)) ! allocate Ritchi-functions' coefficiants for each shell
         endif

         do k = 1, Shl ! each shell
            Target_atoms(j)%KOCS_SHI(k) = 1   ! Mark SHI CDF cross section for this shell
            Target_atoms(j)%KOCS(k) = 1      ! Mark electron CDF cross section for this shell

            ! Special treatement for VB - it may be defined already:
            if ( (j == 1) .and. (k == Shl) ) then ! VB
               if ( NumPar%VB_CDF_defined ) then ! and only then
                  Nsiz = size(Ritchi_temp%A)
                  if (.not. allocated(Target_atoms(j)%Ritchi(k)%E0)) allocate(Target_atoms(j)%Ritchi(k)%E0(Nsiz))
                  if (.not. allocated(Target_atoms(j)%Ritchi(k)%A)) allocate(Target_atoms(j)%Ritchi(k)%A(Nsiz))
                  if (.not. allocated(Target_atoms(j)%Ritchi(k)%Gamma)) allocate(Target_atoms(j)%Ritchi(k)%Gamma(Nsiz))
                  Target_atoms(j)%Ritchi(k)%E0 = Ritchi_temp%E0
                  Target_atoms(j)%Ritchi(k)%A = Ritchi_temp%A
                  Target_atoms(j)%Ritchi(k)%Gamma = Ritchi_temp%Gamma
                  if (NumPar%kind_of_DR == 4) then
                     if (.not. allocated(Target_atoms(j)%Ritchi(k)%alpha)) allocate(Target_atoms(j)%Ritchi(k)%alpha(Nsiz))
                     Target_atoms(j)%Ritchi(k)%alpha = Ritchi_temp%alpha
                  endif

                  ! clean up:
                  deallocate(Ritchi_temp%E0)
                  deallocate(Ritchi_temp%A)
                  deallocate(Ritchi_temp%Gamma)
                  if (allocated(Ritchi_temp%alpha)) deallocate(Ritchi_temp%alpha)
               endif
            endif

            ! Define single CDF oscillator:
            if (.not. allocated(Target_atoms(j)%Ritchi(k)%E0)) allocate(Target_atoms(j)%Ritchi(k)%E0(1))
            if (.not. allocated(Target_atoms(j)%Ritchi(k)%A)) allocate(Target_atoms(j)%Ritchi(k)%A(1))
            if (.not. allocated(Target_atoms(j)%Ritchi(k)%Gamma)) allocate(Target_atoms(j)%Ritchi(k)%Gamma(1))
            if (NumPar%kind_of_DR == 4) then
               if (.not. allocated(Target_atoms(j)%Ritchi(k)%alpha)) allocate(Target_atoms(j)%Ritchi(k)%alpha(1))
            endif
            !print*, 'spcdf', j, k, Target_atoms(j)%KOCS_SHI(k)
         enddo ! k
      enddo ! j

   !-----------
   else SP_CDF ! it is a number, read it below again
      NumPar%kind_of_CDF = 0  ! Full user-provided CDF

      backspace(FN2) ! read this line again
      i = i - 1   ! reread this line

      do j = 1, N  ! read for each element its shells data:
         READ(FN2,*,IOSTAT=Reason) Shl ! number of shells in the first atom kind:
      
         call read_file(Reason, i, read_well) ! reports if everything read well
         if (.not. read_well) goto 2014
         ! Now we know how many shells are in this atom:
         Target_atoms(j)%N_shl = Shl   ! save it
         if (.not. allocated(Target_atoms(j)%Shell_name)) allocate(Target_atoms(j)%Shell_name(Shl)) ! allocate shell-names for each shell
         if (.not. allocated(Target_atoms(j)%Shl_num)) allocate(Target_atoms(j)%Shl_num(Shl)) ! allocate shell disignator for each shell
         if (.not. allocated(Target_atoms(j)%Nel)) then
            allocate(Target_atoms(j)%Nel(Shl)) ! allocate numbers of electrons for each shell
            Target_atoms(j)%Nel(:) = 0.0d0
         endif
         if (.not. allocated(Target_atoms(j)%Ip)) then
            allocate(Target_atoms(j)%Ip(Shl)) ! allocate ionization potentials for each shell
            Target_atoms(j)%Ip(:) = -1.0d15
         endif
         if (.not. allocated(Target_atoms(j)%Ek)) then
            allocate(Target_atoms(j)%Ek(Shl)) ! allocate mean kinetic energies for each shell
            Target_atoms(j)%Ek(:) = -1.0d-15
         endif
         if (.not. allocated(Target_atoms(j)%Auger)) then
            allocate(Target_atoms(j)%Auger(Shl)) ! allocate auger-times for each shell
            Target_atoms(j)%Auger(:) = 1.0d31
         endif
         if (.not. allocated(Target_atoms(j)%Radiat)) then
            allocate(Target_atoms(j)%Radiat(Shl)) ! allocate radiative-times for each shell
            Target_atoms(j)%Radiat(:) = 2.0d31
         endif
         if (.not. allocated(Target_atoms(j)%PQN)) then
            allocate(Target_atoms(j)%PQN(Shl)) ! allocate principle quantum numbers for each shell
            Target_atoms(j)%PQN(:) = 0
         endif
         if (.not. allocated(Target_atoms(j)%KOCS)) then
            allocate(Target_atoms(j)%KOCS(Shl)) ! allocate kind of inelastic cross sections
            Target_atoms(j)%KOCS(:) = 0
         endif
         if (.not. allocated(Target_atoms(j)%KOCS_SHI)) then
            allocate(Target_atoms(j)%KOCS_SHI(Shl)) ! allocate kind of inelastic cross sections
            Target_atoms(j)%KOCS_SHI(:) = 0
         endif
         if (.not. allocated(Target_atoms(j)%Ritchi)) then
            allocate(Target_atoms(j)%Ritchi(Shl)) ! allocate Ritchi-functions' coefficiants for each shell
         endif

         do k = 1, Shl ! read all the data for each shell:
            READ(FN2,*,IOSTAT=Reason) CDF_coef, Shl_num, Target_atoms(j)%Ip(k), Target_atoms(j)%Nel(k), Target_atoms(j)%Auger(k) ! number of CDF functions, shell-designator, ionization potential, number of electrons, Auger-time
            call read_file(Reason, i, read_well) ! reports if everything read well
            if (.not. read_well) goto 2014
            Target_atoms(j)%Shl_num(k) = ABS(Shl_num)  ! save shell designator
            call define_PQN(Target_atoms(j)%Shl_num(k), Shell_name, Target_atoms(j)%PQN(k))  ! from module 'Dealing_with_EADL'
            if (Target_atoms(j)%Shl_num(k) .GE. 63) Matter%N_VB_el = Target_atoms(j)%Nel(k) ! then it is the valence band, save number of VB electrons
            Target_atoms(j)%Shell_name(k) = Shell_name    ! save the name of this shell

            ! If some data are missing in the input file, get it from EADL database:
            call check_atomic_parameters(NumPar, Target_atoms, j, k, shl, Error_message, read_well, MPI_param) ! from module 'Dealing_with_EADL'

            if (Target_atoms(j)%Ip(k) .LT. 1.0d-1) Target_atoms(j)%Ip(k) = 1.0d-1 ! [eV], introduce at least a minimum "band gap"
            if (CDF_coef .GT. 0) then ! there is CDF to be used in a cross section
             Target_atoms(j)%KOCS_SHI(k) = 1 ! CDF cross section for this shell
             if (Shl_num .GT. 0) then ! use normal CDF cross sections for electrons:
                Target_atoms(j)%KOCS(k) = 1 ! CDF cross section for this shell
             else ! use BEB for electrons (but still CDF for the SHI):
                Target_atoms(j)%KOCS(k) = 2 ! BEB cross section for this shell
             endif
             if (.not. allocated(Target_atoms(j)%Ritchi(k)%E0)) allocate(Target_atoms(j)%Ritchi(k)%E0(CDF_coef))   ! that's how many CDF function for this shell of this atom
             if (.not. allocated(Target_atoms(j)%Ritchi(k)%A)) allocate(Target_atoms(j)%Ritchi(k)%A(CDF_coef))   ! that's how many CDF function for this shell of this atom
             if (.not. allocated(Target_atoms(j)%Ritchi(k)%Gamma)) allocate(Target_atoms(j)%Ritchi(k)%Gamma(CDF_coef))   ! that's how many CDF function for this shell of this atom
             if (NumPar%kind_of_DR .EQ. 4) then
                 if (.not. allocated(Target_atoms(j)%Ritchi(k)%alpha)) allocate(Target_atoms(j)%Ritchi(k)%alpha(CDF_coef))   ! that's how many CDF functions
             endif
             do l = 1, CDF_coef    ! read all the CDF parameters:
                READ(FN2,*,IOSTAT=Reason) Target_atoms(j)%Ritchi(k)%E0(l), Target_atoms(j)%Ritchi(k)%A(l), Target_atoms(j)%Ritchi(k)%Gamma(l)
                call read_file(Reason, i, read_well) ! reports if everything read well
                if (.not. read_well) goto 2014
                ! For Delta-CDF:
                if (NumPar%kind_of_DR .EQ. 4) then
                   call define_alpha(Target_atoms(j)%Ritchi(k)%A(l), Target_atoms(j)%Ritchi(k)%Gamma(l), Target_atoms(j)%Ritchi(k)%E0(l), Target_atoms(j)%Ip(k), Target_atoms(j)%Ritchi(k)%alpha(l))  ! below
                endif
             enddo
            else ! there is no CDF, so let's use atomic BEB cross section:
             Target_atoms(j)%KOCS(k) = 2 ! BEB cross section for this shell
             Target_atoms(j)%KOCS_SHI(k) = 2 ! BEB cross section for this shell
            endif
         enddo ! k = 1, Shl
      enddo ! j = 1, N
   
      READ(FN2,*,IOSTAT=Reason) CDF_coef
      i = i + 1
      IF (Reason .GT. 0)  THEN ! ... something wrong ...
         write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', wrong type of variable'
         read_well = .false.
         NumPar%kind_of_CDF_ph = 1 ! use single-pole approximation ofr phonon CDf
      ELSE IF (Reason .LT. 0) THEN ! ... end of file reached ...
         write(*,'(a,i3,a)') 'Line ', i, ', END of file; no CDF data for phonon-peak'
         write(*,'(a,i3,a)') 'Using atomic cross-sections for elastic scattering'
         NumPar%kind_of_CDF_ph = 1 ! use single-pole approximation of phonon CDF
         read_well = .true.
      ELSE   ! normal reading
         read_well = .true.  ! it read well, nothing to report
         if (.not. allocated(CDF_Phonon%A)) allocate(CDF_Phonon%A(CDF_coef))
         if (.not. allocated(CDF_Phonon%E0)) allocate(CDF_Phonon%E0(CDF_coef))
         if (.not. allocated(CDF_Phonon%Gamma)) allocate(CDF_Phonon%Gamma(CDF_coef))
         if (NumPar%kind_of_DR .EQ. 4) then
            if (.not. allocated(CDF_Phonon%alpha)) allocate(CDF_Phonon%alpha(CDF_coef))   ! that's how many CDF functions
         endif
         do l = 1, CDF_coef    ! read all the CDF parameters for phonon peak:
            READ(FN2,*,IOSTAT=Reason) CDF_Phonon%E0(l), CDF_Phonon%A(l), CDF_Phonon%Gamma(l)
            call read_file(Reason, i, read_well) ! reports if everything read well
            if (.not. read_well) goto 2014
            ! For Delta-CDF:
            if (NumPar%kind_of_DR .EQ. 4) then
              call define_alpha(CDF_Phonon%A(l), CDF_Phonon%Gamma(l), CDF_Phonon%E0(l), 0.0d0, CDF_Phonon%alpha(l))  ! below
            endif
         enddo
         NumPar%kind_of_CDF_ph = 0 ! user provided phonon CDF
      END IF
   endif SP_CDF

   ! Check if phonon CDF requires single-pole approximation:
   if (NumPar%kind_of_CDF_ph == 1) then   ! allocate phonon CDF
      if (.not. allocated(CDF_Phonon%A)) allocate(CDF_Phonon%A(1), source = 0.0d0)
      if (.not. allocated(CDF_Phonon%E0)) allocate(CDF_Phonon%E0(1), source = 0.0d0)
      if (.not. allocated(CDF_Phonon%Gamma)) allocate(CDF_Phonon%Gamma(1), source = 0.0d0)
   endif


   ! Allocate differential cross section arrays if required:
   select case (NumPar%CS_method)
   case (1)
      if (.not. allocated(aidCS%EIdCS)) then
         allocate(aidCS%EIdCS(size(Target_atoms))) ! for each atom type
      endif
      ! Each shell:
      do j = 1, N  ! read for each element its shells data:
         if (.not.allocated(aidCS%EIdCS(j)%Int_diff_CS)) then
            allocate(aidCS%EIdCS(j)%Int_diff_CS( size(Target_atoms(j)%Ip) )) ! allocate table of integral of differential cross secion for each shell
         endif
      enddo
   end select

2014 continue   ! if we must skip to the end for some reason
   if (.not. read_well) print*, trim(adjustl(Material_file))
   ! Independent file opens, no need for MPI-specific file-close:
   inquire(unit=FN2,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN2)             ! and if it is, close it
   !call MPI_fileclose_wrapper(MPI_param, FN2, Error_message)   ! module "MPI_subroutines"
end subroutine reading_material_parameters



subroutine broadcast_material_parameters(Target_atoms, NumPar, CDF_Phonon, Matter, aidCS, MPI_param) ! MPI-related routine
   type(Atom), dimension(:), allocatable, intent(inout) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   type(Flag), intent(inout) :: NumPar ! numerical parameters
   type(CDF), intent(inout) :: CDF_Phonon   ! CDF parameters for phonon to be read from a file
   type(Solid), intent(inout) :: Matter   ! all material parameters
   type(All_diff_CS), intent(inout) :: aidCS    ! all integrated differential cross sections
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   !------------------------------------------
   integer :: Nsiz, i, j, k, l, Nat, N1, N2, N3, N4
   character(100) :: error_part
   logical :: do_broadcast, do_broadcast2, do_broadcast1

#ifdef MPI_USED
   ! initialize:
   do_broadcast = .false.
   do_broadcast1 = .false.
   do_broadcast2 = .false.

   error_part = 'ERROR in broadcast_material_parameters'

   Nsiz = LEN(NumPar%CDF_file)
   call mpi_bcast(NumPar%CDF_file, Nsiz, MPI_CHARACTER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {CDF_file}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%VB_CDF_defined, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {VB_CDF_defined}') ! module "MPI_subroutines"

   Nsiz = LEN(Matter%Target_name)
   call mpi_bcast(Matter%Target_name, Nsiz, MPI_CHARACTER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_name}') ! module "MPI_subroutines"

   Nsiz = LEN(Matter%Chem)
   call mpi_bcast(Matter%Chem, Nsiz, MPI_CHARACTER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Chem}') ! module "MPI_subroutines"

   if (MPI_param%process_rank == 0) then
      Nat = size(Target_atoms)   ! define number of atoms in the compound
   else
      Nat = 0  ! unknown yet, get from the master process
   endif
   call mpi_bcast(Nat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Nat}') ! module "MPI_subroutines"
   if (.not.allocated(Target_atoms)) allocate(Target_atoms(Nat))

   !------------------------------------------------------
   ! Synchronize MPI processes: make sure all processes know the size of the Target_atoms and allocated it
   call MPI_barrier_wrapper(MPI_param)  ! module "MPI_subroutines"
   !------------------------------------------------------
   ! Now we know the numer of atoms in the compound, define the parameters of each:
   do i = 1, Nat
      call mpi_bcast(Target_atoms(i)%Zat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Zat}') ! module "MPI_subroutines"

      call mpi_bcast(Target_atoms(i)%Pers, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Pers}') ! module "MPI_subroutines"

      call mpi_bcast(Target_atoms(i)%Mass, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Mass}') ! module "MPI_subroutines"

      Nsiz = LEN(Target_atoms(i)%Name)
      call mpi_bcast(Target_atoms(i)%Name, Nsiz, MPI_CHARACTER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Name}') ! module "MPI_subroutines"

      Nsiz = LEN(Target_atoms(i)%Full_Name)
      call mpi_bcast(Target_atoms(i)%Full_Name, Nsiz, MPI_CHARACTER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Full_Name}') ! module "MPI_subroutines"
      !print*, '[Process #', MPI_param%process_rank, '] B:', i, Target_atoms(i)%Zat, Target_atoms(i)%Pers, Target_atoms(i)%Mass, Target_atoms(i)%Name, Target_atoms(i)%Full_Name

   enddo

   call mpi_bcast(Matter%Dens, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Dens}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%At_Dens, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {At_Dens}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%v_f, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {v_f}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%Vsound, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Vsound}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%E_F, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {E_F}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%Egap, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Egap}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%N_VB_el, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {N_VB_el}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%Layer, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Layer}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%work_function, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {work_function}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%bar_length, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {bar_length}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%bar_height, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {bar_height}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%hole_mass, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {hole_mass}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%El_eff_mass, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {El_eff_mass}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%v_f, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {v_f}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%cut_off, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {cut_off}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%temp, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {temp}') ! module "MPI_subroutines"

   call mpi_bcast(Matter%cut_off, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {cut_off}') ! module "MPI_subroutines"

   ! That's how to deal with an array, if it is unknown whether it is allocated or not:
   do_broadcast = allocated(Matter%form_factor) ! to chech if we need to do the broadcas or not
   call mpi_bcast(do_broadcast, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {do_broadcast #1}') ! module "MPI_subroutines"
   !------------------------------------------------------
   ! Synchronize MPI processes: make sure all processes know the size of the form_factor and allocated it
   call MPI_barrier_wrapper(MPI_param)  ! module "MPI_subroutines"
   !------------------------------------------------------
   if (do_broadcast) then  ! broadcast the form-factors, only if they were allocated
      if (MPI_param%process_rank == 0) then
         N1 = size(Matter%form_factor,1)
         N2 = size(Matter%form_factor,2)
      else
         N1 = 0
         N2 = 0
      endif
      call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {N1}') ! module "MPI_subroutines"
      call mpi_bcast(N2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {N2}') ! module "MPI_subroutines"
      !------------------------------------------------------
      ! Synchronize MPI processes: make sure all processes know the size of the form_factor and allocated it
      call MPI_barrier_wrapper(MPI_param)  ! module "MPI_subroutines"
      !------------------------------------------------------
      if (MPI_param%process_rank /= 0) then  ! knowing the size, allocate the arrays
         if (.not.allocated(Matter%form_factor)) allocate(Matter%form_factor(N1,N2))
      endif
      call mpi_bcast(Matter%form_factor, N1*N2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {form_factor}') ! module "MPI_subroutines"
      !print*, '[Process #', MPI_param%process_rank, '] N:', size(Matter%form_factor,1), size(Matter%form_factor,2), Matter%form_factor(2,3)
   endif

   call mpi_bcast(NumPar%kind_of_CDF, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {kind_of_CDF}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%kind_of_CDF_ph, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {kind_of_CDF_ph}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%VB_CDF_defined, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {VB_CDF_defined}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%phonon_CDF_defined, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {phonon_CDF_defined}') ! module "MPI_subroutines"

      ! Now we know the numer of atoms in the compound, define the parameters of each:
   do i = 1, Nat
      call mpi_bcast(Target_atoms(i)%N_shl, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%N_shl}') ! module "MPI_subroutines"
      !------------------------------------------------------
      ! Synchronize MPI processes: make sure all processes know the size of the form_factor and allocated it
      call MPI_barrier_wrapper(MPI_param)  ! module "MPI_subroutines"
      !------------------------------------------------------
      ! Now we can allocate the target-atoms shells-related arrays:
      if (MPI_param%process_rank /= 0) then  ! knowing the size, allocate the arrays
         if (.not.allocated(Target_atoms(i)%Shell_name)) allocate(Target_atoms(i)%Shell_name(Target_atoms(i)%N_shl))
         if (.not.allocated(Target_atoms(i)%Shl_num)) allocate(Target_atoms(i)%Shl_num(Target_atoms(i)%N_shl))
         if (.not.allocated(Target_atoms(i)%Nel)) allocate(Target_atoms(i)%Nel(Target_atoms(i)%N_shl))
         if (.not.allocated(Target_atoms(i)%Ip)) allocate(Target_atoms(i)%Ip(Target_atoms(i)%N_shl))
         if (.not.allocated(Target_atoms(i)%Ek)) allocate(Target_atoms(i)%Ek(Target_atoms(i)%N_shl))
         if (.not.allocated(Target_atoms(i)%Auger)) allocate(Target_atoms(i)%Auger(Target_atoms(i)%N_shl))
         if (.not.allocated(Target_atoms(i)%Radiat)) allocate(Target_atoms(i)%Radiat(Target_atoms(i)%N_shl))
         if (.not.allocated(Target_atoms(i)%Ritchi)) allocate(Target_atoms(i)%Ritchi(Target_atoms(i)%N_shl))
         if (.not.allocated(Target_atoms(i)%PQN)) allocate(Target_atoms(i)%PQN(Target_atoms(i)%N_shl))
         if (.not.allocated(Target_atoms(i)%KOCS)) allocate(Target_atoms(i)%KOCS(Target_atoms(i)%N_shl))
         if (.not.allocated(Target_atoms(i)%KOCS_SHI)) allocate(Target_atoms(i)%KOCS_SHI(Target_atoms(i)%N_shl))
      endif
      ! Now, broadcast the values into these arrays:
      Nsiz = size(Target_atoms(i)%Shell_name)
      N1 = LEN(Target_atoms(i)%Shell_name(1))
      call mpi_bcast(Target_atoms(i)%Shell_name, Nsiz*N1, MPI_CHARACTER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Shell_name}') ! module "MPI_subroutines"

      call mpi_bcast(Target_atoms(i)%Shl_num, Nsiz, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Shl_num}') ! module "MPI_subroutines"

      call mpi_bcast(Target_atoms(i)%Nel, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Nel}') ! module "MPI_subroutines"

      call mpi_bcast(Target_atoms(i)%Ip, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Ip}') ! module "MPI_subroutines"

      call mpi_bcast(Target_atoms(i)%Ek, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Ek}') ! module "MPI_subroutines"

      call mpi_bcast(Target_atoms(i)%Auger, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Auger}') ! module "MPI_subroutines"

      do_broadcast = allocated(Target_atoms(i)%Radiat) ! to chech if we need to do the broadcas or not
      call mpi_bcast(do_broadcast, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {do_broadcast #2}') ! module "MPI_subroutines"
      if (do_broadcast) then
         call mpi_bcast(Target_atoms(i)%Radiat, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Radiat}') ! module "MPI_subroutines"
      endif

      call mpi_bcast(Target_atoms(i)%PQN, Nsiz, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%PQN}') ! module "MPI_subroutines"

      call mpi_bcast(Target_atoms(i)%KOCS, Nsiz, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%KOCS}') ! module "MPI_subroutines"

      call mpi_bcast(Target_atoms(i)%KOCS_SHI, Nsiz, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%KOCS_SHI}') ! module "MPI_subroutines"

      do_broadcast1 = allocated(Target_atoms(i)%Ritchi)
      call mpi_bcast(do_broadcast1, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {do_broadcast1 #2.5}') ! module "MPI_subroutines"

      if (do_broadcast1) then
       do j = 1, size(Target_atoms(i)%Ritchi) ! all CDF oscillators for each shell
         ! Check if Ritchi oscillators are allocated:
         do_broadcast = allocated(Target_atoms(i)%Ritchi(j)%E0) ! to chech if we need to do the broadcas or not
         call mpi_bcast(do_broadcast, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {do_broadcast #3}') ! module "MPI_subroutines"

         if (do_broadcast) then
            N1 = size(Target_atoms(i)%Ritchi(j)%E0)
            call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
            call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Ritchi(j)%E0:N1}') ! module "MPI_subroutines"
            if (MPI_param%process_rank /= 0) then  ! knowing the size, allocate the arrays
               if (.not.allocated(Target_atoms(i)%Ritchi(j)%E0)) allocate(Target_atoms(i)%Ritchi(j)%E0(N1))
               if (.not.allocated(Target_atoms(i)%Ritchi(j)%A)) allocate(Target_atoms(i)%Ritchi(j)%A(N1))
               if (.not.allocated(Target_atoms(i)%Ritchi(j)%Gamma)) allocate(Target_atoms(i)%Ritchi(j)%Gamma(N1))
               if (.not.allocated(Target_atoms(i)%Ritchi(j)%alpha)) allocate(Target_atoms(i)%Ritchi(j)%alpha(N1))
            endif
            !------------------------------------------------------
            ! Synchronize MPI processes: make sure all processes know the size of the form_factor and allocated it
            call MPI_barrier_wrapper(MPI_param)  ! module "MPI_subroutines"
            !------------------------------------------------------
            call mpi_bcast(Target_atoms(i)%Ritchi(j)%E0, N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
            call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Ritchi(j)%E0}') ! module "MPI_subroutines"

            call mpi_bcast(Target_atoms(i)%Ritchi(j)%A, N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
            call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Ritchi(j)%A}') ! module "MPI_subroutines"

            call mpi_bcast(Target_atoms(i)%Ritchi(j)%Gamma, N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
            call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Ritchi(j)%Gamma}') ! module "MPI_subroutines"

            do_broadcast2 = allocated(Target_atoms(i)%Ritchi(j)%alpha) ! to chech if we need to do the broadcas or not
            call mpi_bcast(do_broadcast2, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
            call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {do_broadcast2 #4}') ! module "MPI_subroutines"
            if (do_broadcast2) then
               call mpi_bcast(Target_atoms(i)%Ritchi(j)%alpha, N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
               call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Target_atoms(i)%Ritchi(j)%alpha}') ! module "MPI_subroutines"
            endif
         endif ! do_broadcast
         !print*, '[Process #', MPI_param%process_rank, '] ShN:', i, j, N1, size(Target_atoms(i)%Ritchi(j)%E0), Target_atoms(i)%Ritchi(j)%E0(:)
       enddo ! j
      endif ! do_broadcast1
   enddo ! i

   !------------------------------------------------------
   ! Synchronize MPI processes
   call MPI_barrier_wrapper(MPI_param)  ! module "MPI_subroutines"
   !------------------------------------------------------
   do_broadcast = allocated(CDF_Phonon%A)
   call mpi_bcast(do_broadcast, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {do_broadcast #5}') ! module "MPI_subroutines"
   if (do_broadcast) then
      Nsiz = size(CDF_Phonon%A)
      call mpi_bcast(Nsiz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {CDF_Phonon:Nsiz}') ! module "MPI_subroutines"

      if (MPI_param%process_rank /= 0) then  ! knowing the size, allocate the arrays
         if (.not. allocated(CDF_Phonon%A)) allocate(CDF_Phonon%A(Nsiz))
         if (.not. allocated(CDF_Phonon%E0)) allocate(CDF_Phonon%E0(Nsiz))
         if (.not. allocated(CDF_Phonon%Gamma)) allocate(CDF_Phonon%Gamma(Nsiz))
      endif

      call mpi_bcast(CDF_Phonon%E0, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {CDF_Phonon%E0}') ! module "MPI_subroutines"

      call mpi_bcast(CDF_Phonon%A, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {CDF_Phonon%A}') ! module "MPI_subroutines"

      call mpi_bcast(CDF_Phonon%Gamma, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {CDF_Phonon%Gamma}') ! module "MPI_subroutines"

      do_broadcast2 = allocated(CDF_Phonon%alpha) ! to chech if we need to do the broadcas or not
      call mpi_bcast(do_broadcast2, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {do_broadcast2 #6}') ! module "MPI_subroutines"
      if (do_broadcast2) then
         if (MPI_param%process_rank /= 0) then  ! knowing the size, allocate the arrays
            if (.not. allocated(CDF_Phonon%alpha)) allocate(CDF_Phonon%alpha(Nsiz))
         endif
         call mpi_bcast(CDF_Phonon%alpha, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {CDF_Phonon%alpha}') ! module "MPI_subroutines"
      endif ! do_broadcast2
      !print*, '[Process #', MPI_param%process_rank, '] Ph:', Nsiz, size(CDF_Phonon%E0), CDF_Phonon%E0(:)
   endif ! do_broadcast

   ! Allocate differential cross section arrays if required:
   call mpi_bcast(NumPar%CS_method, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {CS_method}') ! module "MPI_subroutines"

   select case (NumPar%CS_method)
   case (1)
      do_broadcast = allocated(aidCS%EIdCS)
      call mpi_bcast(do_broadcast, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {do_broadcast #7}') ! module "MPI_subroutines"
      if (do_broadcast) then
         Nsiz = size(Target_atoms)
         if (MPI_param%process_rank /= 0) then  ! knowing the size, allocate the arrays
            if (.not. allocated(aidCS%EIdCS)) allocate(aidCS%EIdCS(Nsiz)) ! for each atom type
            ! Each shell:
            do j = 1, Nsiz  ! read for each element its shells data:
               if (.not.allocated(aidCS%EIdCS(j)%Int_diff_CS)) then
                  allocate(aidCS%EIdCS(j)%Int_diff_CS( size(Target_atoms(j)%Ip) )) ! allocate table of integral of differential cross secion for each shell
               endif
            enddo
         endif
      endif
   end select

   !pause 'broadcast_material_parameters'
#endif
end subroutine broadcast_material_parameters



subroutine make_valence_band(Target_atoms, NumPar, Matter, Error_message, read_well)
   type(Atom), dimension(:), allocatable, intent(inout) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   type(Flag), intent(inout) :: NumPar ! numerical parameters
   type(Solid), intent(inout) :: Matter   ! all material parameters
   type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
   logical, intent(inout) :: read_well  ! did we read the file well?
   !------------------------
   type(Atom), dimension(:), allocatable :: Atoms_temp
   type(atomic_data), dimension(:), allocatable :: Periodic_table ! this is an internal module variable
   character(200) :: Folder_name, File_name
   character(10) :: temp
   real(8) :: N_el, N_e_VB
   integer :: FN, N, Reason, i, j, N_at, N_VB_siz, Nsiz
   logical :: file_exists, file_opened
   logical, dimension(:), allocatable :: valent

   !-----------
   ! 1) Read periodic table:
   read_well = .true.
   Folder_name = trim(adjustl(m_atomic_folder)) ! module "Dealing_with_EADL"
   File_name = trim(adjustl(Folder_name))//trim(adjustl(NumPar%path_sep))//trim(adjustl(m_atomic_data_file)) ! module "Dealing_with_EADL"
   inquire(file=trim(adjustl(File_name)),exist=file_exists) ! check if input file is there
   exists:if (file_exists) then
      open (newunit=FN, file=trim(adjustl(File_name)))
      call Count_lines_in_file(FN, N)
      if (.not.allocated(Periodic_table)) allocate(Periodic_table(N-1))
      read(FN,*,IOSTAT=Reason) ! skip first line with names
      do i = 1, N-1
         read(FN,*,IOSTAT=Reason) Periodic_table(i)%Z, Periodic_table(i)%Full_name, Periodic_table(i)%Name, &
                                    Periodic_table(i)%Mass, Periodic_table(i)%Nvb
         !print*, Periodic_table(i)%Z, Periodic_table(i)%Full_name, Periodic_table(i)%Name, Periodic_table(i)%Mass, Periodic_table(i)%Nvb
         if (Reason /= 0) exit
      enddo
      if (Reason /= 0) then
         write(temp,'(i4)') i
         Error_message%Err_descript = 'Error in make_valence_band: could not read line '//trim(adjustl(temp))//' in file '//trim(adjustl(File_name))
         goto 9999
      endif
   else exists
      Error_message%Err_descript = 'Error in make_valence_band: could not find file '//trim(adjustl(File_name))
      goto 9999
   endif exists

   !-----------
   ! 2) Temporarily save atomic data:
   call copy_atomic_data(Target_atoms, Atoms_temp) ! below

   !-----------
   ! 3) Define which levels are valent:
   N_e_VB = 0.0d0 ! total number of valence electrons
   N_at = size(Target_atoms)  ! how many elements in the target
   do i = 1, N_at
      N_el = 0.0d0   ! to start counting
      j = size(Target_atoms(i)%Nel)
      do while (N_el < Periodic_table( Target_atoms(i)%Zat )%Nvb)
         Atoms_temp(i)%Ek(j) = -2.0d0  ! mark it for valent shells
         N_el = N_el + Target_atoms(i)%Nel(j)
         j = j - 1
      enddo
      N_e_VB = N_e_VB + N_el*Target_atoms(i)%Pers
   enddo

   ! Save number of VB electrons:
   Matter%N_VB_el = N_e_VB

   !-----------
   ! 4) Combine valence levels into valence band:
   do i = 1, N_at
      ! Redefine the number of core shells
      N_VB_siz = COUNT (Atoms_temp(i)%Ek < -1.0d0) ! number of valence shells
      Nsiz = size(Atoms_temp(i)%Ek) - N_VB_siz     ! number of core shells

      ! Redefine the core shell and valence band:
      call copy_atomic_data_back(Target_atoms, Atoms_temp, Matter%Egap, N_e_VB, i, Nsiz) ! below
   enddo

9999   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it
end subroutine make_valence_band



subroutine copy_atomic_data_back(Target_atoms, Atoms_temp, Egap, N_e_VB, i, Nsiz) ! below
   type(Atom), dimension(:), intent(inout) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   type(Atom), dimension(:), allocatable, intent(in) :: Atoms_temp  ! define target atoms as objects, we don't know yet how many they are
   real(8), intent(in) :: Egap, N_e_VB ! bandgap [eV]; number of VB electrons per molecule
   integer, intent(in) :: i, Nsiz
   !------------
   integer :: j, Shl, N_temp
   real(8) :: N_VB, N_at_tot

   if (i == 1) then
      Shl = Nsiz + 1 ! core shells + valence band
   else
      Shl = Nsiz  ! core shells only (VB was accounted for in the first element)
   endif

   Target_atoms(i)%N_shl = Shl

   deallocate(Target_atoms(i)%Shell_name)
   deallocate(Target_atoms(i)%Shl_num)
   deallocate(Target_atoms(i)%Nel)
   deallocate(Target_atoms(i)%Ip)
   deallocate(Target_atoms(i)%Ek)
   deallocate(Target_atoms(i)%Auger)
   deallocate(Target_atoms(i)%Radiat)
   deallocate(Target_atoms(i)%PQN)
   deallocate(Target_atoms(i)%KOCS)
   deallocate(Target_atoms(i)%KOCS_SHI)
   if (allocated(Target_atoms(i)%Ritchi)) deallocate(Target_atoms(i)%Ritchi)

   allocate(Target_atoms(i)%Shell_name(Shl)) ! allocate shell-names for each shell
   allocate(Target_atoms(i)%Shl_num(Shl)) ! allocate shell disignator for each shell
   allocate(Target_atoms(i)%Nel(Shl)) ! allocate numbers of electrons for each shell
   allocate(Target_atoms(i)%Ip(Shl)) ! allocate ionization potentials for each shell
   allocate(Target_atoms(i)%Ek(Shl)) ! allocate mean kinetic energies for each shell
   allocate(Target_atoms(i)%Auger(Shl)) ! allocate auger-times for each shell
   allocate(Target_atoms(i)%Radiat(Shl)) ! allocate radiative-times for each shell
   allocate(Target_atoms(i)%PQN(Shl)) ! allocate principle quantum numbers for each shell
   allocate(Target_atoms(i)%KOCS(Shl)) ! allocate kind of inelastic cross sections
   allocate(Target_atoms(i)%KOCS_SHI(Shl)) ! allocate kind of inelastic cross sections
   allocate(Target_atoms(i)%Ritchi(Shl)) ! allocate Ritchi-functions' coefficiants for each shell

   Target_atoms(i)%Shell_name(1:Shl) = Atoms_temp(i)%Shell_name(1:Shl)
   Target_atoms(i)%Shl_num(1:Shl) = Atoms_temp(i)%Shl_num(1:Shl)
   Target_atoms(i)%Nel(1:Shl) = Atoms_temp(i)%Nel(1:Shl)
   Target_atoms(i)%Ip(1:Shl) = Atoms_temp(i)%Ip(1:Shl)
   Target_atoms(i)%Ek(1:Shl) = Atoms_temp(i)%Ek(1:Shl)
   Target_atoms(i)%Auger(1:Shl) = Atoms_temp(i)%Auger(1:Shl)
   Target_atoms(i)%Radiat(1:Shl) = Atoms_temp(i)%Radiat(1:Shl)
   Target_atoms(i)%PQN(1:Shl) = Atoms_temp(i)%PQN(1:Shl)
   Target_atoms(i)%KOCS(1:Shl) = Atoms_temp(i)%KOCS(1:Shl)
   Target_atoms(i)%KOCS_SHI(1:Shl) = Atoms_temp(i)%KOCS_SHI(1:Shl)
   !Target_atoms(i)%Ritchi(1:Shl) = Atoms_temp(i)%Ritchi(1:Shl)   ! undefined yet

   if (i == 1) then  ! there is valence band, define it:
      Target_atoms(i)%Shell_name(Shl) = 'Valence'
      Target_atoms(i)%Shl_num(Shl) = 63
      Target_atoms(i)%Nel(Shl) = N_e_VB
      Target_atoms(i)%Ip(Shl) = Egap
      Target_atoms(i)%Ek(Shl) = 0.0d0
      Target_atoms(i)%Auger(Shl) = 1.0d26
      Target_atoms(i)%Radiat(Shl) = 1.0d27
      ! The rest we don't need to change:
      !Target_atoms(i)%PQN(Shl)
      !Target_atoms(i)%KOCS(Shl)
      !Target_atoms(i)%KOCS_SHI(Shl)
      !Target_atoms(i)%Ritchi(Shl)
   endif
end subroutine copy_atomic_data_back




subroutine copy_atomic_data(Target_atoms, Atoms_temp) ! below
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   type(Atom), dimension(:), allocatable, intent(inout) :: Atoms_temp  ! define target atoms as objects, we don't know yet how many they are
   !------------
   integer :: i, j, k, N_at, Shl, Nj
   N_at = size(Target_atoms)

   allocate(Atoms_temp(N_at))
   do i = 1, N_at ! all atoms
      Atoms_temp(i)%Zat = Target_atoms(i)%Zat
      Atoms_temp(i)%Mass = Target_atoms(i)%Mass
      Atoms_temp(i)%Pers = Target_atoms(i)%Pers
      Atoms_temp(i)%Name = Target_atoms(i)%Name
      Atoms_temp(i)%Full_Name = Target_atoms(i)%Full_Name
      Atoms_temp(i)%N_shl = Target_atoms(i)%N_shl
      Shl = Atoms_temp(i)%N_shl
      allocate(Atoms_temp(i)%Shell_name(Shl)) ! allocate shell-names for each shell
      allocate(Atoms_temp(i)%Shl_num(Shl)) ! allocate shell disignator for each shell
      allocate(Atoms_temp(i)%Nel(Shl)) ! allocate numbers of electrons for each shell
      allocate(Atoms_temp(i)%Ip(Shl)) ! allocate ionization potentials for each shell
      allocate(Atoms_temp(i)%Ek(Shl)) ! allocate mean kinetic energies for each shell
      allocate(Atoms_temp(i)%Auger(Shl)) ! allocate auger-times for each shell
      allocate(Atoms_temp(i)%Radiat(Shl)) ! allocate radiative-times for each shell
      allocate(Atoms_temp(i)%PQN(Shl)) ! allocate principle quantum numbers for each shell
      allocate(Atoms_temp(i)%KOCS(Shl)) ! allocate kind of inelastic cross sections
      allocate(Atoms_temp(i)%KOCS_SHI(Shl)) ! allocate kind of inelastic cross sections
      allocate(Atoms_temp(i)%Ritchi(Shl)) ! allocate Ritchi-functions' coefficiants for each shell

      Atoms_temp(i)%Shell_name = Target_atoms(i)%Shell_name
      Atoms_temp(i)%Shl_num = Target_atoms(i)%Shl_num
      Atoms_temp(i)%Nel = Target_atoms(i)%Nel
      Atoms_temp(i)%Ip = Target_atoms(i)%Ip
      Atoms_temp(i)%Ek = Target_atoms(i)%Ek
      Atoms_temp(i)%Auger = Target_atoms(i)%Auger
      Atoms_temp(i)%Radiat = Target_atoms(i)%Radiat
      Atoms_temp(i)%PQN = Target_atoms(i)%PQN
      Atoms_temp(i)%KOCS = Target_atoms(i)%KOCS
      Atoms_temp(i)%KOCS_SHI = Target_atoms(i)%KOCS_SHI
      if (allocated(Target_atoms(i)%Ritchi)) then
         do j = 1, size(Atoms_temp(i)%Ritchi)
            Nj = size(Target_atoms(i)%Ritchi(j)%A)
            allocate(Atoms_temp(i)%Ritchi(j)%A(Nj))
            allocate(Atoms_temp(i)%Ritchi(j)%E0(Nj))
            allocate(Atoms_temp(i)%Ritchi(j)%Gamma(Nj))
            allocate(Atoms_temp(i)%Ritchi(j)%alpha(Nj))
            Atoms_temp(i)%Ritchi(j)%A = Target_atoms(i)%Ritchi(j)%A
            Atoms_temp(i)%Ritchi(j)%E0 = Target_atoms(i)%Ritchi(j)%E0
            Atoms_temp(i)%Ritchi(j)%Gamma = Target_atoms(i)%Ritchi(j)%Gamma
            Atoms_temp(i)%Ritchi(j)%alpha = Target_atoms(i)%Ritchi(j)%alpha
         enddo
      endif
   enddo ! i
end subroutine copy_atomic_data



pure subroutine define_alpha(Ai, Gammai, E0i, x_min, alpha)
   real(8), intent(in) :: Ai, Gammai, E0i, x_min   ! [eV] the coefficients of Ritchi's CDF, and x_min=Ip or Egap, the starting point
   real(8), intent(out) :: alpha ! coefficient in delta-CDF to be fitted
   real(8) :: low_lim, high_lim
   low_lim =  Int_Ritchi_x(Ai,E0i,Gammai,x_min)     ! below
   high_lim = Int_Ritchi_x(Ai,E0i,Gammai,1.0d30)   ! below
   alpha = high_lim - low_lim   ! normalized to reproduce K-sum rule
end subroutine define_alpha

pure function Int_Ritchi_x(A,E,Gamma,x) ! analytical integral of the Ritchi*x (k-sum rule), analytical integral of Eq.(8) [1]
    real(8), intent(in) :: A, E, Gamma, x  ! parameters and variable
    real(8) Int_Ritchi_x ! function itself
    real(8) S, sq2, G, s_plus, s_minus,B, g2
    complex(8) Sc, Gc, s_plus_c, s_minus_c, Bc, Ic, arg,arg2, oneI, GammaSc, Ic1
    sq2 = sqrt(2.0d0)
!     oneI = cmplx(0.0d0,1.0d0)
    oneI = g_CI
    g2 = Gamma*Gamma
    if ((-(2.0d0*E)*(2.0d0*E) + g2) .GE. 0.0d0) then
        Sc = cmplx(sqrt(-(2.0d0*E)*(2.0d0*E) + g2),0.0d0)
    else
        Sc = cmplx(0.0d0, sqrt(ABS(-(2.0d0*E)*(2.0d0*E) + g2)))
    endif
    Gc = cmplx(g2 - 2.0d0*E*E,0.0d0)
    GammaSc = Gamma*Sc
    s_plus_c = sq2*sqrt(Gc + GammaSc)
    s_minus_c = sq2*sqrt(Gc - GammaSc)
    Bc=Gc/(GammaSc)
    !Int_Ritchi_x = real(A*Gamma*( ATAN(2.0e0*x/s_minus)/s_minus*(1.0e0-B) + ATAN(2.0e0*x/s_plus)/s_plus*(1.0e0+B) ))
    if (s_minus_c == cmplx(0.0d0,0.0d0)) then
       Ic1 = cmplx(0.0d0,0.0d0)
    else
       arg = 2.0d0*x/(s_minus_c)
       Ic1 = 0.5d0*oneI*(log(1.0d0-oneI*arg)-log(1.0d0+oneI*arg))/s_minus_c*(1.0d0-Bc)
    endif
    arg2 = 2.0d0*x/(s_plus_c)
    !Ic = (A*Gamma*( ATAN2(aimag(arg),real(arg))/s_minus_c*(1.0e0-Bc) + ATAN2(aimag(arg2),real(arg2))/s_plus_c*(1.0e0+Bc) ))
    Ic = Ic1 + (0.5d0*oneI*(log(1.0d0-oneI*arg2)-log(1.0d0+oneI*arg2)))/s_plus_c*(1.0d0+Bc)
    Ic = A*Gamma*Ic
    Int_Ritchi_x = dble(Ic)
end function Int_Ritchi_x




subroutine read_short_scdf(FN2, Target_atoms, NumPar, CDF_Phonon, Matter, Error_message, read_well, MPI_param)
   integer, intent(in) :: FN2 ! file number to read from
   type(Atom), dimension(:), allocatable, intent(inout) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   type(Flag), intent(inout) :: NumPar ! numerical parameters
   type(CDF), intent(inout) :: CDF_Phonon   ! CDF parameters for phonon to be read from a file
   type(Solid), intent(inout) :: Matter   ! all material parameters
   type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
   logical, intent(inout) :: read_well  ! did we read the file well?
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   !----------------------------
   real(8) M
   integer Reason, i, j, k, l, N, Shl, CDF_coef, Shl_num, comment_start
   character(100) Error_descript, read_line
   character(3) Name
   character(11) Shell_name
   character(30) Full_Name

   i = 0
   !READ(FN2,*,IOSTAT=Reason) Matter%Target_name ! first line is the full material name
   READ(FN2,'(a)',IOSTAT=Reason) read_line   ! read the full line
   comment_start = index(read_line, '!')  ! to exclude comment, if any
   ! Save the name into the variable:
   if (comment_start > 1) then   ! trim comment
      Matter%Target_name = read_line(1:comment_start-1)
   else  ! no comment
      Matter%Target_name = trim(adjustl(read_line))
   endif
    ! Also, remove TABs if any:
   comment_start = index(read_line, CHAR(9))  ! to exclude TABs, if any
   if (comment_start > 1) then   ! trim comment
      Matter%Target_name = read_line(1:comment_start-1)
   elseif (comment_start == 1) then   ! omit leading TAB
      Matter%Target_name = read_line(2:)
   else  ! no comment
      Matter%Target_name = trim(adjustl(read_line))
   endif
   Matter%Target_name = trim(adjustl(Matter%Target_name))


   READ(FN2,*,IOSTAT=Reason) N   ! number of elements in this compound
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2014
   if (.not. allocated(Target_atoms)) allocate(Target_atoms(N)) ! that's how many atom kinds we have
   
   do j = 1, N  ! read for each element it's basic data:
      READ(FN2,*,IOSTAT=Reason) Target_atoms(j)%Zat, Target_atoms(j)%Pers ! atmoic number and its persentage in the compound's molecule
      call read_file(Reason, i, read_well) ! reports if everything read well
      if (.not. read_well) goto 2014
      call Find_element_name(Target_atoms(j)%Zat, Name, Full_Name, M) ! from module 'Dealing_with_EADL'
      Target_atoms(j)%Name = Name  ! name of the element
      Target_atoms(j)%Full_Name = Full_Name  ! full name of the element
      Target_atoms(j)%Mass = M  ! mass of the element in the proton-mass units
   enddo

   ! material density [g/cm^3], speed of sound in the material [m/s]; Fermi energy[eV]; bandgap [eV]
   READ(FN2,*,IOSTAT=Reason) Matter%Dens, Matter%Vsound, Matter%E_F, Matter%Egap
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2014
   ! [1/cm^3] atomic density:
   !Matter%At_Dens = 1.0d-3*Matter%Dens/(g_Mp*SUM(Target_atoms(:)%Mass)/size(Target_atoms)) ! [1/cm^3] atomic density
   Matter%At_Dens = 1.0d-3*Matter%Dens/(g_Mp*SUM(Target_atoms(:)%Mass*Target_atoms(:)%Pers)/SUM(Target_atoms(:)%Pers))
   Matter%v_f = sqrt(2.0d0*Matter%E_f/g_me) ! units as used in one-pole approximation (not SI!)
   
   call check_atomic_parameters(NumPar, Target_atoms, Error_message=Error_message, read_well=read_well, MPI_param=MPI_param) ! from module 'Dealing_with_EADL'
   
   do j = 1, N ! now mark all the cross sections as BEB:
      Target_atoms(j)%KOCS(:) = 2 ! BEB cross section for all shells
      Target_atoms(j)%KOCS_SHI(:) = 2 ! BEB cross section for all shells
   enddo
2014 continue
end subroutine read_short_scdf


subroutine reading_material_DOS(DOS_file, Mat_DOS, Matter, Target_atoms, NumPar, Error_message, read_well, MPI_param)    ! read material DOS
    character(100), intent(in) :: DOS_file  ! file with material DOS
    type(Density_of_states), intent(inout) :: Mat_DOS  ! materail DOS
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    type(Flag), intent(inout) :: NumPar ! numerical parameters
    logical, intent(inout) :: read_well  ! did we read the file well?
    type(Solid), intent(in) :: Matter
    type(Atom), dimension(:), intent(in) :: Target_atoms
    type(Used_MPI_parameters), intent(inout) :: MPI_param ! MPI parameters
    !--------------------------------------
    real(8), dimension(:,:), allocatable :: Temp_DOS
    real(8) Sum_DOS, loc_DOS, E, dE, Sum_DOS_inv
    integer FN2, i, N, Reason, M
    character(100) :: Error_descript, Free_DOS
    logical file_opened, file_exist, file_exist2


    if (MPI_param%process_rank /= 0) then   ! only MPI master process does it
      return   ! nothing to do here for non-master processes
    endif


    FN2 = 202
    inquire(file=trim(adjustl(DOS_file)),exist=file_exist)    ! check if input file excists
    !inquire(file='INPUT_DOS/Free_electron_DOS.dos',exist=file_exist2)    ! check if input file excists
    Free_DOS = trim(adjustl(m_INPUT_DOS))//trim(adjustl(NumPar%path_sep))//'Free_electron_DOS.dos'
    inquire(file=trim(adjustl(Free_DOS)),exist=file_exist2)    ! check if input file excists

    if (file_exist) then
       if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         print*, 'DOS file is there: ', trim(adjustl(DOS_file))
       endif
       open(unit = FN2, FILE = trim(adjustl(DOS_file)), status = 'old', readonly)   ! if yes, open it and read
       ! Save DOS file name for output:
       NumPar%DOS_file = trim(adjustl(DOS_file))
    elseif (file_exist2) then   ! Free-electron DOS approximation
       if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         print*, 'DOS file ', trim(adjustl(DOS_file)), ' is not found, '
         print*, 'free-electron DOS approximation is used for valence band holes'
       endif
       open(unit = FN2, FILE = trim(adjustl(Free_DOS)), status = 'old', readonly)   ! if yes, open it and read
       ! Save DOS file name for output:
       NumPar%DOS_file = trim(adjustl(Free_DOS))
    else ! no DOS, try atomic approxiamtion instead...
       if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         print*, 'Neither file ', trim(adjustl(DOS_file)), ' nor Free_electron_DOS.dos'
         print*, 'containing density of states are found.'
         print*, 'The calculations proceed within the atomic approximation for energy levels.'
       endif
       Error_descript = 'Files '//trim(adjustl(DOS_file))//' and '//trim(adjustl(Free_DOS))//' are not found!'    ! description of an error
       call Save_error_details(Error_message, 6, Error_descript, MPI_param) ! write it into the error-log file
       read_well = .false.  ! no file found
       ! Save DOS file name for output:
       NumPar%DOS_file = ''   ! nothing
       goto 2016
    endif

    call Count_lines_in_file(FN2, N) ! count how many line the file contains
    if (.not. allocated(Temp_DOS)) allocate(Temp_DOS(2,N)) ! to read from file
    do i = 1, N  ! read all lines
        read(FN2, *, IOSTAT=Reason) Temp_DOS(1,i), Temp_DOS(2,i)
        call read_file_here(Reason, i, read_well)
        if (.not. read_well) goto 2016
    enddo
    Temp_DOS(1,:) = Temp_DOS(1,:) - Temp_DOS(1,N)    ! shift the topmost level of the VB to 'zero'
    dE = 0.1d0   ! [eV] energy grid
    M = FLOOR((Temp_DOS(1,N) - Temp_DOS(1,1))/dE)
    if (.not. allocated(Mat_DOS%E)) allocate(Mat_DOS%E(M))
    if (.not. allocated(Mat_DOS%DOS)) allocate(Mat_DOS%DOS(M))
    if (.not. allocated(Mat_DOS%int_DOS)) allocate(Mat_DOS%int_DOS(M))
    if (.not. allocated(Mat_DOS%k)) allocate(Mat_DOS%k(M))
    if (.not. allocated(Mat_DOS%Eff_m)) allocate(Mat_DOS%Eff_m(M))
    
    if (.not. allocated(Mat_DOS%DOS_inv)) allocate(Mat_DOS%DOS_inv(M))
    if (.not. allocated(Mat_DOS%int_DOS_inv)) allocate(Mat_DOS%int_DOS_inv(M))
    if (.not. allocated(Mat_DOS%k_inv)) allocate(Mat_DOS%k_inv(M))
    if (.not. allocated(Mat_DOS%Eff_m_inv)) allocate(Mat_DOS%Eff_m_inv(M))
    
    E = Temp_DOS(1,1)    ! start from the start [eV]
    do i = 1, M  ! fill all points in the array:
        E = E + dE   ! energy grid [eV]
        call Linear_approx(Temp_DOS, E, loc_DOS, (Temp_DOS(1,1)-dE), 0.0d0)
        Mat_DOS%E(i) = E ! [eV] energy
        Mat_DOS%DOS(i) = loc_DOS ! [1/eV] DOS
        !print*, 'DOS', i,  Mat_DOS%E(i), Mat_DOS%DOS(i)
    enddo
    Mat_DOS%E = ABS(Mat_DOS%E - Mat_DOS%E(size(Mat_DOS%E))) ! shift it to 'zero', and make it positive [eV]
    Mat_DOS%E = Mat_DOS%E(size(Mat_DOS%E):1:-1) ! make the array increasing
    Mat_DOS%DOS = Mat_DOS%DOS(size(Mat_DOS%DOS):1:-1) ! make the array according to the increasing energy
    Mat_DOS%DOS_inv = Mat_DOS%DOS(size(Mat_DOS%DOS):1:-1) ! invert array for metals

    Sum_DOS_inv = 0.0d0  ! start to sum up:
    
    Sum_DOS = 0.0d0  ! start to sum up:
    do i = 1, M  ! fill all points in the array:
        Sum_DOS = Sum_DOS + Mat_DOS%DOS(i)  ! for integrated DOS
        Mat_DOS%int_DOS(i) = Sum_DOS ! integrated DOS
        
        Sum_DOS_inv = Sum_DOS_inv + Mat_DOS%DOS_inv(i)
        Mat_DOS%int_DOS_inv(i) = Sum_DOS_inv
    enddo
    
    SUM_DOS = Mat_DOS%int_DOS(size(Mat_DOS%E))   ! sum of all electrons in the DOS, to normalize it
    Mat_DOS%DOS = Mat_DOS%DOS/SUM_DOS*Matter%N_VB_el    ! normalize it to the number of VB electron per molecule
    Mat_DOS%int_DOS = Mat_DOS%int_DOS/SUM_DOS*Matter%N_VB_el    ! normalize it to the number of VB electron per molecule
       
    Sum_DOS = 0.0d0  ! start to sum up:
    do i = 1, M  ! fill all points in the array:
        Sum_DOS = Sum_DOS + Mat_DOS%DOS(i)  ! for integrated DOS
!         if (Mat_DOS%DOS(i) .LT. 1.0d-6) then  ! NOT READY, gives error for multiple sub-bands DOS
!             Sum_DOS = 0
!         endif
!         Mat_DOS%k(i) = (3.0d0*2.0d0*g_Pi*g_Pi/2.0d0*Sum_DOS*Matter%At_Dens/size(Target_atoms)*1d6)**(1.0d0/3.0d0)  ! [1/m]
        Mat_DOS%k(i) = (3.0d0*2.0d0*g_Pi*g_Pi/2.0d0*Sum_DOS*Matter%At_Dens/sum(Target_atoms(:)%Pers)*1d6)**(1.0d0/3.0d0)  ! [1/m]
        
        if (Mat_DOS%E(i) < 1.0d-10) then  ! too slow electron/hole
           !Mat_DOS%Eff_m(i) = 1.0d20
           Mat_DOS%Eff_m(i) = 1.0d0 ! too slow to get rom DOS, assume free particle mass
        else    ! regular energy particle
           Mat_DOS%Eff_m(i) = g_h*g_h*Mat_DOS%k(i)*Mat_DOS%k(i)/(2.0*Mat_DOS%E(i)*g_e)/g_me
        endif
        !print*, i, Mat_DOS%E(i), Mat_DOS%DOS(i), Mat_DOS%int_DOS(i), Mat_DOS%int_DOS_inv(i)
    enddo

    SUM_DOS = Mat_DOS%int_DOS_inv(size(Mat_DOS%E))   ! sum of all electrons in the DOS, to normalize it
    Mat_DOS%DOS_inv = Mat_DOS%DOS_inv/SUM_DOS*Matter%N_VB_el    ! normalize it to the number of VB electron per molecule
    Mat_DOS%int_DOS_inv = Mat_DOS%int_DOS_inv/SUM_DOS*Matter%N_VB_el    ! normalize it to the number of VB electron per molecule
    Mat_DOS%k_inv = (3.0d0*2.0d0*g_Pi*g_Pi/2.0d0*Mat_DOS%Int_DOS_inv*Matter%At_Dens/sum(Target_atoms(:)%Pers)*1d6)**(1.0d0/3.0d0)  ! [1/m]
    where(Mat_DOS%E(:) < 1.0d-10) ! too slow electron/hole
       !Mat_DOS%Eff_m_inv(:) = 1.0d20
       Mat_DOS%Eff_m_inv(:) = 1.0d0 ! too slow to get rom DOS, assume free particle mass
    elsewhere   ! regular energy particle
       Mat_DOS%Eff_m_inv(:) = g_h*g_h*Mat_DOS%k_inv(:)*Mat_DOS%k_inv(:)/(2.0*Mat_DOS%E(:)*g_e)/g_me
    endwhere

2016 inquire(unit=FN2,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN2)             ! and if it is, close it
end subroutine reading_material_DOS


subroutine broadcast_material_DOS(Mat_DOS, MPI_param)
   type(Density_of_states), intent(inout) :: Mat_DOS  ! materail DOS
   type(Used_MPI_parameters), intent(inout) :: MPI_param ! MPI parameters
   !--------------------------------------
   integer :: Nsiz
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in broadcast_material_DOS'
   if (MPI_param%process_rank == 0) then
      Nsiz = size(Mat_DOS%E)
   else
      Nsiz = 0
   endif
   call mpi_bcast(Nsiz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Mat_DOS:Nsiz}') ! module "MPI_subroutines"

   if (MPI_param%process_rank /= 0) then
      if (.not.allocated(Mat_DOS%E)) allocate(Mat_DOS%E(Nsiz))
      if (.not.allocated(Mat_DOS%k)) allocate(Mat_DOS%k(Nsiz))
      if (.not.allocated(Mat_DOS%DOS)) allocate(Mat_DOS%DOS(Nsiz))
      if (.not.allocated(Mat_DOS%int_DOS)) allocate(Mat_DOS%int_DOS(Nsiz))
      if (.not.allocated(Mat_DOS%Eff_m)) allocate(Mat_DOS%Eff_m(Nsiz))
      if (.not.allocated(Mat_DOS%DOS_inv)) allocate(Mat_DOS%DOS_inv(Nsiz))
      if (.not.allocated(Mat_DOS%int_DOS_inv)) allocate(Mat_DOS%int_DOS_inv(Nsiz))
      if (.not.allocated(Mat_DOS%k_inv)) allocate(Mat_DOS%k_inv(Nsiz))
      if (.not.allocated(Mat_DOS%Eff_m_inv)) allocate(Mat_DOS%Eff_m_inv(Nsiz))
   endif

   call mpi_bcast(Mat_DOS%E, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Mat_DOS%E}') ! module "MPI_subroutines"

   call mpi_bcast(Mat_DOS%DOS, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Mat_DOS%DOS}') ! module "MPI_subroutines"

   call mpi_bcast(Mat_DOS%int_DOS, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Mat_DOS%int_DOS}') ! module "MPI_subroutines"

   call mpi_bcast(Mat_DOS%k, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Mat_DOS%k}') ! module "MPI_subroutines"

   call mpi_bcast(Mat_DOS%Eff_m, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Mat_DOS%Eff_m}') ! module "MPI_subroutines"

   call mpi_bcast(Mat_DOS%DOS_inv, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Mat_DOS%DOS_inv}') ! module "MPI_subroutines"

   call mpi_bcast(Mat_DOS%int_DOS_inv, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Mat_DOS%int_DOS_inv}') ! module "MPI_subroutines"

   call mpi_bcast(Mat_DOS%k_inv, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Mat_DOS%k_inv}') ! module "MPI_subroutines"

   call mpi_bcast(Mat_DOS%Eff_m_inv, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {Mat_DOS%Eff_m_inv}') ! module "MPI_subroutines"
#endif
end subroutine broadcast_material_DOS


subroutine reading_DSF_cross_sections(DSF_file, DSF_DEMFP, NumPar, Error_message, read_well, MPI_param)    ! read material DOS
    character(100), intent(in) :: DSF_file  ! file with DSF crossections
    type(Differential_MFP), dimension(:), allocatable, intent(inout) :: DSF_DEMFP  ! diffential elastic MFPs calculated with DSF
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    logical, intent(inout) :: read_well  ! did we read the file well?
    type(Flag), intent(inout) :: NumPar
    type(Used_MPI_parameters), intent(inout) :: MPI_param ! MPI parameters
    !-------------------------------
    real(8), dimension(:,:), allocatable :: Temp_EMFP, Temp_array
    type(Differential_MFP), dimension(:), allocatable :: Temp_DEMFP
    real(8), dimension(:), allocatable :: Temp_E
    real(8) Sum_MFP, loc_EMFP, E, dE, Emin, Emax, Sum_MFP_emit
    integer FN2, i, j, N, Reason, M, NEPo, NTEPo, k
    character(100) Error_descript
    logical file_opened, file_exist, file_exist2

    if (MPI_param%process_rank /= 0) then   ! only MPI master process does it
      return   ! nothing to do here for non-master process
    endif

    FN2 = 212
    inquire(file=trim(adjustl(DSF_file)),exist=file_exist)    ! check if input file excists
    if (file_exist) then
       if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         print*, 'File with DSF elastic cross-sections is there: ', trim(adjustl(DSF_file))
       endif
       open(unit = FN2, FILE = trim(adjustl(DSF_file)), status = 'old', readonly)   ! if yes, open it and read
    else ! no DSF file, try CDF phonon peaks or atomic approxiamtion instead...
       if (MPI_param%process_rank == 0) then   ! only MPI master process does it
         print*, 'File ', trim(adjustl(DSF_file)), ' is not found.'
         print*, 'The calculations proceed with Mott atomic cross-sections.'
       endif
       Error_descript = 'File '//trim(adjustl(DSF_file))//' is not found.'    ! description of an error
       call Save_error_details(Error_message, 6, Error_descript, MPI_param) ! write it into the error-log file
       read_well = .true.  ! no file found
       NumPar%kind_of_EMFP = 0
       goto 2017
    endif

    call Count_lines_in_file(FN2, N) ! count how many line the file contains
    M = N - 1

    allocate(Temp_EMFP(2,M)) ! to read from file
    allocate(Temp_E(M)) ! to read from file
    read(FN2, *, IOSTAT=Reason) NTEPo        ! Number of transferred Energy Points for each energy

    do i = 2, M ! read all lines
        read(FN2, *, IOSTAT=Reason) Temp_E(i-1), Temp_EMFP(1,i-1), Temp_EMFP(2,i-1)
        if (Temp_EMFP(2,i-1) .LT. 1.0d-10) Temp_EMFP(2,i-1) = 0.0d0
        call read_file_here(Reason, i-1, read_well)
        if (.not. read_well) goto 2017
    enddo

    NEPo = M/NTEPo
    if(allocated(DSF_DEMFP)) deallocate(DSF_DEMFP)
    allocate(DSF_DEMFP(NEPo))
    allocate(Temp_DEMFP(NEPo))
    allocate(Temp_array(2,NTEPo))
    k = 0
    do i = 1, NEPo
        ! Used variable:
        allocate(DSF_DEMFP(i)%dE(NTEPo), source = 0.0d0)
        allocate(DSF_DEMFP(i)%dL(NTEPo), source = 0.0d0)
        allocate(DSF_DEMFP(i)%dL_absorb(NTEPo), source = 0.0d0)
        allocate(DSF_DEMFP(i)%dL_emit(NTEPo), source = 0.0d0)
        ! Temporary variable:
        allocate(Temp_DEMFP(i)%dE(NTEPo), source = 0.0d0)
        allocate(Temp_DEMFP(i)%dL(NTEPo), source = 0.0d0)
        allocate(Temp_DEMFP(i)%dL_absorb(NTEPo), source = 0.0d0)
        allocate(Temp_DEMFP(i)%dL_emit(NTEPo), source = 0.0d0)

        Temp_DEMFP(i)%E = Temp_E(1+NTEPo*(i-1))

        do j = 1, NTEPo  ! fill all points in the array:
            k = k+1
            Temp_DEMFP(i)%dE(j) = Temp_EMFP(1,k) ! [eV] energy
            Temp_DEMFP(i)%dL(j) = Temp_EMFP(2,k) ! [1/eV] DSF
            if (Temp_DEMFP(i)%dE(j) >= 0.0d0) then ! emission
               Temp_DEMFP(i)%dL_emit(j) = Temp_EMFP(2,k) ! [1/eV] DSF
            else ! absorption
               Temp_DEMFP(i)%dL_absorb(j) = Temp_EMFP(2,k) ! [1/eV] DSF
            endif
        enddo
    enddo
    DSF_DEMFP(:)%E = Temp_DEMFP(:)%E

    do i = 1, NEPo
        do j = 1, NTEPo  ! fill all points in the array:
            DSF_DEMFP(i)%dE(j) = Temp_DEMFP(i)%dE(NTEPo+1-j)
            DSF_DEMFP(i)%dL(j) = Temp_DEMFP(i)%dL(NTEPo+1-j)
        enddo
    enddo

    do i = 1, NEPo  ! fill all points in the array:
        Sum_MFP = 0.0d0
        Sum_MFP_emit = 0.0d0

        Temp_array(1,:) = DSF_DEMFP(i)%dE(:)
        Temp_array(2,:) = DSF_DEMFP(i)%dL(:)

        Emin = -0.2d0 !Temp_array(1,1)    ! start from the first point of transferred energy [eV]
        E = Emin
        Emax = Temp_array(1,NTEPo)
        if (Emax .GT. 0.2d0) Emax = 0.2d0

        dE = (Emax - Emin)/dble(NTEPo)

        do j = 1, NTEPo
            E = E + dE   ! energy grid [eV]

            call Linear_approx(Temp_array, E, loc_EMFP, (Emin-dE), 0.0d0)
            Temp_DEMFP(i)%dE(j) = E                         ! [eV] transferred energy
            Temp_DEMFP(i)%dL(j) = loc_EMFP                  ! [(A*eV)^-1] cross-section

            Sum_MFP = Sum_MFP + Temp_DEMFP(i)%dL(j)*dE      ! for integrated EMFP
            DSF_DEMFP(i)%dE(j) = Temp_DEMFP(i)%dE(j)
            if (abs(Sum_MFP) > 1.0d-10) then
                DSF_DEMFP(i)%dL(j) = 1.0d0/Sum_MFP          ! integrated MFP
            else
                DSF_DEMFP(i)%dL(j) = 1.0d30
            endif

            ! Absorbtion and emittion separately:
            if (DSF_DEMFP(i)%dE(j) >= 0.0d0) then ! emission
               Sum_MFP_emit = Sum_MFP_emit + Temp_DEMFP(i)%dL(j)*dE ! for separated integrated EMFP
               if (abs(Sum_MFP_emit) > 1.0d-10) then
                  DSF_DEMFP(i)%dL_emit(j) = 1.0d0/Sum_MFP_emit ! integrated MFP
               else
                  DSF_DEMFP(i)%dL_emit(j) = 1.0d30
               endif
               ! Absorption does not change:
               if (j == 1) then
                  DSF_DEMFP(i)%dL_absorb(j) = 1.0d30
               else
                  DSF_DEMFP(i)%dL_absorb(j) = DSF_DEMFP(i)%dL_absorb(j-1)
               endif
            else ! absorption
               if (abs(Sum_MFP) > 1.0d-10) then
                  DSF_DEMFP(i)%dL_absorb(j) = 1.0d0/Sum_MFP ! integrated MFP
               else
                  DSF_DEMFP(i)%dL_absorb(j) = 1.0d30
               endif
               ! Emission does not change:
               if (j == 1) then
                  DSF_DEMFP(i)%dL_emit(j) = 1.0d30
               else
                  DSF_DEMFP(i)%dL_emit(j) = DSF_DEMFP(i)%dL_emit(j-1)
               endif
            endif ! (DSF_DEMFP(i)%dE(j) >= 0.0d0)

!             write(*,'(i6, f, f, es, es, es)'), j, DSF_DEMFP(i)%E, DSF_DEMFP(i)%dE(j), DSF_DEMFP(i)%dL(j), DSF_DEMFP(i)%dL_emit(j), DSF_DEMFP(i)%dL_absorb(j)
        enddo
    enddo

    deallocate(Temp_E)
    deallocate(Temp_DEMFP)
    deallocate(Temp_EMFP)
    deallocate(Temp_array)

2017 inquire(unit=FN2,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN2)             ! and if it is, close it
!     pause 'reading_DSF_cross_sections'
end subroutine reading_DSF_cross_sections



subroutine broadcast_DSF_cross_sections(DSF_DEMFP, read_well, NumPar, MPI_param, marker)
   type(Differential_MFP), dimension(:), allocatable, intent(inout) :: DSF_DEMFP  ! diffential elastic MFPs calculated with DSF
   logical, intent(inout) :: read_well  ! did we read the file well?
   type(Flag), intent(inout) :: NumPar
   type(Used_MPI_parameters), intent(inout) :: MPI_param ! MPI parameters
   character(*), intent(in) :: marker
   !------------------------
   integer :: Nsiz, N1, i
   character(100) :: error_part

#ifdef MPI_USED
   error_part = 'ERROR in broadcast_DSF_cross_sections ('//trim(adjustl(marker))//')'

   call mpi_bcast(read_well, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {DSF_DEMFP:read_well}') ! module "MPI_subroutines"

   call mpi_bcast(NumPar%kind_of_EMFP, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {DSF_DEMFP:kind_of_EMFP}') ! module "MPI_subroutines"

   if (MPI_param%process_rank == 0) then
      Nsiz = size(DSF_DEMFP)
   else
      Nsiz = 0
   endif
   call mpi_bcast(Nsiz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {DSF_DEMFP:Nsiz}') ! module "MPI_subroutines"

   if (MPI_param%process_rank == 0) then
      N1 = size(DSF_DEMFP(1)%dE)
   else
      N1 = 0
   endif
   call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {DSF_DEMFP:N1}') ! module "MPI_subroutines"

   if (MPI_param%process_rank /= 0) then
      if (.not.allocated(DSF_DEMFP)) allocate(DSF_DEMFP(Nsiz))
      do i = 1, Nsiz
         if (.not.allocated(DSF_DEMFP(i)%dE)) allocate(DSF_DEMFP(i)%dE(N1))
         if (.not.allocated(DSF_DEMFP(i)%dL)) allocate(DSF_DEMFP(i)%dL(N1))
         if (.not.allocated(DSF_DEMFP(i)%dL_absorb)) allocate(DSF_DEMFP(i)%dL_absorb(N1))
         if (.not.allocated(DSF_DEMFP(i)%dL_emit)) allocate(DSF_DEMFP(i)%dL_emit(N1))
      enddo
   endif

   do i = 1, Nsiz
      call mpi_bcast(DSF_DEMFP(i)%E, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {DSF_DEMFP(i)%E}') ! module "MPI_subroutines"

      call mpi_bcast(DSF_DEMFP(i)%dE, N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {DSF_DEMFP(i)%dE}') ! module "MPI_subroutines"

      call mpi_bcast(DSF_DEMFP(i)%dL, N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {DSF_DEMFP(i)%dL}') ! module "MPI_subroutines"

      call mpi_bcast(DSF_DEMFP(i)%dL_absorb, N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {DSF_DEMFP(i)%dL_absorb}') ! module "MPI_subroutines"

      call mpi_bcast(DSF_DEMFP(i)%dL_emit, N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//' {DSF_DEMFP(i)%dL_emit}') ! module "MPI_subroutines"
   enddo
#endif
end subroutine broadcast_DSF_cross_sections



subroutine reading_DSF_cross_sections_OLD(DSF_file, DSF_DEMFP, NumPar, Error_message, read_well, MPI_param)    ! read material DOS
    character(100), intent(in) :: DSF_file  ! file with DSF crossections
    type(Differential_MFP), dimension(:), allocatable, intent(inout) :: DSF_DEMFP  ! diffential elastic MFPs calculated with DSF
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    logical, intent(inout) :: read_well  ! did we read the file well?
    type(Flag), intent(inout) :: NumPar
    type(Used_MPI_parameters), intent(inout) :: MPI_param
    !--------------------------
    real(8), dimension(:,:), allocatable :: Temp_EMFP, Temp_array
    
    type(Differential_MFP), dimension(:), allocatable :: Temp_DEMFP
    real(8), dimension(:), allocatable :: Temp_E
    real(8) Sum_MFP, loc_EMFP, E, dE, Emin, Emax
    integer FN2, i, j, N, Reason, M, NEPo, NTEPo, k
    character(100) Error_descript
    logical file_opened, file_exist, file_exist2
    
    FN2 = 212
    inquire(file=trim(adjustl(DSF_file)),exist=file_exist)    ! check if input file excists
    if (file_exist) then
       print*, 'File with DSF elastic cross-sections is there: ', trim(adjustl(DSF_file))
       open(unit = FN2, FILE = trim(adjustl(DSF_file)), status = 'old', readonly)   ! if yes, open it and read
    else ! no DSF file, try CDF phonon peaks or atomic approxiamtion instead...
       print*, 'File ', trim(adjustl(DSF_file)), ' is not found.'
       print*, 'The calculations proceed with Mott atomic cross-sections.'
       Error_descript = 'File '//trim(adjustl(DSF_file))//' is not found.'    ! description of an error
       call Save_error_details(Error_message, 6, Error_descript, MPI_param) ! write it into the error-log file
       read_well = .true.  ! no file found
       NumPar%kind_of_EMFP = 0
       goto 2017
    endif

    call Count_lines_in_file(FN2, N) ! count how many line the file contains
    M = N - 1
    if (.not. allocated(Temp_EMFP)) allocate(Temp_EMFP(2,M)) ! to read from file
    if (.not. allocated(Temp_E)) allocate(Temp_E(M)) ! to read from file
    read(FN2, *, IOSTAT=Reason) NTEPo        ! Number of transferred Energy Points for each energy
    do i = 2, M ! read all lines
        read(FN2, *, IOSTAT=Reason) Temp_E(i-1), Temp_EMFP(1,i-1), Temp_EMFP(2,i-1)
        call read_file_here(Reason, i-1, read_well)
        if (.not. read_well) goto 2017
    enddo
    
    NEPo = M/NTEPo
                
    if (.not. allocated(DSF_DEMFP)) allocate(DSF_DEMFP(NEPo))
    if (.not. allocated(Temp_DEMFP)) allocate(Temp_DEMFP(NEPo))
    if (.not. allocated(Temp_array)) allocate(Temp_array(2,NTEPo))
    k = 0
    do i = 1, NEPo
        if (.not. allocated(DSF_DEMFP(i)%dE)) allocate(DSF_DEMFP(i)%dE(NTEPo))
        if (.not. allocated(Temp_DEMFP(i)%dE)) allocate(Temp_DEMFP(i)%dE(NTEPo))
        if (.not. allocated(DSF_DEMFP(i)%dL)) allocate(DSF_DEMFP(i)%dL(NTEPo))
        if (.not. allocated(Temp_DEMFP(i)%dL)) allocate(Temp_DEMFP(i)%dL(NTEPo))
        
        Temp_DEMFP(i)%E = Temp_E(1+NTEPo*(i-1))
        
        do j = 1, NTEPo  ! fill all points in the array:
            k = k+1
            Temp_DEMFP(i)%dE(j) = Temp_EMFP(1,k) ! [eV] energy
            Temp_DEMFP(i)%dL(j) = Temp_EMFP(2,k) ! [1/eV] DSF
        enddo
    enddo
    DSF_DEMFP(:)%E = Temp_DEMFP(:)%E
    
    do i = 1, NEPo
        do j = 1, NTEPo  ! fill all points in the array:
            DSF_DEMFP(i)%dE(j) = Temp_DEMFP(i)%dE(NTEPo+1-j)
            DSF_DEMFP(i)%dL(j) = Temp_DEMFP(i)%dL(NTEPo+1-j)
        enddo
    enddo
       
    do i = 1, NEPo  ! fill all points in the array:
        sum_MFP = 0.0d0
                
        Temp_array(1,:) = DSF_DEMFP(i)%dE(:)
        Temp_array(2,:) = DSF_DEMFP(i)%dL(:)
                
        Emin = Temp_array(1,1)    ! start from the first point of transferred energy [eV]
        E = Emin
        Emax = Temp_array(1,NTEPo)
        dE = (Emax - Emin)/real(NTEPo)
        
        do j = 1, NTEPo
            E = E + dE   ! energy grid [eV]

            call Linear_approx(Temp_array, E, loc_EMFP, (Emin-dE), 0.0d0)
            Temp_DEMFP(i)%dE(j) = E                         ! [eV] transferred energy
            Temp_DEMFP(i)%dL(j) = loc_EMFP                  ! [(A*eV)^-1] cross-section
            
            Sum_MFP = Sum_MFP + Temp_DEMFP(i)%dL(j)*dE      ! for integrated EMFP
            DSF_DEMFP(i)%dL(j) = 1.0d0/Sum_MFP              ! integrated MFP
            DSF_DEMFP(i)%dE(j) = Temp_DEMFP(i)%dE(j)
            
        enddo
    enddo
       
    deallocate(Temp_E)
    deallocate(Temp_DEMFP)
    deallocate(Temp_EMFP)
    deallocate(Temp_array)
        
2017 inquire(unit=FN2,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN2)             ! and if it is, close it
end subroutine reading_DSF_cross_sections_OLD


subroutine read_SHI_MFP(FN, FN2, Nat, Target_atoms, SHI_MFP, read_well)
   integer, intent(in) :: FN, FN2   ! files numbers from where to read
   integer, intent(in) :: Nat   ! number of atoms
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
   type(All_MFP), dimension(:),allocatable, intent(inout) :: SHI_MFP    ! SHI mean free paths and dEdx to read from file
   logical, intent(inout) :: read_well  ! did we read the file well?
   integer j, i,k,N, Reason, Nshl
   logical file_opened, file_exist
   
   call Count_lines_in_file(FN, N) ! count how many lines are in this file
   if (.not. allocated(SHI_MFP)) allocate(SHI_MFP(Nat)) ! number of atoms
   do j = 1, Nat
        Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
        if (.not. allocated(SHI_MFP(j)%ELMFP)) allocate(SHI_MFP(j)%ELMFP(Nshl))
        do k = 1, Nshl
            if (.not. allocated(SHI_MFP(j)%ELMFP(k)%E)) allocate(SHI_MFP(j)%ELMFP(k)%E(N))
            if (.not. allocated(SHI_MFP(j)%ELMFP(k)%L)) allocate(SHI_MFP(j)%ELMFP(k)%L(N))
            if (.not. allocated(SHI_MFP(j)%ELMFP(k)%dEdx)) allocate(SHI_MFP(j)%ELMFP(k)%dEdx(N))
        enddo
   enddo
   ! Now write the output into the file:
   do i = 1, N
      read(FN,'(e)', advance='no', IOSTAT=Reason) SHI_MFP(1)%ELMFP(1)%E(i)
      SHI_MFP(1)%ELMFP(1)%E(i) = SHI_MFP(1)%ELMFP(1)%E(i)*1d6 ! MeV -> eV
      call read_file_here(Reason, i, read_well)
      if (.not. read_well) goto 2019
      read(FN2,'(e)', advance='no', IOSTAT=Reason) !SHI_MFP(1)%ELMFP(1)%E(i)
      call read_file_here(Reason, i, read_well)
      if (.not. read_well) goto 2019
      do j = 1, Nat
         Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
         do k = 1, Nshl
            read(FN,'(e)', advance='no', IOSTAT=Reason) SHI_MFP(j)%ELMFP(k)%L(i)    ! write IMFP for all shells
            call read_file_here(Reason, i, read_well)
            if (.not. read_well) goto 2019
            read(FN2,'(e)', advance='no', IOSTAT=Reason) SHI_MFP(j)%ELMFP(k)%dEdx(i)    ! write IMFP for all shells
            call read_file_here(Reason, i, read_well)
            if (.not. read_well) goto 2019
         enddo
      enddo
      read(FN,'(e)',IOSTAT=Reason)
      call read_file_here(Reason, i, read_well)
      if (.not. read_well) goto 2019
      read(FN2,'(e)',IOSTAT=Reason)
      call read_file_here(Reason, i, read_well)
      if (.not. read_well) goto 2019
   enddo
2019 continue
end subroutine read_SHI_MFP


subroutine broadcast_SHI_MFP(SHI_MFP, MPI_param)
   type(All_MFP), dimension(:),allocatable, intent(inout) :: SHI_MFP    ! SHI mean free paths and dEdx to read from file
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   !------------------------------------------
   integer :: Nsiz, Nshl, i, j, k, N1, N2
   character(100) :: error_part
   logical :: do_broadcast
#ifdef MPI_USED
   ! initialize:
   do_broadcast = .false.
   error_part = 'ERROR in broadcast_SHI_MFP'

   if (MPI_param%process_rank == 0) then
      Nsiz = size(SHI_MFP)   ! define number of atoms in the compound
   else
      Nsiz = 0  ! unknown yet, get from the master process
   endif
   call mpi_bcast(Nsiz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{Nsiz}') ! module "MPI_subroutines"

   if (MPI_param%process_rank /= 0) then
      if (.not. allocated(SHI_MFP)) allocate(SHI_MFP(Nsiz)) ! number of atoms
   endif

   do j = 1, Nsiz
      if (MPI_param%process_rank == 0) then
         Nshl = size(SHI_MFP(j)%ELMFP)    ! how mamy shells
      else
         Nshl = 0  ! unknown yet, get from the master process
      endif
      call mpi_bcast(Nshl, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{Nshl}') ! module "MPI_subroutines"
      !------------------------------------------------------
      ! Synchronize MPI processes to make sure the directory is created before going further
      call MPI_barrier_wrapper(MPI_param)  ! module "MPI_subroutines"
      !------------------------------------------------------
      if (MPI_param%process_rank /= 0) then
         if (.not. allocated(SHI_MFP(j)%ELMFP)) allocate(SHI_MFP(j)%ELMFP(Nshl))
      endif

      do k = 1, Nshl
         if (MPI_param%process_rank == 0) then
            N1 = size(SHI_MFP(j)%ELMFP(k)%E)    ! how mamy grid points
         else
            N1 = 0  ! unknown yet, get from the master process
         endif

         call mpi_bcast(N1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{N1}') ! module "MPI_subroutines"
         !------------------------------------------------------
         ! Synchronize MPI processes to make sure the directory is created before going further
         call MPI_barrier_wrapper(MPI_param)  ! module "MPI_subroutines"
         !------------------------------------------------------
         if (MPI_param%process_rank /= 0) then
            if (.not. allocated(SHI_MFP(j)%ELMFP(k)%E)) allocate(SHI_MFP(j)%ELMFP(k)%E(N1))
            if (.not. allocated(SHI_MFP(j)%ELMFP(k)%L)) allocate(SHI_MFP(j)%ELMFP(k)%L(N1))
            if (.not. allocated(SHI_MFP(j)%ELMFP(k)%dEdx)) allocate(SHI_MFP(j)%ELMFP(k)%dEdx(N1))
         endif

         call mpi_bcast(SHI_MFP(j)%ELMFP(k)%E, N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{SHI_MFP(j)%ELMFP(k)%E}') ! module "MPI_subroutines"

         call mpi_bcast(SHI_MFP(j)%ELMFP(k)%L, N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{SHI_MFP(j)%ELMFP(k)%L}') ! module "MPI_subroutines"

         call mpi_bcast(SHI_MFP(j)%ELMFP(k)%dEdx, N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{SHI_MFP(j)%ELMFP(k)%dEdx}') ! module "MPI_subroutines"
      enddo ! k
   enddo ! j

#endif
end subroutine broadcast_SHI_MFP


subroutine broadcast_el_MFPs(temp_MFP, do_range, Total_el_MFPs, MPI_param)
   real(8), dimension(:,:), intent(inout) :: Temp_MFP
   logical, intent(inout) :: do_range
   type(All_MFP), dimension(:), intent(inout) :: Total_el_MFPs   ! electron mean free paths for all shells
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   !------------------------------------------
   integer :: Nat, Nshl, i, j, k, N1, N2
   character(100) :: error_part

#ifdef MPI_USED
   ! initialize:
   error_part = 'ERROR in broadcast_el_MFPs'

   call mpi_bcast(do_range, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{do_range}') ! module "MPI_subroutines"

   N1 = size(Temp_MFP,1)
   N2 = size(Temp_MFP,2)
   call mpi_bcast(Temp_MFP, N1*N2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{Temp_MFP}') ! module "MPI_subroutines"

   Nat = size(Total_el_MFPs)  ! how many atoms
   do j = 1, Nat   ! for all atoms
      Nshl = size(Total_el_MFPs(j)%ELMFP)    ! how mamy shells
      do k = 1, Nshl
         N1 = size(Total_el_MFPs(j)%ELMFP(k)%E)

         call mpi_bcast(Total_el_MFPs(j)%ELMFP(k)%E, N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{Total_el_MFPs(j)%ELMFP(k)%E}') ! module "MPI_subroutines"

         call mpi_bcast(Total_el_MFPs(j)%ELMFP(k)%L, N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{Total_el_MFPs(j)%ELMFP(k)%L}') ! module "MPI_subroutines"

         call mpi_bcast(Total_el_MFPs(j)%ELMFP(k)%dEdx, N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{Total_el_MFPs(j)%ELMFP(k)%dEdx}') ! module "MPI_subroutines"
      enddo ! k
   enddo ! j

#endif
end subroutine broadcast_el_MFPs


subroutine broadcast_el_aidCS_electrons(aidCS, MPI_param)
   type(All_diff_CS), intent(inout) :: aidCS    ! all integrated differential cross sections
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   !------------------------------------------
   integer :: Nat, Nshl, i, j, k, N1, N2, N_diff_siz
   character(100) :: error_part

#ifdef MPI_USED
   ! initialize:
   error_part = 'ERROR in broadcast_el_aidCS_electrons'

   Nat = size(aidCS%EIdCS)

   do j = 1, Nat   ! for all types of atoms
      Nshl = size(aidCS%EIdCS(j)%Int_diff_CS)    ! how mamy shells
      do k = 1, Nshl  ! for all orbitals

         N1 = size(aidCS%EIdCS(j)%Int_diff_CS(k)%E)

         call mpi_bcast(aidCS%EIdCS(j)%Int_diff_CS(k)%E, N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
         call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{aidCS%EIdCS(j)%Int_diff_CS(k)%E}') ! module "MPI_subroutines"

         do i = 1, N1 ! for all energy grid points
            if (MPI_param%process_rank == 0) then
               N_diff_siz = size(aidCS%EIdCS(j)%Int_diff_CS(k)%diffCS(i)%dsdhw)
            else
               N_diff_siz = 0
            endif
            call mpi_bcast(N_diff_siz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
            call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{N_diff_siz}') ! module "MPI_subroutines"
            !------------------------------------------------------
            ! Synchronize MPI processes to make sure the directory is created before going further
            call MPI_barrier_wrapper(MPI_param)  ! module "MPI_subroutines"
            !------------------------------------------------------

            if (MPI_param%process_rank /= 0) then
               if (.not.allocated(aidCS%EIdCS(j)%Int_diff_CS(k)%diffCS(i)%dsdhw)) then
                  allocate(aidCS%EIdCS(j)%Int_diff_CS(k)%diffCS(i)%dsdhw(N_diff_siz))
               endif
               if (.not.allocated(aidCS%EIdCS(j)%Int_diff_CS(k)%diffCS(i)%hw)) then
                  allocate(aidCS%EIdCS(j)%Int_diff_CS(k)%diffCS(i)%hw(N_diff_siz))
               endif
            endif

            call mpi_bcast(aidCS%EIdCS(j)%Int_diff_CS(k)%diffCS(i)%hw, N_diff_siz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
            call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{aidCS%EIdCS(j)%Int_diff_CS(k)%diffCS(i)%hw}') ! module "MPI_subroutines"

            call mpi_bcast(aidCS%EIdCS(j)%Int_diff_CS(k)%diffCS(i)%dsdhw, N_diff_siz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
            call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{aidCS%EIdCS(j)%Int_diff_CS(k)%diffCS(i)%dsdhw}') ! module "MPI_subroutines"
         enddo ! i
         !print*, '[MPI process #',MPI_param%process_rank, '] aidCS=', j, k, aidCS%EIdCS(j)%Int_diff_CS(k)%E
      enddo ! k
   enddo ! j
#endif
end subroutine broadcast_el_aidCS_electrons


subroutine broadcast_el_aidCS_holes(aidCS, MPI_param)
   type(All_diff_CS), intent(inout) :: aidCS    ! all integrated differential cross sections
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   !------------------------------------------
   integer :: Nsiz, Nshl, i, j, k, N1, N2, N_diff_siz
   character(100) :: error_part

#ifdef MPI_USED
   ! initialize:
   error_part = 'ERROR in broadcast_el_aidCS_holes'

   Nsiz = size(aidCS%HIdCS%E)

   call mpi_bcast(aidCS%HIdCS%E, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{aidCS%HIdCS%E}') ! module "MPI_subroutines"


   do i = 1, Nsiz   ! for all types of atoms
      if (MPI_param%process_rank == 0) then
         N_diff_siz = size(aidCS%HIdCS%diffCS(i)%hw)
      else
         N_diff_siz = 0
      endif
      call mpi_bcast(N_diff_siz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{N_diff_siz}') ! module "MPI_subroutines"
      !------------------------------------------------------
      ! Synchronize MPI processes to make sure the directory is created before going further
      call MPI_barrier_wrapper(MPI_param)  ! module "MPI_subroutines"
      !------------------------------------------------------
      if (MPI_param%process_rank /= 0) then
         if (.not.allocated(aidCS%HIdCS%diffCS(i)%hw)) then
            allocate(aidCS%HIdCS%diffCS(i)%hw(N_diff_siz))
         endif
         if (.not.allocated(aidCS%HIdCS%diffCS(i)%dsdhw)) then
            allocate(aidCS%HIdCS%diffCS(i)%dsdhw(N_diff_siz))
         endif
      endif
      call mpi_bcast(aidCS%HIdCS%diffCS(i)%hw, N_diff_siz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{aidCS%HIdCS%diffCS(i)%hw}') ! module "MPI_subroutines"

      call mpi_bcast(aidCS%HIdCS%diffCS(i)%dsdhw, N_diff_siz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{aidCS%HIdCS%diffCS(i)%dsdhw}') ! module "MPI_subroutines"

   enddo ! i

#endif
end subroutine broadcast_el_aidCS_holes




subroutine broadcast_Elastic_MFP(Elastic_MFP, MPI_param)
   type(MFP_elastic), intent(inout) :: Elastic_MFP         ! elastic mean free path of an electron
   type(Used_MPI_parameters), intent(inout) :: MPI_param
   !------------------------------------------
   integer :: Nsiz, Nshl, i, j, k, N1, N2, N_diff_siz
   character(100) :: error_part

#ifdef MPI_USED
   ! initialize:
   error_part = 'ERROR in broadcast_Elastic_MFP'

   Nsiz = size(Elastic_MFP%Total%E)

   call mpi_bcast(Elastic_MFP%Total%E, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{Elastic_MFP%Total%E}') ! module "MPI_subroutines"

   call mpi_bcast(Elastic_MFP%Total%L, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{Elastic_MFP%Total%L}') ! module "MPI_subroutines"

   call mpi_bcast(Elastic_MFP%Total%dEdx, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
   call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{Elastic_MFP%Total%dEdx}') ! module "MPI_subroutines"

   if (allocated(Elastic_MFP%Emit%E)) then
      call mpi_bcast(Elastic_MFP%Emit%E, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{Elastic_MFP%Emit%E}') ! module "MPI_subroutines"
   endif

   if (allocated(Elastic_MFP%Emit%L)) then
      call mpi_bcast(Elastic_MFP%Emit%L, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{Elastic_MFP%Emit%L}') ! module "MPI_subroutines"
   endif

   if (allocated(Elastic_MFP%Emit%dEdx)) then
      call mpi_bcast(Elastic_MFP%Emit%dEdx, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{Elastic_MFP%Emit%dEdx}') ! module "MPI_subroutines"
   endif

   if (allocated(Elastic_MFP%Absorb%E)) then
      call mpi_bcast(Elastic_MFP%Absorb%E, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{Elastic_MFP%Absorb%E}') ! module "MPI_subroutines"
   endif

   if (allocated(Elastic_MFP%Absorb%E)) then
      call mpi_bcast(Elastic_MFP%Absorb%L, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{Elastic_MFP%Absorb%L}') ! module "MPI_subroutines"
   endif

   if (allocated(Elastic_MFP%Absorb%dEdx)) then
      call mpi_bcast(Elastic_MFP%Absorb%dEdx, Nsiz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPI_param%ierror)  ! module "mpi"
      call MPI_error_wrapper(MPI_param%process_rank, MPI_param%ierror, trim(adjustl(error_part))//'{Elastic_MFP%Absorb%dEdx}') ! module "MPI_subroutines"
   endif
#endif
end subroutine broadcast_Elastic_MFP




subroutine get_num_shells(Target_atoms, Nshtot)
   type(Atom), dimension(:), intent(in), target :: Target_atoms  ! all data for target atoms
   integer, intent(out) :: Nshtot   ! total number of shells on all atoms
   integer Nat, Nshl, j, k
   Nshtot = 0
   Nat = size(Target_atoms)    ! how many atoms
   do j = 1, Nat
      Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
      do k = 1, Nshl
          Nshtot = Nshtot + 1
      enddo
   enddo
end subroutine get_num_shells


subroutine read_file(Reason, i, read_well, no_end, do_silent)
   integer, intent(in) :: Reason    ! file number where to read from
   integer, intent(inout) :: i      ! line number
   logical, intent(inout) :: read_well  ! did we read ok?
   logical, intent(in), optional :: no_end, do_silent ! ignore end of file warning (ocationally useful); don't printout anything
   !----------------------------
   logical :: ignor_end, be_silent

   if (present(no_end)) then
      ignor_end = no_end
   else
      ignor_end = .false.
   endif

   if (present(do_silent)) then
      be_silent = do_silent
   else
      be_silent = .false.
   endif

   i = i + 1    ! it's next line
   IF (Reason .GT. 0)  THEN ! ... something wrong ...
       if (.not.be_silent) write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', wrong type of variable'
       read_well = .false.
   ELSE IF (Reason .LT. 0) THEN ! ... end of file reached ...
       if ((.not.be_silent) .and. (.not.ignor_end)) write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', unexpected END of file'
       read_well = .false.
   ELSE   ! normal reading
       read_well = .true.  ! it read well, nothing to report
   END IF
end subroutine read_file

subroutine read_file_here(Reason, i, read_well)
   integer, intent(in) :: Reason    ! file number where to read from
   integer, intent(in) :: i      ! line number
   logical, intent(inout) :: read_well  ! did we read ok?
   IF (Reason .GT. 0)  THEN ! ... something wrong ...
       write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', wrong type of variable'
       read_well = .false.
   ELSE IF (Reason .LT. 0) THEN ! ... end of file reached ...
       write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', unexpected END of file'
       read_well = .false.
   ELSE   ! normal reading
       read_well = .true.  ! it read well, nothing to report
   END IF
end subroutine read_file_here


!=======================================================
! Useful subroutines to work with arrays:
subroutine Find_VB_numbers(Target_atoms, Lowest_Ip_At, Lowest_Ip_Shl)
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
    integer, intent(out) :: Lowest_Ip_At, Lowest_Ip_Shl
    real(8) temp
    integer i, j
    temp = 1d10
    do i = 1, size(Target_atoms)    ! all atomic spicies
       do j = 1, size(Target_atoms(i)%Shl_num)  ! all shells of each atom
          if ((Target_atoms(i)%Ip(j) .LT. temp) .OR. Target_atoms(i)%Shl_num(j) .GE. 63) then ! find the lowest Ip
             temp = Target_atoms(i)%Ip(j)
             Lowest_Ip_At = i
             Lowest_Ip_Shl = j
          endif
          if (Target_atoms(i)%Shl_num(j) .GE. 63) exit  ! if we found VB, that's the lowest Ip anyway
       enddo
    enddo
    !print*, Lowest_Ip_At, Lowest_Ip_Shl, Target_atoms(Lowest_Ip_At)%Ip(Lowest_Ip_Shl)
    !pause 'Find_VB_numbers'
end subroutine Find_VB_numbers

subroutine Linear_approx_2d(Array, In_val, Value1, El1, El2)
   REAL(8), dimension(:,:), INTENT(in) :: Array ! in this array make an approximation
   real(8), INTENT(in) :: In_val    ! this is the value
   real(8), intent(in), optional :: El1, El2    ! where to start approximation from, if needed
   real(8), intent(out) :: Value1   ! output, approximated value
   real(8) el_one
   integer Number
   
   call Find_in_array_monoton(Array, In_val, 1, Number)    ! find the closest value in the array to a given one
   
   if (Number .EQ. 1) then
       if (present(El1)) then  ! the starting points are known, use them:
          if (present(El2)) then ! the starting points are known, use them:
             Value1 = El2+(Array(2,Number)-El2)/(Array(1,Number)-El1)*(In_val - El1)
          else ! only the X-starting point is known, Y is not, use some assumptions:
             if (size(Array,2) .GE. 2) then ! array is long enough to assume something:
                if (Array(2,1) .GT. Array(2,1)) then ! it is locally decreasing array, assume starting point "infinity"
                    Value1 = 1d21
                else    ! it is locally increasing, assume starting point 'zero'
                    Value1 = Array(2,Number)/(Array(1,Number)-El1)*(In_val - El1)
                endif
             else   ! array is too short, no assumption can be made, just make it equal to the first value:
                Value1 = El1
             endif
          endif
       else ! no starting points are present, nothing to assume, use first point as approximation:
          Value1 = Array(2,1) ! [A] total mean free path
       endif
   else
      if (Array(2,Number-1) .GT. 1d20) then ! if it starts from infinity, approximate as 'infinity'
         Value1 = Array(2,Number-1)
      else  ! if it's normal array, just interpolate:
         Value1 = Array(2,Number-1)+(Array(2,Number)-Array(2,Number-1))/(Array(1,Number)-Array(1,Number-1))*(In_val - Array(1,Number-1))
      endif
   endif
end subroutine Linear_approx_2d


subroutine Linear_approx_2x1d(Array, Array2, In_val, Value1, El1, El2)
   REAL(8), dimension(:), INTENT(in) :: Array       ! Array correcp to In_val 
   REAL(8), dimension(:), INTENT(in) :: Array2      ! in this array make an ouput approximation
   real(8), INTENT(in) :: In_val    ! this is the value
   real(8), intent(in), optional :: El1, El2    ! where to start approximation from, if needed
   real(8), intent(out) :: Value1   ! output, approximated value
   real(8) el_one
   integer Number
   
   call Find_in_array_monoton(Array, In_val, Number)    ! find the closest value in the array to a given one
   
   if (Number .EQ. 1) then
       if (present(El1)) then  ! the starting points are known, use them:
          if (present(El2)) then ! the starting points are known, use them:
             Value1 = El2+(Array2(Number)-El2)/(Array(Number)-El1)*(In_val - El1)
          else ! only the X-starting point is known, Y is not, use some assumptions:
             if (size(Array2) .GE. 2) then ! array is long enough to assume something:
                if (Array2(1) .GT. Array2(2)) then ! it is locally decreasing array, assume starting point "infinity"
                    Value1 = 1d21
                else    ! it is locally increasing, assume starting point 'zero'
                    Value1 = Array2(Number)/(Array(Number)-El1)*(In_val - El1)
                endif
             else   ! array is too short, no assumption can be made, just make it equal to the first value:
                Value1 = El1
             endif
          endif
       else ! no starting points are present, nothing to assume, use first point as approximation:
          Value1 = Array2(1) ! [A] total mean free path
       endif
   else
      if (Array2(Number-1) .GT. 1d20) then ! if it starts from infinity, approximate as 'infinity'
         Value1 = Array2(Number-1)
      else  ! if it's normal array, just interpolate:
         Value1 = Array2(Number-1)+(Array2(Number)-Array2(Number-1))/(Array(Number)-Array(Number-1))*(In_val - Array(Number-1))
      endif
   endif
end subroutine Linear_approx_2x1d

subroutine Find_in_1D_array(Array, Value, Number)
   REAL(8), dimension(:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value   ! to be found in the array as near as possible
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i
   i = 1
   do while (Array(i) .LT. Value)
      i = i + 1
   enddo
   Number = i
end subroutine Find_in_1D_array

subroutine Find_in_2D_array(Array, Value, Indx, Number)
   REAL(8), dimension(:,:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value   ! to be found in the array as near as possible
   integer, INTENT(in) :: Indx    ! index of the array, showing in which colonm we search
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i
   i = 1
   do while (Array(Indx,i) .LT. Value)
      i = i + 1
   enddo
   Number = i
end subroutine Find_in_2D_array




subroutine Find_in_monoton_array_decreasing(Array, Value0, Number)
   REAL(8), dimension(:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value0   ! to be found in the array as near as possible
   integer, INTENT(out) :: Number ! number of the element which we are looking for
   !----------------------
   integer :: Nsiz, i, N, i_cur, i_1, i_2, coun
   real(8) :: temp_val, val_1, val_2

   Nsiz = size(Array)

   i_1 = 1     ! to start with
   i_2 = Nsiz  ! to start with
   val_1 = Array(i_1)   ! to start with
   val_2 = Array(i_2)   ! to start with

   i_cur = FLOOR((i_1+i_2)/2.0)
   temp_val = Array(i_cur)

   if (isnan(Value0)) then
       print*, 'The subroutine Find_in_monoton_array_decreasing'
       print*, 'cannot proceed, the value of Value0 is', Value0
       write(*, '(f,f,f,f)') Value0, Array(i_cur), Array(i_1), Array(i_2)
       pause 'STOPPED WORKING...'
   else
       if (Value0 < Array(Nsiz)) then ! smaller than the last value
           i_cur = Nsiz
       elseif (Value0 > Array(1)) then ! bigger than the first value
           i_cur = 1
       else
           coun = 0
           do while ( abs(i_1 - i_2) > 1) ! until find the interval where the value is in
               if (temp_val > Value0) then
                   i_1 = i_cur
                   val_1 = Array(i_1)
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(i_cur)
                else
                   i_2 = i_cur
                   val_2 = temp_val
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(i_cur)
                endif
                coun = coun + 1
                if (coun > 1e3) then
                    print*, 'PROBLEM WITH CONVERGANCE IN'
                    print*, 'Find_in_monoton_array_decreasing', coun
                    write(*, '(f,f,f,f)') Value0, Array(i_cur), Array(i_1), Array(i_2)
                    pause 'STOPPED WORKING...'
                endif
           enddo
       endif
   endif
   Number = i_cur
end subroutine Find_in_monoton_array_decreasing



subroutine Find_in_monotonous_1D_array(Array, Value0, Number, text_to_print)
   REAL(8), dimension(:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value0   ! to be found in the array as near as possible
   integer, INTENT(out) :: Number ! number of the element which we are looking for
   character(*), intent(in), optional :: text_to_print
   !-----------------------------
   integer i, N, i_cur, i_1, i_2, coun
   real(8) temp_val, val_1, val_2
   character(200) :: error_text

   if (present(text_to_print)) then
      error_text = text_to_print
   else
      error_text = ''
   endif
   
   N = size(Array)
   i_1 = 1
   val_1 = Array(i_1)
   i_2 = N
   val_2 = Array(i_2)
   i_cur = FLOOR((i_1+i_2)/2.0)
   temp_val = Array(i_cur)
   if (isnan(Value0)) then
        print*, 'The subroutine Find_in_monotonous_1D_array cannot proceed, the value of Value0 is', Value0
        if (LEN(trim(adjustl(error_text))) > 0) print*, trim(adjustl(error_text))
        write(*, '(f,f,f,f)') Value0, Array(i_cur), Array(i_1), Array(i_2)
        !pause 'STOPPED WORKING...'
        stop
   else
       if (Value0 .LT. Array(1)) then ! it's the first value, no need to search
           i_cur = 0
       else if (Value0 .GE. Array(N)) then ! it's the last value, no need to search
           i_cur = N-1
       else
           coun = 0
           do ! until the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                !if ((Value0 .GE. Array(i_cur)) .AND. (Value0 .LE. Array(i_cur+1))) exit ! when the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if (i_1 == i_2-1) exit
                if (temp_val .LE. Value0) then
                   i_1 = i_cur
                   val_1 = Array(i_1)
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(i_cur)
                else
                   i_2 = i_cur
                   val_2 = temp_val
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(i_cur)
                endif
                coun = coun + 1
                if (coun .GT. 1e3) then
                    print*, 'PROBLEM WITH CONVERGANCE IN'
                    print*, 'Find_in_monotonous_1D_array', coun
                    write(*, '(f,f,f,f)') Value0, Array(i_cur), Array(i_1), Array(i_2)
                    pause 'STOPPED WORKING...'
                endif
           enddo
       endif
   endif    ! isnan
   Number = i_cur+1
end subroutine Find_in_monotonous_1D_array

subroutine Find_in_monotonous_2D_array(Array, Value0, Indx, Number, text_to_print)
   REAL(8), dimension(:,:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value0   ! to be found in the array as near as possible
   integer, INTENT(in) :: Indx    ! index of the array, showing in which colonm we search
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   character(*), intent(in), optional :: text_to_print
   !-----------------------------
   integer i, N, i_cur, i_1, i_2, coun
   real(8) temp_val, val_1, val_2
   character(200) :: error_text
   
   if (present(text_to_print)) then
      error_text = text_to_print
   else
      error_text = ''
   endif

   N = size(Array,2)
   i_1 = 1
   val_1 = Array(Indx,i_1)
   i_2 = N
   val_2 = Array(Indx,i_2)
   i_cur = FLOOR((i_1+i_2)/2.0)
   temp_val = Array(Indx,i_cur)
   
   if (isnan(Value0)) then
        write(*, '(a,e)') 'The subroutine Find_in_monotonous_2D_array got Value0 =', Value0
        if (LEN(trim(adjustl(error_text))) > 0) write(*, '(a)') trim(adjustl(error_text))
        write(*, '(f,f,f,f)') Value0, Array(Indx,i_cur), Array(Indx,i_1), Array(Indx,i_2)
        !pause 'STOPPED WORKING...'
        stop
   else
       if (Value0 .LT. Array(Indx,1)) then ! it's the first value, no need to search
           i_cur = 0
       else if (Value0 .GE. Array(Indx,N)) then ! it's the last value, no need to search
           i_cur = N-1
       else
           coun = 0
           do ! until the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if ((Value0 .GE. Array(Indx,i_cur)) .AND. (Value0 .LE. Array(Indx,i_cur+1))) exit ! when the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if (temp_val .LE. Value0) then
                   i_1 = i_cur
                   val_1 = temp_val
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(Indx,i_cur)
                else
                   i_2 = i_cur
                   val_2 = temp_val
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(Indx,i_cur)
                endif
                
                coun = coun + 1
                if (coun .GT. 1e3) then
                    print*, 'PROBLEM WITH CONVERGANCE IN'
                    print*, 'Find_in_monotonous_2D_array', coun
                    write(*, '(f,f,f,f)') Value0, Array(Indx,i_cur), Array(Indx,i_1), Array(Indx,i_2)
                    pause 'STOPPED WORKING...'
                endif
           enddo
       endif
   endif    ! isnan
   Number = i_cur+1
end subroutine Find_in_monotonous_2D_array

subroutine Integrate_function_one(Int_type, x, f, x0, xn, res, Error_message, MPI_param)
    integer, intent(in) :: Int_type ! type of integration to be used: 0=trapeziod, 1=Simpson-3/8, 2=...
    real(8), dimension(:), intent(in) :: x  ! grid points
    real(8), dimension(:), intent(in) :: f  ! function
    real(8), intent(in) :: x0   ! starting point of integration
    real(8), intent(in) :: xn   ! ending point of integration
    real(8), intent(out) :: res ! result of integration
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    type(Used_MPI_parameters), intent(inout) :: MPI_param
    !------------------------------------------
    real(8) temp_x0, temp_xn
    integer i, N, Nx, Nf, i_0, i_n
    character(100) Err_data
    Nx = size(x)
    Nf = size(f)
    temp_x0 = x0
    temp_xn = xn
    if (Nx .NE. Nf) then
        Err_data = 'Trapeziod integration failed, size of x is not equal to size of f' ! no input data found
        call Save_error_details(Error_message, 3, Err_data, MPI_param)
    else if (temp_x0 .GT. x(Nx)) then
        Err_data = 'Trapeziod integration failed, starting point is larger than x(size(x))' ! no input data found
        call Save_error_details(Error_message, 4, Err_data, MPI_param)
    else if (temp_xn .GT. x(Nx)) then
        Err_data = 'Trapeziod integration failed, ending point is larger than x(size(x))' ! no input data found
        call Save_error_details(Error_message, 5, Err_data, MPI_param)
    else if (x0 .EQ. xn) then
        res = 0.0d0
    else    ! everything seems to be fine...
        if (temp_x0 .LT. x(1)) then
            temp_x0 = x(1)   ! start from the first point
            i_0 = 1     ! it's number is, obviously, 1
        else
            call Find_in_array_monoton(x, temp_x0, i_0)  ! find from where to start
            i_0 = i_0 - 1
        endif ! -> i_0
        
        if (temp_xn .GT. x(Nx)) then
            temp_xn = x(Nx)   ! end at the last point
            i_n = Nx     ! it's number is, obviously, Nx
        else
            call Find_in_array_monoton(x, temp_xn, i_n)  ! find from where to start
        endif ! -> i_n
        
        !write(*,'(a,i,i)') 'Int:', i_0, i_n
        select case(Int_type)   ! for different types of integration:
        case default
            call Trapeziod(x,f,temp_x0,temp_xn,i_0,i_n,res) ! here is the method of integration
        endselect
        !pause 'TESTING INTEGRATION'
    endif
end subroutine 


subroutine Trapeziod_one(x,f,x0,xn,i_0,i_n,res)
    real(8), dimension(:), intent(in) :: x, f   ! grid and function
    real(8), intent(in) :: x0,xn    ! starting and ending points
    integer, intent(in) :: i_0, i_n ! number of starting and ending points in the array
    real(8), intent(out) :: res     ! result of integration
    real(8) temp
    integer i
    res = 0.0d0
    do i = i_0, i_n-1 ! integration in this limits
        if (i .EQ. i_0) then
            call Linear_approx(x, f, x0, temp)
            res = res + (temp+f(i+1))/2.0d0*(x(i+1)-x0)
        else if (i .EQ. i_n-1) then
            call Linear_approx(x, f, xn, temp)
            res = res + (f(i)+temp)/2.0d0*(xn-x(i))
        else
            res = res + (f(i)+f(i+1))/2.0d0*(x(i+1)-x(i))
        endif
        !print*, i, res
    enddo
end subroutine




subroutine Integrate_function_save(Int_type, x, f, x0, xn, res, Error_message, MPI_param)
    integer, intent(in) :: Int_type ! type of integration to be used: 0=trapeziod, 1=Simpson-3/8, 2=...
    real(8), dimension(:), intent(in) :: x  ! grid points
    real(8), dimension(:), intent(in) :: f  ! function
    real(8), intent(in) :: x0   ! starting point of integration
    real(8), intent(in) :: xn   ! ending point of integration
    real(8), dimension(:), intent(out) :: res ! result of integration saved for each grid-point
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    type(Used_MPI_parameters), intent(inout) :: MPI_param
    !----------------------------------
    real(8) temp_x0, temp_xn
    integer i, N, Nx, Nf, i_0, i_n
    character(100) Err_data
    Nx = size(x)
    Nf = size(f)
    temp_x0 = x0
    temp_xn = xn
    if (Nx .NE. Nf) then
        Err_data = 'Trapeziod integration failed, size of x is not equal to size of f' ! no input data found
        call Save_error_details(Error_message, 3, Err_data, MPI_param)
    else if (temp_x0 .GT. x(Nx)) then
        Err_data = 'Trapeziod integration failed, starting point is larger than x(size(x))' ! no input data found
        call Save_error_details(Error_message, 4, Err_data, MPI_param)
    else if (temp_xn .LT. temp_x0) then
        Err_data = 'Trapeziod integration failed, ending point is the starting one' ! no input data found
        call Save_error_details(Error_message, 5, Err_data, MPI_param)
    else if (x0 .EQ. xn) then
        res = 0.0d0
    else    ! everything seems to be fine...
        if (temp_x0 .LT. x(1)) then
            temp_x0 = x(1)   ! start from the first point
            i_0 = 1     ! it's number is, obviously, 1
        else
            call Find_in_array_monoton(x, temp_x0, i_0)  ! find from where to start
            i_0 = i_0 - 1
        endif ! -> i_0
        
        if (temp_xn .GT. x(Nx)) then
            temp_xn = x(Nx)   ! end at the last point
            i_n = Nx     ! it's number is, obviously, Nx
        else
            call Find_in_array_monoton(x, temp_xn, i_n)  ! find from where to start
        endif ! -> i_n
        
        !write(*,'(a,i,i)') 'Int:', i_0, i_n
        select case(Int_type)   ! for different types of integration:
        case default
            call Trapeziod_save(x,f,temp_x0,temp_xn,i_0,i_n,res) ! here is the method of integration
        endselect
        !pause 'TESTING INTEGRATION'
    endif
end subroutine


subroutine Trapeziod_save(x,f,x0,xn,i_0,i_n,res)
    real(8), dimension(:), intent(in) :: x, f   ! grid and function
    real(8), intent(in) :: x0,xn    ! starting and ending points
    integer, intent(in) :: i_0, i_n ! number of starting and ending points in the array
    real(8), dimension(:), intent(out) :: res     ! result of integration
    real(8) temp
    integer i
    res = 0.0d0
    do i = i_0, i_n-1 ! integration in this limits
        if (i .EQ. i_0) then
            call Linear_approx(x, f, x0, temp)
            res(i) = (temp+f(i+1))/2.0d0*(x(i+1)-x0)
        else if (i .EQ. i_n-1) then
            call Linear_approx(x, f, xn, temp)
            res(i) = res(i-1) + (f(i)+temp)/2.0d0*(xn-x(i))
        else
            res(i) = res(i-1) + (f(i)+f(i+1))/2.0d0*(x(i+1)-x(i))
        endif
        !print*, i, res(i)
    enddo
    res(i_n) = res(i_n-1)
end subroutine


subroutine Linear_approx_2x1d_DSF(Array_dL, Array_dE, In_val, Value1)
   REAL(8), dimension(:), INTENT(in) :: Array_dL ! Ls (for two particle energies on the grid)
   REAL(8), dimension(:), INTENT(in) :: Array_dE ! dE (for two particle energies on the grid)
   real(8), INTENT(in) :: In_val    ! Sampled dL
   real(8), intent(out) :: Value1   ! interpolater output value
   !----------------
   integer :: Number, Nsiz

   Nsiz = size(Array_dL)

   ! The following subroutine works for monotoneusly increasing array,
   ! whereas the Array passed is decreaseing, so use negative array to find the correct index:
   call Find_in_array_monoton(-Array_dL, -In_val, Number)    ! find the closest value in the array to a given one

   if (Number == 1) then
      Value1 = Array_dE(Number)+(Array_dE(Number+1)-Array_dE(Number)) / &
                    (Array_dL(Number+1)-Array_dL(Number))*(In_val - Array_dL(Number))
   elseif ( abs(Array_dL(Number)-Array_dL(Number-1)) < 1.0d-9 ) then
      Value1 = Array_dE(Number-1)
   else
      if (Array_dL(Number-1) > 1d20) then ! if it starts from infinity, approximate as 'infinity'
         Value1 = Array_dE(Number-1)
      else  ! if it's normal array, just interpolate:
         Value1 = Array_dE(Number-1)+(Array_dE(Number)-Array_dE(Number-1)) / &
                    (Array_dL(Number)-Array_dL(Number-1))*(In_val - Array_dL(Number-1))
      endif
   endif
end subroutine Linear_approx_2x1d_DSF



subroutine Linear_approx_2x1d_DSF_OLD(Array, Array2, In_val, Value1, El1, El2)
   REAL(8), dimension(:), INTENT(in) :: Array       ! Array corresp to In_val
   REAL(8), dimension(:), INTENT(in) :: Array2      ! in this array make an ouput approximation
   real(8), INTENT(in) :: In_val    ! this is the value
   real(8), intent(in), optional :: El1, El2    ! where to start approximation from, if needed
   real(8), intent(out) :: Value1   ! output, approximated value
   real(8) el_one
   integer Number, Number1, i, kk

   kk = size(array)
!     i = 1
!     do while (Array(i) .GT. In_val) ! testing
!          if (i .EQ. kk) exit
!          i = i+1
!     enddo
!     Number1 = i
   !print*, Number1

   ! The following subroutine works for monotoneusly increasing array,
   ! where as Array passed is decreaseing, so use negative array to find the correct index:
   call Find_in_array_monoton(-Array, -In_val, Number)    ! find the closest value in the array to a given one
!     if (Number1 /= Number) then ! testing
!         print*, 'Linear_approx_2x1d_DSF', Number1, Number
!         pause
!     endif

   if (Number .EQ. 1) then
       if (present(El1)) then  ! the starting points are known, use them:
          if (present(El2)) then ! the starting points are known, use them:
             Value1 = El2+(Array2(Number)-El2)/(Array(Number)-El1)*(In_val - El1)
          else ! only the X-starting point is known, Y is not, use some assumptions:
             if (size(Array2) .GE. 2) then ! array is long enough to assume something:
                if (Array2(1) .GT. Array2(2)) then ! it is locally decreasing array, assume starting point "infinity"
                    Value1 = 1d21
                else    ! it is locally increasing, assume starting point 'zero'
                    Value1 = Array2(Number)/(Array(Number)-El1)*(In_val - El1)
                endif
             else   ! array is too short, no assumption can be made, just make it equal to the first value:
                Value1 = El1
             endif
          endif
       else ! no starting points are present, nothing to assume, use first point as approximation:
          Value1 = Array2(1) ! [A] total mean free path
       endif
   else if (Number .EQ. kk) then
      Value1 = Array2(Number)
   else
      if (Array2(Number-1) .GT. 1d20) then ! if it starts from infinity, approximate as 'infinity'
         Value1 = Array2(Number-1)
      else  ! if it's normal array, just interpolate:
         Value1 = Array2(Number-1)+(Array2(Number)-Array2(Number-1))/(Array(Number)-Array(Number-1))*(In_val - Array(Number-1))
      endif
   endif
!   if (Value1 .GT. 1.0) then
!        print*, 'Linear approx subr = '
!        print*, Number, kk, Value1
!        print*, array(Number), array2(Number)
!        pause
!   endif
end subroutine Linear_approx_2x1d_DSF_OLD


end module Reading_files_and_parameters

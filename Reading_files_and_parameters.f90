!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains subroutines for reading and interpreting input files:

module Reading_files_and_parameters
  use Universal_Constants   ! let it use universal constants
  use Objects   ! since it uses derived types, it must know about them from module 'Objects'
  use Dealing_with_EADL ! all subroutines dealing with EADL and EPDL97 are there
  
  implicit none

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

!private  ! hides items not listed on public statement 
public :: Find_in_array, Find_in_array_monoton, Linear_approx
public :: Integrate_function, Trapeziod

contains

subroutine Read_input_file(Target_atoms, CDF_Phonon, Matter, Mat_DOS, SHI, Tim, dt, Output_path, Output_path_SHI, &
           Material_name, NMC, Num_th, Error_message, read_well, DSF_DEMFP, DSF_DEMFP_H, NumPar, File_names)
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
   
   real(8) M, SUM_DOS
   integer FN, Reason, i, temp, temp1, temp2(1), temper
   character(100)   Material_file   ! file with material parameters, name MUST coinside with 'Material_name'!!!
   character(100)   Short_material_file ! short version of the 'material parameters file'
   character(100)   File_name   ! file name if needed to use
   character(100)   DOS_file    ! file with material DOS, name MUST coinside with 'Material_name'!!!
   character(100)   DSF_file, DSF_file_h    ! file with DSF differential cross-sections, name MUST coinside with 'Material_name'!!!
   character(100)   Error_descript  ! to write a description of an error, if any
   character(100)   Temp_char, temp_char1
   character(3) Name    ! for reading elements names
   character(30) Full_Name    ! for reading full elements names
   character(100) command   ! to pass to cmd a command
   logical file_exist    ! to check where file to be open exists
   logical file_opened   ! to check if a file is still opened

   FN = 200
   inquire(file='INPUT_PARAMETERS.txt',exist=file_exist)     ! check if input file excists
   if (file_exist) then
      open(unit = FN, FILE = 'INPUT_PARAMETERS.txt', status = 'old', readonly)   ! if yes, open it and read
   else ! if no, save error message about it:
      File_name = 'OUTPUT_Error_log.dat'
      open(unit = 100, FILE = trim(adjustl(File_name)))
      Error_message%File_Num = 100	! file number with error-log
      Error_descript = 'File INPUT_PARAMETERS.txt is not found!'    ! description of an error
      call Save_error_details(Error_message, 1, Error_descript) ! write it into the error-log file
      print*, trim(adjustl(Error_descript)) ! print it also on the sreen
      read_well = .false.   ! it didn't read well the input file...
      goto 2013 ! go to the end of the subroutine, there is nothing else we could do
   endif
   
   i = 0
   READ(FN,*,IOSTAT=Reason) Material_name
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013
   Material_file = 'INPUT_CDF/'//trim(adjustl(Material_name))//'.cdf' ! that's the file where material properties must be storred, or alternatively:
   Short_material_file = 'INPUT_EADL/'//trim(adjustl(Material_name))//'.scdf' ! that's the file where short version of material properties must be storred
   DOS_file = 'INPUT_DOS/'//trim(adjustl(Material_name))//'.dos'      ! that's the file where material DOS must be storred!
   !//trim(adjustl(Material_name))//'_Electron_DSF_Differential_EMFPs.dat'     ! that's the file where DSF differential EMFP for this material must be storred!
   
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
   DSF_file = 'INPUT_DSF/'//trim(adjustl(Material_name))//'/'//trim(adjustl(Material_name))//'_Electron_DSF_Differential_EMFPs_'//trim(adjustl(temp_char1))//'K.dat'     ! that's the file where DSF differential EMFP for this material must be storred!
   DSF_file_h = 'INPUT_DSF/'//trim(adjustl(Material_name))//'/'//trim(adjustl(Material_name))//'_Hole_DSF_Differential_EMFPs_'//trim(adjustl(temp_char1))//'K.dat'     ! that's the file where DSF differential EMFP for this material must be storred!

   
   READ(FN,*,IOSTAT=Reason) SHI%Kind_Zeff, SHI%fixed_Zeff   ! 0=Barkas; 1=Bohr; 2=Nikolaev-Dmitriev; 3=Schiwietz-Grande, 4 - fixed value;
   call read_file(Reason, i, read_well) ! reports if everything read well
   if ((SHI%fixed_Zeff .LE. 0.0) .OR. (SHI%fixed_Zeff .GT. SHI%Zat)) SHI%fixed_Zeff = SHI%Zat
   if (.not. read_well) goto 2013
   
   READ(FN,*,IOSTAT=Reason) SHI%Kind_ion   ! 0=point-like charge; 1=Brandt-Kitagawa ion
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2013
   
   ! kind of elastic cross-section: -1=no elastic, 0=Mott; 1=CDF phonons; 2=DSF
   ! and kind of effective target charge (used in CDF): 0 = effective Barkas-like; 1=fixed value=1
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
   
   ! This option is not ready, so it is excluded from release:
!    READ(FN,'(a)',IOSTAT=Reason) Temp_char  ! path to Gnuplot installed
!    call read_file(Reason, i, read_well) ! reports if everything read well
!    if (.not. read_well) goto 2013
   Temp_char = '0'
   if (trim(adjustl(Temp_char)) .EQ. '0') then ! no gmuplot installed, no need to create scripts
      print*, 'No Gnuplot script will be created'
   else 
      allocate(File_names%F(10))
      File_names%F(1) = Temp_char
      print*, 'Gnuplot scripts will be created, calling Gnuplot from'
      print*, trim(adjustl(File_names%F(1)))
   endif
   
   !------------------------------------------------------
   ! Create an output folder:
   Output_path = 'OUTPUT_'//trim(adjustl(Material_name)) ! that should be a folder with output
   if (SHI%Zat .GT. 0) Output_path_SHI = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//'OUTPUT_'//trim(adjustl(SHI%Name))//'_in_'//trim(adjustl(Material_name)) ! that should be a folder with output
   inquire(DIRECTORY=trim(adjustl(Output_path)),exist=file_exist)    ! check if input file excists
   if (file_exist) then
      write(*,'(a,a,a)') 'Folder ', trim(adjustl(Output_path)), ' already exists, save output files into it'
   else
      write(*,'(a,a,a)') 'Folder ', trim(adjustl(Output_path)), ' was created, save output files into it'
      command='mkdir '//trim(adjustl(Output_path)) ! to create a folder use this command
      CALL system(command)  ! create the folder
   endif

   Error_message%File_Num = 100 ! this is the file where error-log will be storred
   if (SHI%Zat .GT. 0) then
    File_name = trim(adjustl(Output_path))//'/'//trim(adjustl(SHI%Name))//'_in_'//trim(adjustl(Material_name))//'_error_log.txt'
   else
    File_name = trim(adjustl(Output_path))//'/'//trim(adjustl(Material_name))//'_error_log.txt'
   endif
   open(unit = Error_message%File_Num, FILE = trim(adjustl(File_name)))
   
   call reading_material_parameters(Material_file, Short_material_file, Target_atoms, NumPar, CDF_Phonon, Matter, Error_message, read_well)    ! read material parameters
   if (.not. read_well) goto 2015
   
   do i = 1, size(Target_atoms) ! all atoms
       temp2(1) = maxval(Target_atoms(i)%KOCS_SHI(:)) ! check if there is at least one shell for which we use BEB instead of CDF
       if (temp2(1) .GE. 2) then ! BEB cross sections are used, then:
         print*, 'Cannot use Brandt-Kitagawa model with BEB cross section, using point charge instead'
         SHI%Kind_ion = 0
         exit
      endif
   enddo
   
   if (SHI%Zat .GT. 0) then ! if there is an SHI:
       ! If SHI energy is too small:
       if (SHI%E .LE. (SHI%Mass*g_Mp + g_me)*(SHI%Mass*g_Mp + g_me)/(SHI%Mass*g_Mp*g_me)*Target_atoms(1)%Ip(size(Target_atoms(1)%Ip))/4.0d0) then
           write(*,'(a,e,a,e,a)') 'The SHI energy is ', SHI%E, ' [eV] which is smaller than the minimum allowed energy of ', (SHI%Mass*g_Mp + g_me)*(SHI%Mass*g_Mp + g_me)/(SHI%Mass*g_Mp*g_me)*Target_atoms(1)%Ip(size(Target_atoms(1)%Ip))/4.0d0, ' [eV]'
           print*, 'Calculations cannot be performed, sorry. Try to increase the ion energy.'
           read_well = .false.  ! under such conditions, calculations cannot be performed
       endif
       if (SHI%E .LE. (SHI%Mass*g_Mp + g_me)*(SHI%Mass*g_Mp + g_me)/(SHI%Mass*g_Mp*g_me)*Target_atoms(1)%Ip(size(Target_atoms(1)%Ip))/4.0d0) goto 2015
       ! if SHI energy is too high (relativistic):
       if (SHI%E .GE. 175.0d6/2.0d0*SHI%Mass) then
           write(*,'(a,f17.3,a,f17.3,a)') 'The SHI energy is ', SHI%E, ' [eV] which is higher than the maximum allowed energy of ', 175.0d6/2.0d0*SHI%Mass, ' [eV]'
           print*, 'Calculations cannot be performed, sorry. Try to decrease the ion energy.'
           read_well = .false.  ! under such conditions, calculations cannot be performed
       else if (SHI%E .GE. 175.0d6/4.0d0*SHI%Mass) then
           write(*,'(a,f17.3,a,f17.3,a)') 'The SHI energy is ', SHI%E, ' [eV] which is higher than the energy of ', 17.50d6*SHI%Mass, ' [eV]'
           print*, 'Note that the calculations might not be very reliable, since NO relativistic effects are included!'
       else
       endif
       if (SHI%E .GE. 175.0d6/2.0d0*SHI%Mass) goto 2015
   endif

   call reading_material_DOS(DOS_file, Mat_DOS, Matter, Target_atoms, Error_message, read_well)    ! read material DOS
   if (.not. read_well) goto 2015
         
   open(25, file = trim(adjustl(Output_path))//'/'//'OUTPUT_'//trim(adjustl(Material_name))//'_DOS_analysis.dat')
   write(25, *) 'E[eV]   k[1/m]   DOS[a.u.]    Int_DOS[a.u.]   Eff_mass[me]'
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
     
   if (NumPar%kind_of_EMFP .EQ. 2) then
       call reading_DSF_cross_sections(DSF_file, DSF_DEMFP, NumPar, Error_message, read_well)
       call reading_DSF_cross_sections(DSF_file_h, DSF_DEMFP_H, NumPar, Error_message, read_well)
   else
      allocate(DSF_DEMFP(0))
   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Printing out the parameters:
   print*, '-------------------------------------------------------'
   if (SHI%Zat .GT. 0) then
       write(*,'(a,a,a,a)') 'Performing calculations for ', trim(adjustl(SHI%Full_Name)), ' in ', trim(adjustl(Material_name))
       write(*,'(a,a,a,f9.2,a)') 'where ', trim(adjustl(SHI%Name)), '-ion has energy ', SHI%E/1d6, ' [MeV].'
       write(*,'(a,f6.2,a)') 'The total time to be analyzed is ', Tim, ' [fs]'
       select case (NumPar%dt_flag)
        case (:0)
            write(*,'(a,f6.2,a)') 'with the timestep of ', dt, ' [fs], linear time-scale'
        case (1:)
            write(*,'(a,f6.2,a)') 'starting with 0.01 [fs] increasing by ', dt, ' in logarithmic time-scale'
       endselect
   else
       write(*,'(a,a,a,a)') 'Performing calculations of the electron MFP in ', trim(adjustl(Material_name)), ' only.'
   endif
   print*, '-------------------------------------------------------'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
2013 if (.not. read_well) print*, 'Error in INPUT_PARAMETERS.txt file or inrut files. See log!!'
2015 continue   ! if we must skip to the end for some reason
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it
end subroutine Read_input_file


subroutine reading_material_parameters(Material_file, Short_material_file, Target_atoms, NumPar, CDF_Phonon, Matter, Error_message, read_well)
   character(100), intent(in) :: Material_file, Short_material_file  ! files with material parameters (full or short)
   type(Atom), dimension(:), allocatable, intent(inout) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   type(Flag), intent(inout) :: NumPar ! numerical parameters
   type(CDF), intent(inout) :: CDF_Phonon   ! CDF parameters for phonon to be read from a file
   type(Solid), intent(inout) :: Matter   ! all material parameters
   type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
   logical, intent(inout) :: read_well  ! did we read the file well?
   
   real(8) M
   integer FN2, Reason, i, j, k, l, N, Shl, CDF_coef, Shl_num
   character(100) Error_descript, temp
   character(3) Name
   character(11) Shell_name
   character(30) Full_Name
   logical file_opened, file_exist, file_exist2
   
   FN2 = 201
   inquire(file=trim(adjustl(Material_file)),exist=file_exist)    ! check if the full input file exists
   inquire(file=trim(adjustl(Short_material_file)),exist=file_exist2)    ! check if the short input file exists
   if (file_exist) then
      open(unit = FN2, FILE = trim(adjustl(Material_file)), status = 'old', readonly)   ! yes, open full file and read
   else if (file_exist2) then ! if at least short version exists, the rest can be used within atomic approximation
      open(unit=FN2, FILE = trim(adjustl(Short_material_file)), status = 'old', readonly)   ! yes, open full file and read
      call read_short_scdf(FN2, Target_atoms, NumPar, CDF_Phonon, Matter, Error_message, read_well)
      goto 2014 ! short version is done, skip reading the long one below
   else ! if no, save error message about it:
      Error_descript = 'File '//trim(adjustl(Material_file))//' is not found!'    ! description of an error
      call Save_error_details(Error_message, 2, Error_descript) ! write it into the error-log file
      print*, trim(adjustl(Error_descript)) ! print it also on the sreen
      read_well = .false.   ! it didn't read well the input file...
      goto 2014 ! go to the end of the subroutine, there is nothing else we could do
   endif
   
   i = 0
   READ(FN2,'(a100)',IOSTAT=Reason) Matter%Target_name ! first line is the full material name
   
   READ(FN2,*,IOSTAT=Reason) N   ! number of elements in this compound
   IF (Reason .GT. 0)  THEN ! if it's a formula:
      backspace(FN2) ! read this line again:
      READ(FN2,*,IOSTAT=Reason) temp   ! chemical formula of this compound
      call read_file(Reason, i, read_well) ! reports if everything read well
      if (.not. read_well) goto 2014
      call Decompose_compound(temp, Numpar%path_sep, Target_atoms, read_well) ! from module "Dealing_with_EADL"
      N = size(Target_atoms) ! that's how many atoms we have
      print*, 'The material formula ', trim(adjustl(temp)), ' was interpreted as:'
      do j = 1, N  ! read for each element it's basic data:
         write(*,'(a,a,a,i3,a,f5.1,a,f7.2)') 'Atom ', Target_atoms(j)%Name, ' number #', Target_atoms(j)%Zat, ', composition: ', Target_atoms(j)%Pers, ', mass', Target_atoms(j)%Mass
      enddo
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
      enddo
   endif

   READ(FN2,*,IOSTAT=Reason) Matter%Dens, Matter%Vsound, Matter%E_F   ! material density [g/cm^3], speed of sound in the material [m/s], Fermi energy []
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2014
   !Matter%At_Dens = 1.0d-3*Matter%Dens/(g_Mp*SUM(Target_atoms(:)%Mass)/size(Target_atoms)) ! [1/cm^3] atomic density
   Matter%At_Dens = 1.0d-3*Matter%Dens/(g_Mp*SUM(Target_atoms(:)%Mass*Target_atoms(:)%Pers)/SUM(Target_atoms(:)%Pers))
   Matter%v_f = sqrt(2.0d0*Matter%E_f/g_me) ! units as used in one-pole approximation (not SI!)
  
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
          call check_atomic_parameters(NumPar, Target_atoms, j, k, shl, Error_message, read_well) ! from module 'Dealing_with_EADL'

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
      enddo
   enddo
   
   READ(FN2,*,IOSTAT=Reason) CDF_coef
   i = i + 1
   IF (Reason .GT. 0)  THEN ! ... something wrong ...
       write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', wrong type of variable'
       read_well = .false.
   ELSE IF (Reason .LT. 0) THEN ! ... end of file reached ...
       write(*,'(a,i3,a)') 'Line ', i, ', END of file; no CDF data for phonon-peak' 
       write(*,'(a,i3,a)') 'Using atomic cross-sections for elastic scattering'
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
   END IF

2014 continue   ! if we must skip to the end for some reason
   if (.not. read_well) print*, trim(adjustl(Material_file))
   inquire(unit=FN2,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN2)             ! and if it is, close it
end subroutine reading_material_parameters


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




subroutine read_short_scdf(FN2, Target_atoms, NumPar, CDF_Phonon, Matter, Error_message, read_well)
   integer, intent(in) :: FN2 ! file number to read from
   type(Atom), dimension(:), allocatable, intent(inout) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   type(Flag), intent(inout) :: NumPar ! numerical parameters
   type(CDF), intent(inout) :: CDF_Phonon   ! CDF parameters for phonon to be read from a file
   type(Solid), intent(inout) :: Matter   ! all material parameters
   type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
   logical, intent(inout) :: read_well  ! did we read the file well?
   real(8) M
   integer Reason, i, j, k, l, N, Shl, CDF_coef, Shl_num
   character(100) Error_descript
   character(3) Name
   character(11) Shell_name
   character(30) Full_Name

   i = 0
   READ(FN2,'(a100)',IOSTAT=Reason) Matter%Target_name ! first line is the full material name
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

   READ(FN2,*,IOSTAT=Reason) Matter%Dens, Matter%Vsound, Matter%E_F   ! material density [g/cm^3], speed of sound in the material [m/s]
   call read_file(Reason, i, read_well) ! reports if everything read well
   if (.not. read_well) goto 2014
   ! [1/cm^3] atomic density:
   !Matter%At_Dens = 1.0d-3*Matter%Dens/(g_Mp*SUM(Target_atoms(:)%Mass)/size(Target_atoms)) ! [1/cm^3] atomic density
   Matter%At_Dens = 1.0d-3*Matter%Dens/(g_Mp*SUM(Target_atoms(:)%Mass*Target_atoms(:)%Pers)/SUM(Target_atoms(:)%Pers))
   Matter%v_f = sqrt(2.0d0*Matter%E_f/g_me) ! units as used in one-pole approximation (not SI!)
   
   call check_atomic_parameters(NumPar, Target_atoms, Error_message=Error_message, read_well=read_well) ! from module 'Dealing_with_EADL'
   
   do j = 1, N ! now mark all the cross sections as BEB:
      Target_atoms(j)%KOCS(:) = 2 ! BEB cross section for all shells
      Target_atoms(j)%KOCS_SHI(:) = 2 ! BEB cross section for all shells
   enddo
2014 continue
end subroutine read_short_scdf


subroutine reading_material_DOS(DOS_file, Mat_DOS, Matter, Target_atoms, Error_message, read_well)    ! read material DOS
    character(100), intent(in) :: DOS_file  ! file with material DOS
    type(Density_of_states), intent(inout) :: Mat_DOS  ! materail DOS
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    logical, intent(inout) :: read_well  ! did we read the file well?
    type(Solid), intent(in) :: Matter
    type(Atom), dimension(:), intent(in) :: Target_atoms
    
    real(8), dimension(:,:), allocatable :: Temp_DOS
    real(8) Sum_DOS, loc_DOS, E, dE, Sum_DOS_inv
    integer FN2, i, N, Reason, M
    character(100) Error_descript
    logical file_opened, file_exist, file_exist2
    FN2 = 202
    inquire(file=trim(adjustl(DOS_file)),exist=file_exist)    ! check if input file excists
    inquire(file='INPUT_DOS/Free_electron_DOS.dos',exist=file_exist2)    ! check if input file excists
    if (file_exist) then
       print*, 'DOS file is there: ', trim(adjustl(DOS_file))
       open(unit = FN2, FILE = trim(adjustl(DOS_file)), status = 'old', readonly)   ! if yes, open it and read
    elseif (file_exist2) then   ! Free-electron DOS approximation
       print*, 'DOS file ', trim(adjustl(DOS_file)), ' is not found, '
       print*, 'free-electron DOS approximation is used for valence band holes'
       open(unit = FN2, FILE = 'INPUT_DOS/Free_electron_DOS.dos', status = 'old', readonly)   ! if yes, open it and read
    else ! no DOS, try atomic approxiamtion instead...
       print*, 'Neither file ', trim(adjustl(DOS_file)), ' nor Free_electron_DOS.dos'
       print*, 'containing density of states are not found.'
       print*, 'The calculations proceed within the atomic approximation for energy levels.'
       Error_descript = 'Files '//trim(adjustl(DOS_file))//' and INPUT_DOS/Free_electron_DOS.dos are not found!'    ! description of an error
       call Save_error_details(Error_message, 6, Error_descript) ! write it into the error-log file
       read_well = .false.  ! no file found
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



subroutine reading_DSF_cross_sections(DSF_file, DSF_DEMFP, NumPar, Error_message, read_well)    ! read material DOS
    character(100), intent(in) :: DSF_file  ! file with DSF crossections
    type(Differential_MFP), dimension(:), allocatable, intent(inout) :: DSF_DEMFP  ! diffential elastic MFPs calculated with DSF
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    logical, intent(inout) :: read_well  ! did we read the file well?
    type(Flag), intent(inout) :: NumPar

    real(8), dimension(:,:), allocatable :: Temp_EMFP, Temp_array

    type(Differential_MFP), dimension(:), allocatable :: Temp_DEMFP
    real(8), dimension(:), allocatable :: Temp_E
    real(8) Sum_MFP, loc_EMFP, E, dE, Emin, Emax, Sum_MFP_emit
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
       call Save_error_details(Error_message, 6, Error_descript) ! write it into the error-log file
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



subroutine reading_DSF_cross_sections_OLD(DSF_file, DSF_DEMFP, NumPar, Error_message, read_well)    ! read material DOS
    character(100), intent(in) :: DSF_file  ! file with DSF crossections
    type(Differential_MFP), dimension(:), allocatable, intent(inout) :: DSF_DEMFP  ! diffential elastic MFPs calculated with DSF
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    logical, intent(inout) :: read_well  ! did we read the file well?
    type(Flag), intent(inout) :: NumPar
    
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
       call Save_error_details(Error_message, 6, Error_descript) ! write it into the error-log file
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


subroutine read_file(Reason, i, read_well)
   integer, intent(in) :: Reason    ! file number where to read from
   integer, intent(inout) :: i      ! line number
   logical, intent(inout) :: read_well  ! did we read ok?
   i = i + 1    ! it's next line
   IF (Reason .GT. 0)  THEN ! ... something wrong ...
       write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', wrong type of variable'
       read_well = .false.
   ELSE IF (Reason .LT. 0) THEN ! ... end of file reached ...
       write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', unexpected END of file'
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

subroutine Find_in_monotonous_1D_array(Array, Value0, Number)
   REAL(8), dimension(:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value0   ! to be found in the array as near as possible
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i, N, i_cur, i_1, i_2, coun
   real(8) temp_val, val_1, val_2
   
   N = size(Array)
   i_1 = 1
   val_1 = Array(i_1)
   i_2 = N
   val_2 = Array(i_2)
   i_cur = FLOOR((i_1+i_2)/2.0)
   temp_val = Array(i_cur)
   if (isnan(Value0)) then
        print*, 'The subroutine Find_in_monotonous_1D_array'
        print*, 'cannot proceed, the value of Value0 is', Value0
        write(*, '(f,f,f,f)') Value0, Array(i_cur), Array(i_1), Array(i_2)
        pause 'STOPPED WORKING...'
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

subroutine Find_in_monotonous_2D_array(Array, Value0, Indx, Number)
   REAL(8), dimension(:,:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value0   ! to be found in the array as near as possible
   integer, INTENT(in) :: Indx    ! index of the array, showing in which colonm we search
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i, N, i_cur, i_1, i_2, coun
   real(8) temp_val, val_1, val_2
   
   N = size(Array,2)
   i_1 = 1
   val_1 = Array(Indx,i_1)
   i_2 = N
   val_2 = Array(Indx,i_2)
   i_cur = FLOOR((i_1+i_2)/2.0)
   temp_val = Array(Indx,i_cur)
   
   if (isnan(Value0)) then
        print*, 'The subroutine Find_in_monotonous_2D_array'
        print*, 'cannot proceed, the value of Value0 is', Value0
        write(*, '(f,f,f,f)') Value0, Array(Indx,i_cur), Array(Indx,i_1), Array(Indx,i_2)
        pause 'STOPPED WORKING...'
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

subroutine Integrate_function_one(Int_type, x, f, x0, xn, res, Error_message)
    integer, intent(in) :: Int_type ! type of integration to be used: 0=trapeziod, 1=Simpson-3/8, 2=...
    real(8), dimension(:), intent(in) :: x  ! grid points
    real(8), dimension(:), intent(in) :: f  ! function
    real(8), intent(in) :: x0   ! starting point of integration
    real(8), intent(in) :: xn   ! ending point of integration
    real(8), intent(out) :: res ! result of integration
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    real(8) temp_x0, temp_xn
    integer i, N, Nx, Nf, i_0, i_n
    character(100) Err_data
    Nx = size(x)
    Nf = size(f)
    temp_x0 = x0
    temp_xn = xn
    if (Nx .NE. Nf) then
        Err_data = 'Trapeziod integration failed, size of x is not equal to size of f' ! no input data found
        call Save_error_details(Error_message, 3, Err_data)
    else if (temp_x0 .GT. x(Nx)) then
        Err_data = 'Trapeziod integration failed, starting point is larger than x(size(x))' ! no input data found
        call Save_error_details(Error_message, 4, Err_data)
    else if (temp_xn .GT. x(Nx)) then
        Err_data = 'Trapeziod integration failed, ending point is larger than x(size(x))' ! no input data found
        call Save_error_details(Error_message, 5, Err_data)
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




subroutine Integrate_function_save(Int_type, x, f, x0, xn, res, Error_message)
    integer, intent(in) :: Int_type ! type of integration to be used: 0=trapeziod, 1=Simpson-3/8, 2=...
    real(8), dimension(:), intent(in) :: x  ! grid points
    real(8), dimension(:), intent(in) :: f  ! function
    real(8), intent(in) :: x0   ! starting point of integration
    real(8), intent(in) :: xn   ! ending point of integration
    real(8), dimension(:), intent(out) :: res ! result of integration saved for each grid-point
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    real(8) temp_x0, temp_xn
    integer i, N, Nx, Nf, i_0, i_n
    character(100) Err_data
    Nx = size(x)
    Nf = size(f)
    temp_x0 = x0
    temp_xn = xn
    if (Nx .NE. Nf) then
        Err_data = 'Trapeziod integration failed, size of x is not equal to size of f' ! no input data found
        call Save_error_details(Error_message, 3, Err_data)
    else if (temp_x0 .GT. x(Nx)) then
        Err_data = 'Trapeziod integration failed, starting point is larger than x(size(x))' ! no input data found
        call Save_error_details(Error_message, 4, Err_data)
    else if (temp_xn .LT. temp_x0) then
        Err_data = 'Trapeziod integration failed, ending point is the starting one' ! no input data found
        call Save_error_details(Error_message, 5, Err_data)
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
   ! where as Array passed is decreaseing, so use negative array to find the correct index:
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

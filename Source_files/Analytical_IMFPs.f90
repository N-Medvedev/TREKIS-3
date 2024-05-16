!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains all subroutines needed to calculate cross-sections:

module Analytical_IMFPs
  use Universal_Constants   ! let it use universal constants
  use Objects   ! since it uses derived types, it must know about them from module 'Objects'
  use Reading_files_and_parameters, only : get_file_stat, Find_in_array_monoton, read_file_here, Linear_approx, read_SHI_MFP, &
                                           print_time_step
  use Cross_sections, only : Elastic_cross_section, TotIMFP, Tot_Phot_IMFP, SHI_Total_IMFP, construct_CDF, Total_copmlex_CDF, &
                             allocate_diff_CS_tables
  use Dealing_with_EADL, only : Count_lines_in_file
implicit none
PRIVATE

! Interface to automatically chose from the bubble array-sorting subroutines
interface sort_array
   module procedure sort_array_r ! for real arrays
   module procedure sort_array_c ! for complex arrays
end interface sort_array

interface find_order_of_number
   module procedure find_order_of_number_real
   module procedure find_order_of_number_int
end interface find_order_of_number


public :: Analytical_electron_dEdx, Analytical_ion_dEdx, Interpolate, printout_optical_CDF



contains


! Printout optical CDF reconstructed from Ritchie-Howie loss-function:
subroutine printout_optical_CDF(Output_path, Target_atoms, Matter, NumPar, Mat_DOS)
   character(100), intent(in) :: Output_path   ! path to the folder where the file is/will be stored
   type(Atom), dimension(:), intent(in), target :: Target_atoms  ! all data for target atoms
   type(Solid), intent(in) :: Matter   ! material properties
   type(Flag), intent(in) :: NumPar ! numerical parameters
   type(Density_of_states), intent(in) :: Mat_DOS
   !------------------------
   real(8), dimension(:), allocatable :: Temp_grid
   complex(8) :: complex_CDF ! constructed CDF
   real(8) :: Ele
   integer :: i, N, FN
   character(250) :: Output_file, temp_char1, KCS


    ! Make energy grid for printing out optical CDF:
    call get_grid_4CS(N, Temp_grid, Target_atoms, 0.001d0, 100.0d0-0.1d0, rescale_dE=0.1d0)  ! below

    ! if user requested, construct and printout optical CDF from the Ritchie-Howie fitted loss function:
    if (NumPar%print_CDF_optical) then
        KCS = ''
        select case (NumPar%kind_of_CDF)
        case (0)    ! Ritchie-Howie
            KCS = 'Ritchie_CDF'
        case (1)    ! single-pole CDF
            KCS = 'single_pole_CDF'
        endselect

        write(temp_char1, '(a)') 'OUTPUT_Optical_CDF_from_'//trim(adjustl(KCS))//'.dat'
        Output_file = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(temp_char1))

        ! Prepare output file:
        FN = 201
        open(FN, file=trim(adjustl(Output_file)))
        write(FN, '(a)') 'Energy(eV)    Re(CDF) Im(CDF) n   k'

        ! Write the total CDF into it:
        do i = 1, N
            Ele = Temp_grid(i)  ! energy grid [eV]
            ! Reconstruct CDF:
            call Total_copmlex_CDF(Target_atoms, Mat_DOS, Matter, NumPar, Ele, 0.0d0, &
                                        complex_CDF, photon=.true.) ! module "Cross_sections"
            ! Printout into the file:
            write(FN, '(f16.5, es25.6, es25.6, es25.6, es25.6)') Ele, dble(complex_CDF), aimag(complex_CDF), &
                                        dble(sqrt(complex_CDF)), aimag(sqrt(complex_CDF))
        enddo

        ! Clean up:
        close(FN)
    endif ! (NumPar%print_CDF_optical)
end subroutine printout_optical_CDF



!EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
! Calculates electron mean free paths with parallelization via openmp:
subroutine Analytical_electron_dEdx(Output_path, Material_name, Target_atoms, CDF_Phonon, Matter, Total_el_MFPs, &
                        Elastic_MFP, Error_message, read_well, DSF_DEMFP, Mat_DOS, NumPar, kind_of_particle, File_names)
    character(100), intent(in) :: Output_path   ! path to the folder where the file is/will be stored
    character(100), intent(in) :: Material_name ! name of the material
    type(Atom), dimension(:), intent(inout), target :: Target_atoms  ! all data for target atoms
    type(CDF), intent(in), target, optional :: CDF_Phonon ! CDF parameters of a phonon peak if excist
    type(Solid), intent(inout) :: Matter ! all material parameters
    type(All_MFP), dimension(:), allocatable, intent(inout) :: Total_el_MFPs   ! electron mean free paths for all shells
    type(MFP_elastic), intent(inout) :: Elastic_MFP ! elastic mean free path
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    logical, intent(out) :: read_well    ! did we read the file without an error?
    type(Density_of_states), intent(in) :: Mat_DOS
    type(Differential_MFP), dimension(:), intent(in) :: DSF_DEMFP
    type(All_names), intent(inout) :: File_names    ! file names to use later for gnuplot printing
    type(Flag), intent(inout) :: NumPar    
    character(8), intent(in) :: kind_of_particle    
    !--------------------------
    integer :: FN, FN1, FN2, FN3, FN4, FN_diff     ! file numbers where to save the output
    integer :: N, Nelast, Nsiz, N_diff_siz
    integer temp(1)
    real(8) Ele, IMFP_calc, dEdx, dEdx1, dEdx0, dE, Emin, Emax, L_tot, vel, Mass, InelMFP, ElasMFP, e_range, InelMFP_inv, ElasMFP_inv
    real(8), dimension(:,:), allocatable :: Temp_MFP
    real(8), dimension(:), allocatable :: Temp_grid
    integer Num_th, my_id, OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS, IMFP_last_modified
    integer i, j, k, Nat, Nshl, Reason, Va, Ord, Mnum, MFPnum, i_diff_CS
    character(200) Input_files, Input_elastic_file, File_el_range, File_hole_range, command
    character(200) temp_char, temp_char1, temp_ch, File_name, temp_char2, folder_diff_CS, diff_CS_file, full_CS_file
    !character(3) KCS
    character(10) KCS
    logical :: file_exist, file_opened, file_opened2, do_range, redo_MFP_default, diff_CS_file_exists

    read_well = .true.  ! so far so good
    do_range = .false.  ! don't recalculate electron range by default
    redo_MFP_default = NumPar%redo_IMFP ! save what the user defined
    Nat = size(Target_atoms)    ! how many atoms
    
    Emin = 1.0d0    ! [eV] we start with this minimum
    if (kind_of_particle .EQ. 'Electron') then ! it's an electrons
        Emax = 175.6d6*2.0d0/1836.0d0   ! [eV] ~maximum energy where relativism still can be neglected
    else if (kind_of_particle .EQ. 'Hole') then ! it's a hole
        Emax = Maxval(Mat_DOS%E(:)) + 1.0d0 ! [eV] only within the width of VB (or a bit more, just in case...)
    else ! it's a photon
        Emax = 175.6d6*2.0d0/1836.0d0   ! [eV] excluding Bremsstrahlung, photon can't be more energetic than this
    endif

    ! Make an energy grid for inelastic CS, taking into account finer resolution around ionization potentials:
    ! 1) For inelastic scttering grid:
    call get_grid_4CS(N, Temp_grid, Target_atoms, 0.1d0, Emax)  ! below
    if (.not. allocated(Total_el_MFPs)) allocate(Total_el_MFPs(Nat)) ! how many atoms    
    do j = 1, Nat   ! declair variables if they are not yet
        Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
        if (.not. allocated(Total_el_MFPs(j)%ELMFP)) allocate(Total_el_MFPs(j)%ELMFP(Nshl)) ! how many shells
        do k = 1, Nshl
            if (.not. allocated(Total_el_MFPs(j)%ELMFP(k)%E)) then
                allocate(Total_el_MFPs(j)%ELMFP(k)%E(N), source = Temp_grid)    ! we already know the grid, reuse it
                !Total_el_MFPs(j)%ELMFP(k)%E = 0.0d0
            endif
            if (.not. allocated(Total_el_MFPs(j)%ELMFP(k)%L)) then
                allocate(Total_el_MFPs(j)%ELMFP(k)%L(N))
                Total_el_MFPs(j)%ELMFP(k)%L = 1.0d24
            endif
            if (.not. allocated(Total_el_MFPs(j)%ELMFP(k)%dEdx)) then
                allocate(Total_el_MFPs(j)%ELMFP(k)%dEdx(N))
                Total_el_MFPs(j)%ELMFP(k)%dEdx = 0.0d0
            endif

            ! For diff.CS save tables:
            select case (NumPar%CS_method)
            case (1)
                if (.not. allocated(Target_atoms(j)%Int_diff_CS(k)%diffCS)) then
                    ! for each energy grid point there will be its own table:
                    allocate(Target_atoms(j)%Int_diff_CS(k)%E(N))
                    allocate(Target_atoms(j)%Int_diff_CS(k)%diffCS(N))
                endif
            endselect

        enddo
    enddo
    allocate(temp_MFP(2,N))

    ! 2) For elastic scttering grid:
    if (kind_of_particle .NE. 'Photon') then ! elastic only for massive particles:
        ! Set the grid for electron energies in ELASTIC cross-section:
        call get_grid_4CS(Nelast, Elastic_MFP%Total%E, Target_atoms, 0.01d0, Emax)  ! below

        if (.not. allocated(Elastic_MFP%Total%L)) then
            allocate(Elastic_MFP%Total%L(Nelast))  ! [A] MFP itself
            Elastic_MFP%Total%L = 1.0d24
        endif
        if (.not. allocated(Elastic_MFP%Total%dEdx)) then
            allocate(Elastic_MFP%Total%dEdx(Nelast))  ! [eV/A] mean energy loss
            Elastic_MFP%Total%dEdx = 0.0d0
        endif
    endif
    

    ! Which cross sections we use - CDF vs BEB:
    temp = 0

    KCS = ''
    select case (NumPar%kind_of_CDF)
    case (0)    ! Ritchie-Howie
        KCS = '_CDF'
    case (1)    ! single-pole CDF
        KCS = '_spCDF'
    endselect

    do i = 1, size(Target_atoms) ! all atoms
       temp(1) = maxval(Target_atoms(i)%KOCS(:)) ! check if there is at least one shell for which we use BEB instead of CDF
       if (temp(1) .GE. 2) then ! all CDF
          KCS = 'BEB'
          exit ! if there is at least one BEB cross section name the file 'BEB'
       endif
    enddo
    
    !==============================
    kind_of_particle1:if (kind_of_particle .EQ. 'Electron') then

        write(temp_ch, '(f9.2)') Matter%temp
        if(Matter%El_eff_mass .EQ. 0) then
            write(temp_char, '(a, a, a)') 'DOS_', trim(adjustl(temp_ch)), '_K'
        else if ((Matter%El_eff_mass .LT. 0) .OR. (Matter%El_eff_mass .EQ. 1.0d0)) then
            write(temp_char, '(a, a, a)') 'Me_', trim(adjustl(temp_ch)), '_K'
        else
            write(temp_char, '(f5.2,a,a,a)') Matter%El_eff_mass, '_me_', trim(adjustl(temp_ch)), '_K'
        endif
        temp_char = trim(adjustl(KCS(2:)))//'_'//trim(adjustl(temp_char))

        if (KCS .EQ. 'BEB') then ! BEB vs CDF cross section:
          temp_char1 = 'OUTPUT_Electron_IMFPs_'//trim(adjustl(KCS))//'.dat'
          !Input_files = trim(adjustl(Output_path))//'/OUTPUT_Electron_IMFPs_'//trim(adjustl(KCS))//'.dat'
          !if (allocated(File_names%F)) File_names%F(2) = '/OUTPUT_Electron_IMFPs_'//trim(adjustl(KCS))//'.dat' ! save for later use
          File_el_range = trim(adjustl(Output_path))//'/OUTPUT_Electron_range_'//trim(adjustl(KCS))//'.dat'
        else ! CDF:
          if (NumPar%kind_of_DR .EQ. 4) then    ! Delta-CDF
            temp_char1 = 'OUTPUT_Electron_IMFPs_Delta_'//trim(adjustl(temp_char))//'.dat'
            !Input_files = trim(adjustl(Output_path))//'/OUTPUT_Electron_IMFPs_Delta_'//trim(adjustl(temp_char))//'.dat'
            !if (allocated(File_names%F)) File_names%F(2) = '/OUTPUT_Electron_IMFPs_Delta_'//trim(adjustl(temp_char))//'.dat' ! save for later use
            File_el_range = trim(adjustl(Output_path))//'/OUTPUT_Electron_Delta_range_'//trim(adjustl(temp_char))//'.dat'
          else if (NumPar%kind_of_DR .EQ. 3) then
            temp_char1 = 'OUTPUT_Electron_IMFPs_Ritchie_'//trim(adjustl(temp_char))//'.dat'
            !Input_files = trim(adjustl(Output_path))//'/OUTPUT_Electron_IMFPs_Ritchie_'//trim(adjustl(temp_char))//'.dat'
            !if (allocated(File_names%F)) File_names%F(2) = '/OUTPUT_Electron_IMFPs_Ritchie_'//trim(adjustl(temp_char))//'.dat' ! save for later use
            File_el_range = trim(adjustl(Output_path))//'/OUTPUT_Electron_Ritchie_range_'//trim(adjustl(temp_char))//'.dat'
          else if (NumPar%kind_of_DR .EQ. 2) then
            temp_char1 = 'OUTPUT_Electron_IMFPs_Plasmon_pole_'//trim(adjustl(temp_char))//'.dat'
            !Input_files = trim(adjustl(Output_path))//'/OUTPUT_Electron_IMFPs_Plasmon_pole_'//trim(adjustl(temp_char))//'.dat'
            !if (allocated(File_names%F)) File_names%F(2) = '/OUTPUT_Electron_IMFPs_Plasmon_pole_'//trim(adjustl(temp_char))//'.dat' ! save for later use
            File_el_range = trim(adjustl(Output_path))//'/OUTPUT_Electron_Plasmon_pole_range_'//trim(adjustl(temp_char))//'.dat'
          else
            temp_char1 = 'OUTPUT_Electron_IMFPs_Free_'//trim(adjustl(temp_char))//'.dat'
            !Input_files = trim(adjustl(Output_path))//'/OUTPUT_Electron_IMFPs_Free_'//trim(adjustl(temp_char))//'.dat'
            !if (allocated(File_names%F)) File_names%F(2) = '/OUTPUT_Electron_IMFPs_Free_'//trim(adjustl(temp_char))//'.dat' ! save for later use
            File_el_range = trim(adjustl(Output_path))//'/OUTPUT_Electron_range_Free_'//trim(adjustl(temp_char))//'.dat'
          endif
        endif ! which name
        if (allocated(File_names%F)) File_names%F(2) = trim(adjustl(temp_char1)) ! save for later use
        Input_files = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(temp_char1))

        ! IMFP files for electron:
        FN = 201
        inquire(file=trim(adjustl(Input_files)),exist=file_exist)    ! check if input file excists

        ! Check, if file with MFP was created with paramters in actual CDF file, find out when this file was last modified:
        if (file_exist) then
            call get_file_stat(trim(adjustl(Input_files)), Last_modification_time=IMFP_last_modified) ! above
            !print*, 'IMFP file last modified on:', IMFP_last_modified
            if (IMFP_last_modified < NumPar%Last_mod_time_CDF) then
                NumPar%redo_IMFP = .true. ! Material parameters changed, recalculate IMFPs
                print*, 'File with CDF was modified more recently than the electron MFP => recalculating MFP'
            endif

            ! Check if the file is consistent with the grid set:
            open(FN, file=trim(adjustl(Input_files)), action='read')
            call count_lines_in_file(FN, Nsiz) ! module "Dealing_with_EADL"
            if (Nsiz /= N) then
                NumPar%redo_IMFP = .true. ! Grid mismatch, recalculate IMFPs
                print*, 'Energy grid mismatch in electron MFP file => recalculating MFP'
            endif
            close(FN)
        endif


        ! Check if precalculated diff.CS is required:
        select case (NumPar%CS_method)
        case (1)    ! save files are required
            full_CS_file = trim(adjustl(temp_char1))    ! to be reused below for diff.CS file naming
            folder_diff_CS = 'diff_CS'  ! folder name
            ! Add the path to it:
            folder_diff_CS = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(folder_diff_CS))
#ifndef __GFORTRAN__
            ! for intel fortran compiler:
            inquire(DIRECTORY=trim(adjustl(folder_diff_CS)),exist=diff_CS_file_exists)    ! check if folder excists
#else
            ! for gfortran compiler:
            inquire(FILE=trim(adjustl(folder_diff_CS)),exist=diff_CS_file_exists)    ! check if folder excists
#endif
            if (.not.diff_CS_file_exists) then
                print*, 'Files with diff.CS for electrons are required => recalculating MFP'
                NumPar%redo_IMFP = .true.   ! diff.CS need to be calculated
                command='mkdir '//trim(adjustl(folder_diff_CS)) ! to create a folder use this command
                CALL system(command)  ! create the folder
                write(*,'(a)') ' Folder '//trim(adjustl(folder_diff_CS))//' created for diff.CS files'
            endif
        case default  ! no files, calculate on the fly
            ! Nothing to do
        end select

        if (file_exist .and. .not.NumPar%redo_IMFP) then    ! read from the file:
            write(*,'(a,a,a)') 'IMFPs of an electron in ', trim(adjustl(Material_name)), ' are already in the file:'
            write(*, '(a)') trim(adjustl(Input_files))
            write(*, '(a)') ' '
            open(FN, file=trim(adjustl(Input_files)), action='read')
        else    ! create and write to the file:
            call All_shells_Electron_MFP(N, Target_atoms, Total_el_MFPs, Mat_DOS, Matter, NumPar, kind_of_particle) ! calculate all IMFPs
            open(FN, file=trim(adjustl(Input_files)))
            write(*,'(a,a,a)') 'Calculated inelastic mean free paths of an electron in ', trim(adjustl(Material_name)), ' are stored in the file'
            write(*, '(a)') trim(adjustl(Input_files))
            write(*, '(a)') ' '
        endif

    !==============================
    else if (kind_of_particle .EQ. 'Hole') then
        FN = 202
        if (KCS .EQ. 'BEB') then ! BEB vs CDF cross section:
           temp_char1 = 'OUTPUT_Hole_IMFPs_BEB.dat'
           !Input_files = trim(adjustl(Output_path))//'/OUTPUT_Hole_IMFPs_BEB.dat'
           !if (allocated(File_names%F)) File_names%F(3) = 'OUTPUT_Hole_IMFPs_BEB.dat' ! save for later use
           File_hole_range = trim(adjustl(Output_path))//'/OUTPUT_Hole_range_BEB.dat'
        else
           KCS = ''
           select case (NumPar%kind_of_CDF)
           case (0)    ! Ritchie-Howie
            KCS = '_CDF'
           case (1)    ! single-pole CDF
            KCS = '_spCDF'
           endselect
           write(temp_char, '(f7.2, a)') Matter%temp, '_K'
           temp_char1 = 'OUTPUT_Hole_IMFPs_CDF'//trim(adjustl(KCS))//'_'//trim(adjustl(temp_char))//'.dat'
           !Input_files = trim(adjustl(Output_path))//'/OUTPUT_Hole_IMFPs_CDF'//trim(adjustl(KCS))//'_'//trim(adjustl(temp_char))//'.dat'
           !if (allocated(File_names%F)) File_names%F(3) = 'OUTPUT_Hole_IMFPs_CDF'//trim(adjustl(KCS))//'_'
                                            !//trim(adjustl(temp_char))//'.dat' ! save for later use
           File_hole_range = trim(adjustl(Output_path))//'/OUTPUT_Hole_range_CDF'//trim(adjustl(KCS))//'_'//trim(adjustl(temp_char))//'.dat'
        endif
        Input_files = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(temp_char1))
        if (allocated(File_names%F)) File_names%F(3) = trim(adjustl(temp_char1)) ! save for later use

        !file_exist = .false.                !This file must be overwriten before each calculation.
        inquire(file=trim(adjustl(Input_files)),exist=file_exist)    ! check if input file excists
        ! Check, if file with MFP was created with paramters in actual CDF file, find out when this file was last modified:
        if (file_exist) then
            call get_file_stat(trim(adjustl(Input_files)), Last_modification_time=IMFP_last_modified) ! above
            !print*, 'IMFP file last modified on:', IMFP_last_modified
            if (IMFP_last_modified < NumPar%Last_mod_time_CDF) then
                NumPar%redo_IMFP = .true. ! Material parameters changed, recalculate IMFPs
                print*, 'File with CDF was modified more recently than the hole MFP => recalculating MFP'
            endif

            ! Check if the file is consistent with the grid set:
            open(FN, file=trim(adjustl(Input_files)), action='read')
            call count_lines_in_file(FN, Nsiz) ! module "Dealing_with_EADL"
            if (Nsiz /= N) then
                NumPar%redo_IMFP = .true. ! Grid mismatch, recalculate IMFPs
                print*, 'Energy grid mismatch in hole MFP file => recalculating MFP'
            endif
            close(FN)
        endif

        if (file_exist .and. .not.NumPar%redo_IMFP) then    ! read from the file:
            write(*,'(a,a,a)') 'IMFPs of an electron in ', trim(adjustl(Material_name)), ' are already in the file:'
            write(*, '(a)') trim(adjustl(Input_files))
            write(*, '(a)') ' '
            open(FN, file=trim(adjustl(Input_files)), action='read')
        else    ! create and write to the file:
            call All_shells_Electron_MFP(N, Target_atoms, Total_el_MFPs, Mat_DOS, Matter, NumPar, kind_of_particle) ! calculate all IMFPs
            open(FN, file=trim(adjustl(Input_files)))
            write(*,'(a,a,a)') 'Calculated inelastic mean free paths of a hole in ', trim(adjustl(Material_name)), ' are stored in the file'
            write(*, '(a)') trim(adjustl(Input_files))
            write(*, '(a)') ' '
        endif

    !==============================
    else kind_of_particle1 ! photon
        !if (KCS .EQ. 'BEB') then ! BEB vs CDF cross section:
        if (.not.allocated(Target_atoms(1)%Ritchi(size(Target_atoms))%E0)) then ! no CDF known, use EPDL97 database:
           temp_char1 = 'OUTPUT_Photon_IMFPs_EPDL.dat'
           !Input_files = trim(adjustl(Output_path))//'/OUTPUT_Photon_IMFPs_EPDL.dat'
           !if (allocated(File_names%F)) File_names%F(7) = 'OUTPUT_Photon_IMFPs_EPDL.dat' ! save for later use
        else
           KCS = ''
           select case (NumPar%kind_of_CDF)
           case (0)    ! Ritchie-Howie
            KCS = '_CDF'
           case (1)    ! single-pole CDF
            KCS = '_spCDF'
           endselect
           write(temp_char, '(f7.2, a)') Matter%temp, '_K'
           temp_char1 = 'OUTPUT_Photon_IMFPs'//trim(adjustl(KCS))//'_'//trim(adjustl(temp_char))//'.dat'
           !Input_files = trim(adjustl(Output_path))//'/OUTPUT_Photon_IMFPs'//trim(adjustl(KCS))//'_'//trim(adjustl(temp_char))//'.dat'
           !if (allocated(File_names%F)) then
           !     File_names%F(7) = 'OUTPUT_Photon_IMFPs'//trim(adjustl(KCS))//'_'//trim(adjustl(temp_char))//'.dat' ! save for later use
           !endif
        endif
        Input_files = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(temp_char1))
        if (allocated(File_names%F)) File_names%F(7) = trim(adjustl(temp_char1)) ! save for later use

        inquire(file=trim(adjustl(Input_files)),exist=file_exist)    ! check if input file excists

        ! Check, if file with MFP was created with paramters in actual CDF file, find out when this file was last modified:
        if (file_exist) then
            call get_file_stat(trim(adjustl(Input_files)), Last_modification_time=IMFP_last_modified) ! above
            !print*, 'IMFP file last modified on:', IMFP_last_modified
            if (IMFP_last_modified < NumPar%Last_mod_time_CDF) then
                NumPar%redo_IMFP = .true. ! Material parameters changed, recalculate IMFPs
                print*, 'File with CDF was modified more recently than the photon MFP => recalculating MFP'
            endif
            ! Check if the file is consistent with the grid set:
            open(newunit=FN, file=trim(adjustl(Input_files)), ACTION='READ')
            call count_lines_in_file(FN, Nsiz) ! module "Dealing_with_EADL"
            if (Nsiz /= N) then
                NumPar%redo_IMFP = .true. ! Grid mismatch, recalculate IMFPs
                print*, 'Energy grid mismatch in photon MFP file => recalculating MFP'
            endif
            close(FN)
        endif


        if (file_exist .and. .not.NumPar%redo_IMFP) then    ! read from the file:
            write(*,'(a,a,a)') 'IMFPs of a photon in ', trim(adjustl(Material_name)), ' are already in the file:'
            write(*, '(a)') trim(adjustl(Input_files))
            write(*, '(a)') ' '
            open(newunit=FN, file=trim(adjustl(Input_files)), ACTION='READ')
        else    ! create and write to the file:
            call All_shells_Photon_MFP(N, Target_atoms, Total_el_MFPs, Matter, NumPar, Mat_DOS) ! calculate all IMFPs
            open(newunit=FN, file=trim(adjustl(Input_files)))
            write(*,'(a,a,a)') 'Calculated attenuation lengths (IMFPs) of a photon in ', trim(adjustl(Material_name)), ' are stored in the file'
            write(*, '(a)') trim(adjustl(Input_files))
            write(*, '(a)') ' '
        endif
    endif kind_of_particle1


   !===================================
   !Now write the output into the file:
   e_range = 0.0d0 ! start with the range calculations
   dEdx1 = 0.0d0
   dEdx0 = 1.0d30
   do i = 1, N
      if (.not. file_exist .or. NumPar%redo_IMFP) then  ! if file didn't exist and we just created it:
        write(FN,'(f)', advance='no') Total_el_MFPs(1)%ELMFP(1)%E(i)
        IMFP_calc = 0.0d0 ! to sum up for a total IMFP
      else
        read(FN,'(f)', advance='no', IOSTAT=Reason) Total_el_MFPs(1)%ELMFP(1)%E(i)
        call read_file_here(Reason, i, read_well)
        if (.not. read_well) goto 2015
      endif
      do j = 1, Nat
         Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
         do k = 1, Nshl
            if (.not. file_exist .or. NumPar%redo_IMFP) then  ! if file didn't exist and we just created it:
                write(FN,'(e)', advance='no') Total_el_MFPs(j)%ELMFP(k)%L(i)    ! write IMFP for all shells
                IMFP_calc = IMFP_calc + 1.0d0/Total_el_MFPs(j)%ELMFP(k)%L(i)    ! sum them all up to get the total value
            else
                read(FN,'(e)', advance='no') Total_el_MFPs(j)%ELMFP(k)%L(i)    ! write IMFP for all shells
                call read_file_here(Reason, i, read_well)
                if (.not. read_well) goto 2015
                if ((j .NE. 1) .OR. (k .NE. 1)) then
                    Total_el_MFPs(j)%ELMFP(k)%E(i) = Total_el_MFPs(1)%ELMFP(1)%E(i) ! save it for all elements and shells
                endif
            endif
         enddo
      enddo
      temp_MFP(1,i) = Total_el_MFPs(1)%ELMFP(1)%E(i)
      if (.not. file_exist .or. NumPar%redo_IMFP) then  ! if file didn't exist and we just created it:
        write(FN,'(e)') 1.0d0/IMFP_calc ! write total IMFP
        temp_MFP(2,i) = 1.0d0/IMFP_calc
        do_range = .true.  ! recalculate electron range only if dEdx is calculated
      else
        read(FN,'(e)') IMFP_calc    ! read total IMFP
        temp_MFP(2,i) = IMFP_calc
        call read_file_here(Reason, i, read_well)
        if (.not. read_well) goto 2015
      endif
   enddo
   flush(FN)

   ! Diff.CS tables, if required:
   if (kind_of_particle .EQ. 'Electron') then
    select case (NumPar%CS_method)
    case (1)    ! save files are required
        if (NumPar%verbose) call print_time_step('Starting dealing with diff. CS tables and files:', msec=.true.)

        FN_diff = 333   ! file number to be reused

        do j = 1, Nat   ! for all types of atoms
            Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
            do k = 1, Nshl  ! for all orbitals
                do i = 1, N ! for all energy grid points
                    ! Construct the file name:
                    Target_atoms(j)%Int_diff_CS(k)%E(i) = Total_el_MFPs(j)%ELMFP(k)%E(i)

                    ! Name of the atom and shell:
                    temp_char = trim(adjustl(full_CS_file( 8 : LEN(trim(adjustl(full_CS_file)))-4 ))) //'_'// &
                                trim(adjustl(Target_atoms(j)%Name))//'_'//trim(adjustl(Target_atoms(j)%Shell_name(k)))
                    ! Add energy grid point:
                    write(temp_char1,'(f14.3)') Target_atoms(j)%Int_diff_CS(k)%E(i)
                    ! Combine the info into file name:
                    temp_char = trim(adjustl(temp_char))//'_'//trim(adjustl(temp_char1))//'.dat'
                    ! Full name of the file with diff.CS for this energy grid point, element and shell:
                    diff_CS_file = trim(adjustl(folder_diff_CS))//trim(adjustl(NumPar%path_sep))//trim(adjustl(temp_char))

                    if ((.not. file_exist) .or. (.not.diff_CS_file_exists) .or. NumPar%redo_IMFP) then  ! if file doesn't exist
                        open(FN_diff, file=trim(adjustl(diff_CS_file)))   ! create it
                        ! Write the data into this file:
                        N_diff_siz = size(Target_atoms(j)%Int_diff_CS(k)%diffCS(i)%dsdhw)
                        do i_diff_CS = 1, N_diff_siz
                            write(FN_diff, '(es,es)') Target_atoms(j)%Int_diff_CS(k)%diffCS(i)%hw(i_diff_CS), &
                                                      Target_atoms(j)%Int_diff_CS(k)%diffCS(i)%dsdhw(i_diff_CS)
                        enddo ! i_diff_CS
                    else    ! if file exist, read from it:
                        open(FN_diff, file=trim(adjustl(diff_CS_file)), action='read')
                        ! Read the data from this file:
                        call Count_lines_in_file(FN_diff, N_diff_siz) ! count how many line the file contains
                        ! Allocate the diff.CS tables:
                        if (.not.allocated(Target_atoms(j)%Int_diff_CS(k)%diffCS(i)%dsdhw)) then
                            allocate(Target_atoms(j)%Int_diff_CS(k)%diffCS(i)%dsdhw(N_diff_siz))
                        endif
                        if (.not.allocated(Target_atoms(j)%Int_diff_CS(k)%diffCS(i)%hw)) then
                            allocate(Target_atoms(j)%Int_diff_CS(k)%diffCS(i)%hw(N_diff_siz))
                        endif

                        do i_diff_CS = 1, N_diff_siz
                            read(FN_diff,*) Target_atoms(j)%Int_diff_CS(k)%diffCS(i)%hw(i_diff_CS), &
                                            Target_atoms(j)%Int_diff_CS(k)%diffCS(i)%dsdhw(i_diff_CS)
                        enddo ! i_diff_CS
                    endif
                    !print*, j,k, i, Target_atoms(j)%Int_diff_CS(k)%E(i), trim(adjustl(diff_CS_file))
                    close (FN_diff)
                enddo ! k
            enddo ! k
        enddo ! j

        if (NumPar%verbose) call print_time_step('Done with diff. CS tables and files:', msec=.true.)
    case default  ! no files, calculate on the fly
        ! Nothing to do
    end select
   endif ! 'Electron'


   NumPar%redo_IMFP = redo_MFP_default ! defualt it for the next kind of particle


   !######################### Now do the same for elastic mean free path of an electron and hole:
   kind_of_part2:if (kind_of_particle .EQ. 'Electron') then
         0001 continue
         select case (NumPar%kind_of_EMFP) ! el_elastic_CS
         case (2) ! Read DSF elastic MFP
            write(temp_char1, '(f7.2, a)') Matter%temp, '_K'
            temp_char2 = 'OUTPUT_Electron_DSF_EMFPs_'//trim(adjustl(temp_char1))//'.dat'
            !Input_elastic_file = trim(adjustl(Output_path))//'/OUTPUT_Electron_DSF_EMFPs_'//trim(adjustl(temp_char1))//'.dat'
            Input_elastic_file = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(temp_char2))
            if (allocated(File_names%F)) File_names%F(4) = trim(adjustl(temp_char2)) ! save for later
            FN2 = 2032
            inquire(file=trim(adjustl(Input_elastic_file)),exist=file_exist)    ! check if input file excists

            ! Check, if file with MFP was created with paramters in actual CDF file, find out when this file was last modified:
            if (file_exist) then
                call get_file_stat(trim(adjustl(Input_elastic_file)), Last_modification_time=IMFP_last_modified) ! above
                !print*, 'IMFP file last modified on:', IMFP_last_modified
                if (IMFP_last_modified < NumPar%Last_mod_time_DSF) then
                    NumPar%redo_EMFP = .true. ! Material parameters changed, recalculate EMFPs
                    print*, 'File with DSF was modified more recently than the MFP => recalculating MFP'
                endif
            endif
            
            if (file_exist .and. .not.NumPar%redo_EMFP) then    ! read from the file:
                write(*,'(a,a,a)') 'DSF EMFPs of an electron in ', trim(adjustl(Material_name)), ' are in the file:'
                write(*, '(a)') trim(adjustl(Input_elastic_file))
                write(*, '(a)') ' '
                open(FN2, file=trim(adjustl(Input_elastic_file)), ACTION='READ')
                call count_lines_in_file(FN2, Nelast)
                deallocate(Elastic_MFP%Total%E)
                deallocate(Elastic_MFP%Total%L)
                deallocate(Elastic_MFP%Total%dEdx)
                allocate(Elastic_MFP%Total%E(Nelast))
                allocate(Elastic_MFP%Total%L(Nelast))
                allocate(Elastic_MFP%Total%dEdx(Nelast))
                ! And separately absorption and emission of energy:
                if (allocated(Elastic_MFP%Emit%E)) deallocate(Elastic_MFP%Emit%E)
                if (allocated(Elastic_MFP%Emit%L)) deallocate(Elastic_MFP%Emit%L)
                if (allocated(Elastic_MFP%Emit%dEdx)) deallocate(Elastic_MFP%Emit%dEdx)
                if (allocated(Elastic_MFP%Absorb%E)) deallocate(Elastic_MFP%Absorb%E)
                if (allocated(Elastic_MFP%Absorb%L)) deallocate(Elastic_MFP%Absorb%L)
                if (allocated(Elastic_MFP%Absorb%dEdx)) deallocate(Elastic_MFP%Absorb%dEdx)
                allocate(Elastic_MFP%Emit%E(Nelast), source = 0.0d0)
                allocate(Elastic_MFP%Emit%L(Nelast), source = 0.0d0)
                allocate(Elastic_MFP%Emit%dEdx(Nelast), source = 0.0d0)
                allocate(Elastic_MFP%Absorb%E(Nelast), source = 0.0d0)
                allocate(Elastic_MFP%Absorb%L(Nelast), source = 0.0d0)
                allocate(Elastic_MFP%Absorb%dEdx(Nelast), source = 0.0d0)
            else    ! create and write to the file:
                write(*,'(a,a,a)') 'Calculated elastic mean free paths of an electron in ', trim(adjustl(Material_name)), ' are stored in the file:'
                write(*, '(a)') trim(adjustl(Input_elastic_file))
                write(*, '(a)') ' '
                open(FN2, file=trim(adjustl(Input_elastic_file)))
                deallocate(Elastic_MFP%Total%E)
                deallocate(Elastic_MFP%Total%L)
                deallocate(Elastic_MFP%Total%dEdx)
                Nelast = size(DSF_DEMFP%E)
                allocate(Elastic_MFP%Total%E(Nelast))
                allocate(Elastic_MFP%Total%L(Nelast))
                allocate(Elastic_MFP%Total%dEdx(Nelast))
                ! And separately absorption and emission of energy:
                if (allocated(Elastic_MFP%Emit%E)) deallocate(Elastic_MFP%Emit%E)
                if (allocated(Elastic_MFP%Emit%L)) deallocate(Elastic_MFP%Emit%L)
                if (allocated(Elastic_MFP%Emit%dEdx)) deallocate(Elastic_MFP%Emit%dEdx)
                if (allocated(Elastic_MFP%Absorb%E)) deallocate(Elastic_MFP%Absorb%E)
                if (allocated(Elastic_MFP%Absorb%L)) deallocate(Elastic_MFP%Absorb%L)
                if (allocated(Elastic_MFP%Absorb%dEdx)) deallocate(Elastic_MFP%Absorb%dEdx)
                allocate(Elastic_MFP%Emit%E(Nelast), source = 0.0d0)
                allocate(Elastic_MFP%Emit%L(Nelast), source = 0.0d0)
                allocate(Elastic_MFP%Emit%dEdx(Nelast), source = 0.0d0)
                allocate(Elastic_MFP%Absorb%E(Nelast), source = 0.0d0)
                allocate(Elastic_MFP%Absorb%L(Nelast), source = 0.0d0)
                allocate(Elastic_MFP%Absorb%dEdx(Nelast), source = 0.0d0)

                do i = 1, Nelast
                    Elastic_MFP%Total%E(i) = DSF_DEMFP(i)%E
                    Elastic_MFP%Total%L(i) = DSF_DEMFP(i)%dL(size(DSF_DEMFP(1)%dL))
                    Elastic_MFP%Emit%L(i)   = DSF_DEMFP(i)%dL_emit(size(DSF_DEMFP(1)%dL_emit)) ! Emission
                    Elastic_MFP%Absorb%L(i) = DSF_DEMFP(i)%dL_absorb(size(DSF_DEMFP(1)%dL_absorb)) ! Absorption
                    !print*, i, Elastic_MFP%Total%E(i), Elastic_MFP%Total%L(i), 1.0d0/(1.0d0/Elastic_MFP%Emit%L(i) + 1.0d0/Elastic_MFP%Absorb%L(i))
                enddo
                !pause 'DSF Test electron'
           endif
         case (1) ! Calculate or read CDF elastic mean free path

            !print*, 'Phonon CDF:', allocated(CDF_Phonon%A)
            if (allocated(CDF_Phonon%A)) then
                ! Info about cross section model:
                KCS = ''
                select case (NumPar%kind_of_CDF_ph)
                case (0)    ! Ritchie-Howie
                    KCS = '_CDF'
                case (1)    ! single-pole CDF
                    KCS = '_spCDF'
                endselect

                ! Add info about tempreature:
                write(temp_char1, '(f7.2, a)') Matter%temp, '_K'

                ! Add info about effective charge:
                select case(NumPar%CDF_elast_Zeff)
                case default
                    temp_char1 = 'Zeff_'//trim(adjustl(temp_char1))
                case (1)
                    temp_char1 = 'Z=1_'//trim(adjustl(temp_char1))
                case (2)
                    temp_char1 = 'Z_CDFe_FF_'//trim(adjustl(temp_char1))
                case (3)
                    temp_char1 = 'Z_CDFe_'//trim(adjustl(temp_char1))
                endselect

                temp_char2 = 'OUTPUT_Electron_EMFPs'//trim(adjustl(KCS))//'_'//trim(adjustl(temp_char1))//'.dat'
                Input_elastic_file = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(temp_char2))
                !Input_elastic_file = trim(adjustl(Output_path))//'/OUTPUT_Electron_EMFPs'//trim(adjustl(KCS))//'_'// &
                !    trim(adjustl(temp_char1))//'.dat'
                if (allocated(File_names%F)) then
                    !File_names%F(4) = 'OUTPUT_Electron_EMFPs'//trim(adjustl(KCS))//'_'//trim(adjustl(temp_char1))//'.dat' ! save for later use
                    File_names%F(4) = trim(adjustl(temp_char2)) ! save for later use
                endif
                FN2 = 203
                inquire(file=trim(adjustl(Input_elastic_file)),exist=file_exist)    ! check if input file excists

                ! Check, if file with MFP was created with paramters in actual CDF file, find out when this file was last modified:
                if (file_exist) then
                    call get_file_stat(trim(adjustl(Input_elastic_file)), Last_modification_time=IMFP_last_modified) ! above
                    !print*, 'IMFP file last modified on:', IMFP_last_modified
                    if (IMFP_last_modified < NumPar%Last_mod_time_CDF) then
                        NumPar%redo_EMFP = .true. ! Material parameters changed, recalculate EMFPs
                        print*, 'File with CDF was modified more recently than the MFP => recalculating MFP'
                    endif
                    ! Check if the file is consistent with the grid set:
                    open(FN2, file=trim(adjustl(Input_elastic_file)), ACTION='READ')
                    call count_lines_in_file(FN2, Nsiz) ! module "Dealing_with_EADL"
                    if (Nsiz /= Nelast) then
                        NumPar%redo_EMFP = .true. ! Grid mismatch, recalculate IMFPs
                        print*, 'Energy grid mismatch in elastic electron MFP file => recalculating MFP'
                    endif
                    close(FN2)
                endif

                if (file_exist .and. .not.NumPar%redo_EMFP) then    ! read from the file:
                    write(*,'(a,a,a)') 'Calculated with '//trim(adjustl(KCS(2:)))//' EMFPs of an electron in ', &
                        trim(adjustl(Material_name)), ' are already in the file:'
                    write(*, '(a)') trim(adjustl(Input_elastic_file))
                    write(*, '(a)') ' '
                    open(FN2, file=trim(adjustl(Input_elastic_file)), ACTION='READ')
                else    ! create and write to the file:
                    call All_elastic_scattering(Nelast, Target_atoms, CDF_Phonon, Matter, Elastic_MFP%Total, &
                                                NumPar, Mat_DOS, kind_of_particle)
                    open(FN2, file=trim(adjustl(Input_elastic_file)))
                    write(*,'(a,a,a)') 'Elastic mean free paths of an electron calculated with '//trim(adjustl(KCS(2:)))// &
                        ' phonon peaks in ', trim(adjustl(Material_name)), ' are stored in the file'
                    write(*, '(a)') trim(adjustl(Input_elastic_file))
                    write(*, '(a)') ' '
                endif
            else
                write(*,'(a,a,a)') 'Coefficients of CDF phonon peaks for', trim(adjustl(Material_name)), 'is not specified.'
                write(*, '(a)') 'Calculation will proceed with Mott atomic cross-sections.'
                write(*, '(a)') ' '
                NumPar%kind_of_EMFP = 0
                go to 0001
            endif
         case (0) ! Calculate or read Mott elastic MFP
            write(temp_char1, '(f7.2, a)') Matter%temp, '_K'

            temp_char2 = 'OUTPUT_Electron_EMFPs_Mott.dat'

            Input_elastic_file = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(temp_char2))
            if (allocated(File_names%F)) File_names%F(4) = trim(adjustl(temp_char2))
            FN2 = 2031
            inquire(file=trim(adjustl(Input_elastic_file)),exist=file_exist)    ! check if input file excists

            ! Check, if file with MFP was created with paramters in actual CDF file, find out when this file was last modified:
            if (file_exist) then
                call get_file_stat(trim(adjustl(Input_elastic_file)), Last_modification_time=IMFP_last_modified) ! above
                !print*, 'IMFP file last modified on:', IMFP_last_modified
                if (IMFP_last_modified < NumPar%Last_mod_time_CDF) then
                    NumPar%redo_EMFP = .true. ! Material parameters changed, recalculate EMFPs
                    print*, 'File with CDF was modified more recently than the MFP => recalculating MFP'
                endif
                ! Check if the file is consistent with the grid set:
                open(FN2, file=trim(adjustl(Input_elastic_file)), ACTION='READ')
                call count_lines_in_file(FN2, Nsiz) ! module "Dealing_with_EADL"
                !print*, 'Nelast', Nelast, Nsiz

                if (Nsiz /= Nelast) then
                    NumPar%redo_EMFP = .true. ! Grid mismatch, recalculate IMFPs
                    print*, 'Energy grid mismatch in MFP file => recalculating MFP'
                endif
                close(FN2)
            endif

            if (file_exist .and. .not.NumPar%redo_EMFP) then    ! read from the file:
                write(*,'(a,a,a)') 'Mott EMFPs of an electron in ', trim(adjustl(Material_name)), ' are already in the file:'
                write(*, '(a)') trim(adjustl(Input_elastic_file))
                write(*, '(a)') ' '
                open(FN2, file=trim(adjustl(Input_elastic_file)), ACTION='READ')
            else    ! create and write to the file:
                call All_elastic_scattering(Nelast, Target_atoms, CDF_Phonon, Matter, Elastic_MFP%Total, NumPar, Mat_DOS, kind_of_particle)
                open(FN2, file=trim(adjustl(Input_elastic_file)))
                write(*,'(a,a,a)') 'Elastic mean free paths of an electron calculated using Mott formula in ', trim(adjustl(Material_name)), ' are stored in the file'
                write(*, '(a)') trim(adjustl(Input_elastic_file))
                write(*, '(a)') ' '
            endif
         case default ! el_elastic_CS             ! Elastic scattering is disabled
            write(temp_char1, '(f7.2, a)') Matter%temp, '_K'
            temp_char2 = 'OUTPUT_Electron_No_elas_EMFPs.dat'
            Input_elastic_file = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(temp_char2))
            if (allocated(File_names%F)) File_names%F(4) = trim(adjustl(temp_char2))    ! save for later use
            FN2 = 2031
            open(FN2, file=trim(adjustl(Input_elastic_file)))
            write(*,'(a)') 'Electron kinetics will be traces without elastic scatterings on target atoms'
            file_exist = .false.
            call All_elastic_scattering(Nelast, Target_atoms, CDF_Phonon, Matter, Elastic_MFP%Total, NumPar, Mat_DOS, kind_of_particle)
            Elastic_MFP%Total%L(:) = 1.0d30
         endselect ! el_elastic_CS

    !==============================
    else if (kind_of_particle .EQ. 'Hole') then
!         write(temp_char, '(f7.2, a)') Matter%temp, '_K'
!         Input_elastic_file = trim(adjustl(Output_path))//'/OUTPUT_Hole_EMFPs_'//trim(adjustl(temp_char))//'.dat'
!         if (allocated(File_names%F)) File_names%F(5) = 'OUTPUT_Hole_EMFPs_'//trim(adjustl(temp_char))//'.dat' ! save for later use
!         FN2 = 204
!         file_exist = .false.
!         call All_elastic_scattering(Nelast, Target_atoms, CDF_Phonon, Matter, Elastic_MFP%Total, NumPar, Mat_DOS, kind_of_particle)
!         open(FN2, file=trim(adjustl(Input_elastic_file)))
!         write(*,'(a,a,a)') 'Calculated elastic mean free paths of a hole in ', trim(adjustl(Material_name)), ' are stored in the file'
!         write(*, '(a)') trim(adjustl(Input_elastic_file))
        select case (NumPar%kind_of_EMFP)
        case (2)    ! Read DSF elastic MFP
            write(temp_char1, '(f7.2, a)') Matter%temp, '_K'

            temp_char2 = 'OUTPUT_Hole_DSF_EMFPs_'//trim(adjustl(temp_char1))//'.dat'

            Input_elastic_file = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(temp_char2))
            if (allocated(File_names%F)) File_names%F(5) = trim(adjustl(temp_char2)) ! save for later use
            FN2 = 2032

            file_exist = .false.
            write(*,'(a,a,a)') 'Calculated elastic mean free paths of a hole in ', trim(adjustl(Material_name)), ' are stored in the file:'
            write(*, '(a)') trim(adjustl(Input_elastic_file))
            write(*, '(a)') ' '
            open(FN2, file=trim(adjustl(Input_elastic_file)))
            deallocate(Elastic_MFP%Total%E)
            deallocate(Elastic_MFP%Total%L)
            deallocate(Elastic_MFP%Total%dEdx)
            Nelast = size(DSF_DEMFP%E)
            allocate(Elastic_MFP%Total%E(Nelast))
            allocate(Elastic_MFP%Total%L(Nelast))
            allocate(Elastic_MFP%Total%dEdx(Nelast))
            ! And separately absorption and emission of energy:
            if (allocated(Elastic_MFP%Emit%E)) deallocate(Elastic_MFP%Emit%E)
            if (allocated(Elastic_MFP%Emit%L)) deallocate(Elastic_MFP%Emit%L)
            if (allocated(Elastic_MFP%Emit%dEdx)) deallocate(Elastic_MFP%Emit%dEdx)
            if (allocated(Elastic_MFP%Absorb%E)) deallocate(Elastic_MFP%Absorb%E)
            if (allocated(Elastic_MFP%Absorb%L)) deallocate(Elastic_MFP%Absorb%L)
            if (allocated(Elastic_MFP%Absorb%dEdx)) deallocate(Elastic_MFP%Absorb%dEdx)
            allocate(Elastic_MFP%Emit%E(Nelast), source = 0.0d0)
            allocate(Elastic_MFP%Emit%L(Nelast), source = 0.0d0)
            allocate(Elastic_MFP%Emit%dEdx(Nelast), source = 0.0d0)
            allocate(Elastic_MFP%Absorb%E(Nelast), source = 0.0d0)
            allocate(Elastic_MFP%Absorb%L(Nelast), source = 0.0d0)
            allocate(Elastic_MFP%Absorb%dEdx(Nelast), source = 0.0d0)

            do i = 1, Nelast
                Elastic_MFP%Total%E(i) = DSF_DEMFP(i)%E
                Elastic_MFP%Total%L(i) = DSF_DEMFP(i)%dL(size(DSF_DEMFP(1)%dL)) ! Total
                Elastic_MFP%Emit%L(i)   = DSF_DEMFP(i)%dL_emit(size(DSF_DEMFP(1)%dL_emit)) ! Emission
                Elastic_MFP%Absorb%L(i) = DSF_DEMFP(i)%dL_absorb(size(DSF_DEMFP(1)%dL_absorb)) ! Absorption
                !print*, i, Elastic_MFP%Total%E(i), Elastic_MFP%Total%L(i), &
                !Elastic_MFP%Emit%L(i), Elastic_MFP%Absorb%L(i) &
                !1.0d0/(1.0d0/Elastic_MFP%Emit%L(i) + 1.0d0/Elastic_MFP%Absorb%L(i))
            enddo
            !pause 'DSF Test VB hole'
         case (1) ! Calculate elastic MFP
            KCS = ''
            select case (NumPar%kind_of_CDF_ph)
            case (0)    ! Ritchie-Howie
                KCS = '_CDF'
            case (1)    ! single-pole CDF
                KCS = '_spCDF'
            endselect

            write(temp_char, '(f7.2, a)') Matter%temp, '_K'

            temp_char2 = 'OUTPUT_Hole'//trim(adjustl(KCS))//'_EMFPs_'//trim(adjustl(temp_char))//'.dat'

            Input_elastic_file = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(temp_char2))
            if (allocated(File_names%F)) then
                File_names%F(5) = trim(adjustl(temp_char2)) ! save for later use
            endif
            FN2 = 204
            !file_exist = .false.
            inquire(file=trim(adjustl(Input_elastic_file)),exist=file_exist)    ! check if input file excists

            ! Check, if file with MFP was created with paramters in actual CDF file, find out when this file was last modified:
            if (file_exist) then
                call get_file_stat(trim(adjustl(Input_elastic_file)), Last_modification_time=IMFP_last_modified) ! above
                !print*, 'IMFP file last modified on:', IMFP_last_modified
                if (IMFP_last_modified < NumPar%Last_mod_time_CDF) then
                    NumPar%redo_EMFP = .true. ! Material parameters changed, recalculate EMFPs
                    print*, 'File with CDF was modified more recently than the MFP => recalculating MFP'
                endif
                ! Check if the file is consistent with the grid set:
                open(FN2, file=trim(adjustl(Input_elastic_file)), ACTION='READ')
                call count_lines_in_file(FN2, Nsiz) ! module "Dealing_with_EADL"
                !print*, 'Nelast', Nelast, Nsiz

                if (Nsiz /= Nelast) then
                    NumPar%redo_EMFP = .true. ! Grid mismatch, recalculate IMFPs
                    print*, 'Energy grid mismatch in hole elastic MFP file => recalculating MFP'
                endif
                close(FN2)
            endif

            if (file_exist .and. .not.NumPar%redo_EMFP) then    ! read from the file:
                write(*,'(a,a,a)') 'Calculated with '//trim(adjustl(KCS(2:)))//' EMFPs of a hole in ', &
                    trim(adjustl(Material_name)), ' are already in the file:'
                write(*, '(a)') trim(adjustl(Input_elastic_file))
                write(*, '(a)') ' '
                open(FN2, file=trim(adjustl(Input_elastic_file)), ACTION='READ')
            else    ! create and write to the file:
                call All_elastic_scattering(Nelast, Target_atoms, CDF_Phonon, Matter, Elastic_MFP%Total, NumPar, Mat_DOS, kind_of_particle)
                open(FN2, file=trim(adjustl(Input_elastic_file)))
                write(*,'(a,a,a)') 'Calculated elastic mean free paths of a hole in ', trim(adjustl(Material_name)), ' are stored in the file'
                write(*, '(a)') trim(adjustl(Input_elastic_file))
                write(*, '(a)') ' '
            endif

         case (0) ! Mott cross-sections
            write(temp_char1, '(f7.2, a)') Matter%temp, '_K'

            temp_char2 = 'OUTPUT_Hole_Mott_EMFPs_'//trim(adjustl(temp_char1))//'.dat'

            Input_elastic_file = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(temp_char2))
            if (allocated(File_names%F)) File_names%F(5) = trim(adjustl(temp_char2)) ! save for later use
            FN2 = 2043
            inquire(file=trim(adjustl(Input_elastic_file)),exist=file_exist)    ! check if input file excists

            ! Check, if file with MFP was created with paramters in actual CDF file, find out when this file was last modified:
            if (file_exist) then
                call get_file_stat(trim(adjustl(Input_elastic_file)), Last_modification_time=IMFP_last_modified) ! above
                !print*, 'IMFP file last modified on:', IMFP_last_modified
                if (IMFP_last_modified < NumPar%Last_mod_time_CDF) then
                    NumPar%redo_EMFP = .true. ! Material parameters changed, recalculate EMFPs
                    print*, 'File with CDF was modified more recently than the MFP => recalculating MFP'
                endif
                ! Check if the file is consistent with the grid set:
                open(FN2, file=trim(adjustl(Input_elastic_file)), ACTION='READ')
                call count_lines_in_file(FN2, Nsiz) ! module "Dealing_with_EADL"
                if (Nsiz /= Nelast) then
                    NumPar%redo_EMFP = .true. ! Grid mismatch, recalculate IMFPs
                    print*, 'Energy grid mismatch hole elastic  MFP file => recalculating MFP'
                endif
                close(FN2)
            endif


            if (file_exist .and. .not.NumPar%redo_EMFP) then    ! read from the file:
                write(*,'(a,a,a)') 'Mott EMFPs of a hole in ', trim(adjustl(Material_name)), ' are already in the file:'
                write(*, '(a)') trim(adjustl(Input_elastic_file))
                write(*, '(a)') ' '
                open(FN2, file=trim(adjustl(Input_elastic_file)), ACTION='READ')
            else    ! create and write to the file:
                call All_elastic_scattering(Nelast, Target_atoms, CDF_Phonon, Matter, Elastic_MFP%Total, NumPar, Mat_DOS, kind_of_particle)
                open(FN2, file=trim(adjustl(Input_elastic_file)))
                write(*,'(a,a,a)') 'Elastic mean free paths of an hole calculated using Mott formulae is in ', trim(adjustl(Material_name)), ' are stored in the file'
                write(*, '(a)') trim(adjustl(Input_elastic_file))
                write(*, '(a)') ' '
            endif
         case default ! disabled
            write(temp_char1, '(f7.2, a)') Matter%temp, '_K'

            temp_char2 = 'OUTPUT_Hole_No_elas_EMFPs.dat'

            Input_elastic_file = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(temp_char2))
            if (allocated(File_names%F)) File_names%F(4) = trim(adjustl(temp_char2)) ! save for later use
            FN2 = 2031
            open(FN2, file=trim(adjustl(Input_elastic_file)))
            write(*,'(a)') 'Valence hole kinetics will be traces without elastic scatterings on target atoms'
            file_exist = .false.
            call All_elastic_scattering(Nelast, Target_atoms, CDF_Phonon, Matter, Elastic_MFP%Total, NumPar, Mat_DOS, kind_of_particle)
            Elastic_MFP%Total%L(:) = 1.0d30
         endselect
    endif kind_of_part2
    

    ! Now write the elastic output into (or read from) the file:
    if (kind_of_particle .NE. 'Photon') then

        select case (NumPar%kind_of_EMFP)

        case (2)    ! DSF: resolved emission vs absorption
         do i = 1, Nelast
           if (.not. file_exist .or. NumPar%redo_EMFP) then  ! if file didn't exist and we just created it:
             write(FN2,'(f,es,es,es)') Elastic_MFP%Total%E(i), Elastic_MFP%Total%L(i), Elastic_MFP%Emit%L(i), Elastic_MFP%Absorb%L(i)
           else
             read(FN2,*, IOSTAT=Reason) Elastic_MFP%Total%E(i), Elastic_MFP%Total%L(i), Elastic_MFP%Emit%L(i), Elastic_MFP%Absorb%L(i)
             call read_file_here(Reason, i, read_well)
             if (.not. read_well) print*, trim(adjustl(Input_elastic_file))
             if (.not. read_well) goto 2016
           endif
         enddo

        case default ! no resolution of emission vs absorption
         do i = 1, Nelast
           if (.not. file_exist .or. NumPar%redo_EMFP) then  ! if file didn't exist and we just created it:
             write(FN2,'(f,e)') Elastic_MFP%Total%E(i), Elastic_MFP%Total%L(i)
           else
             read(FN2,*, IOSTAT=Reason) Elastic_MFP%Total%E(i), Elastic_MFP%Total%L(i)
             call read_file_here(Reason, i, read_well)
             if (.not. read_well) print*, trim(adjustl(Input_elastic_file))
             if (.not. read_well) goto 2016
           endif
         enddo
        end select
    endif
    NumPar%redo_EMFP = .false. ! default it for the next kind of particles
    

    !rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
    ! And only in this case calculate the electron (or hole) range:
    if (do_range) then
       if ((kind_of_particle .EQ. 'Electron') .OR. (kind_of_particle .EQ. 'Hole')) then
         select case (trim(adjustl(kind_of_particle)))
         case ('Electron', 'electron', 'ELECTRON')
            FN3 = 303 ! electron range
            open(FN3, file=trim(adjustl(File_el_range)))
         case ('Hole', 'hole', 'HOLE')
            FN4 = 304 ! hole range
            open(FN4, file=trim(adjustl(File_hole_range)))
         end select
         !print*, trim(adjustl(File_el_range)), ' '//trim(adjustl(File_hole_range))
         do i = 1, N ! all grid-points
          dEdx1 = 0.0d0
          if (allocated(Total_el_MFPs(1)%ELMFP(1)%dEdx)) then
           do j = 1, Nat ! all atoms
              Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
              do k = 1, Nshl
                 if (Total_el_MFPs(j)%ELMFP(k)%dEdx(i) > 0.0d0) then
                    dEdx1 = dEdx1 + Total_el_MFPs(j)%ELMFP(k)%dEdx(i) ! sum over all shell contributions
                 endif
              enddo
           enddo
           !Egap = Target_atoms(1)%Ip(size(Target_atoms(1)%Ip))   ! band gap [eV]
           if ((dEdx1 > 0.0d0) .AND. (Total_el_MFPs(1)%ELMFP(1)%E(i) > Target_atoms(1)%Ip(size(Target_atoms(1)%Ip)) + 10.0d0)) then
           !if ((dEdx1 > 0.0d0) .AND. (Total_el_MFPs(1)%ELMFP(1)%E(i) > 10.0d0)) then
              if (dEdx0 > 0.0d0) then
                 e_range = e_range + 0.5d0*(1.0d0/dEdx1 + 1.0d0/dEdx0)*(Total_el_MFPs(1)%ELMFP(1)%E(i) - Total_el_MFPs(1)%ELMFP(1)%E(i-1)) ! range as integral of inverse losses
              else
                 e_range = e_range + 0.5d0*(1.0d0/dEdx1)*(Total_el_MFPs(1)%ELMFP(1)%E(i) - Total_el_MFPs(1)%ELMFP(1)%E(i-1)) ! range as integral of inverse losses
              endif
           endif
           !print*, Total_el_MFPs(1)%ELMFP(1)%E(i), dEdx1, dEdx0, e_range
           dEdx0 = dEdx1 ! save for integration
           select case (trim(adjustl(kind_of_particle)))
           case ('Electron', 'electron', 'ELECTRON')
              if (i == 1) write(FN3,'(a,a,a)') 'Energy(eV)	', 'Loss(eV/A)	', 'Range(A)'
              write(FN3,'(e,e,e)') Total_el_MFPs(1)%ELMFP(1)%E(i), dEdx1, e_range 
           case ('Hole', 'hole', 'HOLE')
              if (i == 1) write(FN4,'(a,a,a)') 'Energy(eV)	', 'Loss(eV/A)	', 'Range(A)'
              write(FN4,'(e,e,e)') Total_el_MFPs(1)%ELMFP(1)%E(i), dEdx1, e_range
           end select
          else
           print*, 'Cannot compute ', trim(adjustl(kind_of_particle)), ' range, dEdx-array is not allocated'
          endif
         enddo
         select case (trim(adjustl(kind_of_particle)))
           case ('Electron', 'electron', 'ELECTRON')
              inquire(unit=FN3,opened=file_opened2)    ! check if this file is opened
              if (file_opened2) close(FN3) ! electron range
            case ('Hole', 'hole', 'HOLE')
              inquire(unit=FN4,opened=file_opened2)    ! check if this file is opened
              if (file_opened2) close(FN4) ! hole range
         end select
!          inquire(unit=FN3,opened=file_opened2)    ! check if this file is opened
!          if (file_opened2) close(FN3) ! electron range
!          inquire(unit=FN4,opened=file_opened2)    ! check if this file is opened
!          if (file_opened2) close(FN4) ! hole range
       endif
    endif
    !rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
    

    ! Calculating hole diffusion coefficients:
    if ((kind_of_particle .EQ. 'Hole') .AND. (Matter%Hole_mass .LT. 1.0d5)) then
        open(333, File = trim(adjustl(Output_path))//'/Hole_selfdiffusion_coeff_in_'//trim(adjustl(Material_name))//'.dat')
        write(333, '(A)') 'E[eV]   Eff_Mass[me]   DOS   Selfdiffusion_coeff[cm^2/s]   L_tot[A]'

        do i = 1, size(Mat_DOS%E)
            call linear_approx(Temp_MFP, Mat_DOS%E(i), InelMFP)             ! find inelastic MFP
            call linear_approx(Elastic_MFP%Total%E, Elastic_MFP%Total%L, Mat_DOS%E(i), ElasMFP)
            if (Matter%Hole_mass .LT. 0) then
                Mass = Mat_DOS%Eff_m(i)
            else
                Mass = Matter%Hole_Mass
            endif

            if (InelMFP > 0.0d0) then
                InelMFP_inv = 1.0d0/InelMFP
            else
                InelMFP_inv = 1.0d30
            endif

            if (ElasMFP > 0.0d0) then
                ElasMFP_inv = 1.0d0/ElasMFP
            else
                ElasMFP_inv = 1.0d30
            endif

            L_tot = 1.0d0/(ElasMFP_inv + InelMFP_inv)               ! total MFP [A]
            if (Mass > 1.0d-10) then
                Vel = sqrt(2*Mat_DOS%E(i)*g_e/(Mass*g_me))*1.0d10       ! Velosity [A/s]
            else
                Vel = 0.0d0       ! Velosity [A/s]
            endif
            write(333, '(e,e,e,e,e)') Mat_DOS%E(i), Mass, Mat_DOS%DOS(i), 1.0d0/3.0d0*L_tot*vel/1.0d16, L_tot
        enddo
        close(333)
    endif


    ! Printout CDF file if requested or for sp-approximation:
    if ((kind_of_particle .EQ. 'Electron') .and. &
        ( (NumPar%kind_of_CDF == 1) .or. (NumPar%print_CDF) ) ) then
        FN1 = 3121
        if (NumPar%kind_of_CDF == 1) then
            KCS = '_spCDF'
        else
            KCS = '_CDF'
        endif

        File_name = trim(adjustl(Output_path))//trim(adjustl(NumPar%path_sep))//trim(adjustl(Matter%Target_name)) &
                    //trim(adjustl(KCS))//'.cdf'
        open(unit = FN1, FILE = trim(adjustl(File_name)))

        call write_CDF_file(FN1, Matter%Target_name, Matter%Chem, Matter%Dens, Matter%E_F, &
                            Matter%Egap, Matter%Vsound, Target_atoms, CDF_Phonon) ! below

        inquire(unit=FN1,opened=file_opened)    ! check if this file is opened
        if (file_opened) close(FN1)             ! and if it is, close it
    endif

    
2015    continue
   if (.not. read_well) print*, trim(adjustl(Input_files))
2016   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)     
   if (kind_of_particle .NE. 'Photon') then
      inquire(unit=FN2,opened=file_opened)    ! check if this file is opened
      if (file_opened) close(FN2) 
   endif
end subroutine Analytical_electron_dEdx




subroutine write_CDF_file(FN, Material, Chem, dens, Efermi, Egap, c_sound, Atoms, CDF_Phonon)
   integer, intent(in) :: FN  ! file to write CDF data into
   character(*), intent(in) :: Material   ! material name
   character(*), intent(in) :: Chem       ! chemical formula
   real(8), intent(in) :: dens   ! density of the target [g/cm^3]
   real(8), intent(in) :: Efermi, Egap ! [eV] fermi energy; bandgap [ev]
   real(8), intent(in) :: c_sound  ! speed of sound [m/s]
   type(Atom), dimension(:), intent(in) :: Atoms !Target_atoms
   type(CDF), intent(in) :: CDF_Phonon   ! CDF parameters for phonon to be read from a file
   !-------------------------
   integer :: N_KOA, i, j, k
   character(200) :: line_to_write

   N_KOA = size(Atoms)  ! kinds of atoms

   write(FN,'(a)') trim(adjustl(Material)) ! material name
   write(FN,'(a)') trim(adjustl(Chem))     ! chemical formula
   write(line_to_write,'(f15.6, f15.6, f15.6, f15.6, a)') dens, c_sound, Efermi, Egap, &
                        '   ! density [g/cm^3], speed of sound [m/s], Fermi level [eV], bandgap [eV]'
   write(FN,'(a)') trim(adjustl(line_to_write))

   AT_NUM:do i = 1, N_KOA ! for each kind of atoms:
      write(line_to_write,'(i4,a)') Atoms(i)%N_shl, '   ! number of shells in element: '//trim(adjustl(Atoms(i)%Name))
      write(FN,'(a)') trim(adjustl(line_to_write))

      if (allocated(Atoms(i)%Ritchi)) then
         SH_NUM:do j = 1, Atoms(i)%N_shl  ! for all shells
            write(line_to_write,'(i3, i4, f15.6, f15.6, es15.6E3, a)') size(Atoms(i)%Ritchi(j)%E0), Atoms(i)%Shl_num(j), &
               Atoms(i)%Ip(j), Atoms(i)%Nel(j), Atoms(i)%Auger(j), &
               '  ! Oscillators, designator:'//trim(adjustl(Atoms(i)%Shell_name(j)))//', Ip [eV], Ne, Auger [fs]'
            write(FN,'(a)') trim(adjustl(line_to_write))

            CDF_NUM:do k = 1, size(Atoms(i)%Ritchi(j)%A)  ! for all CDF-functions for this shell
               write(line_to_write,'(f15.6, f15.6, f15.6, a)') Atoms(i)%Ritchi(j)%E0(k), Atoms(i)%Ritchi(j)%A(k), &
               Atoms(i)%Ritchi(j)%Gamma(k), '   ! E0, A, G'
               write(FN,'(a)') trim(adjustl(line_to_write))
            enddo CDF_NUM
         enddo SH_NUM

         ! Phononic part, if used:
         if (allocated(CDF_Phonon%A) .and. (i == N_KOA)) then   ! only after last element
            write(line_to_write,'(i4,a)') size(CDF_Phonon%A), '   ! number of phonon peaks'
            write(FN,'(a)') trim(adjustl(line_to_write))
            do k = 1, size(CDF_Phonon%A)
                write(line_to_write,'(f15.6, f15.6, f15.6, a)') CDF_Phonon%E0(k), CDF_Phonon%A(k), &
                CDF_Phonon%Gamma(k), '   ! E0, A, G'
                write(FN,'(a)') trim(adjustl(line_to_write))
            enddo
        endif

      else
         write(FN,'(a)') '************************************************'
         write(FN,'(a)') 'The Ritchie-Howie CDF-coefficients are undefined, cannot print them out!'
         write(FN,'(a)') 'Probably, the cdf-file was not specified or specified incorrectly in INPUT.txt'
      endif
   enddo AT_NUM
end subroutine write_CDF_file





subroutine All_elastic_scattering(Nelast, Target_atoms, CDF_Phonon, Matter, Elastic_MFP, NumPar, Mat_DOS, kind_of_particle, dont_do)
    integer, intent(in) :: Nelast   ! number of grid points
    type(Atom), dimension(:), intent(in), target :: Target_atoms  ! all data for target atoms
    type(CDF), intent(in), target :: CDF_Phonon ! CDF parameters of a phonon peak if excist
    type(Solid), intent(inout) :: Matter   ! all material parameters
    type(MFP), intent(inout) :: Elastic_MFP ! elastic mean free path
    type(Flag), intent(in) :: NumPar
    type(Density_of_states), intent(in) :: Mat_DOS
    logical, optional, intent(in) :: dont_do
    character(8), intent(in) :: kind_of_particle
        
    integer i, Va, Ord
    real(8) Ele, EMFP, dEdx, Emin, dE(Nelast)

!     Emin = 1.0d0    ! [eV] we start with this minimum
!     Ord = -2 ! start with 1
!     Va = int(Emin)
!     dE(1) = 10.0d0**Ord
!     do i = 1, Nelast-1
!         if (dE(i) .LE. 1.0d0) then
!             if (Va .GE. 10) then
!                 Va = Va - 9
!                 Ord = Ord + 1
!             endif
!         else
!             if (Va .GE. 100) then
!                 Va = Va - 90
!                 Ord = Ord + 1
!             endif
!         endif
!         dE(i+1) = dE(i) + 10.0d0**Ord
!         Va = Va + 1
!     enddo   ! while (dE .LT. Emax)


if (.not.present(dont_do)) then ! only do it when we have the CDF
!$omp parallel &
!$omp private (i, Ele, EMFP, dEdx)
!$omp do schedule(dynamic)
    do i = 1, Nelast
        !Ele = dE(i)
        Ele = Elastic_MFP%E(i)
        !print*, 'All_elastic_scattering', i, Ele
        call Elastic_cross_section(Ele, CDF_Phonon, Target_atoms, Matter, EMFP, dEdx, NumPar, Mat_DOS, kind_of_particle) ! from module Cross_sections
        !Elastic_MFP%E(i) = Ele      ! [eV] energy; Grid was already preset, reuse it!
        Elastic_MFP%L(i) = EMFP     ! [A] elastic mean free path
        Elastic_MFP%dEdx(i) = dEdx  ! [eV/A] energy loss
        call progress(' Progress of calculation: ', i, Nelast)
    enddo
!$omp end do
!$omp end parallel
endif
end subroutine All_elastic_scattering


! Forms all IMFPs for all shells of all atoms, parallelized with openmp
subroutine All_shells_Electron_MFP(N, Target_atoms, Total_el_MFPs, Mat_DOS, Matter, NumPar, kind_of_particle)
    integer, intent(in) :: N  ! number of energy points
    type(Atom), dimension(:), intent(inout), target :: Target_atoms  ! all data for target atoms
    type(All_MFP), dimension(:), allocatable, intent(inout) :: Total_el_MFPs   ! electron mean free paths for all shells
    type(Density_of_states), intent(in) :: Mat_DOS
    type(Solid), intent(in) :: Matter ! properties of material
    type(Flag), intent(in) :: NumPar
    character(8), intent(in) :: kind_of_particle
    !-----------------------
    real(8) IMFP_calc, dEdx, Ele, dE(N), Emin
    integer i, j, k, Nshl, Nat, Va, Ord
    integer Num_th, my_id, OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
    !complex(8) :: complex_CDF ! constructed CDF
    
    Nat = size(Target_atoms)    ! number of atoms
    
    ! If diff.CS tables are required, they need to be allocated:
    if (trim(adjustl(kind_of_particle)) .EQ. 'Electron') then
      select case (NumPar%CS_method)
      case (1)
        do i = 1, N
            Ele = Total_el_MFPs(1)%ELMFP(1)%E(i)
            do j = 1, Nat  ! for all atoms:
                Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
                do k = 1, Nshl  ! for all shells of each atom:

                    ! Save energy point for diff.CS tables:
                    Target_atoms(j)%Int_diff_CS(k)%E(i) = Ele  ! [eV] electron energy on the grid

                    call allocate_diff_CS_tables(Ele, i, Target_atoms, j, k, Matter, Mat_DOS, NumPar, kind_of_particle) ! from module "Cross_sections"
                enddo ! k
            enddo ! k = 1, size(Target_atoms(j)%Ip)  ! for all shells of each atom:
        enddo ! j = 1,size(Target_atoms)  ! for all atoms:
      endselect
    endif ! 'Electron'


!$omp parallel &
!$omp private (i, j, k, Ele, IMFP_calc, dEdx)
!$omp do schedule(dynamic)
    do i = 1, N
        !Ele = dE(i)
        Ele = Total_el_MFPs(1)%ELMFP(1)%E(i)
        !my_id = 1 + OMP_GET_THREAD_NUM() ! identify which thread it is
        do j = 1, Nat  ! for all atoms:

          Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
          do k = 1, Nshl  ! for all shells of each atom:

             call TotIMFP(Ele, i, Target_atoms, j, k, IMFP_calc, dEdx, Matter, Mat_DOS, NumPar, kind_of_particle) ! from module "Cross_sections"
             !Total_el_MFPs(j)%ELMFP(k)%E(i) = Ele
             Total_el_MFPs(j)%ELMFP(k)%L(i) = IMFP_calc
             Total_el_MFPs(j)%ELMFP(k)%dEdx(i) = dEdx

             !print*, j, k, i, Ele, Total_el_MFPs(j)%ELMFP(k)%L(i), Target_atoms(j)%Int_diff_CS(k)%diffCS(i)%dsdhw( size(Target_atoms(j)%Int_diff_CS(k)%diffCS(i)%dsdhw) )
          enddo ! k = 1, size(Target_atoms(j)%Ip)  ! for all shells of each atom:
        enddo ! j = 1,size(Target_atoms)  ! for all atoms:
        call progress(' Progress of calculation: ', i, N)
    enddo
!$omp end do
!$omp end parallel
end subroutine All_shells_Electron_MFP


subroutine All_shells_Photon_MFP(N, Target_atoms, Total_el_MFPs, Matter, NumPar, Mat_DOS) ! calculate all IMFPs
    integer, intent(in) :: N  ! number of energy points
    type(Atom), dimension(:), intent(in), target :: Target_atoms  ! all data for target atoms
    type(All_MFP), dimension(:), allocatable, intent(inout) :: Total_el_MFPs   ! electron mean free paths for all shells
    type(Solid), intent(in) :: Matter ! properties of material
    type(Flag), intent(in) :: NumPar
    type(Density_of_states), intent(in) :: Mat_DOS ! unused here, but needs to be passed to further subroutines
    !----------------------
    real(8) IMFP_calc, dEdx, Ele, dE(N), Emin
    integer i, j, k, Nshl, Nat, Va, Ord
    integer Num_th, my_id, OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
    complex(8) :: complex_CDF ! constructed CDF
    
    Nat = size(Target_atoms)    ! number of atoms
    
!     Emin = 1.0d0    ! [eV] we start with this minimum
!     Ord = 0 ! start with 1
!     Va = int(Emin)
!     dE(1) = Emin
!     do i = 1, N-1
!         if (Va .GE. 100) then
!             Va = Va - 90
!             Ord = Ord + 1
!         endif
!         dE(i+1) = dE(i) + 10.0d0**Ord
!         Va = Va + 1
!     enddo


!$omp parallel &
!$omp private (i, j, k, Ele, IMFP_calc, dEdx, complex_CDF)
!$omp do schedule(dynamic)
    do i = 1, N
        !Ele = dE(i)
        Ele = Total_el_MFPs(1)%ELMFP(1)%E(i)
        !my_id = 1 + OMP_GET_THREAD_NUM() ! identify which thread it is
        do j = 1, Nat  ! for all atoms:
          Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
          do k = 1, Nshl  ! for all shells of each atom:
             call Tot_Phot_IMFP(Ele, Target_atoms, j, k, IMFP_calc, dEdx, Matter, NumPar, Mat_DOS) ! from module "Cross_sections"
             !Total_el_MFPs(j)%ELMFP(k)%E(i) = Ele
             Total_el_MFPs(j)%ELMFP(k)%L(i) = IMFP_calc
             Total_el_MFPs(j)%ELMFP(k)%dEdx(i) = dEdx

             !print*, 'All_shells_Photon_MFP', j, k
             ! test of CDF:
              if ((j == 1) .and. (k==Nshl)) then ! VB
                 call construct_CDF(complex_CDF, Target_atoms, j, k, Mat_DOS, Matter, NumPar, Ele, 0.0d0, photon=.true.) ! module "Cross_sections"
              endif

          enddo ! k = 1, size(Target_atoms(j)%Ip)  ! for all shells of each atom:
        enddo ! j = 1,size(Target_atoms)  ! for all atoms:
        !call progress(' Progress of calculation: ', i, N)
    enddo
!$omp end do
!$omp end parallel
end subroutine All_shells_Photon_MFP


subroutine Analytical_ion_dEdx(Output_path_SHI, Material_name, Target_atoms, SHI, SHI_MFP, Error_message, read_well, NumPar, Matter, Mat_DOS, File_names)
    character(100), intent(in) :: Output_path_SHI   ! path to the folder where the file is/will be stored
    character(100), intent(in) :: Material_name ! name of the material
    type(Atom), dimension(:), intent(in), target :: Target_atoms  ! all data for target atoms
    class(Ion), intent (inout), target :: SHI  ! all the data for the SHI
    type(All_MFP), dimension(:), allocatable, intent(inout) :: SHI_MFP         ! SHI mean free paths for all shells
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    logical, intent(out) :: read_well    ! did we read the file without an error?
    type(All_names), intent(inout) :: File_names
    type(Solid), intent(in) :: Matter
    type(Density_of_states), intent(in) :: Mat_DOS
    type(Flag), intent(inout) :: NumPar
    !------------------------------
    real(8), dimension(:), allocatable :: dEdx_tot, Temp_grid
    real(8) SHI_E, Emin, Emax, dE
    integer N, Ord, Va, IMFP_last_modified, Nsiz
    integer i, j, k, Nat, Nshl, FN, FN2, FN3
    character(100) Input_files, Input_files11, Input_files2, Input_files3, Path_name, command, charge_name, charge_kind, CS_name, ch_temp
    logical file_exist, file_exist2, file_exist3, file_exist4, all_files_exist
    !------------------------------

    read_well = .true.  ! so far so good
    Nat = size(Target_atoms)    ! how many atoms
    if (.not. allocated(SHI_MFP)) allocate(SHI_MFP(size(Target_atoms))) ! that's how many atoms
    ! How many energy points will be here:
    Emin = real(CEILING(((SHI%Mass*g_Mp + g_me)*(SHI%Mass*g_Mp + g_me)/(SHI%Mass*g_Mp*g_me)*Target_atoms(1)%Ip(size(Target_atoms(1)%Ip))/4.0d0)))  ! [eV]
    Emax = 175.6d6/2.0d0*SHI%Mass ! [eV]  maximal energy that still has no relativism
!     dE = Emin
!     call Find_order_of_magn(Emin, i, j)    ! order of magnitude
!     dE = j*10**i  ! start with this value
!     Va = j
!     Ord = i
!     N = 0
!     do while (dE .LT. Emax)
!         N = N + 1
!         if (Va .GE. 10) then
!             Va = Va - 9
!             Ord = Ord + 1
!         endif
!         dE = dE + 10.0d0**Ord
!         Va = Va + 1
!     enddo   ! while (dE .LT. Emax)
!     N = N + 2
    call get_grid_4CS(N, Temp_grid, Target_atoms, Emin, Emax)  ! below
    
    do j = 1, Nat   ! declair variables if they are not yet
        Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
        if (.not. allocated(SHI_MFP(j)%ELMFP)) allocate(SHI_MFP(j)%ELMFP(Nshl)) ! how many shells
        do k = 1, Nshl
            if (.not. allocated(SHI_MFP(j)%ELMFP(k)%E)) then    ! energy array
                allocate(SHI_MFP(j)%ELMFP(k)%E(N), source = Temp_grid)  ! grid already constructed, reuse it
                !SHI_MFP(j)%ELMFP(k)%E = 0.0d0
            endif
            if (.not. allocated(SHI_MFP(j)%ELMFP(k)%L)) then    ! mean free path array
                allocate(SHI_MFP(j)%ELMFP(k)%L(N))
                SHI_MFP(j)%ELMFP(k)%L = 1.0d24
            endif
            if (.not. allocated(SHI_MFP(j)%ELMFP(k)%dEdx)) then ! energy loss
                allocate(SHI_MFP(j)%ELMFP(k)%dEdx(N))
                SHI_MFP(j)%ELMFP(k)%dEdx = 0.0d0
            endif
        enddo
    enddo
    
    Path_name = trim(adjustl(Output_path_SHI))  ! here it should be

#ifndef __GFORTRAN__
   ! for intel fortran compiler:
   inquire(DIRECTORY=trim(adjustl(Path_name)),exist=file_exist)    ! check if input file excists
#else
   ! for gfortran compiler:
   inquire(FILE=trim(adjustl(Path_name)),exist=file_exist)    ! check if input file excists
#endif


    if (.not. file_exist) then  ! create the directory
        command='mkdir '//trim(adjustl(Path_name)) ! to create a folder use this command
        CALL system(command)  ! create the folder
    endif
    
    CS_name = ''
    select case (NumPar%kind_of_CDF)
    case (0)    ! Ritchie-Howie
        CS_name = '_CDF'
    case (1)    ! single-pole CDF
        CS_name = '_spCDF'
    endselect

    SELECT CASE (SHI%Kind_Zeff) ! 0=Barkas; 1=Bohr; 2=Nikolaev-Dmitriev; 3=Schiwietz-Grande; 4 - Fixed value;
    CASE (1)
        charge_name = '_Bohr'
    CASE (2)
        charge_name = '_ND'
    CASE (3)
        charge_name = '_SG'
    CASE (4)
        charge_name = '_fixed'
    CASE DEFAULT
        charge_name = '_Barkas'
    END SELECT
    
    SELECT CASE (SHI%Kind_ion) ! 0=Point-like charge; 1=Brandt-Kitagawa
    CASE (1)
       charge_kind = '_BK'
    CASE DEFAULT
       charge_kind = '_P'
    END SELECT
    
    ch_temp = trim(adjustl(Path_name))//trim(adjustl(NumPar%path_sep))//'OUTPUT_'//trim(adjustl(SHI%Name))//trim(adjustl(CS_name))// &
        trim(adjustl(charge_name))//trim(adjustl(charge_kind))

    Input_files = trim(adjustl(ch_temp))//'_IMFP.dat'
    Input_files2 = trim(adjustl(ch_temp))//'_dEdx.dat'
    Input_files11 = trim(adjustl(ch_temp))//'_effective_charges.dat'
    Input_files3 = trim(adjustl(ch_temp))//'_Range.dat'

    if (allocated(File_names%F)) then   ! name of the file without the path
        File_names%F(6) = 'OUTPUT_'//trim(adjustl(SHI%Name))//trim(adjustl(CS_name))// &
        trim(adjustl(charge_name))//trim(adjustl(charge_kind))
    endif

    inquire(file=trim(adjustl(Input_files)),exist=file_exist)    ! check if input file excists
    inquire(file=trim(adjustl(Input_files2)),exist=file_exist2)    ! check if input file excists
    inquire(file=trim(adjustl(Input_files11)),exist=file_exist3)    ! check if input file excists
    inquire(file=trim(adjustl(Input_files3)),exist=file_exist4)    ! check if input file excists
    if (file_exist .and. file_exist2 .and. file_exist3 .and. file_exist4) then ! all needed files exist
        all_files_exist = .true.
    else ! not all needed files exist
        all_files_exist = .false.
    endif

    ! Check, if file with MFP was created with paramters in actual CDF file, find out when this file was last modified:
    if (file_exist) then
        call get_file_stat(trim(adjustl(Input_files)), Last_modification_time=IMFP_last_modified) ! above
        !print*, 'IMFP file last modified on:', IMFP_last_modified
        if (IMFP_last_modified < NumPar%Last_mod_time_CDF) then
            NumPar%redo_IMFP_SHI = .true. ! Material parameters changed, recalculate IMFPs
            print*, 'File with CDF was modified more recently than the SHI MFP => recalculating MFP'
        endif
        ! Check if the file is consistent with the grid set:
        FN = 200
        open(FN, file=trim(adjustl(Input_files)))   ! just to check the energy grid inside
        call count_lines_in_file(FN, Nsiz) ! module "Dealing_with_EADL"
        if (Nsiz /= N) then
            NumPar%redo_IMFP_SHI = .true. ! Grid mismatch, recalculate IMFPs
            print*, 'Energy grid mismatch in SHI MFP file => recalculating MFP'
        endif
        close(FN)
    endif

    if (all_files_exist .and. .not.NumPar%redo_IMFP_SHI) then
        write(*,'(a,a,a,a,a)') 'IMFP and dEdx of ', SHI%Name ,' in ', trim(adjustl(Material_name)), ' are already in the files:'
        write(*, '(a,a,a)') trim(adjustl(Input_files)), ' and ', trim(adjustl(Input_files2))
        write(*, '(a)') ' ' 
        FN = 200
        open(FN, file=trim(adjustl(Input_files)), status='old', readonly)
        FN2 = 201
        open(FN2, file=trim(adjustl(Input_files2)), status='old', readonly)
        !call read_SHI_MFP(FN, SHI_MFP, read_well)
        call read_SHI_MFP(FN, FN2, Nat, Target_atoms, SHI_MFP, read_well)
        close(FN)
        close(FN2)
        if (.not. read_well) then
            print*, 'Something was wrong in the file, recalculating it...'
            SHI_E = SHI%E   ! just save it for future
            write(*,'(a,a,a,a,a)') 'IMFP and dEdx of ', SHI%Name ,' in ', trim(adjustl(Material_name)), ' will be stored in the files:'
            write(*, '(a,a,a)') trim(adjustl(Input_files)), ' and ', trim(adjustl(Input_files2))
            write(*, '(a)') ' ' 
            call Analytical_SHI_dEdx(Input_files, Input_files2, Input_files11, N, Emin, Emax, SHI, SHI_MFP, Target_atoms, dEdx_tot, Matter, Mat_DOS, NumPar) ! from module Analytical_IMPS / openmp parallelization
            SHI%E = SHI_E ! restore the original value
        endif
    else
        SHI_E = SHI%E   ! just save it for future
        write(*,'(a,a,a,a,a)') 'IMFP and dEdx of ', SHI%Name ,' in ', trim(adjustl(Material_name)), ' will be stored in the files:'
        write(*, '(a,a,a)') trim(adjustl(Input_files)), ' and ', trim(adjustl(Input_files2))
        write(*, '(a)') ' ' 
        call Analytical_SHI_dEdx(Input_files, Input_files2, Input_files11, N, Emin, Emax, SHI, SHI_MFP, Target_atoms, dEdx_tot, Matter, Mat_DOS, NumPar) ! from module Analytical_IMPS / openmp parallelization
        SHI%E = SHI_E ! restore the original value
    endif
    inquire(file=trim(adjustl(Input_files3)),exist=file_exist2)    ! check if file with Ranges excists
    if (.not.file_exist2 .or. NumPar%redo_IMFP_SHI) then  ! if not, create it
        write(*,'(a,a,a,a,a)') 'Ranges of ', SHI%Name ,' in ', trim(adjustl(Material_name)), ' will be stored in the file:'
        write(*, '(a)') trim(adjustl(Input_files3))
        write(*, '(a)') ' '
        call Get_ion_range(Input_files3,N,SHI_MFP,Target_atoms,dEdx_tot) ! calculate ion range out of its energy-loss function
    endif
    !NumPar%redo_IMFP = .false. ! default it for the next kind of particle
end subroutine Analytical_ion_dEdx


subroutine Get_ion_range(Input_files3,N,SHI_MFP,Target_atoms,dEdx_tot) ! calculate ion range out of its energy-loss function
    integer, intent(in) :: N
    character(100), intent(in) :: Input_files3
    type(All_MFP), dimension(:), intent(in) :: SHI_MFP         ! SHI mean free paths for all shells
    real(8), dimension(:), allocatable :: SHI_range  ! total ion range [A]
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    real(8), dimension(:), allocatable, intent(inout) :: dEdx_tot ! total energy loss [eV/A]
    integer i, j, k, Nat, Nshl, FN3
    Nat = size(Target_atoms) ! how many atoms
    allocate(SHI_range(N))
    if (.not.allocated(dEdx_tot)) then
        allocate(dEdx_tot(N))
        do i = 1, N
            dEdx_tot(i) = 0.0d0 ! to sum up for a total IMFP
            do j = 1, Nat
                Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
                do k = 1, Nshl
                    dEdx_tot(i) = dEdx_tot(i) + SHI_MFP(j)%ELMFP(k)%dEdx(i)
                enddo
            enddo
        enddo
    endif
    SHI_range = 0.0d0
    do i = 2,N  ! for all calculated energies
        do j = 2,i  ! coarse integration...
            if (dEdx_tot(j-1) .EQ. 0.0d0) then
                SHI_range(i) = SHI_range(i) + 0.5d0/dEdx_tot(j)*(SHI_MFP(1)%ELMFP(1)%E(j)-SHI_MFP(1)%ELMFP(1)%E(j-1))
            else
                SHI_range(i) = SHI_range(i) + 0.5d0*(1.0d0/dEdx_tot(j)+1.0d0/dEdx_tot(j-1))*(SHI_MFP(1)%ELMFP(1)%E(j)-SHI_MFP(1)%ELMFP(1)%E(j-1))
            endif
        enddo
        call progress(' Progress of calculation: ', i, N)
    enddo
    
    FN3 = 202
    open(FN3, file=trim(adjustl(Input_files3)))
    write(FN3,'(A)')    '# Energy dEdx    Range'
    write(FN3,'(A)')    '# [eV] [eV/A]    [A]'
    do i = 1,N  ! save into file
        write(FN3,'(es,es,es)') SHI_MFP(1)%ELMFP(1)%E(i), dEdx_tot(i), SHI_range(i)
    enddo
    deallocate(SHI_range)
    deallocate(dEdx_tot)
    close(FN3)
end subroutine Get_ion_range


! This version is with parallelization via openmp:
subroutine Analytical_SHI_dEdx(Input_files, Input_files2, Input_files11, N, Emin, Emax, SHI, SHI_MFP, Target_atoms, dEdx_tot, Matter, Mat_DOS, NumPar)   ! calculates dEdx for range of SHI energies
   character(100), intent(in) :: Input_files ! path to the folder where the file is/will be stored
   character(100), intent(in) :: Input_files2
   character(100), intent(in) :: Input_files11
   integer, intent(in) :: N ! number of grid point in the SHI energy
   real(8), intent(in) :: Emin, Emax    ! [eV] limiting energies for SHI dE/dx calculations
   type (Solid), intent(in) :: Matter
   class(Ion), intent (inout), target :: SHI  ! all the data for the SHI
   type(All_MFP), dimension(:), allocatable, intent(inout) :: SHI_MFP         ! SHI mean free paths for all shells
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
   real(8), dimension(:), allocatable, optional, intent(out) :: dEdx_tot    ! total ion energy loss
   type(Density_of_states), intent(in) :: Mat_DOS
   type(Flag), intent(in) :: NumPar
   !-------------------------------
   integer Nat, Nshl, j, i, k, Va, Ord, FN, FN2, FN3, N_grid
   integer Num_th, my_id, OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
   real(8), dimension(N, size(Target_atoms), size(Target_atoms(1)%Ip)) :: SHI_IMFP, SHI_dEdx
   real(8), dimension(N) :: E, Zeff  ! SHI inverse mean free path [1/A], and dEdx [eV/A], Zeff
   real(8), dimension(size(Target_atoms), size(Target_atoms(1)%Ip)) :: IMFP, dEdx
   real(8) IMFP_calc, dEdx_calc
   real(8), pointer :: Z
   type(Ion) :: SHI_1
   real(8), dimension(N) :: dE  ! energy grid [eV]
   real(8), dimension(:), allocatable :: Temp_grid

   Nat = size(Target_atoms) ! how many atoms
   
!    dE = 0.0d0
!    call Find_order_of_magn(Emin, i, j)    ! order of magnitude
!    dE(1) = j*10**i  ! start with this value
!    Va = j
!    Ord = i
!    do k = 1, N-1
!        if (Va .GE. 10) then
!            Va = Va - 9
!            Ord = Ord + 1
!        endif
!        dE(k+1) = dE(k) + 10.0d0**Ord
!        Va = Va + 1
!    enddo   ! while (dE .LT. Emax)
   call get_grid_4CS(N_grid, Temp_grid, Target_atoms, Emin, Emax)  ! below
   if (size(dE) /= size(Temp_grid)) print*, 'Error: grid mismatch in Analytical_SHI_dEdx'
   dE = Temp_grid

   
   SHI_dEdx = 0.0d0
   E = 0.0d0
   SHI_IMFP = 0.0d0
!$omp parallel &
!$omp private (j, SHI_1, IMFP, dEdx)
   SHI_1 = SHI  ! this is a current value defined for each thread, multiple copy of SHI to use
!$omp do schedule(dynamic)  !!reduction( + : E, SHI_dEdx, SHI_IMFP)  
   do j = 1,N
       !my_id = 1 + OMP_GET_THREAD_NUM() ! identify which thread it is
       SHI_1%E = dE(j)
       E(j) =  SHI_1%E  ! save energy to write into file later
       call All_shells_SHI_dEdx(SHI_1, Target_atoms, IMFP, dEdx, Matter, Mat_DOS, NumPar)  ! get dEdx for all shells summed up
       SHI_dEdx(j, :, :) = dEdx(:,:)   ! [eV/A]
       !SHI_IMFP(j, :, :) = 1.0d0/IMFP(:,:) ! [A]
       where (IMFP(:,:) > 1.0d-10)  ! calculated inverse MFP
          SHI_IMFP(j, :, :) = 1.0d0/IMFP(:,:) ! [A]
       elsewhere    ! infinity
          SHI_IMFP(j, :, :) = 1.0d28     ! [A]
       endwhere
       Zeff(j) = SHI_1%Zeff   ! equilibrium charge
       call progress(' Progress of calculation: ', j, N)
   enddo
!$omp end do
!$omp end parallel
   
   FN = 200
   open(FN, file=trim(adjustl(Input_files)))
   FN2 = 201
   open(FN2, file=trim(adjustl(Input_files2)))
   FN3 = 202
   open(FN3, file=trim(adjustl(Input_files11)))
   
   if (present(dEdx_tot)) allocate(dEdx_tot(N)) ! array of total ion energy loss [eV/A]
   SHI_MFP(1)%ELMFP(1)%E(:) = E(:)
   ! Now write the output into the file:
   do i = 1, N
      write(FN,'(e)', advance='no') SHI_MFP(1)%ELMFP(1)%E(i)/1.0d6                      ! MeV
      write(FN2,'(e)', advance='no') SHI_MFP(1)%ELMFP(1)%E(i)/1.0d6
      write(FN3,'(e,e)') SHI_MFP(1)%ELMFP(1)%E(i)/1.0d6, Zeff(i)                        ! Save effective charge
      IMFP_calc = 0.0d0 ! to sum up for a total IMFP
      dEdx_calc = 0.0d0 ! to sum up for a total IMFP
      do j = 1, Nat
         Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
         do k = 1, Nshl
            SHI_MFP(j)%ELMFP(k)%E(i) = E(i)
            SHI_MFP(j)%ELMFP(k)%L(i) = SHI_IMFP(i,j,k)
            if (SHI_MFP(j)%ELMFP(k)%L(i) > 1d30) SHI_MFP(j)%ELMFP(k)%L(i) = 1d30 ! get rid of infunities
            SHI_MFP(j)%ELMFP(k)%dEdx(i) = SHI_dEdx(i,j,k)
            write(FN,'(e)', advance='no') SHI_MFP(j)%ELMFP(k)%L(i)    ! write IMFP for all shells
            write(FN2,'(e)', advance='no') SHI_MFP(j)%ELMFP(k)%dEdx(i)    ! write IMFP for all shells
            IMFP_calc = IMFP_calc + 1.0d0/SHI_MFP(j)%ELMFP(k)%L(i)    ! sum them all up to get the total value
            dEdx_calc = dEdx_calc + SHI_MFP(j)%ELMFP(k)%dEdx(i)    ! sum them all up to get the total value
         enddo
      enddo
      write(FN,'(e)') 1.0d0/IMFP_calc   ! write total IMFP
      write(FN2,'(e)') dEdx_calc        ! write total dEdx
      if (present(dEdx_tot)) dEdx_tot(i) = dEdx_calc ! array of total ion energy loss [eV/A]
   enddo
   close(FN)
   close(FN2)
   close(FN3)
end subroutine Analytical_SHI_dEdx


subroutine Find_order_of_magn(Val, N, Rem)    ! order of magnitude
    real(8), intent(in) :: Val  ! value order of which needs to be found
    integer, intent(out) :: N, Rem   ! order of magnitude found, and the value
    real(8) Val0
    integer i
    Val0 = Val
    N = 0    ! order
    Rem = 0  ! value
    i = MOD(Val0,10.0d0)
    Rem = i  ! value
    do while (Val0 .GT. 10)
       N = N + 1  ! order
       Val0 = Val0/10.0d0
       i = MOD(Val0,10.0d0)
       Rem = i  ! value
    enddo
end subroutine Find_order_of_magn



! This subroutine calculates total SHI dEdx for one given energy:
subroutine All_shells_SHI_dEdx(SHI, Target_atoms, SHI_IMFP, SHI_dEdx, Matter, Mat_DOS, NumPar)
   class(Ion), intent (inout), target :: SHI  ! all the data for the SHI
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
   real(8), dimension(:,:), intent(out) :: SHI_IMFP, SHI_dEdx ! calculated inverse mean free path (cross-section) [1/A], and the energy losses [eV/A]
   type(Solid), intent(in) :: Matter
   type(Density_of_states), intent(in) :: Mat_DOS
   type(Flag), intent(in) :: NumPar
       
   real(8) IMFP, dEdx
   integer Nat, Nshl
   SHI_dEdx = 0.0e0
   SHI_IMFP = 0.0e0
   do Nat = 1, size(Target_atoms)
      do Nshl = 1, size(Target_atoms(Nat)%Ip)
         call SHI_Total_IMFP(SHI, Target_atoms, Nat, Nshl, IMFP, dEdx, Matter, Mat_DOS, NumPar) ! from module Cross_sections
         SHI_dEdx(Nat, Nshl) = SHI_dEdx(Nat, Nshl) + dEdx
         SHI_IMFP(Nat, Nshl) = SHI_IMFP(Nat, Nshl) + IMFP
      enddo
   enddo
end subroutine All_shells_SHI_dEdx


! Interpolation between two values with different methods:
subroutine Interpolate(Iflag, E1, E2, Sigma1, Sigma2, E_needed, OUT_value)
   integer, intent(in) :: Iflag ! what kind of interpolation to use
   real(8), intent(in) :: E1, E2, Sigma1, Sigma2, E_needed  ! input data: X and Y points
   real(8), intent(out) :: OUT_value    ! interpolated value
   real(8) E2log, E1log, E_needed_log, Sigma1log, Sigma2log
   select case(Iflag) ! what interpolation to use:
      case(3)	! logarithmic x, linear y
         E2log = log(E2)
         E1log = log(E1)
         E_needed_log = log(E_needed)
         OUT_value = Sigma1 + (Sigma2 - Sigma1)/(E2log - E1log)*(E_needed_log - E1log)
      case(4)	! linear x, logarithmic y
         Sigma1log = log(Sigma1)
         Sigma2log = log(Sigma2)
         OUT_value = Sigma1log + (Sigma2log - Sigma1log)/(E2 - E1)*(E_needed - E1)
         OUT_value = exp(OUT_value)
      case(5)	! logarithmic x and y
         E2log = log(E2)
         E1log = log(E1)
         E_needed_log = log(E_needed)
         Sigma1log = log(Sigma1)
         Sigma2log = log(Sigma2)
         OUT_value = Sigma1log + (Sigma2log - Sigma1log)/(E2log - E1log)*(E_needed_log - E1log)
         OUT_value = exp(OUT_value)
      case default ! linear x and y
         OUT_value = Sigma1 + (Sigma2 - Sigma1)/(E2 - E1)*(E_needed - E1) 
   end select
end subroutine Interpolate

subroutine progress(string,ndone,ntotal)
    implicit none
    character*(*) string
    character*255 prog,oldprog
    integer ndone,ntotal,i, z
    !save prog
    !z = 10
    if (100.0*ndone/ntotal .GE. 100.0d0) then
        write(0,'(a,$)') '                                                                                   ',char(13)
    else
        write(prog,'(a25,1x,''['')') string
        do i=1,40
            prog(27+i:27+i)=' '
        enddo
        write(prog(44:51),'(f7.2,''%'')') 100.0*ndone/ntotal
        do i=1,40
            if ((1.0*ndone/ntotal).gt.(1.0*i/40)) then
                if (prog(27+i:27+i).eq.' ') prog(27+i:27+i)='-'
            endif
        enddo
        prog(67:67)=']'
        write(0,'(a,a,$)') prog(1:77),char(13)
        return
    endif
end subroutine progress



!GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
! Define a fine grid for MFPs:
subroutine get_grid_4CS(N, grid_array, Target_atoms, Emin_in, Emax_in, rescale_dE)
   integer, intent(out) :: N ! number of grid points (array size to be defined)
   real(8), intent(out), dimension(:), allocatable :: grid_array ! array of these grid points
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
   real(8), intent(in), optional :: Emin_in, Emax_in    ! start and end of the grid [eV]
   real(8), optional :: rescale_dE   ! scaling factor for dE
   !-----------------------------
   real(8), dimension(:), allocatable :: special_point
   real(8) :: Emin, Emax, E_sp_eps, scale_dE
   integer :: i, NP

   ! Initial definitions:
   if (present(Emin_in)) then
      Emin = Emin_in   ! [eV] we start with this given minimum
   else ! default value
      Emin = 0.01d0
   endif

   if (present(Emax_in)) then
      Emax = Emax_in   ! [eV] defaul value, may be changed if needed
   else ! default value
      Emax = 2.0d5   ! [eV] defaul value, may be changed if needed
   endif

   if (present(rescale_dE)) then    ! rescale (refine) grid by this value
      scale_dE = rescale_dE
   else ! no rescaling, default grid
      scale_dE = 1.0d0
   endif

   E_sp_eps = 1.0d-3 ! how close a grid point should be around a special point

   ! Get the special points associated with the ionization potentials of all shells:
   call define_special_points(Target_atoms, special_point)  ! below

   ! Count how many points, to allocate the grid:
   call go_thru_grid(Emin, Emax, E_sp_eps, special_point, scale_dE, Ngrid=N)   ! below

   ! Save the grid:
   allocate(grid_array(N), source = 0.0d0)
   call go_thru_grid(Emin, Emax, E_sp_eps, special_point, scale_dE, array=grid_array)   ! below

!    do i = 1, size(grid_array)
!       print*, i, grid_array(i)
!    enddo
!    pause 'get_grid_4CS'

    if (allocated(special_point)) deallocate(special_point)
end subroutine get_grid_4CS


subroutine go_thru_grid(Emin, Emax, E_sp_eps, special_point, scale_dE, Ngrid, array)
   real(8), intent(in) :: Emin, Emax, E_sp_eps   ! grid start and end; precision around a special point
   real(8), dimension(:), intent(in) :: special_point
   real(8), intent(in) :: scale_dE   ! scaling factor for dE
   integer, intent(inout), optional :: Ngrid ! number of grid points
   real(8), dimension(:), intent(inout), optional :: array ! save grid
   !--------------
   integer :: SP_count, N
   real(8) :: E_cur, dE, dE_min
   logical :: point_is_here

   SP_count = 1
   N = 0
   dE_min = 0.01d0 * scale_dE
   E_cur = Emin - dE_min  ! start from min
   do while (E_cur < Emax)
      N = N + 1   ! count points
      if (E_cur < 0.1d0-(dE_min)*0.5d0) then
         dE = dE_min
      elseif (E_cur < 1.0d0-dE_min*0.5d0) then
         !dE = 0.1d0
         dE = dE_min * 10.0d0
      else if (E_cur < 100.0d0-dE_min*0.5d0) then
         !dE = 1.0d0
         dE = dE_min * 100.0d0
      else
         dE = 10.0d0**(find_order_of_number(E_cur)-2) ! below
      endif
      E_cur = E_cur + dE

      ! save grid points:
      if (present(array)) then
         if (N <= size(array)) then
            array(N) = E_cur
         else
            print*, 'Mismatch #1 of array size in go_thru_grid:', N, size(array)
            return
         endif
      endif

      ! Check if the special points are still there:
      point_is_here = .false. ! by default
      if (SP_count <= size(special_point)) then ! there may be a specila point:
         if ( (E_cur >= special_point(SP_count)) .and. ((E_cur-dE) < special_point(SP_count)) ) then ! special point inside interval
            SP_count = SP_count + 1 ! this special point is done, do the next one
            point_is_here = .true.
         endif
      endif

      ! Save grid point:
      if (point_is_here) then
         if (present(array)) then
            if (N+2 <= size(array)) then
               ! add two points:
               ! 1) below the special point:
               N = N + 1
               array(N) = special_point(SP_count-1) - E_sp_eps
               ! 2) above the special point:
               N = N + 1
               array(N) = special_point(SP_count-1) + E_sp_eps
            else
               print*, 'Mismatch #2 of array size in go_thru_grid', N, size(array)
               return
            endif
         else
            ! add 2 extra points around the special point
            N = N + 2
         endif
      endif

   enddo ! while (E_cur .LT. Emax)

   ! prepare output:
   if (present(Ngrid)) Ngrid = N ! save grid size

   if (present(array)) then   ! make sure the array is sorted increasing
      call sort_array(array)  ! below
   endif
end subroutine go_thru_grid



subroutine define_special_points(Target_atoms, special_point)
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
   real(8), dimension(:), allocatable :: special_point
   !-----------------------------
   integer :: Nsiz, i, k, sh, coun

   ! Count how many Ip's are there:
   Nsiz = 0
   do i = 1, size(Target_atoms)  ! for all elements
      Nsiz = Nsiz + size(Target_atoms(i)%Ip)
   enddo

   ! allocate and set special_points:
   allocate(special_point(Nsiz), source = 0.0d0)

   ! Copy the special points (ionization potentials):
   coun = 0
   do i = 1, size(Target_atoms)  ! for all elements
      sh = size(Target_atoms(i)%Ip)
      do k = 1, sh
         coun = coun + 1
         special_point(coun) = Target_atoms(i)%Ip(k)
      enddo
   enddo

   ! Sort the array increasing:
   call sort_array(special_point)   ! below

   ! Exclude copies of points, if any:
   call exclude_doubles(special_point) ! below

!    do i = 1, size(special_point)  ! for all elements
!       print*, i, special_point(i)
!    enddo
end subroutine define_special_points


pure subroutine exclude_doubles(array)
   real(8), dimension(:), allocatable, intent(inout) :: array
   !-----------
   real(8), dimension(:), allocatable :: array_copy
   integer, dimension(:), allocatable :: ind
   integer :: i, k, Nsiz, Ndub, coun, coun_ind

   Nsiz = size(array)
   allocate(array_copy(Nsiz), source = array)   ! copy array
   allocate(ind(Nsiz), source = 0)  ! indices of the points with doubles, to exclude

   ! Mark doubles:
   coun_ind = 0
   do i = 1, Nsiz-1
      do k = i+1, Nsiz
         if (array(i) == array(k)) then ! exclude this point
            coun_ind = coun_ind + 1
            ind(coun_ind) = min(i,k) ! this element to be excluded
         endif
      enddo
   enddo

   ! count doubles:
   Ndub = count(ind > 0)

   ! exclude doubles:
   deallocate(array)
   allocate(array(Nsiz-Ndub), source = 0.0d0)   ! make it smaller, excluding doubles
   coun = 0
   coun_ind = 1
   do i = 1, Nsiz
      if (i == ind(coun_ind)) then  ! it's a double, exclude
         coun_ind = coun_ind + 1
      else  ! it's not a double, save
         coun = coun + 1
         array(coun) = array_copy(i)
      endif
   enddo

   ! Clean up:
   deallocate(array_copy, ind)
end subroutine exclude_doubles



subroutine sort_array_r(array_in)  ! bubble sorting algorithm for real 1d array
   real(8), dimension(:), intent(inout) :: array_in
   real(8) :: temp
   integer N,i,j
   logical :: swapped
   N = size(array_in)
   do j = N-1, 1, -1
      swapped = .false. ! nothing swapped at the start
      do i = 1, j
         if (array_in(i) > array_in(i+1)) then ! swap elements
            temp = array_in(i)
            array_in(i) = array_in(i+1)
            array_in(i+1) = temp
            swapped = .true.  ! at least one pair of elements needed swapping
         end if
      enddo
      if (.not. swapped) exit
   enddo
end subroutine sort_array_r

subroutine sort_array_c(array_in)  ! bubble sorting algorithm for complex 1d array (sorted by real part)
   complex, dimension(:), intent(inout) :: array_in
   complex :: temp
   integer N,i,j
   logical :: swapped
   N = size(array_in)
   do j = N-1, 1, -1
      swapped = .false. ! nothing swapped at the start
      do i = 1, j
         if (real(array_in(i)) > real(array_in(i+1))) then ! swap elements
            temp = array_in(i)
            array_in(i) = array_in(i+1)
            array_in(i+1) = temp
            swapped = .true.  ! at least one pair of elements needed swapping
         end if
      enddo
      if (.not. swapped) exit
   enddo
end subroutine sort_array_c



pure function find_order_of_number_real(num)
   integer find_order_of_number_real
   real(8), intent(in) :: num
   character(64) :: temp
   !--------------------
   real(8) :: r_tmp
   !print*, num
   if (num > 1e9) then
      r_tmp = num/1.0e8 ! order is too high, split it
      write(temp,'(i16)') ceiling(r_tmp) ! make it a string
      find_order_of_number_real = 8 + LEN(TRIM(adjustl(temp))) ! find how many characters in this string
   else
      write(temp,'(i12)') ceiling(num) ! make it a string
      find_order_of_number_real = LEN(TRIM(adjustl(temp))) ! find how many characters in this string
   endif
   !print*, num, find_order_of_number_real
end function find_order_of_number_real

pure function find_order_of_number_int(num)
   integer find_order_of_number_int
   integer, intent(in) :: num
   character(64) :: temp
   write(temp,'(i16)') num ! make it a string
   find_order_of_number_int = LEN(TRIM(adjustl(temp))) ! find how many characters in this string
end function find_order_of_number_int



end module Analytical_IMFPs

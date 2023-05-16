!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains all subroutines needed to calculate cross-sections:

module Analytical_IMFPs
  use Universal_Constants   ! let it use universal constants
  use Objects   ! since it uses derived types, it must know about them from module 'Objects'
  use Reading_files_and_parameters, only : get_file_stat, Find_in_array_monoton, read_file_here, Linear_approx, read_SHI_MFP
  use Cross_sections, only : Elastic_cross_section, TotIMFP, Tot_Phot_IMFP, SHI_Total_IMFP
  use Dealing_with_EADL, only : Count_lines_in_file
implicit none
PRIVATE


public :: Analytical_electron_dEdx, Analytical_ion_dEdx, Interpolate



contains

! Calculates electron mean free paths with parallelization via openmp:
subroutine Analytical_electron_dEdx(Output_path, Material_name, Target_atoms, CDF_Phonon, Matter, Total_el_MFPs, &
                        Elastic_MFP, Error_message, read_well, DSF_DEMFP, Mat_DOS, NumPar, kind_of_particle, File_names)
    character(100), intent(in) :: Output_path   ! path to the folder where the file is/will be storred
    character(100), intent(in) :: Material_name ! name of the material
    type(Atom), dimension(:), intent(in), target :: Target_atoms  ! all data for target atoms
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
        
    integer :: FN      ! file number where to save the output
    integer :: FN2     ! file number where to save the output
    integer :: FN3, FN4     ! file numbers where to save the output
    integer :: N ! = 300       ! how many points for IMFP
    integer :: Nelast ! = 318  ! how many point for EMFP
    integer temp(1)
    real(8) Ele, IMFP_calc, dEdx, dEdx1, dEdx0, dE, Emin, Emax, L_tot, vel, Mass, InelMFP, ElasMFP, e_range
    real(8), dimension(:,:), allocatable :: Temp_MFP
    integer Num_th, my_id, OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS, IMFP_last_modified
    integer i, j, k, Nat, Nshl, Reason, Va, Ord, Mnum, MFPnum
    character(100) Input_files, Input_elastic_file, File_el_range, File_hole_range
    character(100) temp_char, temp_char1, temp_ch
    character(3) KCS
    logical file_exist, file_opened, file_opened2, do_range
       
    read_well = .true.  ! so far so good
    do_range = .false.  ! don't recalculate electron range by default
    Nat = size(Target_atoms)    ! how many atoms
    
    Emin = 1.0d0    ! [eV] we start with this minimum
    if (kind_of_particle .EQ. 'Electron') then ! it's an electrons
        Emax = 175.6d6*2.0d0/1836.0d0   ! [eV] ~maximum energy where relativism still can be neglected
    else if (kind_of_particle .EQ. 'Hole') then ! it's a hole
        Emax = Maxval(Mat_DOS%E(:)) + 10.0d0 ! [eV] only within the width of VB (or a bit more, just in case...)
    else ! it's a photon
        Emax = 175.6d6*2.0d0/1836.0d0   ! [eV] excluding Bremsstrahlung, photon can't be more energetic than this
    endif   
    N = 0
    Ord = 0 ! start with 1
    Va = int(Emin)
    dE = Emin
    do while (dE .LT. Emax)
        N = N + 1
        if (Va .GE. 100) then
            Va = Va - 90
            Ord = Ord + 1
        endif
        dE = dE + 10.0d0**Ord
        Va = Va + 1
    enddo   ! while (dE .LT. Emax)
    N = N + 1
    
    elast:if (kind_of_particle .NE. 'Photon') then ! elastic only for massive particles:
        Emin = 1.0d0    ! [eV] we start with this minimum
        Nelast = 0
        Ord = -2 ! start with 1
        Va = int(Emin)
        dE = 10.0d0**Ord
        do while (dE .LT. Emax)
            Nelast = Nelast + 1
            if (dE .LE. 1.0d0) then
                if (Va .GE. 10) then
                    Va = Va - 9
                    Ord = Ord + 1
                endif
            else
                if (Va .GE. 100) then
                    Va = Va - 90
                    Ord = Ord + 1
                endif
            endif
            dE = dE + 10.0d0**Ord
            Va = Va + 1
        enddo   ! while (dE .LT. Emax)
        Nelast = Nelast + 1 
        if (.not. allocated(Elastic_MFP%Total%E)) then
            allocate(Elastic_MFP%Total%E(Nelast))  ! [eV] energies for MFP
            Elastic_MFP%Total%E = 0.0d0   ! just to start
        endif
        if (.not. allocated(Elastic_MFP%Total%L)) then
            allocate(Elastic_MFP%Total%L(Nelast))  ! [A] MFP itself
            Elastic_MFP%Total%L = 1.0d24
        endif
        if (.not. allocated(Elastic_MFP%Total%dEdx)) then
            allocate(Elastic_MFP%Total%dEdx(Nelast))  ! [eV/A] mean energy loss
            Elastic_MFP%Total%dEdx = 0.0d0
        endif
    endif elast
        
    if (.not. allocated(Total_el_MFPs)) allocate(Total_el_MFPs(Nat)) ! how many atoms    
    do j = 1, Nat   ! declair variables if they are not yet
        Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
        if (.not. allocated(Total_el_MFPs(j)%ELMFP)) allocate(Total_el_MFPs(j)%ELMFP(Nshl)) ! how many shells
        do k = 1, Nshl
            if (.not. allocated(Total_el_MFPs(j)%ELMFP(k)%E)) then
                allocate(Total_el_MFPs(j)%ELMFP(k)%E(N))    !
                Total_el_MFPs(j)%ELMFP(k)%E = 0.0d0
            endif
            if (.not. allocated(Total_el_MFPs(j)%ELMFP(k)%L)) then
                allocate(Total_el_MFPs(j)%ELMFP(k)%L(N))
                Total_el_MFPs(j)%ELMFP(k)%L = 1.0d24
            endif
            if (.not. allocated(Total_el_MFPs(j)%ELMFP(k)%dEdx)) then
                allocate(Total_el_MFPs(j)%ELMFP(k)%dEdx(N))
                Total_el_MFPs(j)%ELMFP(k)%dEdx = 0.0d0
            endif
        enddo
    enddo
    
    allocate(temp_MFP(2,N))
    
    ! Which cross sections we use - CDF vs BEB:
    temp = 0
    KCS = 'CDF'
    do i = 1, size(Target_atoms) ! all atoms
       temp(1) = maxval(Target_atoms(i)%KOCS(:)) ! check if there is at least one shell for which we use BEB instead of CDF
       if (temp(1) .GE. 2) then ! all CDF
          KCS = 'BEB'
          exit ! if there is at least one BEB cross section name the file 'BEB'
       endif
    enddo
    
    kind_of_particle1:if (kind_of_particle .EQ. 'Electron') then
        write(temp_ch, '(f9.2)') Matter%temp
        if(Matter%El_eff_mass .EQ. 0) then
            write(temp_char, '(a, a, a)') 'DOS_', trim(adjustl(temp_ch)), '_K'
        else if ((Matter%El_eff_mass .LT. 0) .OR. (Matter%El_eff_mass .EQ. 1.0d0)) then
            write(temp_char, '(a, a, a)') 'Me_', trim(adjustl(temp_ch)), '_K'
        else
            write(temp_char, '(f5.2,a,a,a)') Matter%El_eff_mass, '_me_', trim(adjustl(temp_ch)), '_K'
        endif

        if (KCS .EQ. 'BEB') then ! BEB vs CDF cross section:
          Input_files = trim(adjustl(Output_path))//'/OUTPUT_Electron_IMFPs_'//trim(adjustl(KCS))//'.dat'
          if (allocated(File_names%F)) File_names%F(2) = '/OUTPUT_Electron_IMFPs_'//trim(adjustl(KCS))//'.dat' ! save for later use
          File_el_range = trim(adjustl(Output_path))//'/OUTPUT_Electron_range_'//trim(adjustl(KCS))//'.dat'
        else ! CDF:
          if (NumPar%kind_of_DR .EQ. 4) then    ! Delta-CDF
            Input_files = trim(adjustl(Output_path))//'/OUTPUT_Electron_IMFPs_Delta_'//trim(adjustl(temp_char))//'.dat'
            if (allocated(File_names%F)) File_names%F(2) = '/OUTPUT_Electron_IMFPs_Delta_'//trim(adjustl(temp_char))//'.dat' ! save for later use
            File_el_range = trim(adjustl(Output_path))//'/OUTPUT_Electron_Delta_range_'//trim(adjustl(temp_char))//'.dat'
          else if (NumPar%kind_of_DR .EQ. 3) then
            Input_files = trim(adjustl(Output_path))//'/OUTPUT_Electron_IMFPs_Ritchie_'//trim(adjustl(temp_char))//'.dat'
            if (allocated(File_names%F)) File_names%F(2) = '/OUTPUT_Electron_IMFPs_Ritchie_'//trim(adjustl(temp_char))//'.dat' ! save for later use
            File_el_range = trim(adjustl(Output_path))//'/OUTPUT_Electron_Ritchie_range_'//trim(adjustl(temp_char))//'.dat'
          else if (NumPar%kind_of_DR .EQ. 2) then
            Input_files = trim(adjustl(Output_path))//'/OUTPUT_Electron_IMFPs_Plasmon_pole_'//trim(adjustl(temp_char))//'.dat'
            if (allocated(File_names%F)) File_names%F(2) = '/OUTPUT_Electron_IMFPs_Plasmon_pole_'//trim(adjustl(temp_char))//'.dat' ! save for later use
            File_el_range = trim(adjustl(Output_path))//'/OUTPUT_Electron_Plasmon_pole_range_'//trim(adjustl(temp_char))//'.dat'
          else
            Input_files = trim(adjustl(Output_path))//'/OUTPUT_Electron_IMFPs_Free_'//trim(adjustl(temp_char))//'.dat'
            if (allocated(File_names%F)) File_names%F(2) = '/OUTPUT_Electron_IMFPs_Free_'//trim(adjustl(temp_char))//'.dat' ! save for later use
            File_el_range = trim(adjustl(Output_path))//'/OUTPUT_Electron_range_Free_'//trim(adjustl(temp_char))//'.dat'
          endif
        endif ! which name

        ! IMFP files:
        FN = 201
        inquire(file=trim(adjustl(Input_files)),exist=file_exist)    ! check if input file excists

        ! Check, if file with MFP was created with paramters in actual CDF file, find out when this file was last modified:
        if (file_exist) then
            call get_file_stat(trim(adjustl(Input_files)), Last_modification_time=IMFP_last_modified) ! above
            !print*, 'IMFP file last modified on:', IMFP_last_modified
            if (IMFP_last_modified < NumPar%Last_mod_time_CDF) NumPar%redo_IMFP = .true. ! Material parameters changed, recalculate IMFPs
        endif

        if (file_exist .and. .not.NumPar%redo_IMFP) then    ! read from the file:
            write(*,'(a,a,a)') 'IMFPs of an electron in ', trim(adjustl(Material_name)), ' are already in the file:'
            write(*, '(a)') trim(adjustl(Input_files))
            write(*, '(a)') ' '
            open(FN, file=trim(adjustl(Input_files)), ACTION='READ')
        else    ! create and write to the file:
            call All_shells_Electron_MFP(N, Target_atoms, Total_el_MFPs, Mat_DOS, Matter, NumPar, kind_of_particle) ! calculate all IMFPs
            open(FN, file=trim(adjustl(Input_files)))
            write(*,'(a,a,a)') 'Calculated inelastic mean free paths of an electron in ', trim(adjustl(Material_name)), ' are storred in the file'
            write(*, '(a)') trim(adjustl(Input_files))
            write(*, '(a)') ' '
        endif
    else if (kind_of_particle .EQ. 'Hole') then
        FN = 202
        if (KCS .EQ. 'BEB') then ! BEB vs CDF cross section:
           Input_files = trim(adjustl(Output_path))//'/OUTPUT_Hole_IMFPs_BEB.dat'
           if (allocated(File_names%F)) File_names%F(3) = 'OUTPUT_Hole_IMFPs_BEB.dat' ! save for later use
           File_hole_range = trim(adjustl(Output_path))//'/OUTPUT_Hole_range_BEB.dat'
        else
           write(temp_char, '(f7.2, a)') Matter%temp, '_K'
           Input_files = trim(adjustl(Output_path))//'/OUTPUT_Hole_IMFPs_CDF_'//trim(adjustl(temp_char))//'.dat'
           if (allocated(File_names%F)) File_names%F(3) = 'OUTPUT_Hole_IMFPs_CDF_'//trim(adjustl(temp_char))//'.dat' ! save for later use
           File_hole_range = trim(adjustl(Output_path))//'/OUTPUT_Hole_range_CDF_'//trim(adjustl(temp_char))//'.dat'
        endif
        file_exist = .false.                !This file must be overwriten before each calculation.
        call All_shells_Electron_MFP(N, Target_atoms, Total_el_MFPs, Mat_DOS, Matter, NumPar, kind_of_particle) ! calculate all IMFPs
        open(FN, file=trim(adjustl(Input_files)))
        write(*,'(a,a,a)') 'Calculated inelastic mean free paths of a hole in ', trim(adjustl(Material_name)), ' are storred in the file'
        write(*, '(a)') trim(adjustl(Input_files))
        write(*, '(a)') ' '
    else kind_of_particle1 ! photon
        !if (KCS .EQ. 'BEB') then ! BEB vs CDF cross section:
        if (.not.allocated(Target_atoms(1)%Ritchi(size(Target_atoms))%E0)) then ! no CDF known, use EPDL97 database:
           Input_files = trim(adjustl(Output_path))//'/OUTPUT_Photon_IMFPs_EPDL.dat'
           if (allocated(File_names%F)) File_names%F(7) = 'OUTPUT_Photon_IMFPs_EPDL.dat' ! save for later use
        else
           write(temp_char, '(f7.2, a)') Matter%temp, '_K'
           Input_files = trim(adjustl(Output_path))//'/OUTPUT_Photon_IMFPs_CDF_'//trim(adjustl(temp_char))//'.dat'
           if (allocated(File_names%F)) File_names%F(7) = 'OUTPUT_Photon_IMFPs_CDF_'//trim(adjustl(temp_char))//'.dat' ! save for later use
        endif
        inquire(file=trim(adjustl(Input_files)),exist=file_exist)    ! check if input file excists

        ! Check, if file with MFP was created with paramters in actual CDF file, find out when this file was last modified:
        if (file_exist) then
            call get_file_stat(trim(adjustl(Input_files)), Last_modification_time=IMFP_last_modified) ! above
            !print*, 'IMFP file last modified on:', IMFP_last_modified
            if (IMFP_last_modified < NumPar%Last_mod_time_CDF) NumPar%redo_IMFP = .true. ! Material parameters changed, recalculate IMFPs
        endif

        if (file_exist .and. .not.NumPar%redo_IMFP) then    ! read from the file:
            write(*,'(a,a,a)') 'IMFPs of a photon in ', trim(adjustl(Material_name)), ' are already in the file:'
            write(*, '(a)') trim(adjustl(Input_files))
            write(*, '(a)') ' '
            open(newunit=FN, file=trim(adjustl(Input_files)), ACTION='READ')
        else    ! create and write to the file:
            call All_shells_Photon_MFP(N, Target_atoms, Total_el_MFPs, Matter, NumPar) ! calculate all IMFPs
            open(newunit=FN, file=trim(adjustl(Input_files)))
            write(*,'(a,a,a)') 'Calculated attenuation lengths (IMFPs) of a photon in ', trim(adjustl(Material_name)), ' are storred in the file'
            write(*, '(a)') trim(adjustl(Input_files))
            write(*, '(a)') ' '
        endif
    endif kind_of_particle1

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
   NumPar%redo_IMFP = .false. ! defualt it for the next kind of particle

   !######################### Now do the same for elastic mean free path of an electron and hole:
   kind_of_part2:if (kind_of_particle .EQ. 'Electron') then
         0001 continue
         select case (NumPar%kind_of_EMFP) ! el_elastic_CS
         case (2) ! Read DSF elastic MFP
            write(temp_char1, '(f7.2, a)') Matter%temp, '_K'
            Input_elastic_file = trim(adjustl(Output_path))//'/OUTPUT_Electron_DSF_EMFPs_'//trim(adjustl(temp_char1))//'.dat'
            if (allocated(File_names%F)) File_names%F(4) = 'OUTPUT_Electron_DSF_EMFP_'//trim(adjustl(temp_char1))//'.dat' ! save for later use
            FN2 = 2032
            inquire(file=trim(adjustl(Input_elastic_file)),exist=file_exist)    ! check if input file excists

            ! Check, if file with MFP was created with paramters in actual CDF file, find out when this file was last modified:
            if (file_exist) then
                call get_file_stat(trim(adjustl(Input_elastic_file)), Last_modification_time=IMFP_last_modified) ! above
                !print*, 'IMFP file last modified on:', IMFP_last_modified
                if (IMFP_last_modified < NumPar%Last_mod_time_DSF) NumPar%redo_EMFP = .true. ! Material parameters changed, recalculate EMFPs
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
                write(*,'(a,a,a)') 'Calculated elastic mean free paths of an electron in ', trim(adjustl(Material_name)), ' are storred in the file:'
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
            if (allocated(CDF_Phonon%A)) then
                write(temp_char1, '(f7.2, a)') Matter%temp, '_K'
                Input_elastic_file = trim(adjustl(Output_path))//'/OUTPUT_Electron_CDF_EMFPs_'//trim(adjustl(temp_char1))//'.dat'
                if (allocated(File_names%F)) File_names%F(4) = 'OUTPUT_Electron_CDF_EMFPs_'//trim(adjustl(temp_char1))//'.dat' ! save for later use
                FN2 = 203
                inquire(file=trim(adjustl(Input_elastic_file)),exist=file_exist)    ! check if input file excists

                ! Check, if file with MFP was created with paramters in actual CDF file, find out when this file was last modified:
                if (file_exist) then
                    call get_file_stat(trim(adjustl(Input_elastic_file)), Last_modification_time=IMFP_last_modified) ! above
                    !print*, 'IMFP file last modified on:', IMFP_last_modified
                    if (IMFP_last_modified < NumPar%Last_mod_time_CDF) NumPar%redo_EMFP = .true. ! Material parameters changed, recalculate EMFPs
                endif

                if (file_exist .and. .not.NumPar%redo_EMFP) then    ! read from the file:
                    write(*,'(a,a,a)') 'Calculated with CDF EMFPs of an electron in ', trim(adjustl(Material_name)), ' are already in the file:'
                    write(*, '(a)') trim(adjustl(Input_elastic_file))
                    write(*, '(a)') ' '
                    open(FN2, file=trim(adjustl(Input_elastic_file)), ACTION='READ')
                else    ! create and write to the file:
                    call All_elastic_scattering(Nelast, Target_atoms, CDF_Phonon, Matter, Elastic_MFP%Total, NumPar, Mat_DOS, kind_of_particle)
                    open(FN2, file=trim(adjustl(Input_elastic_file)))
                    write(*,'(a,a,a)') 'Elastic mean free paths of an electron calculated using CDF phonon peaks in ', trim(adjustl(Material_name)), ' are storred in the file'
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
            Input_elastic_file = trim(adjustl(Output_path))//'/OUTPUT_Electron_Mott_EMFPs.dat'
            if (allocated(File_names%F)) File_names%F(4) = 'OUTPUT_Electron_Mott_EMFPs.dat' ! save for later use
            FN2 = 2031
            inquire(file=trim(adjustl(Input_elastic_file)),exist=file_exist)    ! check if input file excists

            ! Check, if file with MFP was created with paramters in actual CDF file, find out when this file was last modified:
            if (file_exist) then
                call get_file_stat(trim(adjustl(Input_elastic_file)), Last_modification_time=IMFP_last_modified) ! above
                !print*, 'IMFP file last modified on:', IMFP_last_modified
                if (IMFP_last_modified < NumPar%Last_mod_time_CDF) NumPar%redo_EMFP = .true. ! Material parameters changed, recalculate EMFPs
            endif

            if (file_exist .and. .not.NumPar%redo_EMFP) then    ! read from the file:
                write(*,'(a,a,a)') 'Mott EMFPs of an electron in ', trim(adjustl(Material_name)), ' are already in the file:'
                write(*, '(a)') trim(adjustl(Input_elastic_file))
                write(*, '(a)') ' '
                open(FN2, file=trim(adjustl(Input_elastic_file)), ACTION='READ')
            else    ! create and write to the file:
                call All_elastic_scattering(Nelast, Target_atoms, CDF_Phonon, Matter, Elastic_MFP%Total, NumPar, Mat_DOS, kind_of_particle)
                open(FN2, file=trim(adjustl(Input_elastic_file)))
                write(*,'(a,a,a)') 'Elastic mean free paths of an electron calculated using Mott formulae in ', trim(adjustl(Material_name)), ' are storred in the file'
                write(*, '(a)') trim(adjustl(Input_elastic_file))
                write(*, '(a)') ' '
            endif
         case default ! el_elastic_CS             ! Elastic scattering is disabled
            write(temp_char1, '(f7.2, a)') Matter%temp, '_K'
            Input_elastic_file = trim(adjustl(Output_path))//'/OUTPUT_Electron_No_elas_EMFPs.dat'
            if (allocated(File_names%F)) File_names%F(4) = 'OUTPUT_Electron_No_elas_EMFPs.dat' ! save for later use
            FN2 = 2031
            open(FN2, file=trim(adjustl(Input_elastic_file)))
            write(*,'(a)') 'Electron kinetics will be traces without elastic scatterings on target atoms'
            file_exist = .false.
            call All_elastic_scattering(Nelast, Target_atoms, CDF_Phonon, Matter, Elastic_MFP%Total, NumPar, Mat_DOS, kind_of_particle)
            Elastic_MFP%Total%L(:) = 1.0d30
         endselect ! el_elastic_CS
    else if (kind_of_particle .EQ. 'Hole') then
!         write(temp_char, '(f7.2, a)') Matter%temp, '_K'
!         Input_elastic_file = trim(adjustl(Output_path))//'/OUTPUT_Hole_EMFPs_'//trim(adjustl(temp_char))//'.dat'
!         if (allocated(File_names%F)) File_names%F(5) = 'OUTPUT_Hole_EMFPs_'//trim(adjustl(temp_char))//'.dat' ! save for later use
!         FN2 = 204
!         file_exist = .false.
!         call All_elastic_scattering(Nelast, Target_atoms, CDF_Phonon, Matter, Elastic_MFP%Total, NumPar, Mat_DOS, kind_of_particle)
!         open(FN2, file=trim(adjustl(Input_elastic_file)))
!         write(*,'(a,a,a)') 'Calculated elastic mean free paths of a hole in ', trim(adjustl(Material_name)), ' are storred in the file'
!         write(*, '(a)') trim(adjustl(Input_elastic_file))
        select case (NumPar%kind_of_EMFP)
        case (2)    ! Read DSF elastic MFP
            write(temp_char1, '(f7.2, a)') Matter%temp, '_K'
            Input_elastic_file = trim(adjustl(Output_path))//'/OUTPUT_Hole_DSF_EMFPs_'//trim(adjustl(temp_char1))//'.dat'
            if (allocated(File_names%F)) File_names%F(5) = 'OUTPUT_Hole_DSF_EMFPs.dat'//trim(adjustl(temp_char1))//'.dat' ! save for later use
            FN2 = 2032

            file_exist = .false.
            write(*,'(a,a,a)') 'Calculated elastic mean free paths of a hole in ', trim(adjustl(Material_name)), ' are storred in the file:'
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
         case (1) ! Calculate or read CDF elastic MFP
            write(temp_char, '(f7.2, a)') Matter%temp, '_K'
            Input_elastic_file = trim(adjustl(Output_path))//'/OUTPUT_Hole_CDF_EMFPs_'//trim(adjustl(temp_char))//'.dat'
            if (allocated(File_names%F)) File_names%F(5) = 'OUTPUT_Hole_CDF_EMFPs_'//trim(adjustl(temp_char))//'.dat' ! save for later use
            FN2 = 204
            file_exist = .false.
            call All_elastic_scattering(Nelast, Target_atoms, CDF_Phonon, Matter, Elastic_MFP%Total, NumPar, Mat_DOS, kind_of_particle)
            open(FN2, file=trim(adjustl(Input_elastic_file)))
            write(*,'(a,a,a)') 'Calculated elastic mean free paths of a hole in ', trim(adjustl(Material_name)), ' are storred in the file'
            write(*, '(a)') trim(adjustl(Input_elastic_file))
         case (0) ! Mott cross-sections
            write(temp_char1, '(f7.2, a)') Matter%temp, '_K'
            Input_elastic_file = trim(adjustl(Output_path))//'/OUTPUT_Hole_Mott_EMFPs'//trim(adjustl(temp_char1))//'.dat'
            if (allocated(File_names%F)) File_names%F(4) = 'OUTPUT_Hole_Mott_EMFPs'//trim(adjustl(temp_char1))//'.dat' ! save for later use
            FN2 = 2043
            inquire(file=trim(adjustl(Input_elastic_file)),exist=file_exist)    ! check if input file excists

            ! Check, if file with MFP was created with paramters in actual CDF file, find out when this file was last modified:
            if (file_exist) then
                call get_file_stat(trim(adjustl(Input_elastic_file)), Last_modification_time=IMFP_last_modified) ! above
                !print*, 'IMFP file last modified on:', IMFP_last_modified
                if (IMFP_last_modified < NumPar%Last_mod_time_CDF) NumPar%redo_EMFP = .true. ! Material parameters changed, recalculate EMFPs
            endif

            if (file_exist .and. .not.NumPar%redo_EMFP) then    ! read from the file:
                write(*,'(a,a,a)') 'Mott EMFPs of a hole in ', trim(adjustl(Material_name)), ' are already in the file:'
                write(*, '(a)') trim(adjustl(Input_elastic_file))
                write(*, '(a)') ' '
                open(FN2, file=trim(adjustl(Input_elastic_file)), ACTION='READ')
            else    ! create and write to the file:
                call All_elastic_scattering(Nelast, Target_atoms, CDF_Phonon, Matter, Elastic_MFP%Total, NumPar, Mat_DOS, kind_of_particle)
                open(FN2, file=trim(adjustl(Input_elastic_file)))
                write(*,'(a,a,a)') 'Elastic mean free paths of an hole calculated using Mott formulae is in ', trim(adjustl(Material_name)), ' are storred in the file'
                write(*, '(a)') trim(adjustl(Input_elastic_file))
                write(*, '(a)') ' '
            endif
         case default ! disabled
            write(temp_char1, '(f7.2, a)') Matter%temp, '_K'
            Input_elastic_file = trim(adjustl(Output_path))//'/OUTPUT_Hole_No_elas_EMFPs.dat'
            if (allocated(File_names%F)) File_names%F(4) = 'OUTPUT_Hole_No_elas_EMFPs.dat' ! save for later use
            FN2 = 2031
            open(FN2, file=trim(adjustl(Input_elastic_file)))
            write(*,'(a)') 'Valence hole kinetics will be traces without elastic scatterings on target atoms'
            file_exist = .false.
            call All_elastic_scattering(Nelast, Target_atoms, CDF_Phonon, Matter, Elastic_MFP%Total, NumPar, Mat_DOS, kind_of_particle)
            Elastic_MFP%Total%L(:) = 1.0d30
         endselect
    endif kind_of_part2
    
    ! Now write the elastic output into the file:
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
            L_tot = 1.0d0/(1.0d0/ElasMFP + 1.0d0/InelMFP)               ! total MFP [A]
            Vel = sqrt(2*Mat_DOS%E(i)*g_e/(Mass*g_me))*1.0d10       ! Velosity [A/s]
            write(333, '(e,e,e,e,e)') Mat_DOS%E(i), Mass, Mat_DOS%DOS(i), 1.0d0/3.0d0*L_tot*vel/1.0d16, L_tot
        enddo
        close(333)
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

    Emin = 1.0d0    ! [eV] we start with this minimum
    Ord = -2 ! start with 1
    Va = int(Emin)
    dE(1) = 10.0d0**Ord
    do i = 1, Nelast-1    
        if (dE(i) .LE. 1.0d0) then
            if (Va .GE. 10) then
                Va = Va - 9
                Ord = Ord + 1
            endif
        else
            if (Va .GE. 100) then
                Va = Va - 90
                Ord = Ord + 1
            endif
        endif
        dE(i+1) = dE(i) + 10.0d0**Ord
        Va = Va + 1
    enddo   ! while (dE .LT. Emax)

if (.not.present(dont_do)) then ! only do it when we have the CDF
!$omp parallel &
!$omp private (i, Ele, EMFP, dEdx)
!$omp do schedule(dynamic)
    do i = 1, Nelast
        Ele = dE(i)
        call Elastic_cross_section(Ele, CDF_Phonon, Target_atoms, Matter, EMFP, dEdx, NumPar, Mat_DOS, kind_of_particle) ! from module Cross_sections
        Elastic_MFP%E(i) = Ele      ! [eV] energy
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
    type(Atom), dimension(:), intent(in), target :: Target_atoms  ! all data for target atoms
    type(All_MFP), dimension(:), allocatable, intent(inout) :: Total_el_MFPs   ! electron mean free paths for all shells
    type(Density_of_states), intent(in) :: Mat_DOS
    type(Solid), intent(in) :: Matter ! properties of material
    type(Flag), intent(in) :: NumPar
    character(8), intent(in) :: kind_of_particle

    real(8) IMFP_calc, dEdx, Ele, dE(N), Emin
    integer i, j, k, Nshl, Nat, Va, Ord
    integer Num_th, my_id, OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
    
    Nat = size(Target_atoms)    ! number of atoms
    
    Emin = 1.0d0    ! [eV] we start with this minimum
    Ord = 0 ! start with 1
    Va = int(Emin)
    dE(1) = Emin
    do i = 1, N-1
        if (Va .GE. 100) then
            Va = Va - 90
            Ord = Ord + 1
        endif
        dE(i+1) = dE(i) + 10.0d0**Ord
        Va = Va + 1
    enddo

!$omp parallel &
!$omp private (i, j, k, Ele, IMFP_calc, dEdx)
!$omp do schedule(dynamic)
    do i = 1, N
        Ele = dE(i)
        !my_id = 1 + OMP_GET_THREAD_NUM() ! identify which thread it is
        do j = 1, Nat  ! for all atoms:
          Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
          do k = 1, Nshl  ! for all shells of each atom:
             call TotIMFP(Ele, Target_atoms, j, k, IMFP_calc, dEdx, Matter, Mat_DOS, NumPar, kind_of_particle) ! from module "Cross_sections"
             Total_el_MFPs(j)%ELMFP(k)%E(i) = Ele
             Total_el_MFPs(j)%ELMFP(k)%L(i) = IMFP_calc
             Total_el_MFPs(j)%ELMFP(k)%dEdx(i) = dEdx
          enddo ! k = 1, size(Target_atoms(j)%Ip)  ! for all shells of each atom:
        enddo ! j = 1,size(Target_atoms)  ! for all atoms:
        call progress(' Progress of calculation: ', i, N)
    enddo
!$omp end do
!$omp end parallel
end subroutine All_shells_Electron_MFP


subroutine All_shells_Photon_MFP(N, Target_atoms, Total_el_MFPs, Matter, NumPar) ! calculate all IMFPs
    integer, intent(in) :: N  ! number of energy points
    type(Atom), dimension(:), intent(in), target :: Target_atoms  ! all data for target atoms
    type(All_MFP), dimension(:), allocatable, intent(inout) :: Total_el_MFPs   ! electron mean free paths for all shells
    type(Solid), intent(in) :: Matter ! properties of material
    type(Flag), intent(in) :: NumPar

    real(8) IMFP_calc, dEdx, Ele, dE(N), Emin
    integer i, j, k, Nshl, Nat, Va, Ord
    integer Num_th, my_id, OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
    
    Nat = size(Target_atoms)    ! number of atoms
    
    Emin = 1.0d0    ! [eV] we start with this minimum
    Ord = 0 ! start with 1
    Va = int(Emin)
    dE(1) = Emin
    do i = 1, N-1
        if (Va .GE. 100) then
            Va = Va - 90
            Ord = Ord + 1
        endif
        dE(i+1) = dE(i) + 10.0d0**Ord
        Va = Va + 1
    enddo


!$omp parallel &
!$omp private (i, j, k, Ele, IMFP_calc, dEdx)
!$omp do schedule(dynamic)
    do i = 1, N
        Ele = dE(i)
        !my_id = 1 + OMP_GET_THREAD_NUM() ! identify which thread it is
        do j = 1, Nat  ! for all atoms:
          Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
          do k = 1, Nshl  ! for all shells of each atom:
             call Tot_Phot_IMFP(Ele, Target_atoms, j, k, IMFP_calc, dEdx, Matter, NumPar) ! from module "Cross_sections"
             Total_el_MFPs(j)%ELMFP(k)%E(i) = Ele
             Total_el_MFPs(j)%ELMFP(k)%L(i) = IMFP_calc
             Total_el_MFPs(j)%ELMFP(k)%dEdx(i) = dEdx
          enddo ! k = 1, size(Target_atoms(j)%Ip)  ! for all shells of each atom:
        enddo ! j = 1,size(Target_atoms)  ! for all atoms:
        call progress(' Progress of calculation: ', i, N)
    enddo
!$omp end do
!$omp end parallel
end subroutine All_shells_Photon_MFP


subroutine Analytical_ion_dEdx(Output_path_SHI, Material_name, Target_atoms, SHI, SHI_MFP, Error_message, read_well, NumPar, Matter, Mat_DOS, File_names)
    character(100), intent(in) :: Output_path_SHI   ! path to the folder where the file is/will be storred
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

    real(8), dimension(:), allocatable :: dEdx_tot
    real(8) SHI_E, Emin, Emax, dE
    integer N, Ord, Va, IMFP_last_modified
    integer i, j, k, Nat, Nshl, FN, FN2, FN3
    character(100) Input_files, Input_files11, Input_files2, Input_files3, Path_name, command, charge_name, charge_kind
    logical file_exist, file_exist2
    read_well = .true.  ! so far so good
    Nat = size(Target_atoms)    ! how many atoms
    if (.not. allocated(SHI_MFP)) allocate(SHI_MFP(size(Target_atoms))) ! that's how many atoms
    ! How many energy points will be here:
    Emin = real(CEILING(((SHI%Mass*g_Mp + g_me)*(SHI%Mass*g_Mp + g_me)/(SHI%Mass*g_Mp*g_me)*Target_atoms(1)%Ip(size(Target_atoms(1)%Ip))/4.0d0)))  ! [eV]
    Emax = 175.6d6/2.0d0*SHI%Mass ! [eV]  maximal energy that still has no relativism
    dE = Emin
    call Find_order_of_magn(Emin, i, j)    ! order of magnitude
    dE = j*10**i  ! start with this value
    Va = j
    Ord = i
    N = 0
    do while (dE .LT. Emax)
        N = N + 1
        if (Va .GE. 10) then 
            Va = Va - 9
            Ord = Ord + 1
        endif
        dE = dE + 10.0d0**Ord
        Va = Va + 1
    enddo   ! while (dE .LT. Emax)
    N = N + 2
    
    do j = 1, Nat   ! declair variables if they are not yet
        Nshl = size(Target_atoms(j)%Ip)    ! how mamy shells
        if (.not. allocated(SHI_MFP(j)%ELMFP)) allocate(SHI_MFP(j)%ELMFP(Nshl)) ! how many shells
        do k = 1, Nshl
            if (.not. allocated(SHI_MFP(j)%ELMFP(k)%E)) then    ! energy array
                allocate(SHI_MFP(j)%ELMFP(k)%E(N))
                SHI_MFP(j)%ELMFP(k)%E = 0.0d0
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
    inquire(DIRECTORY=trim(adjustl(Path_name)),exist=file_exist)    ! check if input file excists
    if (.not. file_exist) then  ! create the directory
        command='mkdir '//trim(adjustl(Path_name)) ! to create a folder use this command
        CALL system(command)  ! create the folder
    endif
    
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
    
    Input_files = trim(adjustl(Path_name))//'/OUTPUT_'//trim(adjustl(SHI%Name))//trim(adjustl(charge_name))//trim(adjustl(charge_kind))//'_IMFP.dat'
    Input_files2 = trim(adjustl(Path_name))//'/OUTPUT_'//trim(adjustl(SHI%Name))//trim(adjustl(charge_name))//trim(adjustl(charge_kind))//'_dEdx.dat'
    Input_files11 = trim(adjustl(Path_name))//'/OUTPUT_'//trim(adjustl(SHI%Name))//trim(adjustl(charge_name))//trim(adjustl(charge_kind))//'_effective_charges.dat'
    Input_files3 = trim(adjustl(Path_name))//'/OUTPUT_'//trim(adjustl(SHI%Name))//trim(adjustl(charge_name))//trim(adjustl(charge_kind))//'_Range.dat'
    if (allocated(File_names%F)) File_names%F(6) = 'OUTPUT_'//trim(adjustl(SHI%Name))//trim(adjustl(charge_name))//trim(adjustl(charge_kind))
    inquire(file=trim(adjustl(Input_files)),exist=file_exist)    ! check if input file excists

    ! Check, if file with MFP was created with paramters in actual CDF file, find out when this file was last modified:
    if (file_exist) then
        call get_file_stat(trim(adjustl(Input_files)), Last_modification_time=IMFP_last_modified) ! above
        !print*, 'IMFP file last modified on:', IMFP_last_modified
        if (IMFP_last_modified < NumPar%Last_mod_time_CDF) NumPar%redo_IMFP = .true. ! Material parameters changed, recalculate IMFPs
    endif

    if (file_exist .and. .not.NumPar%redo_IMFP) then
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
            write(*,'(a,a,a,a,a)') 'IMFP and dEdx of ', SHI%Name ,' in ', trim(adjustl(Material_name)), ' will be storred in the files:'
            write(*, '(a,a,a)') trim(adjustl(Input_files)), ' and ', trim(adjustl(Input_files2))
            write(*, '(a)') ' ' 
            call Analytical_SHI_dEdx(Input_files, Input_files2, Input_files11, N, Emin, Emax, SHI, SHI_MFP, Target_atoms, dEdx_tot, Matter, Mat_DOS, NumPar) ! from module Analytical_IMPS / openmp parallelization
            SHI%E = SHI_E ! restore the original value
        endif
    else
        SHI_E = SHI%E   ! just save it for future
        write(*,'(a,a,a,a,a)') 'IMFP and dEdx of ', SHI%Name ,' in ', trim(adjustl(Material_name)), ' will be storred in the files:'
        write(*, '(a,a,a)') trim(adjustl(Input_files)), ' and ', trim(adjustl(Input_files2))
        write(*, '(a)') ' ' 
        call Analytical_SHI_dEdx(Input_files, Input_files2, Input_files11, N, Emin, Emax, SHI, SHI_MFP, Target_atoms, dEdx_tot, Matter, Mat_DOS, NumPar) ! from module Analytical_IMPS / openmp parallelization
        SHI%E = SHI_E ! restore the original value
    endif
    inquire(file=trim(adjustl(Input_files3)),exist=file_exist2)    ! check if file with Ranges excists
    if (.not.file_exist2 .or. NumPar%redo_IMFP) then  ! if not, create it
        write(*,'(a,a,a,a,a)') 'Ranges of ', SHI%Name ,' in ', trim(adjustl(Material_name)), ' will be storred in the file:'
        write(*, '(a)') trim(adjustl(Input_files3))
        write(*, '(a)') ' '
        call Get_ion_range(Input_files3,N,SHI_MFP,Target_atoms,dEdx_tot) ! calculate ion range out of its energy-loss function
    endif
    NumPar%redo_IMFP = .false. ! default it for the next kind of particle
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
    write(FN3,'(A)')    'Energy dEdx    Range'
    write(FN3,'(A)')    '[eV] [eV/A]    [A]'
    do i = 1,N  ! save into file
        write(FN3,'(es,es,es)') SHI_MFP(1)%ELMFP(1)%E(i), dEdx_tot(i), SHI_range(i)
    enddo
    deallocate(SHI_range)
    deallocate(dEdx_tot)
    close(FN3)
end subroutine Get_ion_range


! This version is with parallelization via openmp:
subroutine Analytical_SHI_dEdx(Input_files, Input_files2, Input_files11, N, Emin, Emax, SHI, SHI_MFP, Target_atoms, dEdx_tot, Matter, Mat_DOS, NumPar)   ! calculates dEdx for range of SHI energies
   character(100), intent(in) :: Input_files ! path to the folder where the file is/will be storred
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
   
   integer Nat, Nshl, j, i, k, Va, Ord, FN, FN2, FN3
   integer Num_th, my_id, OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
   real(8), dimension(N, size(Target_atoms), size(Target_atoms(1)%Ip)) :: SHI_IMFP, SHI_dEdx
   real(8), dimension(N) :: E, Zeff  ! SHI inverse mean free path [1/A], and dEdx [eV/A], Zeff
   real(8), dimension(size(Target_atoms), size(Target_atoms(1)%Ip)) :: IMFP, dEdx
   real(8) IMFP_calc, dEdx_calc
   real(8), pointer :: Z
   type(Ion) :: SHI_1
   real(8), dimension(N) :: dE  ! energy grid [eV]

   Nat = size(Target_atoms) ! how many atoms
   
   dE = 0.0d0
   call Find_order_of_magn(Emin, i, j)    ! order of magnitude
   dE(1) = j*10**i  ! start with this value
   Va = j
   Ord = i
   do k = 1, N-1
       if (Va .GE. 10) then 
           Va = Va - 9
           Ord = Ord + 1
       endif
       dE(k+1) = dE(k) + 10.0d0**Ord
       Va = Va + 1
   enddo   ! while (dE .LT. Emax)   
   
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

end module Analytical_IMFPs

!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains subroutines that help to sort out output distributions

MODULE Sorting_output_data
  use Universal_Constants   ! let it use universal constants
  use Objects   ! since it uses derived types, it must know about them from module 'Objects'
  use Reading_files_and_parameters
implicit none

contains


subroutine Save_output(Output_path, ctim, NMC, Num_th, Tim, dt, Material_name, Matter, Target_atoms, &
                       SHI, Out_R, Out_tot_Ne, Out_tot_Nphot, Out_tot_E, Out_E_e, Out_E_phot, Out_nphot, Out_Ephot, Out_Ee_vs_E, Out_E_at, Out_E_h, Out_Eat_dens, & 
                       Out_Distr, Out_Elat, Out_theta, Out_field_all, Out_Ne_Em, Out_E_Em, Out_Ee_vs_E_Em, NumPar, &
                       Out_E_field, Out_diff_coeff)
    character(100), intent(in) :: Output_path, Material_name
    integer, dimension(:), intent(in) :: ctim
    integer, intent(in) :: NMC, Num_th  ! number of MC iterations, and number of threads used for openmp
    real(8), intent(in) :: Tim, dt ! [fs]
    type(Solid), intent(in) :: Matter   ! all material parameters
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! target atoms as objects
    type(Ion), intent(in) :: SHI   ! declare SHI as an object with atributes "Ion"
    real(8), dimension(:), allocatable, intent(in) :: Out_R
    real(8), dimension(:), allocatable, intent(in) :: Out_tot_Ne
    real(8), dimension(:), allocatable, intent(in) :: Out_tot_Nphot
    real(8), dimension(:), allocatable, intent(in) :: Out_tot_E
    real(8), dimension(:), allocatable, intent(in) :: Out_E_e
    real(8), dimension(:), allocatable, intent(in) :: Out_E_phot
    real(8), dimension(:,:), intent(in) :: Out_nphot 
    real(8), dimension(:,:), intent(in) :: Out_Ephot
    real(8), dimension(:,:), allocatable, intent(in) :: Out_Ee_vs_E
    real(8), dimension(:), allocatable, intent(in) :: Out_E_at
    real(8), dimension(:,:,:), allocatable, intent(in) :: Out_E_h
    real(8), dimension(:,:), allocatable, intent(in) :: Out_Eat_dens
    real(8), dimension(:,:), allocatable, intent(in) :: Out_Elat
    real(8), dimension(:,:), allocatable, intent(in) :: Out_theta
    real(8), dimension(:,:), allocatable, intent(in) :: Out_field_all
    
    real(8), dimension(:), allocatable, intent(in) :: Out_diff_coeff
    
    real(8), dimension(:), intent(in) :: Out_E_field
    
    real(8), dimension(:,:), intent(in) :: Out_Ee_vs_E_Em
    real(8), dimension(:), intent(in) :: Out_Ne_Em
    real(8), dimension(:), intent(in) :: Out_E_Em
    
    type(Cylinder_distr), intent(in) :: Out_Distr   ! OUTPUT radial distributions
    type(Flag), intent(in) :: NumPar
        
    real(8) t, as1, tim_glob
    integer c1(8), i, j,k,l,N, Nat, N_R, FN, FN1, FN2, FN3, FN31, FN4, Lowest_Ip_At, Lowest_Ip_Shl !, NOTP
    character(200) command, charge_name, charge_kind, File_name, File_name1, File_name2, File_name3, File_name4, C_time
    character(30) ch1, ch2, ch3
    character(LEN=25) FMT
    logical file_exist, file_opened
    real(8), dimension(:), allocatable :: time_grid ! grid_points in time
    
    Nat = size(Target_atoms)    ! number of atoms
    N = size(Out_tot_Ne)        ! number of time-steps
    N_R = size(Out_R)           ! number of grid-points for radius
    
    select case (NumPar%dt_flag)   ! what kind of time-grid to use:
        case (:0)   ! linear time-grid
            allocate(time_grid(CEILING(Tim/dt)+1))
            time_grid = 0.0d0
            time_grid(1) = dt
            do i = 1, size(time_grid)-1
                time_grid(i+1) = min(time_grid(i) + dt,Tim)
            enddo
        case (1:)   ! logarithmic time-grid
            i = 0
            tim_glob = 0.01d0
            do while (tim_glob .LE. Tim)
                i = i + 1
                tim_glob = tim_glob*dt	! [fs]
            enddo
            allocate(time_grid(i+1))
            time_grid(1) = 0.01d0   ! [fs] first point
            do i = 1, size(time_grid)-1
                time_grid(i+1) = time_grid(i)*dt	! [fs]
            enddo
            time_grid(size(time_grid)) = min(time_grid(size(time_grid)),Tim)
    endselect
    
    SELECT CASE (SHI%Kind_Zeff) ! 0=Barkas; 1=Bohr; 2=Nikolaev-Dmitriev; 3=Schiwietz-Grande;
    CASE (1)
        charge_name = 'Bohr'
    CASE (2)
        charge_name = 'Nikolaev-Dmitriev'
    CASE (3)
        charge_name = 'Schiwietz-Grande'
    CASE (4)
        charge_name = 'fixed'
    CASE DEFAULT
        charge_name = 'Barkas'
    END SELECT
    
    SELECT CASE (SHI%Kind_ion) ! 0=Point-like charge; 1=Brandt-Kitagawa
    CASE (1)
       charge_kind = 'Brandt-Kitagawa'
    CASE DEFAULT
       charge_kind = 'point-like charge'
    END SELECT
    
    
    write(ch1,'(f8.2)') SHI%E/1d6
    write(ch2,'(f6.2)') Tim
    write(File_name,'(a,a,a,a,a,a,a,a)') trim(adjustl(Output_path)),trim(adjustl(NumPar%path_sep)),trim(adjustl(SHI%Name)),'_E_',trim(adjustl(ch1)),'_MeV_',trim(adjustl(ch2)),'_fs'
    File_name2 = File_name
    i = 0
    inquire(DIRECTORY=trim(adjustl(File_name2)),exist=file_exist)    ! check if input file excists
    do while (file_exist)
        i = i + 1
        write(ch1,'(i6)') i
        write(File_name2,'(a,a,a)') trim(adjustl(File_name)), '_', trim(adjustl(ch1))
        inquire(DIRECTORY=trim(adjustl(File_name2)),exist=file_exist)    ! check if input file excists
    enddo
    command='mkdir '//trim(adjustl(File_name2)) ! to create a folder use this command
    CALL system(command)  ! create the folder
    write(*,'(a)') '--------------------------------'
    write(*,'(a,a)') 'The outputs with MC results are storred in the folder:', trim(adjustl(File_name2))
    
    !========================================================
    !Parameters of this calculation:
    FN1 = 299
    File_name = trim(adjustl(File_name2))//'/!Parameters.txt'
    open(unit = FN1, FILE = trim(adjustl(File_name)))
    write(FN1, '(a,a,a)') 'Material: ', trim(adjustl(Material_name)), ' ('//trim(adjustl(Matter%Target_name))//')'
    write(FN1, '(a,f10.2,a,e16.7,a)') 'Material density: ', Matter%Dens, ' [g/cm^3] or ', Matter%At_dens, ' [1/cm^3]'
    if (Matter%hole_mass .GT. 0) then
        write(FN1, '(a,e16.7,a)') 'Effective mass of valence holes: ', Matter%hole_mass, ' [me]'
    else
        write(FN1, '(a)') 'Effective mass of valence holes is calculated from DOS'
    endif    
    if (NumPar%kind_of_EMFP .EQ. 2) then
        write(FN1, '(a)') 'DSF cross-sections were used for electron-atom elastic collisions'
    else if (NumPar%kind_of_EMFP .EQ. 1) then
        write(FN1, '(a)') 'Cross-sections calculated from CDF phonon peaks were used for electron-atom elastic collisions'
    else if (NumPar%kind_of_EMFP .EQ. 0) then
        write(FN1, '(a)') 'Mott atomic cross-sections were used for electron-atom elastic collisions'
    else
        write(FN1, '(a)') 'Elastic scatterings are disabled'
    endif
    
    if (matter%El_eff_mass .EQ. 0) then
        write(FN1, '(a)') 'Effective mass of valence electrons is calculated from DOS.'
    else
        write(FN1, '(a,f6.2,a)') 'Effective mass of valence electrons: ', Matter%El_eff_mass, ' [me]'
    endif
    
    if(NumPar%Plasmon_Emax) then
        write(FN1, '(a)') 'Plasmon maximal energy is used as upper integration limit during cross-section calculation'
    endif
    
    if (maxval(target_atoms(1)%KOCS(:)) .EQ. 1) then ! to make better later...
        if (NumPar%kind_of_DR .EQ. 3) then
            write(FN1, '(a)') 'Electron inelastic mean free paths were calculated using dispersion relation from Ritchie and Howie (Phil. Mag. 36 (1977) 463).'
        else if (NumPar%kind_of_DR .EQ. 2) then
            write(FN1, '(a)') 'Electron inelastic mean free paths were calculated using plasmon-pole dispersion relation of scattering centers.'
        else if (NumPar%kind_of_DR .EQ. 4) then
            write(FN1, '(a)') 'Electron inelastic mean free paths were calculated using Delta-function CDF.'
        else
            write(FN1, '(a)') 'Electron inelastic mean free paths were calculated with free electron dispersion relation of scattering centers.'
        endif
    else ! BEB cross sections = 2
        write(FN1, '(a)') 'Electron inelastic mean free paths were calculated using BEB atomic cross sections.'
    endif
        
    if (NumPar%field_include .GT. 0.9) then
        write(FN1, '(a)') 'Produced electric fields are included in the modelling'
    else
        write(FN1, '(a)') 'Produced electric fields are NOT included in the modelling'
    endif
    if (Matter%work_function .GT. 0) then
        write(FN1, '(a, f10.2, a)') 'Electron emission is calculated with parameters: Work function = ', Matter%work_function, ' [eV]'
        write(FN1, '(a, f10.2, a, a, f10.2, a)') 'Potential barrier length =  ', Matter%bar_length, ' [A]', 'Barrier height = ', Matter%bar_height, ' [eV]'
    else
        write(FN1, '(a)') 'Electron emission is NOT included in the modelling'
    endif
    if (NumPar%include_photons) then
        write(FN1, '(a)') 'Radiative decays of deep-shell holes are included'
    else
        write(FN1, '(a)') 'Radiative decays of deep-shell holes are NOT included, only Auger decays are permitted'
    endif
    write(FN1, '(a,a,a)') 'Ion: ', trim(adjustl(SHI%Name)), ' ('//trim(adjustl(SHI%Full_Name))//')'
    write(FN1, '(a,f10.2,a)') 'Ion energy: ', SHI%E/1d6, ' [MeV]'
    
    if (SHI%Kind_Zeff .EQ. 4) then
        write(FN1, '(a,f6.3,a)') 'Ion fixed charge: ', SHI%fixed_Zeff, ' [electron charge] is used'
        write(FN1, '(a,a,a)') 'with ', trim(adjustl(charge_kind)), ' model for SHI.'
    else    
        write(FN1, '(a,f6.3,a)') 'Ion equilibrium charge: ', SHI%Zeff, ' [electron charge]'
        write(FN1, '(a,a,a,a,a)') 'calculated using ', trim(adjustl(charge_name)), ' formula with ', trim(adjustl(charge_kind)), ' model for SHI.'
    endif
    write(FN1, '(a,f9.3,a)') 'Ion mass: ', SHI%Mass, ' [proton mass]'
    write(FN1, '(a,f9.2,a)') 'MC calculated energy loss: ', Out_tot_E(N)/Matter%Layer, ' [eV/A]'
    select case (NumPar%dt_flag)
        case (:0)
            write(FN1, '(a,f7.2,a,f7.3,a)') 'Calculated up to: ', Tim, ' [fs] with time step of ', dt, ' [fs]'
        case (1:)
            write(FN1, '(a,f7.2,a,f7.3,a)') 'Calculated from 0.01 [fs] up to: ', Tim, ' increasing by ', dt, ' in logarithm scale'
       endselect
    write(FN1, '(a,f9.2,a)') 'Electron cut-off energy: ', Matter%cut_off, ' [eV]'
    write(FN1, '(a,f9.2,a)') 'Thickness of the analyzed layer: ', Matter%Layer, ' [A]'
    write(FN1, '(a,f9.2,a)') 'Temperature of the target: ', Matter%temp, ' [K]'
    write(FN1, '(a,i5)') 'Number of MC iterations: ', NMC
    write(FN1, '(a,i5)') 'Number of threads used for openmp:', Num_th
    call date_and_time(values=c1)	    ! For calculation of the time of execution of the program
    as1=REAL(24*60*60*(c1(3)-ctim(3))+3600*(c1(5)-ctim(5))+60*(c1(6)-ctim(6))+(c1(7)-ctim(7))+(c1(8)-ctim(8))*0.001)	! sec
    call parse_time(as1,C_time) ! module "Sorting_output_data.f90"
    write(FN1, '(a,a)') 'Duration of calculation: ', trim(adjustl(C_time))
    write(FN1, '(a)') '-------------------------------------------------------------'
    write(FN1, '(a)') 'The following atomic parameters of the target have been used:'
    do j = 1, size(Target_atoms)  ! read for each element its shells data:
      write(FN1, '(a)') trim(adjustl(Target_atoms(j)%Name))//'-atom:'
      write(FN1, '(a)') 'Shell  Quantum_n  Ne    Ip[eV]  Ekin[eV]    t(Auger)[fs]    t(Radiative)[fs]'
      do k = 1, Target_atoms(j)%N_shl ! read all the data for each shell:
         if ((j .EQ. 1) .AND. (k .EQ. Target_atoms(j)%N_shl)) then
            write(FN1,'(a,a,f8.2,f9.2,f9.2,es12.2,es12.2)') Target_atoms(j)%Shell_name(k), '    ', Target_atoms(j)%Nel(k), Target_atoms(j)%Ip(k), Target_atoms(j)%Ek(k), Target_atoms(j)%Auger(k), Target_atoms(j)%Radiat(k)
         else
            write(FN1,'(a,i3,f8.2,f9.2,f9.2,es12.2,es12.2)') Target_atoms(j)%Shell_name(k), Target_atoms(j)%PQN(k), Target_atoms(j)%Nel(k), Target_atoms(j)%Ip(k), Target_atoms(j)%Ek(k), Target_atoms(j)%Auger(k), Target_atoms(j)%Radiat(k)
         endif
      enddo
    enddo
    inquire(unit=FN1,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN1)             ! and if it is, close it


    ! Electron field vs radius:
    if (NumPar%field_include .GT. 0.9) then
        FN3 = 312
        File_name = trim(adjustl(File_name2))//'/Electric_fields_produced[V_m-1].txt'
        open(unit = FN3, FILE = trim(adjustl(File_name)))
        write(FN3, '(a)', advance='no') 'Radius[A] '
        t = 0.0d0
        do i = 1, N     ! timesteps
            t = time_grid(i)
            write(FN3, '(f6.2,a)', advance='no') t, '[fs]   '
        enddo
        write(FN3, '(a)') ' '
        do i = 1, N_R   ! radii
            write(FN3, '(f9.1)', advance='no') Out_R(i)
            t = 0.0d0
            do k = 1, N ! time-steps
                t = time_grid(k)
                write(FN3, '(e)', advance='no') Out_field_all(k,i)
            enddo
            write(FN3, '(a)') ' '
        enddo
        inquire(unit=FN3,opened=file_opened)    ! check if this file is opened
        if (file_opened) close(FN3)             ! and if it is, close it
    endif
    
    !Angular distribution of electrons:
    FN3 = 500
    File_name = trim(adjustl(File_name2))//'/Theta_distribution.txt'
    open(unit = FN3, FILE = trim(adjustl(File_name)))
    write(FN3, '(a)', advance='no') 'Angle[deg] '
    t = 0.0d0
    do i = 1, N     ! timesteps
        t = time_grid(i)
        write(FN3, '(f6.2,a)', advance='no') t, '[fs]   '
    enddo
    write(FN3, '(a)') ' '
    do i = 1, 180
        write(FN3, '(e)', advance='no') REAL(i)
        t = 0.0d0
        do k = 1, N ! time-steps
            t = time_grid(k)
            write(FN3, '(e)', advance='no') Out_theta(k,i)
        enddo
        write(FN3, '(a)') ' '
    enddo
    inquire(unit=FN3,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN3)             ! and if it is, close it


  !Emitted Electron distribution in energy space (vs E) vs time:
  if (Matter%work_function .GT. 0.0d0) then      
    FN3 = 501
    File_name = trim(adjustl(File_name2))//'/Emitted_electron_distribution_vs_E[1_eV].txt'
    open(unit = FN3, FILE = trim(adjustl(File_name)))
    write(FN3, '(a)', advance='no') 'Energy[eV] '
    t = 0.0d0
    do i = 1, N     ! timesteps
        t = time_grid(i)
        write(FN3, '(f6.2,a)', advance='no') t, '[fs]   '
    enddo
    write(FN3, '(a)') ' '
    do i = 1, N_R   ! radii
        write(FN3, '(f9.1)', advance='no') Out_R(i)/10.0
        t = 0.0d0
        do k = 1, N ! time-steps
            t = time_grid(k)
            write(FN3, '(e)', advance='no') Out_Ee_vs_E_Em(k,i)
        enddo
        write(FN3, '(a)') ' '
    enddo
    inquire(unit=FN3,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN3)             ! and if it is, close it      
  endif


    !Total number of electrons and total energy to check the energy conservation:
    FN = 300
    File_name = trim(adjustl(File_name2))//'/Total_numbers.txt'
    open(unit = FN, FILE = trim(adjustl(File_name)))
    t = 0.0d0
    if (NumPar%include_photons) then ! only if we include photons:
       write(FN, '(a)') 'Time[fs]    Ne    Ne_Emitted    Energy[eV]     Energy_Emitted[eV] N_photons'
       do i = 1, N
            t = time_grid(i)
            write(FN, '(e,e,e,e,e,e)') t, Out_tot_Ne(i), Out_Ne_Em(i), Out_tot_E(i), Out_E_Em(i), Out_tot_Nphot(i)
       enddo
    else
       write(FN, '(a)') 'Time[fs]    Ne    Ne_Emitted    Energy[eV]     Energy_Emitted[eV]'
       do i = 1, N
            t = time_grid(i)
            write(FN, '(e,e,e,e,e)') t, Out_tot_Ne(i), Out_Ne_Em(i), Out_tot_E(i), Out_E_Em(i)
       enddo
    endif
    inquire(unit=FN,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN)             ! and if it is, close it
    
    
    
    !Hole mean diffusion coefficient
    FN = 3001
    File_name = trim(adjustl(File_name2))//'/Hole_mean_diffusion_coefficient.txt'
    open(unit = FN, FILE = trim(adjustl(File_name)))
    t = 0.0d0
    if (Matter%hole_mass .LT. 1.0d3) then ! only if we calculate hole propagation
       write(FN, '(a)') 'Time[fs]    Diffusion_coeff[cm^2/s]'
       do i = 1, N
            t = time_grid(i)
            write(FN, '(e,e)') t, Out_diff_coeff(i)
       enddo
    endif
    inquire(unit=FN,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN)             ! and if it is, close it
    
    
    
    !Total energy of electrons, atoms, and all holes:
    FN2 = 301
    File_name = trim(adjustl(File_name2))//'/Total_energies.txt'
    open(unit = FN2, FILE = trim(adjustl(File_name)))
    if (NumPar%include_photons) then ! only if we include photons:
        write(FN2, '(a)', advance='no') 'Time[fs]    Electrons[eV]   Atoms[eV]   Field[eV]  Photons[eV] '
    else
        write(FN2, '(a)', advance='no') 'Time[fs]    Electrons[eV]   Atoms[eV]   Field[eV]  '
    endif
    do i = 1, Nat   ! all atoms
        do j = 1, size(Target_atoms(i)%Ip)
            if (Target_atoms(i)%Shl_num(j) .LT. 63) then    ! shell of an atom:
                write(FN2, '(a,a,a,a)', advance='no') trim(adjustl(Target_atoms(i)%Name)),'_',trim(adjustl(Target_atoms(i)%Shell_name(j))),'[eV]    '
            else    ! valence band:
                write(FN2, '(a,a)', advance='no') trim(adjustl(Target_atoms(i)%Shell_name(j))),'[eV]    '
            endif
        enddo
    enddo
    write(FN2, '(a)') ' '
    t = 0.0d0
    do k = 1, N
        t = time_grid(k)
        if (NumPar%include_photons) then ! only if we include photons:
            write(FN2, '(e,e,e,e,e)', advance='no') t, Out_E_e(k), Out_E_at(k), Out_E_field(k), Out_E_phot(k)
        else
            write(FN2, '(e,e,e,e)', advance='no') t, Out_E_e(k), Out_E_at(k), Out_E_field(k)
        endif
        do i = 1, Nat   ! all atoms
            do j = 1, size(Target_atoms(i)%Ip)
                write(FN2, '(e)', advance='no') Out_E_h(k,i,j)
            enddo
        enddo
        write(FN2, '(a)') ' '
    enddo
    inquire(unit=FN2,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN2)             ! and if it is, close it
    
    
    !Electron density vs R vs time:
    FN3 = 302
    File_name = trim(adjustl(File_name2))//'/Radial_electron_density[1_cm^-3].txt'
    open(unit = FN3, FILE = trim(adjustl(File_name)))
    write(FN3, '(a)', advance='no') 'Radius[A] '
    t = 0.0d0
    do i = 1, N     ! timesteps
        t = time_grid(i)
        write(FN3, '(f6.2,a)', advance='no') t, '[fs]   '
    enddo
    write(FN3, '(a)') ' '
    do i = 1, N_R   ! radii
        write(FN3, '(f9.1)', advance='no') Out_R(i)
        t = 0.0d0
        do k = 1, N ! time-steps
            t = time_grid(k)
            write(FN3, '(e)', advance='no') Out_Distr%ne(k,i)*1.0d24
        enddo
        write(FN3, '(a)') ' '
    enddo
    inquire(unit=FN3,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN3)             ! and if it is, close it
    
    !Electron energy density vs R vs time:
    FN3 = 303
    File_name = trim(adjustl(File_name2))//'/Radial_electron_energy[eV_A^-3].txt'
    open(unit = FN3, FILE = trim(adjustl(File_name)))
    write(FN3, '(a)', advance='no') 'Radius[A] '
    t = 0.0d0
    do i = 1, N     ! timesteps
        t = time_grid(i)
        write(FN3, '(f6.2,a)', advance='no') t, '[fs]   '
    enddo
    write(FN3, '(a)') ' '
    do i = 1, N_R   ! radii
        write(FN3, '(f9.1)', advance='no') Out_R(i)
        t = 0.0d0
        do k = 1, N ! time-steps
            t = time_grid(k)
            write(FN3, '(e)', advance='no') Out_Distr%Ee(k,i)
        enddo
        write(FN3, '(a)') ' '
    enddo
    inquire(unit=FN3,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN3)             ! and if it is, close it
    
    
    photons_here:if (NumPar%include_photons) then ! only if we include photons:
        !Photon density vs R vs time:
        FN3 = 306
        File_name = trim(adjustl(File_name2))//'/Radial_photon_density[1_cm^-3].txt'
        open(unit = FN3, FILE = trim(adjustl(File_name)))
        write(FN3, '(a)', advance='no') 'Radius[A] '
        t = 0.0d0
        do i = 1, N     ! timesteps
            t = time_grid(i)
            write(FN3, '(f6.2,a)', advance='no') t, '[fs]   '
        enddo
        write(FN3, '(a)') ' '
        do i = 1, N_R   ! radii
            write(FN3, '(f9.1)', advance='no') Out_R(i)
            t = 0.0d0
            do k = 1, N ! time-steps
                t = time_grid(k)
                write(FN3, '(e)', advance='no') Out_nphot(k,i)*1.0d24
            enddo
            write(FN3, '(a)') ' '
        enddo
        inquire(unit=FN3,opened=file_opened)    ! check if this file is opened
        if (file_opened) close(FN3)             ! and if it is, close it
        
        !Photon energy density vs R vs time:
        FN3 = 306
        File_name = trim(adjustl(File_name2))//'/Radial_photon_energy[eV_A^-3].txt'
        open(unit = FN3, FILE = trim(adjustl(File_name)))
        write(FN3, '(a)', advance='no') 'Radius[A] '
        t = 0.0d0
        do i = 1, N     ! timesteps
            t = time_grid(i)
            write(FN3, '(f6.2,a)', advance='no') t, '[fs]   '
        enddo
        write(FN3, '(a)') ' '
        do i = 1, N_R   ! radii
            write(FN3, '(f9.1)', advance='no') Out_R(i)
            t = 0.0d0
            do k = 1, N ! time-steps
                t = time_grid(k)
                write(FN3, '(e)', advance='no') Out_Ephot(k,i)
            enddo
            write(FN3, '(a)') ' '
        enddo
        inquire(unit=FN3,opened=file_opened)    ! check if this file is opened
        if (file_opened) close(FN3)             ! and if it is, close it
    endif photons_here
    
    
    
    
    !Electron distribution in energy space (vs E) vs time:
    FN3 = 303
    File_name = trim(adjustl(File_name2))//'/Electron_distribution_vs_E[1_eV].txt'
    open(unit = FN3, FILE = trim(adjustl(File_name)))
    write(FN3, '(a)', advance='no') 'Energy[eV] '
    t = 0.0d0
    do i = 1, N     ! timesteps
        t = time_grid(i)
        write(FN3, '(f6.2,a)', advance='no') t, '[fs]   '
    enddo
    write(FN3, '(a)') ' '
    do i = 1, N_R   ! radii
        write(FN3, '(f9.1)', advance='no') Out_R(i)
        t = 0.0d0
        do k = 1, N ! time-steps
            t = time_grid(k)
            write(FN3, '(e)', advance='no') Out_Ee_vs_E(k,i)
        enddo
        write(FN3, '(a)') ' '
    enddo
    inquire(unit=FN3,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN3)             ! and if it is, close it
    

    !Lattice energy density vs R vs time:
    FN3 = 303
    File_name = trim(adjustl(File_name2))//'/Radial_Lattice_energy[eV_A^-3].txt'
    open(unit = FN3, FILE = trim(adjustl(File_name)))
    write(FN3, '(a)', advance='no') 'Radius[A] '
    t = 0.0d0
    do i = 1, N     ! timesteps
        t = time_grid(i)
        write(FN3, '(f6.2,a)', advance='no') t, '[fs]   '
    enddo
    write(FN3, '(a)') ' '
    do i = 1, N_R   ! radii
        write(FN3, '(f9.1)', advance='no') Out_R(i)
        t = 0.0d0
        do k = 1, N ! time-steps
            t = time_grid(k)
            write(FN3, '(e)', advance='no') SUM(Out_Elat(1:k,i))
        enddo
        write(FN3, '(a)') ' '
    enddo
    inquire(unit=FN3,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN3)             ! and if it is, close it
    
    
    !Lattice energy density vs R vs time:
    call Find_VB_numbers(Target_atoms, Lowest_Ip_At, Lowest_Ip_Shl)
    FN3 = 304
    File_name = trim(adjustl(File_name2))//'/Radial_Track_energy[eV_A^-3].txt'
    open(unit = FN3, FILE = trim(adjustl(File_name)))
    write(FN3, '(a)', advance='no') 'Radius[A] '
    t = 0.0d0
    do i = 1, N     ! timesteps
        t = time_grid(i)
        write(FN3, '(f6.2,a)', advance='no') t, '[fs]   '
    enddo
    write(FN3, '(a)') ' '
    do i = 1, N_R   ! radii
        write(FN3, '(f9.1)', advance='no') Out_R(i)
        t = 0.0d0
        do k = 1, N ! time-steps
            t = time_grid(k)
            write(FN3, '(e)', advance='no') SUM(Out_Elat(1:k,i)) + Out_Distr%Atom(Lowest_Ip_At)%Shl(Lowest_Ip_Shl)%Eh(k,i) + Out_Distr%Atom(Lowest_Ip_At)%Shl(Lowest_Ip_Shl)%Ehkin(k,i)
        enddo
        write(FN3, '(a)') ' '
    enddo
    inquire(unit=FN3,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN3)             ! and if it is, close it
    
    
    !Holes densities vs R vs time:
    do j = 1, Nat   ! all atoms
        do l = 1, size(Target_atoms(j)%Ip)  ! all shells
            FN3 = 304
            if (Target_atoms(j)%Shl_num(l) .LT. 63) then    ! shell of an atom:
                write(File_name4, '(a,a,a)') trim(adjustl(Target_atoms(j)%Name)),'_',trim(adjustl(Target_atoms(j)%Shell_name(l)))
            else    ! valence band:
                write(File_name4, '(a)') trim(adjustl(Target_atoms(j)%Shell_name(l)))
            endif
            File_name = trim(adjustl(File_name2))//'/Radial_'//trim(adjustl(File_name4))//'_holes_density[1_cm^-3].txt'
            open(unit = FN3, FILE = trim(adjustl(File_name)))
            write(FN3, '(a)', advance='no') 'Radius[A] '
            t = 0.0d0
            do i = 1, N     ! timesteps
                t = time_grid(i)
                write(FN3, '(f6.2,a)', advance='no') t, '[fs]   '
            enddo
            write(FN3, '(a)') ' '
            do i = 1, N_R   ! radii
                write(FN3, '(f9.1)', advance='no') Out_R(i)
                t = 0.0d0
                do k = 1, N ! time-steps
                    t = time_grid(k)
                    write(FN3, '(e)', advance='no') Out_Distr%Atom(j)%Shl(l)%nh(k,i)*1.0d24
                enddo
                write(FN3, '(a)') ' '
            enddo
            inquire(unit=FN3,opened=file_opened)    ! check if this file is opened
            if (file_opened) close(FN3)             ! and if it is, close it
        enddo   ! l
    enddo   ! j
    
    
    !Holes energies density vs R vs time:
    do j = 1, Nat   ! all atoms
        do l = 1, size(Target_atoms(j)%Ip)  ! all shells
            if (Target_atoms(j)%Shl_num(l) .LT. 63) then    ! shell of an atom:
                FN3 = 305
                write(File_name4, '(a,a,a)') trim(adjustl(Target_atoms(j)%Name)),'_',trim(adjustl(Target_atoms(j)%Shell_name(l)))
                File_name = trim(adjustl(File_name2))//'/Radial_'//trim(adjustl(File_name4))//'_holes_energy[eV_A^-3].txt'
                open(unit = FN3, FILE = trim(adjustl(File_name)))
                write(FN3, '(a)', advance='no') 'Radius[A] '
                t = 0.0d0
                do i = 1, N     ! timesteps
                    t = time_grid(i)
                    write(FN3, '(f6.2,a)', advance='no') t, '[fs]   '
                enddo
                write(FN3, '(a)') ' '
                do i = 1, N_R   ! radii
                    write(FN3, '(f9.1)', advance='no') Out_R(i)
                    t = 0.0d0
                    do k = 1, N ! time-steps
                        t = time_grid(k)
                        write(FN3, '(e)', advance='no') Out_Distr%Atom(j)%Shl(l)%Eh(k,i)
                    enddo
                    write(FN3, '(a)') ' '
                enddo
            else    ! valence band:
                write(File_name4, '(a)') trim(adjustl(Target_atoms(j)%Shell_name(l)))
                FN3 = 305
                File_name = trim(adjustl(File_name2))//'/Radial_'//trim(adjustl(File_name4))//'_holes_pot_energy[eV_A^-3].txt'
                FN31 = 3051
                File_name1 = trim(adjustl(File_name2))//'/Radial_'//trim(adjustl(File_name4))//'_holes_kin_energy[eV_A^-3].txt'
                open(unit = FN3, FILE = trim(adjustl(File_name)))
                open(unit = FN31, FILE = trim(adjustl(File_name1)))
                write(FN3, '(a)', advance='no') 'Radius[A] '
                write(FN31, '(a)', advance='no') 'Radius[A] '
                t = 0.0d0
                do i = 1, N     ! timesteps
                    t = time_grid(i)
                    write(FN3, '(f6.2,a)', advance='no') t, '[fs]   '
                    write(FN31, '(f6.2,a)', advance='no') t, '[fs]   '
                enddo
                write(FN3, '(a)') ' '
                write(FN31, '(a)') ' '
                do i = 1, N_R   ! radii
                    write(FN3, '(f9.1)', advance='no') Out_R(i)
                    write(FN31, '(f9.1)', advance='no') Out_R(i)
                    t = 0.0d0
                    do k = 1, N ! time-steps
                        t = time_grid(k)
                        write(FN3, '(e)', advance='no') Out_Distr%Atom(j)%Shl(l)%Eh(k,i)
                        write(FN31, '(e)', advance='no') Out_Distr%Atom(j)%Shl(l)%Ehkin(k,i)
                    enddo
                    write(FN3, '(a)') ' '
                    write(FN31, '(a)') ' '
                enddo
                inquire(unit=FN31,opened=file_opened)    ! check if this file is opened
                if (file_opened) close(FN31)
            endif
            inquire(unit=FN3,opened=file_opened)    ! check if this file is opened
            if (file_opened) close(FN3)             ! and if it is, close it
          
        enddo   ! l
    enddo   ! j
end subroutine Save_output



subroutine Allocate_out_arrays(target_atoms, Out_Distr, Out_tot_Ne, Out_tot_Nphot, Out_tot_E, Out_E_e, Out_E_phot, Out_E_at, Out_E_h, &
                         Out_R, Out_V, Out_ne, Out_Ee, Out_nphot, Out_Ephot, Out_Ee_vs_E, Out_Elat, Out_nh, Out_Eh, Out_Ehkin, Out_Eat_dens, &
                         Out_theta, Out_theta1, Out_field_all, Out_Ne_Em, Out_E_Em, Out_Ee_vs_E_Em, Out_E_field, Out_diff_coeff)
    type(Atom), dimension(:), allocatable, intent(in) :: Target_atoms  ! target atoms as objects
    type(Cylinder_distr), intent(in) :: Out_Distr   ! OUTPUT radial distributions
    real(8), dimension(:), allocatable, intent(inout) :: Out_tot_Ne
    real(8), dimension(:), allocatable, intent(inout) :: Out_tot_Nphot
    real(8), dimension(:), allocatable, intent(inout) :: Out_tot_E
    real(8), dimension(:), allocatable, intent(inout) :: Out_V
    real(8), dimension(:), allocatable, intent(inout) :: Out_E_e
    real(8), dimension(:), allocatable, intent(inout) :: Out_E_phot
    real(8), dimension(:), allocatable, intent(inout) :: Out_E_at
    real(8), dimension(:,:,:), allocatable, intent(inout) :: Out_E_h
    real(8), dimension(:), allocatable, intent(inout) :: Out_R
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_ne
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Ee
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_nphot
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Ephot
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Ee_vs_E
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Elat
    real(8), dimension(:,:,:,:), allocatable, intent(inout) :: Out_nh
    real(8), dimension(:,:,:,:), allocatable, intent(inout) :: Out_Eh
    real(8), dimension(:,:,:,:), allocatable, intent(inout) :: Out_Ehkin
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Eat_dens  ! [eV/A^3] atom's energy energy   
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_theta
    real(8), dimension(:), allocatable, intent(inout) :: Out_theta1 
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_field_all
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Ee_vs_E_Em
    real(8), dimension(:), allocatable, intent(inout) :: Out_Ne_Em
    real(8), dimension(:), allocatable, intent(inout) :: Out_E_Em 
    real(8), dimension(:), allocatable, intent(inout) :: Out_E_field
    real(8), dimension(:), allocatable, intent(inout) :: Out_diff_coeff
    
    
    integer i, Nne1, Nne2, Nee1, Nee2, Nr
    
    Nne1 = size(Out_Distr%ne,1)
    Nne2 = size(Out_Distr%ne,2)
    Nr =  size(Out_Distr%R)
    Nee1 = size(Out_Distr%Ee,1)
    Nee2 = size(Out_Distr%Ee,2)
        
    ! Allocate and set temporary arrays:
    allocate(Out_theta(1+Nee1,180))
    allocate(Out_theta1(180))
    Out_theta = 0.0d0
    do i = 1,180
        Out_theta1 (i) = i
    enddo
    allocate(Out_tot_Ne(Nne1))  ! total number of electrons in time
    allocate(Out_tot_Nphot(Nne1))
    Out_tot_Ne = 0.0d0
    Out_tot_Nphot = 0.0d0
    
    allocate(Out_diff_coeff(Nne1))
    Out_diff_coeff = 0.0d0
    
    allocate(Out_tot_E(Nne1))   ! total energy in time
    Out_tot_E = 0.0d0

    if (.not. allocated(Out_V)) allocate(Out_V(Nr))
    allocate(Out_E_e(Nne1))
    Out_E_e = 0.0d0
    allocate(Out_E_phot(Nne1))
    Out_E_phot = 0.0d0
    allocate(Out_E_at(Nne1))
    Out_E_at = 0.0d0
    allocate(Out_E_h(Nne1, size(target_atoms), size(target_atoms(1)%Ip)))
    Out_E_h = 0.0d0
    allocate(Out_R(Nr))
    Out_R = Out_Distr%R
    allocate(Out_ne(Nne1, Nne2))
    Out_ne = 0.0d0
    allocate(Out_nphot(Nne1, Nne2))
    Out_nphot = 0.0d0
    allocate(Out_Ee(Nee1, Nee2))
    Out_Ee = 0.0d0
    allocate(Out_Ephot(Nee1, Nee2))
    Out_Ephot = 0.0d0
    allocate(Out_Ee_vs_E(Nee1, Nee2))
    Out_Ee_vs_E = 0.0d0
    allocate(Out_Elat(Nee1, Nee2))
    Out_Elat = 0.0d0
    allocate(Out_Eat_dens(Nee1, Nee2))
    Out_Eat_dens = 0.0d0
    allocate(Out_nh(Nne1, Nne2, size(target_atoms), size(target_atoms(1)%Ip)))
    Out_nh = 0.0d0
    allocate(Out_Eh(Nne1, Nne2, size(target_atoms), size(target_atoms(1)%Ip)))
    Out_Eh = 0.0d0
    allocate(Out_Ehkin(Nne1, Nne2, size(target_atoms), size(target_atoms(1)%Ip)))
    Out_Ehkin = 0.0d0
    allocate(Out_field_all(Nne1, Nee2))
    Out_field_all = 0.0d0    
    allocate(Out_E_Em(Nne1))
    Out_E_Em = 0.0d0
    allocate(Out_Ne_Em(Nne1))
    Out_Ne_Em = 0.0d0
    allocate(Out_Ee_vs_E_Em(Nee1,Nee2))
    Out_Ee_vs_E_Em = 0.0d0
    allocate(Out_E_field(Nne1))
    Out_E_field = 0.0d0
end subroutine Allocate_out_arrays



subroutine Allocate_out_arrays_old(target_atoms, Out_Distr, Out_tot_Ne, Out_tot_E, Out_E_e, Out_E_at, Out_E_h, &
                         Out_R, Out_V, Out_ne, Out_Ee, Out_Ee_vs_E, Out_Elat, Out_nh, Out_Eh, Out_Ehkin, Out_Eat_dens, &
                         Out_theta, Out_theta1, Out_field_all, Out_Ne_Em, Out_E_Em, Out_Ee_vs_E_Em, Out_E_field)
    type(Atom), dimension(:), allocatable, intent(in) :: Target_atoms  ! target atoms as objects
    type(Cylinder_distr), intent(in) :: Out_Distr   ! OUTPUT radial distributions
    real(8), dimension(:), allocatable, intent(inout) :: Out_tot_Ne
    real(8), dimension(:), allocatable, intent(inout) :: Out_tot_E
    real(8), dimension(:), allocatable, intent(inout) :: Out_V
    real(8), dimension(:), allocatable, intent(inout) :: Out_E_e
    real(8), dimension(:), allocatable, intent(inout) :: Out_E_at
    real(8), dimension(:,:,:), allocatable, intent(inout) :: Out_E_h
    real(8), dimension(:), allocatable, intent(inout) :: Out_R
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_ne
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Ee
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Ee_vs_E
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Elat
    real(8), dimension(:,:,:,:), allocatable, intent(inout) :: Out_nh
    real(8), dimension(:,:,:,:), allocatable, intent(inout) :: Out_Eh
    real(8), dimension(:,:,:,:), allocatable, intent(inout) :: Out_Ehkin
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Eat_dens  ! [eV/A^3] atom's energy energy   
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_theta
    real(8), dimension(:), allocatable, intent(inout) :: Out_theta1 
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_field_all
    
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Ee_vs_E_Em
    real(8), dimension(:), allocatable, intent(inout) :: Out_Ne_Em
    real(8), dimension(:), allocatable, intent(inout) :: Out_E_Em 
    real(8), dimension(:), allocatable, intent(inout) :: Out_E_field
    
    integer i, Nne1, Nne2, Nee1, Nee2, Nr
    
    Nne1 = size(Out_Distr%ne,1)
    Nne2 = size(Out_Distr%ne,2)
    Nr =  size(Out_Distr%R)
    Nee1 = size(Out_Distr%Ee,1)
    Nee2 = size(Out_Distr%Ee,2)
    
    ! Allocate and set temporary arrays:
    allocate(Out_theta(1+size(Out_Distr%Ee,1),180))
    allocate(Out_theta1(180))
    Out_theta = 0.0d0
    do i = 1,180
        Out_theta1 (i) = i
    enddo
    
    allocate(Out_E_Em(size(Out_Distr%ne,1)))
    Out_E_Em = 0.0d0
    allocate(Out_Ne_Em(size(Out_Distr%ne,1)))
    Out_Ne_Em = 0.0d0
    allocate(Out_Ee_vs_E_Em(size(Out_Distr%Ee,1), size(Out_Distr%Ee,2)))
    Out_Ee_vs_E_Em = 0.0d0
    
    allocate(Out_E_field(size(Out_Distr%ne,1)))
    Out_E_field = 0.0d0
    
    allocate(Out_tot_Ne(size(Out_Distr%ne,1)))  ! total number of electrons in time
    Out_tot_Ne = 0.0d0
    allocate(Out_tot_E(size(Out_Distr%ne,1)))   ! total energy in time
    Out_tot_E = 0.0d0
    if (.not. allocated(Out_V)) allocate(Out_V(size(Out_Distr%R)))
    allocate(Out_E_e(size(Out_Distr%ne,1)))
    Out_E_e = 0.0d0
    allocate(Out_E_at(size(Out_Distr%ne,1)))
    Out_E_at = 0.0d0
    allocate(Out_E_h(size(Out_Distr%ne,1), size(target_atoms), size(target_atoms(1)%Ip)))
    Out_E_h = 0.0d0
    allocate(Out_R(size(Out_Distr%R)))
    Out_R = Out_Distr%R
    allocate(Out_ne(size(Out_Distr%ne,1), size(Out_Distr%ne,2)))
    Out_ne = 0.0d0
    allocate(Out_Ee(size(Out_Distr%Ee,1), size(Out_Distr%Ee,2)))
    Out_Ee = 0.0d0
    allocate(Out_Ee_vs_E(size(Out_Distr%Ee,1), size(Out_Distr%Ee,2)))
    Out_Ee_vs_E = 0.0d0
    allocate(Out_Elat(size(Out_Distr%Ee,1), size(Out_Distr%Ee,2)))
    Out_Elat = 0.0d0
    allocate(Out_Eat_dens(size(Out_Distr%Ee,1), size(Out_Distr%Ee,2)))
    Out_Eat_dens = 0.0d0
    allocate(Out_nh(size(Out_Distr%ne,1), size(Out_Distr%ne,2), size(target_atoms), size(target_atoms(1)%Ip)))
    Out_nh = 0.0d0
    allocate(Out_Eh(size(Out_Distr%ne,1), size(Out_Distr%ne,2), size(target_atoms), size(target_atoms(1)%Ip)))
    Out_Eh = 0.0d0
    allocate(Out_Ehkin(size(Out_Distr%ne,1), size(Out_Distr%ne,2), size(target_atoms), size(target_atoms(1)%Ip)))
    Out_Ehkin = 0.0d0
    allocate(Out_field_all(Nne1, Nee2))
    Out_field_all = 0.0d0
end subroutine Allocate_out_arrays_old



subroutine Deallocate_out_arrays(Out_tot_Ne, Out_tot_Nphot, Out_tot_E, Out_E_e, Out_E_phot, Out_E_at, Out_E_h, Out_R, Out_V, Out_ne, Out_Ee, Out_nphot, Out_Ephot, Out_Ee_vs_E, &
                                 Out_Elat, Out_nh, Out_Eh, Out_Ehkin, Out_Eat_dens, Out_theta, Out_field_all, Out_Ne_Em, Out_E_Em, &
                                 Out_Ee_vs_E_Em, Out_E_field)
    real(8), dimension(:), allocatable, intent(inout) :: Out_tot_Ne
    real(8), dimension(:), allocatable, intent(inout) :: Out_tot_Nphot
    real(8), dimension(:), allocatable, intent(inout) :: Out_tot_E
    real(8), dimension(:), allocatable, intent(inout) :: Out_V
    real(8), dimension(:), allocatable, intent(inout) :: Out_E_e
    real(8), dimension(:), allocatable, intent(inout) :: Out_E_phot
    real(8), dimension(:), allocatable, intent(inout) :: Out_E_at
    real(8), dimension(:,:,:), allocatable, intent(inout) :: Out_E_h
    real(8), dimension(:), allocatable, intent(inout) :: Out_R
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_ne
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Ee
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_nphot
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Ephot
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Ee_vs_E
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Elat
    real(8), dimension(:,:,:,:), allocatable, intent(inout) :: Out_nh
    real(8), dimension(:,:,:,:), allocatable, intent(inout) :: Out_Eh
    real(8), dimension(:,:,:,:), allocatable, intent(inout) :: Out_Ehkin
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Eat_dens
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_theta
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_field_all
    
    real(8), dimension(:,:), allocatable, intent(inout) :: Out_Ee_vs_E_Em
    real(8), dimension(:), allocatable, intent(inout) :: Out_Ne_Em
    real(8), dimension(:), allocatable, intent(inout) :: Out_E_Em 
    
    real(8), dimension(:), allocatable, intent(inout) :: Out_E_field
    
    ! Dellocate temporary arrays:
    deallocate(Out_tot_Ne)
    deallocate(Out_tot_Nphot)
    deallocate(Out_tot_E)
    deallocate(Out_R)
    deallocate(Out_ne)
    deallocate(Out_Ee)
    deallocate(Out_E_phot)
    deallocate(Out_Ee_vs_E)
    deallocate(Out_Elat)
    deallocate(Out_nh)
    deallocate(Out_Eh)
    deallocate(Out_Ehkin)
    deallocate(Out_V)
    deallocate(Out_E_e)
    deallocate(Out_E_at)
    deallocate(Out_E_h)
    deallocate(Out_nphot)
    deallocate(Out_Ephot)
    deallocate(Out_Eat_dens)
    deallocate(Out_theta)
    deallocate(Out_field_all)
    deallocate(Out_Ee_vs_E_Em)
    deallocate(Out_Ne_Em)
    deallocate(Out_E_Em)
    deallocate(Out_E_Field)
end subroutine Deallocate_out_arrays

subroutine Radius_for_distributions(Target_atoms, Out_Distr, Out_V, Tim, dt, dt_flag, L)
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects
    real(8), intent(in) :: Tim, dt  ! total time, timestep [fs]
    integer, intent(in) :: dt_flag  ! kind of time-grid
    real(8), intent(in) :: L    ! [A] layer thickness that we analyse
    type(Cylinder_distr), intent(inout) :: Out_Distr   ! OUTPUT radial distributions
    real(8), dimension(:), allocatable, intent(inout) :: Out_V
    real(8) tim_glob
    real(8) R, R0
    integer i, j, N, Nt
    N = 50
    
    select case (dt_flag)   ! what kind of time-grid to use:
        case (:0)   ! linear time-grid
            Nt = CEILING(Tim/dt) !+1    ! that's how many timesteps will be there
        case (1:)   ! logarithmic time-grid
            i = 0
            tim_glob = 0.01d0
            do while (tim_glob .LT. Tim)
                i = i + 1
                tim_glob = tim_glob*dt	! [fs]
            enddo
            Nt = i + 1
    endselect

    ! set the size of all output cylindrical distribution:
    if (.not. allocated(Out_Distr%R)) allocate(Out_Distr%R(N))
    if (.not. allocated(Out_Distr%ne)) allocate(Out_Distr%ne(Nt, N))
    if (.not. allocated(Out_Distr%Ee)) allocate(Out_Distr%Ee(Nt, N))
    if (.not. allocated(Out_Distr%Atom)) allocate(Out_Distr%Atom(size(Target_atoms)))   ! how many atoms
    if (.not. allocated(Out_V)) allocate(Out_V(N))
    do i = 1, size(Out_Distr%Atom)    ! for all atoms
        if (.not. allocated(Out_Distr%Atom(i)%Shl)) allocate(Out_Distr%Atom(i)%Shl(size(Target_atoms(i)%Ip)))   ! how many atoms
        do j = 1, size(Out_Distr%Atom(i)%Shl)    ! for all shells of atom
            if (.not. allocated(Out_Distr%Atom(i)%Shl(j)%nh)) allocate(Out_Distr%Atom(i)%Shl(j)%nh(Nt, N))
            if (.not. allocated(Out_Distr%Atom(i)%Shl(j)%Eh)) allocate(Out_Distr%Atom(i)%Shl(j)%Eh(Nt, N))
            if (.not. allocated(Out_Distr%Atom(i)%Shl(j)%Ehkin)) allocate(Out_Distr%Atom(i)%Shl(j)%Ehkin(Nt, N))
        enddo ! j
    enddo   ! i
    
    R = 0.0d0
    R0 = R
    do i = 1, N
        if (i .LT. 11) then
            R = real(i)
        else if (i .LT. 21) then
            if (i .EQ. 11) then
                R = real((i-10)*15)
            else
                R = real((i-10)*10)
            endif
        else if (i .LT. 31) then
            if (i .EQ. 21) then
                R = real((i-20)*150)
            else
                R = real((i-20)*100)
            endif
        else if (i .LT. 41) then
            if (i .EQ. 31) then
                R = real((i-30)*1500)
            else
                R = real((i-30)*1000)
            endif
        else if (i .EQ. 41) then
            R = real((i-40)*15000)
        else
            R = real((i-40)*10000)
        endif
        Out_Distr%R(i) = R  ! [A] radii of cylindical layers
        Out_V(i) = 1.0d0/((R*R - R0*R0)*g_Pi*L) ! [1/A^3]
        R0 = R
    enddo ! i = 1, 50
end subroutine Radius_for_distributions



subroutine parse_time(sec,chtest)
   real(8), intent(inout) :: sec ! time interval in [sec]
   character(*), intent(out) :: chtest ! split it into miuns, hours, days...
   character(100) temp
   real(8) days, hours, mins
   days = 0.0d0
   hours = 0.0d0
   mins = 0.0d0
   if (sec .GE. 60.0d0) then
      mins = FLOOR(sec/60.0d0)
      sec = sec - mins*60.0d0
      if (mins .GT. 60.0d0) then
         hours = FLOOR(mins/60.0d0)
         mins = mins - hours*60.0d0
         if (hours .GT. 24.0d0) then
            days = FLOOR(hours/24.0d0)
            hours = hours - days*24.0d0
         endif
      endif
   endif
   chtest = ''
   temp = ''
   if (days .GT. 1.0d0) then
      write(temp, '(i9)') int(days)
      write(chtest, '(a,a)') trim(adjustl(temp)), ' days'
   else if (days .GT. 0.0d0) then
      write(temp, '(i9)') int(days)
      write(chtest, '(a,a)') trim(adjustl(temp)), ' day'
   endif
   if (hours .GT. 1.0d0) then
      write(temp, '(i9)') int(hours)
      write(chtest, '(a,a,a)') trim(adjustl(chtest)), ' '//trim(adjustl(temp)), ' hours'
   else if (hours .GT. 0.0d0) then
      write(temp, '(i9)') int(hours)
      write(chtest, '(a,a,a)') trim(adjustl(chtest)), ' '//trim(adjustl(temp)), ' hour'
   endif
   if (mins .GT. 1.0d0) then
      write(temp, '(i9)') int(mins)
      write(chtest, '(a,a,a)') trim(adjustl(chtest)), ' '//trim(adjustl(temp)), ' mins'
   else if (mins .GT. 0.0d0) then
      write(temp, '(i9)') int(mins)
      write(chtest, '(a,a,a)') trim(adjustl(chtest)), ' '//trim(adjustl(temp)), ' min'
   endif
   write(temp, '(f7.3)') sec
   write(chtest, '(a,a,a)') trim(adjustl(chtest)), ' '//trim(adjustl(temp)), ' sec'
end subroutine parse_time



end MODULE Sorting_output_data

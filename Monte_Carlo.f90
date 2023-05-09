!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains all monte-carlo subroutines

MODULE Monte_Carlo
  use Universal_Constants   ! let it use universal constants
  use Objects   ! since it uses derived types, it must know about them from module 'Objects'
  use Cross_sections    ! it's using cross-sections from that module
  use Analytical_IMFPs  ! some convenient analytical expressions are there
  use Reading_files_and_parameters  ! use some subroutines from the reading file module
implicit none

! this interface finds by itself which of the two subroutine to use depending on the parameters passed:
interface Update_electron_angles ! for updating electron angles
    module procedure Update_electron_angles_SHI
    module procedure Update_electron_angles_El
end interface Update_electron_angles
interface Next_free_path
    module procedure Next_free_path_1d
    module procedure Next_free_path_2d
end interface Next_free_path
!private  ! hides items not listed on public statement 
public :: Update_electron_angles, Next_free_path

contains    ! the MC code itself is all here:


subroutine Monte_Carlo_modelling(SHI, SHI_MFP, diff_SHI_MFP, Target_atoms, Lowest_Ip_At, Lowest_Ip_Shl, CDF_Phonon, &
            Total_el_MFPs, Elastic_MFP, Total_Hole_MFPs, Elastic_Hole_MFP, Total_Photon_MFPs, Mat_DOS, Tim, dt, Matter, NumPar, &
            Out_R, Out_V, Out_ne, Out_Ee, Out_nphot, Out_Ephot, Out_Ee_vs_E, Out_Elat, Out_nh, Out_Eh, Out_Ehkin, Out_tot_Ne, Out_tot_Nphot, Out_tot_E, &
            Out_E_e, Out_E_phot, Out_E_at, Out_E_h, Out_Eat_dens, Out_theta, Out_theta1, &
            Out_Ne_Em, Out_E_Em, Out_Ee_vs_E_Em, Error_message, DSF_DEMFP, DSF_DEMFP_H, Out_field_all, Out_E_field, &
            Out_diff_coeff)
    type(Ion), intent(in) :: SHI   ! declare SHI as an object with atributes "Ion"
    type(All_MFP), dimension(:), allocatable, intent(in) :: SHI_MFP         ! SHI mean free paths for all shells
    type(All_MFP), dimension(:), allocatable, intent(in) :: diff_SHI_MFP    ! SHI differential mean free paths for all shells
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
    integer, intent(in) :: Lowest_Ip_At, Lowest_Ip_Shl ! number of atom and of shell which correspond to the lowest ionization potential
    type(CDF), intent(in) :: CDF_Phonon ! declare CDF for phonons
    real(8), intent(in) :: Tim ! [fs] total time
    real(8), intent(in) :: dt  ! [fs] timestep
    type(All_MFP), dimension(:), intent(in), target :: Total_el_MFPs     ! electron mean free paths for all shells
    type(MFP_elastic), intent(in), target :: Elastic_MFP                 ! elastic mean free path
    type(All_MFP), dimension(:), intent(in), target :: Total_Hole_MFPs   ! hole mean free paths for all shells
    type(MFP_elastic), intent(in), target :: Elastic_Hole_MFP            ! elastic mean free path
    type(All_MFP), dimension(:), intent(in), target :: Total_Photon_MFPs ! photon MFPs for all shells
    type(Solid), intent(in) :: Matter   ! all material parameters
    type(Density_of_states), intent(in) :: Mat_DOS  ! material DOS
    type(Flag), intent(inout) :: NumPar ! numerical parameters
    real(8), dimension(:), intent(inout) :: Out_R   ! [A] radius for distributions
    real(8), dimension(:,:), intent(inout) :: Out_ne  ! [1/cm^3] electron density
    real(8), dimension(:,:), intent(inout) :: Out_Ee  ! [eV/A^3] electron energy density
    real(8), dimension(:,:), intent(inout) :: Out_nphot ! [1/cm^3] photon density
    real(8), dimension(:,:), intent(inout) :: Out_Ephot ! [eV/A^3] photon energy density
    real(8), dimension(:,:), intent(inout) :: Out_Ee_vs_E   ! [1/eV] electron distribution in energy space vs time
    real(8), dimension(:,:), intent(inout) :: Out_Elat  ! [eV/A^3] lattuce energy density vs time vs R
    real(8), dimension(:,:), intent(inout) :: Out_Eat_dens  ! [eV/A^3] atom's energy energy
    real(8), dimension(:,:,:,:), intent(inout) :: Out_nh    ! [1/cm^3] holes densities
    real(8), dimension(:,:,:,:), intent(inout) :: Out_Eh    ! [eV/A^3] holes enegies
    real(8), dimension(:,:,:,:), intent(inout) :: Out_Ehkin    ! [eV/A^3] holes enegies
    real(8), dimension(:), intent(inout) :: Out_V  ! inverse volume of cilinder layers [1/A^3]
    real(8), dimension(:), intent(inout) :: Out_tot_Ne      ! total number of electrons
    real(8), dimension(:), intent(inout) :: Out_tot_Nphot   ! total number of photons
    real(8), dimension(:), intent(inout) :: Out_tot_E       ! total energy of the system
    real(8), dimension(:), intent(inout) :: Out_E_e    ! total energy of electrons in time [eV] vs [fs]
    real(8), dimension(:), intent(inout) :: Out_E_phot ! total energy of photons in time [eV] vs [fs]
    real(8), dimension(:), intent(inout) :: Out_E_at   ! total energy of atoms in time [eV] vs [fs]
    real(8), dimension(:,:,:), intent(inout) :: Out_E_h    ! total energy of holes in each shell in time [eV] vs [fs]
    real(8), dimension(:,:), intent(inout) :: Out_theta
    real(8), dimension(:), intent(inout) :: Out_theta1
    real(8), dimension(:,:), intent(inout) :: Out_field_all ! [V/m] electrical fields vs time vs R 
    real(8), dimension(:,:), intent(inout) :: Out_Ee_vs_E_Em    ! [1/eV] emitted electrons distribution in energy space vs time
    real(8), dimension(:), intent(inout) :: Out_Ne_Em           ! number of emitted electrons
    real(8), dimension(:), intent(inout) :: Out_E_Em            ! total enegry of emitted electrons in time
    real(8), dimension(:), intent(inout) :: Out_E_field         ! total enegry of field
    type(Error_handling), intent(inout) :: Error_message	! error messages are dealed with as objects
    type(Differential_MFP), dimension(:), intent(in) :: DSF_DEMFP, DSF_DEMFP_H
    real(8), dimension(:), intent(inout) :: Out_diff_coeff
    !------------------------------------------
    
    ! Internal subroutine's variables:
    type(Ion) :: SHI_loc ! define SHI for local use
    type(Electron), dimension(:), allocatable, target :: All_electrons ! define array of electrons
    real(8), dimension(:), allocatable :: Em_electrons                 ! define array of energies of emitted electrons 
    type(Hole), dimension(:), allocatable, target :: All_holes         ! define array of holes
    type(Photon), dimension(:), allocatable, target :: All_photons     ! define array of photons
    real(8), dimension(:,:), allocatable :: SHI_path    ! to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), allocatable :: SHI_loss    ! total IMFP for SHI , to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), allocatable :: El_IMFP     ! total IMFP for electron, to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), allocatable :: El_EMFP     ! to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), allocatable :: Hole_IMFP    ! total IMFP for electron, to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), allocatable :: Hole_EMFP    ! to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), allocatable :: Phot_IMFP    ! total IMFP for photons, to use in a subroutine, we need this shape of an array
    real(8), dimension(:), allocatable :: Out_field
    real(8) Tot_field
    real(8) t_cur   ! current time of particle [fs]
    real(8) tim_glob ! global current time [fs]
    real(8) At_NRG  ! energy lost to the lattice [eV]
    real(8) RN, dEdx, MFP_tot, IMFP, EMFP, SHI_IMFP, SHI_dEdx, temp, dE, dE_cur, mh, HIMFP, HEMFP, Egap, Eloc, Ehole
    real(8) L, X, Y, Z, t
    real(8) :: Eel
    real(8) theta, phi, theta0, phi0, theta1, phi1, theta2, psi, psi2   ! for calculations of electron scattering angles
    real(8) htheta, hphi, htheta1, hphi1, htheta2, hphi2                ! for calculations of hole scattering angles
    real(8) Em_gamma, Em_E1                                             ! For calculation of electron emission
    real(8) E_new1, E_new2, R, field_time, tttt
    real(8), dimension(:), allocatable :: time_grid ! grid_points in time
    integer Sh1, KOA1, Sh2, KOA2, KOA, Sh
    integer Tot_Nel, Em_Nel ! number of electrons
    integer Tot_Nphot       ! number of photons (nonequal to number of electrons)
    integer KOP ! kind of particle that is propagating now
    integer NOP ! number of particle in an array
    integer NOA ! number of atom
    integer Nshl    ! number of shell of an atom
    integer Nat     ! total number of atoms
    integer Nat_cur, Nshl_cur
    integer Nel, N, i, j, k, N_temmp, mhole
    character(200) :: Error_descript

    ! Eckart-type surface barrier parameters
    call barrier_parameters (Matter, Em_E1, Em_gamma)
    
    Nat = size(target_atoms)    ! number of atom
    SHI_loc = SHI   ! to use it further, here you can modify it
    !eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
    ! Find how many electrons (approximately) can there be and allocate the arrays for
    ! SHI energy loss, SHI MFP, electron IMFP, and electron EMFP:
    call How_many_electrons(Nat, El_IMFP, El_EMFP, Hole_IMFP, Hole_EMFP, Phot_IMFP, SHI_path, SHI_loss, target_atoms, Matter, &
                            Lowest_Ip_At, Lowest_Ip_Shl, SHI, SHI_MFP, Total_el_MFPs, Elastic_MFP%Total, Total_Hole_MFPs, Elastic_Hole_MFP%Total, Total_Photon_MFPs, &
                            All_electrons, Em_electrons, All_holes, All_photons, Nel, NumPar)
    if (.not.allocated(All_photons)) allocate(All_photons(0))
    if (.not.allocated(Phot_IMFP)) allocate(Phot_IMFP(0,0))
    allocate(Out_field(size(Out_V)))    ! temporary array for fields
    Out_field = 0.0d0   ! no fields at the start
    Tot_field = 0.0d0
    !eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
    ! Set the time-grid:
    call set_time_grid (Tim, dt, NumPar%dt_flag, time_grid)        !from Monte_Carlo_subroutines.f90
    
    !AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    KOP = 1     ! So far we start only with SHI:
    Tot_Nel = 0 ! no excited electrons at the beginning
    Tot_Nphot = 0   ! no photons at the beginning
    Em_Nel = 0      ! no emitted electrons either
    
    i = 0   ! number of propogation timestep
    t_cur = 0.0d0       ! [fs] to start
    tim_glob = 0.0d0    ! [fs] to start
    if (NumPar%field_include .GT. 0.5) then
        field_time = 0.01d0  ! [fs] to start
    else
        field_time = tim
    endif
    call Next_free_path(SHI_loc%E, SHI_path, SHI_IMFP) ! SHI MFP [A]
    call random_number(RN)
    SHI_IMFP = -SHI_IMFP*log(RN)    ! sample free path
    call Get_time_of_next_event(SHI_loc, MFP=SHI_IMFP)  ! SHI_loc%tn is updated [fs]
    call Particle_event(SHI_loc, L=SHI_IMFP)  ! new mean free path [A]
    call Find_min_time_particle(SHI_loc, All_electrons, All_holes, All_photons, KOP, NOP, t_cur, NumPar) ! finds what particle collides next
    At_NRG = 0.0d0  ! [eV] lattice energy
    
    time_do:do while (tim_glob .LE. Tim-1e-6) ! time propagation
        i = i + 1
        tim_glob = min(time_grid(i),Tim) ! [fs]
        tttt = 1.0d0
        ! Propagate particles until next time-grid point:
        grid_do:do while (t_cur .LT. min(tim_glob,Tim))         ! propagate particles until the time grid point

                select case (KOP)
                case (1)    ! SHI
                    call SHI_Monte_Carlo(SHI_MFP, SHI_path, SHI_loc, diff_SHI_MFP, Target_atoms, All_electrons, All_holes, Tot_Nel, &
                        Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, Error_message, El_IMFP, El_EMFP, &
                        Hole_IMFP, Hole_EMFP, Matter, NumPar)
                    ! Now if SHI is out of the layer, we let it go...
                    if (SHI_loc%Z .GE. Matter%Layer) call Particle_event(SHI_loc, tn=1d16) ! SHI is out of the analyzed layer
                case (2)    ! electron
                    call Electron_Monte_Carlo(All_electrons, All_holes, El_IMFP, El_EMFP, Hole_IMFP, Hole_EMFP, &
                        CDF_Phonon, Matter, target_atoms, &
                        Total_el_MFPs, Elastic_MFP%Total, Tot_Nel, NOP, Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, &
                        Error_message, Em_electrons, &
                        Em_Nel, Em_gamma, Em_E1, At_NRG, Out_R, Out_Elat, Out_V, i, DSF_DEMFP, NumPar)
                case (3)    ! hole
                    call Hole_Monte_Carlo(All_electrons, All_holes, All_photons, El_IMFP, El_EMFP, Hole_IMFP, Hole_EMFP, &
                        Phot_IMFP, CDF_Phonon, Matter, target_atoms, &
                        Total_Hole_MFPs, Elastic_Hole_MFP%Total, Tot_Nel, Tot_Nphot, NOP, &
                        Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, Error_message, &
                        At_NRG, Out_R, Out_Elat, Out_V, i, t_cur, DSF_DEMFP_H, NumPar)
                case (4)    ! photon
                    call Photon_Monte_Carlo(All_electrons, All_holes, All_photons, El_IMFP, El_EMFP, &
                        Hole_IMFP, Hole_EMFP, Phot_IMFP, CDF_Phonon, Matter, target_atoms, &
                        Total_Photon_MFPs, Tot_Nel, Tot_Nphot, NOP, Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, Error_message, &
                        t_cur, NumPar)
                endselect
                call Find_min_time_particle(SHI_loc, All_electrons, All_holes, All_photons, KOP, NOP, t_cur, NumPar) ! finds what particle collides next
        enddo grid_do ! propagate particles until the time grid point     
        !dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
        ! Save distributions for this time-step:
        call Calculated_statistics(i, tim_glob, Tot_Nel, Tot_Nphot, At_NRG, All_electrons, Em_electrons, &
                All_holes, All_photons, Target_atoms, Out_R, Out_V, &
                Out_tot_Ne, Out_tot_Nphot, Out_tot_E, Out_E_e, Out_E_phot, Out_E_at, Out_E_h, Out_ne, Out_Ee, &
                Out_nphot, Out_Ephot, Out_Ee_vs_E, Out_Eat_dens, Out_nh, Out_Eh, Out_Ehkin, &
                Out_theta, Out_theta1, Out_Ne_Em, Out_E_Em, Out_Ee_vs_E_Em, Em_Nel, Matter, &
                Out_field, Out_field_all, Tot_field, Out_E_field, Out_diff_coeff, NumPar)
        !dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
        call Find_min_time_particle(SHI_loc, All_electrons, All_holes, All_photons, KOP, NOP, t_cur, NumPar) ! finds what particle collides next
        
    enddo time_do ! time propagation
    !AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    ! Cleal up the mess after you made it...
    deallocate(SHI_path)
    deallocate(SHI_loss)
    deallocate(El_IMFP)
    deallocate(El_EMFP)
    deallocate(Hole_IMFP)
    deallocate(Hole_EMFP)
    deallocate(All_electrons)
    deallocate(All_holes)
    if (allocated(All_photons)) deallocate(All_photons)
    if (allocated(Phot_IMFP)) deallocate(Phot_IMFP)
    if (allocated(Out_field)) deallocate(Out_field)
end subroutine Monte_Carlo_modelling


subroutine check_hole_parameters(Mat_DOS, Eel, dE, Ehole, All_electrons)
    type(Density_of_states), intent(in) :: Mat_DOS  ! material DOS
    real(8), intent(in) :: Eel
    real(8), intent(inout) :: dE
    real(8), intent(out) :: Ehole
    type(Electron), intent(inout), optional :: All_electrons
    real (8) Eloc
    integer mhole
    
    !This subroutine checks new parameters of the incident hole   
    Ehole = Eel-dE
    
    call Find_in_array_monoton(Mat_DOS%E, Ehole, mhole)
    
    if (Mat_DOS%DOS(mhole) .LT. 1.0d-4) then                !DOS is zero at this state
        
         if (present(All_electrons)) then                   !hole-electron scattering
             do while (Mat_DOS%DOS(mhole) .LT. 1.0d-4)      ! Seek the closest nonzero state
                mhole = mhole-1
             enddo
             
             Eloc = Mat_DOS%E(mhole)                        
             All_electrons%E = All_electrons%E + (Ehole-Eloc)     !Excess energy transferred to ionized electron
             Ehole = Eloc                                       ! Shift to new state
         else                                       !hole-lattice interaction
             do while (Mat_DOS%DOS(mhole) .LT. 1.0d-4)
                mhole = mhole+1
             enddo
             Eloc = Mat_DOS%E(mhole)
             dE = dE + (Ehole-Eloc)
             
             if (dE .LT. 1.0d-5) then       ! Check if lattice got negative energy
                dE = 0.0d0                  ! Set to zero lattice energy transfer
                Ehole = Eel                 ! Incident hole energy does not change
             else                           !
                Ehole = Eloc                ! Hole shifted to the energy level Eloc
             endif
         endif   
    endif    
end subroutine check_hole_parameters

!Define specific hole parameters such as effective mass and energy
subroutine Hole_parameters(All_holes, Matter, Mat_DOS, target_atoms, Hole_IMFP, Hole_EMFP, Eh, Lowest_Ip_At, Lowest_Ip_Shl)
    type(Hole), intent(inout) :: All_holes ! define array of holes
    type(Atom), dimension(:), intent(in) :: target_atoms  ! define target atoms as objects, we don't know yet how many they are
    type(Solid), intent(in) :: Matter
    type(Density_of_states), intent(in) :: Mat_DOS
    real(8), intent(in) :: Eh       ! [eV] total energy of this hole (dE-dE_cur)
    real(8), dimension(:,:), intent(in) :: Hole_IMFP  ! array of inelastic scattering MFPs
    real(8), dimension(:,:), intent(in) :: Hole_EMFP  ! array of elastic scattering MFPs
    integer, intent(in) :: Lowest_Ip_At, Lowest_Ip_Shl    ! number of atom, number of shell
    real(8) HIMFP, HEMFP, MFP_tot, RN, Eh_kin, Etemp
        
    if ((All_holes%KOA .EQ. Lowest_Ip_At) .AND. (All_holes%Shl .EQ. Lowest_Ip_Shl)) then   ! it's VB:

        Etemp = All_holes%Ehkin
        All_holes%Ehkin = Eh - Target_atoms(1)%Ip(size(Target_atoms(1)%Ip)) ! [eV] hole kinetic energy (total - Egap)
        All_holes%E = Target_atoms(1)%Ip(size(Target_atoms(1)%Ip)) ! [eV] hole potential energy (for VB = Egap)
        
        call Assign_holes_mass(Matter, All_holes, Mat_DOS, Lowest_Ip_At, Lowest_Ip_Shl)
               
        if (All_holes%Mass .LT. 1.0d3) then
            call Next_free_path(All_holes%Ehkin, Hole_IMFP, HIMFP) ! => IMFP of hole [A]
            if (Etemp .EQ. (Eh - Target_atoms(1)%Ip(size(Target_atoms(1)%Ip)))) then
                HEMFP = 1.0d30
            else
                call Next_free_path(All_holes%Ehkin, Hole_EMFP, HEMFP) ! => EMFP of hole [A]
            endif
            
            call random_number(RN)
            MFP_tot = -log(RN)/(1.0d0/HIMFP + 1.0d0/HEMFP)
            call Get_time_of_next_event(All_holes, MFP=MFP_tot, Target_atoms=Target_atoms, &
                KOA=All_holes%KOA, Shl=All_holes%Shl, Lowest_Ip_At=Lowest_Ip_At, Lowest_Ip_Shl=Lowest_Ip_Shl)
            All_holes%L = MFP_tot
        else
            All_holes%L = 1.0d30
            All_holes%tn = 1.0d30
        endif 
    else    ! it's deep shell:
        All_holes%Mass = 1.0d29
        call random_number(RN)
        ! => All_holes(Tot_Nel)%tn is updated [fs] Auger:
        call Get_time_of_next_event(All_holes, RN=RN, Target_atoms=Target_atoms, &
            KOA=All_holes%KOA, Shl=All_holes%Shl, Lowest_Ip_At=Lowest_Ip_At, Lowest_Ip_Shl=Lowest_Ip_Shl)
        All_holes%L = 1.0d30
        All_holes%E = Eh
        All_holes%Ehkin = 0.0d0
    endif 
end subroutine Hole_parameters


subroutine Assign_holes_mass(Matter, All_holes, Mat_DOS, Lowest_Ip_At, Lowest_Ip_Shl)
    type(Solid), intent(in) :: Matter   ! all material parameters
    type(Hole), intent(inout) :: All_holes ! define array of holes
    integer, intent(in) :: Lowest_Ip_At, Lowest_Ip_Shl ! VB atomic num, VB shell num
    type(Density_of_states), intent(in) :: Mat_DOS  ! material DOS
    integer :: Mnum
    
    if ((All_holes%KOA .EQ. Lowest_Ip_At) .AND. (All_holes%Shl .EQ. Lowest_Ip_Shl)) then   ! it's VB:
        if (Matter%hole_mass .GT. 0) then       ! Certain value of hole mass
            All_holes%Mass = Matter%hole_mass  ! [me]
        else                                    ! Define mass of a hole from DOS
            call find_in_array_monoton(Mat_DOS%E, All_holes%Ehkin, Mnum)
            All_holes%Mass = Mat_DOS%Eff_m(Mnum)
        endif
    else    ! it's deep shell:
        All_holes%Mass = 1.0d29  ! [me], infinite mass, immobile hole
    endif
end subroutine Assign_holes_mass


subroutine Get_time_of_next_event(Particle, RN, MFP, Target_atoms, KOA, Shl, Lowest_Ip_At, Lowest_Ip_Shl)
    class(Basic_particle), intent(inout) :: Particle ! can be Ion, or Electron, or Hole
    real(8), intent(in), optional :: RN   ! random number
    real(8), intent(in), optional :: MFP    ! [A] mean free path
    type(Atom), dimension(:), intent(in), optional :: Target_atoms  ! target atoms as objects
    integer, intent(in), optional :: KOA   ! kind of atom which this hole is sitting in
    integer, intent(in), optional :: Lowest_Ip_At
    integer, intent(in), optional :: Shl   ! number of shell which this hole is in
    integer, intent(in), optional :: Lowest_Ip_Shl
    real(8) V, nu
    
    select type(Particle)
    class is (Ion)
        call Get_velosity(Particle, V)
        if (present(MFP)) Particle%tn = Particle%t0 + MFP/V*1d5   ! [fs]
    class is (Electron)
        call Get_velosity(Particle, V)
        if (present(MFP)) Particle%tn = Particle%t0 + MFP/V*1d5   ! [fs]
    class is (Hole)
        if ((Particle%KOA .EQ. Lowest_Ip_At) .AND. (Particle%Shl .EQ. Lowest_Ip_Shl)) then   ! it's VB:
            call Get_velosity(Particle, V, Target_atoms)
            Particle%tn = Particle%t0 + MFP/V*1d5   ! [fs]
        else    ! it's deep shell:
            !Particle%tn = Particle%t0 - log(RN)*Target_atoms(KOA)%Auger(Shl)
            nu = 1.0d0/Target_atoms(KOA)%Auger(Shl) + 1.0d0/Target_atoms(KOA)%Radiat(Shl) ! Auger ot Radiative decay
            Particle%tn = Particle%t0 - log(RN)/nu
        endif
    class is (Photon) 
        call Get_velosity(Particle, V)
        if (present(MFP)) Particle%tn = Particle%t0 + MFP/V*1d5   ! [fs]
   end select
end subroutine Get_time_of_next_event

subroutine Get_velosity(Particle, V, Target_atoms)
    class(Basic_particle), intent(in) :: Particle ! can be Ion, or Electron, or Hole
    type(Atom), dimension(:), intent(in), optional :: Target_atoms  ! target atoms as objects
    real(8), intent(out) :: V   ! [m/s] velosity corresponding to this energy
    real(8) Egap
    select type(Particle)
    class is (Ion)
        V = SQRT(2.0d0*Particle%E*g_e/(Particle%Mass*g_Mp)) ! [m/s] free-flight time of an ion
    class is (Electron)
        V = SQRT(2.0d0*Particle%E*g_e/g_me) ! [m/s] free-flight time of an electron
    class is (Hole)
        if (Particle%Mass .LT. 1.0e6) then
            !Egap = Target_atoms(1)%Ip(size(Target_atoms(1)%Ip)) ! [eV] bandgap
            !V = SQRT(2.0d0*(Particle%E - Egap)*g_e/(Particle%Mass*g_me)) ! [m/s] free-flight time of a hole
            V = SQRT(2.0d0*(Particle%Ehkin)*g_e/(Particle%Mass*g_me))
        else
            V = 0.0d0   ! [m/s] immobile hole
        endif
    class is (Photon)
        V = g_cvel ! [m/s] those guys always move with the speed of light
    end select
end subroutine Get_velosity


subroutine Calculated_statistics(i, tim, Tot_Nel, Tot_Nphot, At_NRG, All_electrons, Em_electrons, &
        All_holes, All_photons, target_atoms, Out_R, Out_V, &
        Out_tot_Ne, Out_tot_Nphot, Out_tot_E, Out_E_e, Out_E_phot, Out_E_at, Out_E_h, &
        Out_ne, Out_Ee, Out_nphot, Out_Ephot, Out_Ee_vs_E, Out_Eat_dens, Out_nh, Out_Eh, Out_Ehkin, &
        Out_theta, Out_theta1, Out_Ne_Em, Out_E_Em, Out_Ee_vs_E_Em, Em_Nel, Matter, &
        Out_field, Out_field_all, Tot_field, Out_E_field, Out_diff_coeff, NumPar)
    integer, intent(in) :: i    ! number of the time-step
    real(8), intent(in) :: tim
    integer, intent(in) :: Tot_Nel  ! number of electrons
    integer, intent(in) :: Tot_Nphot  ! number of photons
    real(8), intent(in) :: At_NRG   ! [eV] atoms' energy
    type(Electron), dimension(:), intent(in) :: All_electrons ! define array of electrons
    real(8), dimension(:), intent(in) :: Em_electrons
    type(Hole), dimension(:), intent(in) :: All_holes   ! array of all holes
    type(Photon), dimension(:), intent(in) :: All_photons ! array of all photons
    type(Atom), dimension(:), intent(in) :: target_atoms  ! define target atoms as objects, we don't know yet how many they are
    real(8), dimension(:), intent(in) :: Out_R   ! [A] radius for distributions
    real(8), dimension(:), intent(in) :: Out_V  ! inverse volume of cilinder layers [1/A^3]
    real(8), dimension(:), intent(inout) :: Out_tot_Ne
    real(8), dimension(:), intent(inout) :: Out_tot_Nphot
    real(8), dimension(:), intent(inout) :: Out_tot_E
    real(8), dimension(:), intent(inout) :: Out_E_e
    real(8), dimension(:), intent(inout) :: Out_E_phot
    real(8), dimension(:), intent(inout) :: Out_E_at
    real(8), dimension(:,:), intent(inout) :: Out_ne  ! [1/cm^3] electron density
    real(8), dimension(:,:), intent(inout) :: Out_Ee  ! [eV/A^3] electron energy density
    real(8), dimension(:,:), intent(inout) :: Out_nphot ! [1/cm^3] photon density
    real(8), dimension(:,:), intent(inout) :: Out_Ephot ! [eV/A^3] photon energy density
    real(8), dimension(:,:), intent(inout) :: Out_Ee_vs_E   ! [1/eV] electron distribution in energy space
    real(8), dimension(:,:), intent(inout) :: Out_Eat_dens  ! [eV/A^3] atom's energy energy
    real(8), dimension(:,:,:), intent(inout) :: Out_E_h
    real(8), dimension(:,:,:,:), intent(inout) :: Out_nh    ! [1/cm^3] holes densities
    real(8), dimension(:,:,:,:), intent(inout) :: Out_Eh    ! [eV/A^3] holes enegies
    real(8), dimension(:,:,:,:), intent(inout) :: Out_Ehkin
    real(8), dimension(:,:), intent(inout) :: Out_theta
    real(8), dimension(:), intent(inout) :: Out_theta1
    type(Solid), intent(in) :: Matter
    type(Flag), intent(in) :: NumPar
    
    real(8), dimension(:), intent(inout) :: Out_diff_coeff
    
    real(8), dimension(:), allocatable :: Out_E
    real(8), dimension(:,:), intent(inout) :: Out_Ee_vs_E_Em
    real(8), dimension(:), intent(inout) :: Out_Ne_Em
    real(8), dimension(:), intent(inout) :: Out_E_Em
    integer, intent(in) :: Em_Nel 
    real(8), dimension(:), intent(in) :: Out_field      ! current fields
    real(8), dimension(:,:), intent(inout) :: Out_field_all ! fields at all times to save
    real(8), dimension(:), intent(inout) :: Out_E_field
    real(8), intent(in) :: Tot_field
    
    real(8) R, X, Y, Z, L0, theta0, phi0, V,xx
    real(8) Xh, Yh
    integer j, k, l, m, Nat, N0Ee, N_VB_h
    
    Nat = size(target_atoms)    ! number of atom
    
    allocate(Out_E(size(Out_R,1)))
    Out_E(:) = Out_R(:)/10.0
    
    Out_tot_Ne(i) = Out_tot_Ne(i) + real(Tot_Nel) ! total the number of electrons, not normalized
    if (NumPar%include_photons) then ! count photons in:
        Out_tot_Nphot(i) = Out_tot_Nphot(i) + real(Tot_Nphot) ! total the number of photons, not normalized
        Out_E_phot(i) = Out_E_phot(i) + SUM(All_photons(:)%E)   ! [eV] total photon energy
        Out_tot_E(i) = Out_tot_E(i) + SUM(All_electrons(:)%E) + SUM(All_holes(:)%E) + &
            SUM(All_holes(:)%Ehkin) + SUM(All_photons(:)%E) + At_NRG + Tot_field    ! total energy, not normalized
    else    ! no photons:
        Out_tot_E(i) = Out_tot_E(i) + SUM(All_electrons(:)%E) + SUM(All_holes(:)%E) + &
            SUM(All_holes(:)%Ehkin) + At_NRG + Tot_field    ! total energy, not normalized
    endif
    Out_E_e(i) = Out_E_e(i) + SUM(All_electrons(:)%E)   ! [eV] total electron energy
    Out_E_at(i) = Out_E_at(i) + At_NRG                  ! [eV] total lattice energy
    Out_E_field(i) = Out_E_field(i) + Tot_field         ! [eV] total field energy
    Out_field_all(i,:) = Out_field_all(i,:) + Out_field(:)    ! [V/m] save field vs R for the time t(i)
    
    Out_Ne_Em(i) = Out_Ne_Em(i) + real(Em_Nel)
    Out_E_Em(i) = Out_E_Em(i) + SUM(Em_electrons(:))
    
    do k = 1, Nat != size(target_atoms)    ! number of atom
        do j = 1, size(target_atoms(k)%Ip)  ! all shells
            Out_E_h(i,k,j) = Out_E_h(i,k,j) + SUM(All_holes(:)%E, &
                MASK = ((All_holes(:)%KOA .EQ. k) .AND. (All_holes(:)%Shl) .EQ. j)) + &
                SUM(All_holes(:)%Ehkin, MASK = ((All_holes(:)%KOA .EQ. k) .AND. (All_holes(:)%Shl) .EQ. j)) ! energy of this shell [eV]
        enddo
    enddo
   
    if (NumPar%include_photons) then !only if we included photons:
        if (Tot_Nphot .GT. 0) then   !and only if they are there
            do k = 1, Tot_Nphot   ! for all photons
                if (All_photons(k)%E .GT. 0.0d0) then
                    call Get_velosity(All_photons(k), V)
                    L0 = V*(tim - All_photons(k)%t0)*1.0d-5 ! [A] mean free path travelled until this time instance
                    if (L0 .LT. 0.0d0) L0 = 0.0d0  ! just fofr case...
                    theta0 = All_photons(k)%theta   ! old angle
                    phi0 = All_photons(k)%phi       ! old angle

                    X = All_photons(k)%X + L0*sin(theta0)*sin(phi0)    ! [A] new X coordinate
                    Y = All_photons(k)%Y + L0*sin(theta0)*cos(phi0)    ! [A] new Y coordinate
                    R = SQRT(X*X + Y*Y) ! [A] radius of this el
                    if (isnan(R)) then
                        print*, 'For some reason, there is NaN occured in photon statistics:'
                        print*, 'R', R, X, Y, L0, theta0, phi0
                        print*, All_photons(k)%X, All_photons(k)%Y, All_photons(k)%E, All_photons(k)%t0, All_photons(k)%tn
                    endif
                    call Find_in_array_monoton(Out_R, R, j) ! find where in the distribution array it is
                    Out_nphot(i,j) = Out_nphot(i,j) + Out_V(j)   ! [1/A^3] here is the photon
                    Out_Ephot(i,j) = Out_Ephot(i,j) + All_photons(k)%E*Out_V(j)   ! [eV/A^3] here is the photon
                endif
            enddo
        endif
    endif
   
    N_VB_h = 0
    do k = 1, Tot_Nel   ! for all electrons
        if (All_electrons(k)%E .GT. max(0.0d0, Matter%cut_off)) then
            call Get_velosity(All_electrons(k), V)
            L0 = V*(tim - All_electrons(k)%t0)*1.0d-5 ! [A] mean free path travelled until this time instance
            if (L0 .LT. 0.0d0) L0 = 0.0d0  ! just fofr case...
            theta0 = All_electrons(k)%theta   ! old angle
            phi0 = All_electrons(k)%phi       ! old angle
        else
            L0 = 0.0d0
            theta0 = 0.0d0   ! old angle
            phi0 = 0.0d0       ! old angle
        endif
                                        
        X = All_electrons(k)%X + L0*sin(theta0)*sin(phi0)    ! [A] new X coordinate
        Y = All_electrons(k)%Y + L0*sin(theta0)*cos(phi0)    ! [A] new Y coordinate
        R = SQRT(X*X + Y*Y) ! [A] radius of this el
        if (isnan(R)) then
            print*, 'For some reason, there is NaN occured in electron statistics:'
            print*, 'R', R, X, Y, L0, theta0, phi0
            print*, All_electrons(k)%X, All_electrons(k)%Y, All_electrons(k)%E, All_electrons(k)%t0, All_electrons(k)%tn
        endif
        call Find_in_array_monoton(Out_R, R, j) ! find where in the distribution array it is
        Out_ne(i,j) = Out_ne(i,j) + Out_V(j)   ! [1/A^3] here is the electron
        Out_Ee(i,j) = Out_Ee(i,j) + All_electrons(k)%E*Out_V(j)   ! [eV/A^3] here is the electron
        
        if (isnan(All_electrons(k)%E)) print*, 'Negative electron energy in stat ', 'E = ', All_electrons(k)%E
        call Find_in_array_monoton(Out_R, All_electrons(k)%E, j) ! find where in the distribution array it is
        if (j .GT. 1) then
            Out_Ee_vs_E(i,j) = Out_Ee_vs_E(i,j) + 1.0d0/(Out_R(j)-Out_R(j-1))/real(Tot_Nel) ! number of electrons in this ENERGY interval
        else    ! j = 1, energy interval is 1
            Out_Ee_vs_E(i,j) = Out_Ee_vs_E(i,j) + 1.0d0/Out_R(j)/real(Tot_Nel) ! number of electrons in this ENERGY interval
        endif
        
        xx = theta0/g_Pi*180.0d0 ! Electron angles array 1=phi, 2=theta
        if (isnan(xx)) then
            print*, 'The problem is with angles in stst...', xx, All_electrons(k)%theta, All_electrons(k)%E, k
        endif
        if (All_electrons(k)%E .GT. 0.0d0) then ! exclude zero-energy electrons
            call Find_in_array_monoton(Out_theta1, xx, j) ! find where in the distribution array it is
            Out_theta(i,j) = Out_theta(i,j) + 1.0d0/real(Tot_Nel)
        endif
        
        ! DO the same for holes
        !if (All_holes(k)%Mass .LT. 1.0d3) then
        if ((All_holes(k)%Mass .LT. 1.0d3) .and. (All_holes(k)%Ehkin .GT. max(0.0d0, Matter%cut_off))) then
            call Get_velosity(All_holes(k), V)
            L0 = V*(tim - All_holes(k)%t0)*1.0d-5 ! [A] mean free path travelled until this time instance
            if (L0 .LT. 0.0d0) L0 = 0.0d0  ! just for case...
            theta0 = All_holes(k)%theta   ! old angle
            phi0 = All_holes(k)%phi       ! old angle
            Xh = All_holes(k)%X + L0*sin(theta0)*sin(phi0)    ! [A] new X coordinate
            Yh = All_holes(k)%Y + L0*sin(theta0)*cos(phi0)    ! [A] new Y coordinate
        else
            Xh = All_holes(k)%X
            Yh = All_holes(k)%Y
        endif
        R = SQRT(Xh*Xh + Yh*Yh) ! [A] radius of this hole
        
        if (isnan(R)) then
            print*, 'For some reason, there is NaN occured in hole statistics: ' , R
            print*, k, All_holes(k)%X, All_holes(k)%Y
        endif
        call Find_in_array_monoton(Out_R, R, j) ! find where in the distribution array it is
        l = All_holes(k)%KOA
        m = All_holes(k)%Shl
        Out_nh(i,j,l,m) = Out_nh(i,j,l,m) + Out_V(j)   ! [1/A^3] here is the hole
                
        Out_Eh(i,j,l,m) = Out_Eh(i,j,l,m) + All_holes(k)%E*Out_V(j)   ! [eV/A^3] here is the hole
        Out_Ehkin(i,j,l,m) = Out_Ehkin(i,j,l,m) + All_holes(k)%Ehkin*Out_V(j)   ! [eV/A^3] here is the hole
        
        
        !Average hole diffusion coefficient over all VB holes
        if (All_holes(k)%Mass .LT. 1.0d3 ) then
            !if ((All_Holes(k)%Ehkin .LT. 4.0d0) .AND. (All_holes(k)%L .LT. 1.0d3)) then
            if (All_holes(k)%L .LT. 1.0d3) then
                N_VB_h = N_VB_h + 1
                call Get_velosity(All_holes(k), V)
                Out_diff_coeff(i) = Out_diff_coeff(i) + 1.0d0/3.0d0*V*All_holes(k)%L*1.0d-6        ![cm^2/s]
            endif
        endif
    enddo
    if (N_VB_h > 0) Out_diff_coeff(i) = Out_diff_coeff(i)/N_VB_h
    
    ! Emitted electrons:
    do k = 1, Em_Nel
        if (isnan(Em_electrons(k))) print*, 'Negative emitted electron energy in stat ', 'E = ', Em_electrons(k)
        call Find_in_array_monoton(Out_E, Em_electrons(k), j) ! find where in the distribution array it is
        if (j .GT. 1) then
            Out_Ee_vs_E_Em(i,j) = Out_Ee_vs_E_Em(i,j) + 1.0d0/(Out_E(j)-Out_E(j-1))/dble(Em_Nel)     !/Em_Nel ! number of electrons in this ENERGY interval
        else    ! j = 1, energy interval is 1
            Out_Ee_vs_E_Em(i,j) = Out_Ee_vs_E_Em(i,j) + 1.0d0/Out_E(j)/dble(Em_Nel)    !/Em_Nel ! number of electrons in this ENERGY interval
        endif
    end do
end subroutine Calculated_statistics    


function Impact_parameter(SHI_loc, dE)
    real(8) Impact_parameter
    type(Ion), intent(in), target :: SHI_loc ! define SHI for local use
    real(8), intent(in) :: dE
    real(8), pointer :: ZSHI, ESHI
    real(8) MSHI, A
    MSHI = g_Mp*SHI_loc%Mass    ! [kg] mass of SHI
    ESHI => SHI_loc%E           ! [eV] energy
    ZSHI => SHI_loc%Zeff        ! SHI charge
    A = 1.0d0+MSHI/g_me
    Impact_parameter = g_a0*ZSHI*g_Ry/ESHI*sqrt(4.0d0*ESHI/dE*MSHI/g_me - A*A)  ! [A]
    nullify(ZSHI)
    nullify(ESHI)
end function Impact_parameter


! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Dealing with angles:

subroutine Update_holes_angles_SHI(dE, theta, phi)
   REAL(8), INTENT(in) :: dE     	! Transferred electron energy [eV]
   REAL(8), INTENT(out) :: theta, phi   ! scattering angles for this transferred energy
   real(8) RN, RN2
   call random_number(RN)
   theta = g_Pi*RN ! hole theta-angle
   call random_number(RN2)
   phi = 2.0d0*g_Pi*RN2      ! hole phi-angle
end subroutine Update_holes_angles_SHI

subroutine Update_holes_angles_el(All_holes, Eel, dE, theta, phi, theta1, phi1)
   type(Hole), intent(in), target :: All_holes ! define SHI for local use
   REAL(8), INTENT(in) :: dE         	! Transferred energy [eV]
   real(8), intent(in) :: Eel           ! Energy of a hole
   REAL(8), INTENT(out) :: theta, phi   ! scattering angles for this transferred energy
   REAL(8), INTENT(out) :: theta1, phi1 ! scattering angles of incident particle
   REAL(8) Mh     ! hole mass [kg]
   real(8) RN, RN2, E11
   E11 = Eel - dE
   Mh = All_holes%Mass*g_me
   theta = ACOS(SQRT((Mh + g_me)*(Mh + g_me)/(4.0e0*Mh*g_me)*dE/Eel))   ! scattering angle
   call random_number(RN2)
   phi = 2.0d0*g_Pi*RN2      ! electron phi-angle
   
   if (isnan(theta)) then
        call random_number (RN)
        theta = g_Pi*RN
   endif
   
   theta1 = acos((Eel*(Mh-g_me)+E11*(Mh+g_me))/(2*Mh*sqrt(Eel*E11)))           ! Zenith angle of incident hole
   phi1 = phi + g_Pi                                            ! phi angle of incident hole
   
   if (isnan(theta1)) then
        call random_number (RN)
        theta1 = g_Pi*RN
   endif
end subroutine Update_holes_angles_el

subroutine Update_electron_angles_SHI(SHI_loc, dE, theta, phi)
   type(Ion), intent(in), target :: SHI_loc ! define SHI for local use
   REAL(8), INTENT(in) :: dE     	! Transferred electron energy [eV]
   REAL(8), INTENT(out) :: theta, phi   ! scattering angles for this transferred energy
   REAL(8), pointer :: E     	! SHI energy [eV]
   REAL(8) MSHI     ! SHI mass [kg]
   real(8) RN
   MSHI = SHI_loc%Mass*g_Mp     !   [kg] mass of SHI
   E => SHI_loc%E    ! [eV] SHI energy
   if (E .LE. 0.0d0) then
        theta = g_Pi/2.0d0  ! exactly, when energy is zero
   else
        theta = ACOS(SQRT((MSHI + g_me)*(MSHI + g_me)/(4.0e0*MSHI*g_me)*dE/E))   ! scattering angle
   endif
   call random_number(RN)
   phi = 2.0d0*g_Pi*RN      ! Electron phi-angle
   nullify(E)
end subroutine Update_electron_angles_SHI

subroutine Update_electron_angles_El(E, dE, theta, phi)
   REAL(8), INTENT(in) :: E     	! Initial electron energy [eV]
   REAL(8), INTENT(in) :: dE     	! Transferred electron energy [eV]
   REAL(8), INTENT(out) :: theta, phi ! scattering angle for this transferred energy
   real(8) RN
   theta = ACOS((E - dE)/sqrt(E*(E-dE)))   ! electron theta angle
   if (isnan(theta)) then
	call random_number(RN)
	theta = RN*g_pi
   endif
   
   call random_number(RN)
   phi = 2.0d0*g_Pi*RN      ! Electron phi-angle
end subroutine Update_electron_angles_El

subroutine Update_particle_angles_lat(target_atoms, E, dE, theta, phi, Mass)    ! Angle of an incident particle after scattering on lattice atoms
   type(Atom), dimension(:), intent(in) :: target_atoms  ! define target atoms as objects, we don't know yet how many they are
   REAL(8), INTENT(in) :: E     	! Initial electron energy [eV]
   REAL(8), INTENT(in) :: dE     	! Transferred electron energy [eV]
   REAL(8), INTENT(out) :: theta, phi ! scattering angle for this transferred energy
   real(8), intent(in), optional :: Mass    ! Mass of incident particle. Should be present for hole.
   real(8) RN, RN2, Mtarget, E11, arg
   
   !E11 = E-dE
   E11 = abs(E-dE)
   Mtarget = g_Mp*SUM(target_atoms(:)%Mass*dble(target_atoms(:)%Pers))/dble(SUM(target_atoms(:)%Pers))
   if (present(Mass)) then          ! Hole
        arg = (E*(Mass*g_me - Mtarget)+E11*(Mass*g_me + Mtarget))/(2*Mass*g_me*sqrt(E*E11))
   else                             ! Electron
        arg = (E*(g_me - Mtarget)+E11*(g_me + Mtarget))/(2*g_me*sqrt(E*E11))
   endif
   if (abs(arg) > 1.0d0) then
      theta = 0.0d0
   else
      theta = acos(arg)
   endif

   ! In case there is any problem:
   if (isnan(theta)) then
        call random_number (RN)
        theta = g_Pi*RN
   endif     
    
   call random_number(RN2)
   phi = 2.0d0*g_Pi*RN2      ! Electron phi-angle
end subroutine Update_particle_angles_lat

subroutine Check_Angles_both(phi1, theta1) ! new phi1 and theta1 angles
   real(8), INTENT(inout) :: phi1, theta1

   if (theta1 .LT. 0.0d0) then  ! Faster periodic boundaries
        print*, '(theta1 .LT. 0.0e0)', theta1 , theta1+CEILING(ABS(theta1)/g_Pi)*g_Pi
   endif
   if (theta1 .GT. g_Pi) then  ! Faster periodic boundaries
        print*, '(theta1 .GT. g_Pi)', theta1 , theta1-FLOOR(theta1/g_Pi)*g_Pi
   endif
   
   do while (theta1 .LT. 0.0d0)
      theta1 = ABS(theta1)
      phi1 = phi1 + g_Pi
   enddo
   do while (theta1 .GT. g_Pi)
      theta1 = 2.0e0*g_Pi - theta1
      phi1 = phi1 - g_Pi
   enddo
   if (phi1 .LT. 0.0d0) then  ! Faster periodic boundaries
        phi1 = phi1+CEILING(ABS(phi1)/(2.0e0*g_Pi))*2.0e0*g_Pi
   endif
   if (phi1 .GT. 2.0e0*g_Pi)then  ! Faster periodic boundaries
        phi1 = phi1-FLOOR(phi1/(2.0e0*g_Pi))*2.0e0*g_Pi
   endif
end subroutine Check_Angles_both

subroutine New_Angles_both(phi0, theta0, theta, psi, phi1, theta1) ! new phi1 and theta1 angles
   real(8), INTENT(in) :: phi0, theta0, theta, psi
   real(8), INTENT(out) :: phi1, theta1

   call New_phi(phi0, theta0, theta, psi, phi1) ! new phi1 angle after scattering
   call New_theta(theta0, theta, psi, theta1)
   do while (theta1 .LT. 0.0e0)
      theta1 = ABS(theta1)
      phi1 = phi1 + g_Pi
   enddo
   do while (theta1 .GT. g_Pi)
      theta1 = 2.0e0*g_Pi - theta1
      phi1 = phi1 - g_Pi
   enddo
   if (phi1 .GT. 2.0e0*g_Pi)then  ! Faster periodic boundaries
        phi1 = phi1-FLOOR(phi1/(2.0e0*g_Pi))*2.0e0*g_Pi
   endif
   if (phi1 .LT. 0.0e0) then  ! Faster periodic boundaries
        phi1 = phi1+CEILING(ABS(phi1)/(2.0e0*g_Pi))*2.0e0*g_Pi
   endif
end subroutine New_Angles_both

subroutine New_phi(phi0, theta0, theta, psi, phi1) ! new phi1 angle after scattering
   real(8), INTENT(in) :: phi0, theta0, theta, psi
   real(8), INTENT(out) :: phi1
   phi1 = phi0 + theta*cos(theta0)*sin(psi) ! transferring from relative angles (theta, psi) to absolute
end subroutine New_phi

subroutine New_theta(theta0, theta, psi, theta1) ! new theta1 angle after scattering
   real(8), INTENT(in) :: theta0, theta, psi
   real(8), INTENT(out) :: theta1
   theta1 = theta0 + theta*cos(psi) ! transferring from relative angles (theta, psi) to absolute
end subroutine New_theta


! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Dealing with Auger-electrons:

subroutine Auger_decay(KOA, SHL, Target_atoms, Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, &
                       Sh1, KOA1, Sh2, KOA2, Ee, E_new1, E_new2, Error_message)
   integer, intent(in) :: KOA, SHL ! atomic species and shell numbers
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects
   integer, intent(in) :: Lowest_Ip_At, Lowest_Ip_Shl ! number of atom and of shell which correspond to the lowest ionization potential
   type(Density_of_states), intent(in) :: Mat_DOS   ! material DOS
   integer, intent(out) :: KOA1, Sh1	! shell and atom to which deeper hole jumps up
   integer, intent(out) :: KOA2, Sh2	! shell and atom in which the hole is produced
   real(8), intent(out) :: Ee, E_new1, E_new2	! [eV] energies of 1) ejected electron, 2) decayed hole, 3) new hole
   type(Error_handling), optional, intent(inout) :: Error_message ! deals with errors, if any
   character(200) Writing_var
   real(8) RN, Energy_diff, dE_cur
   integer i, j, coun, N, M, iter, Shel, coun_sh, SHL1, SHL2

   N = size(Target_atoms) ! number of elements
   Sh1 = SHL
   Sh2 = 0
   KOA1 = KOA
   KOA2 = 0
   Ee = -1.0d-10
   iter = 0 ! count nmber of iterations
   coun = 0
  !1111111111111111111111111111111111111111111111111111111111111111111111111111
  ! Hole one:
  ! count how many electrons on the shells that can participate in Auger-decay:
  call count_for_Auger_shells(Target_atoms, Target_atoms(KOA)%Ip(SHL), coun)
  call random_number(RN)
  Shel = 1 + nint(RN*(coun-1)) ! this shell is where the hole jumps up
  call Choose_for_Auger_shell(Target_atoms, Target_atoms(KOA)%Ip(SHL), Shel, Sh1, KOA1)
  dE_cur = 0.0d0
  if (allocated(Mat_DOS%E)) then    ! it is VB instead of just a shell
     if ((KOA1 .EQ. Lowest_Ip_At) .AND. (Sh1 .EQ. Lowest_Ip_Shl)) then ! it is VB
        call From_where_in_VB(Mat_DOS, dE_cur)
     endif         
  endif
  E_new1 = dE_cur + Target_atoms(KOA1)%Ip(Sh1)  ! new energy of the hole [eV]
  Energy_diff = Target_atoms(KOA)%Ip(SHL) - E_new1   ! [eV] energy difference between the 2 levels
  !2222222222222222222222222222222222222222222222222222222222222222222222222222
  ! Hole two:
  ! count how many electrons on the shells that can participate in Auger-decay:
  call count_for_Auger_shells(Target_atoms, Energy_diff, coun)

  if (coun .GT. 0) then
     call random_number(RN)
     Shel = 1 + nint(RN*(coun-1)) ! this shell is where the hole jumps up
     call Choose_for_Auger_shell(Target_atoms, Energy_diff, Shel, Sh2, KOA2)
     dE_cur = 0.0d0
     if (allocated(Mat_DOS%E)) then    ! it is VB instead of just a shell
        if ((KOA2 .EQ. Lowest_Ip_At) .AND. (Sh2 .EQ. Lowest_Ip_Shl)) then ! it is VB
           call From_where_in_VB(Mat_DOS, dE_cur, Energy_diff-Target_atoms(KOA2)%Ip(Sh2))
        endif   ! it is deep shell
     endif
     E_new2 = dE_cur + Target_atoms(KOA2)%Ip(Sh2)  ! energy of the new hole [eV]
     Ee = Energy_diff - E_new2  ! [eV] energy of the emmited electron
  else
     print*, 'Now, second shell of Auger is not defined...'
  endif ! (coun .GT. 0)

  if ((Ee .LT. 0.0d0) .OR. (E_new1 .LT. 0.0d0) .OR. (E_new2 .LT. 0.0d0)) then
     write(Writing_var,'(a,e,e,i,i,i,i,i,i)') 'Impossible Auger-decay ', Target_atoms(KOA1)%Ip(Sh1), &
        Target_atoms(KOA1)%Ip(Sh1), KOA, SHL, KOA1, Sh1, KOA2, Sh2
     write(*,'(a)') trim(adjustl(Writing_var))
     call Save_error_details(Error_message, 25, Writing_var)
  endif
end subroutine Auger_decay



subroutine Auger_decay_old(KOA, SHL, Target_atoms, Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, &
                       Sh1, KOA1, Sh2, KOA2, Ee, E_new1, E_new2, Error_message)
   integer, intent(in) :: KOA, SHL ! atomic species and shell numbers
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects
   integer, intent(in) :: Lowest_Ip_At, Lowest_Ip_Shl ! number of atom and of shell which correspond to the lowest ionization potential
   type(Density_of_states), intent(in) :: Mat_DOS   ! material DOS
   integer, intent(out) :: KOA1, Sh1	! shell and atom to which deeper hole jumps up
   integer, intent(out) :: KOA2, Sh2	! shell and atom in which the hole is produced
   real(8), intent(out) :: Ee, E_new1, E_new2	! [eV] energies of 1) ejected electron, 2) decayed hole, 3) new hole
   type(Error_handling), optional, intent(inout) :: Error_message ! deals with errors, if any
   character(200) Writing_var
   real(8) RN, Energy_diff, dE_cur
   integer i, j, coun, N, M, iter, Shel, coun_sh, SHL1, SHL2

   N = size(Target_atoms) ! number of elements
   Sh1 = SHL
   Sh2 = 0
   KOA1 = KOA
   KOA2 = 0
   Ee = -1.0d-10
   iter = 0 ! count nmber of iterations
   coun = 0
   do while ((coun .LE. 0) .OR. (Ee .LT. 0.0d0))
      iter = iter + 1 ! count nmber of iterations
      !1111111111111111111111111111111111111111111111111111111111111111111111111111
      ! Hole one:
      ! count how many electrons on the shells that can participate in Auger-decay:
      call count_for_Auger_shells(Target_atoms, Target_atoms(KOA)%Ip(SHL), coun)
      call random_number(RN)
      Shel = 1 + nint(RN*(coun-1)) ! this shell is where the hole jumps up
      call Choose_for_Auger_shell(Target_atoms, Target_atoms(KOA)%Ip(SHL), Shel, Sh1, KOA1)
      dE_cur = 0.0d0
      if (allocated(Mat_DOS%E)) then    ! it is VB instead of just a shell
         if ((KOA1 .EQ. Lowest_Ip_At) .AND. (Sh1 .EQ. Lowest_Ip_Shl)) then ! it is VB
            call From_where_in_VB(Mat_DOS, dE_cur)
         endif         
      endif
      E_new1 = dE_cur + Target_atoms(KOA1)%Ip(Sh1)  ! new energy of the hole [eV]
      Energy_diff = Target_atoms(KOA)%Ip(SHL) - E_new1   ! [eV] energy difference between the 2 levels
      !2222222222222222222222222222222222222222222222222222222222222222222222222222
      ! Hole two:
      ! count how many electrons on the shells that can participate in Auger-decay:
      call count_for_Auger_shells(Target_atoms, Energy_diff, coun)

      if (coun .GT. 0) then
         call random_number(RN)
         Shel = 1 + nint(RN*(coun-1)) ! this shell is where the hole jumps up
         call Choose_for_Auger_shell(Target_atoms, Energy_diff, Shel, Sh2, KOA2)
         dE_cur = 0.0d0
         if (allocated(Mat_DOS%E)) then    ! it is VB instead of just a shell
            if ((KOA2 .EQ. Lowest_Ip_At) .AND. (Sh2 .EQ. Lowest_Ip_Shl)) then ! it is VB
               call From_where_in_VB(Mat_DOS, dE_cur, Energy_diff)
            endif   ! it is deep shell
         endif
         E_new2 = dE_cur + Target_atoms(KOA2)%Ip(Sh2)  ! energy of the new hole [eV]
         Ee = Energy_diff - E_new2  ! [eV] energy of the emmited electron
      endif ! (coun .GT. 0)

      if (iter .GE. 100) then
         write(Writing_var,'(a,e,e,i,i,i,i,i,i)') 'Impossible Auger-decay ', Target_atoms(KOA1)%Ip(Sh1), &
            Target_atoms(KOA1)%Ip(Sh1), KOA, SHL, KOA1, Sh1, KOA2, Sh2
         write(*,'(a)') trim(adjustl(Writing_var))
         call Save_error_details(Error_message, 25, Writing_var)
      endif                
      if (iter .GE. 100) exit
   enddo ! while
end subroutine Auger_decay_old


subroutine From_where_in_VB(Mat_DOS, dE_cur, E)
    type(Density_of_states), intent(in) :: Mat_DOS   ! material DOS
    real(8), intent(out) :: dE_cur  ! [eV] energy
    real(8), intent(in), optional :: E  ! [eV] given energy, maximum for the search within VB
    real(8) RN, Sum_DOS, Tot_N
    integer N, M_temp, N_temmp
    
    if (.not. present(E)) then  ! search the whole VB randomly :
        if (.not. allocated(Mat_DOS%E)) then  ! VB is treated as atomic energy level:
            dE_cur = 0.0d0  ! add nothing to the band-gap [eV]
        else    ! it is a band, there is DOS, use it:
            N = size(Mat_DOS%int_DOS)   ! the last element tells the number of electrons in VB
            Sum_DOS = Mat_DOS%int_DOS(N)   ! only these electrons of VB are available
            call random_number(RN)
            Tot_N = RN*Sum_DOS   ! chose this electron
            call Find_in_array_monoton(Mat_DOS%int_DOS, Tot_N, N_temmp) ! make increasing array
            dE_cur = Mat_DOS%E(N_temmp)    ! add this energy to the band gap [eV]
        endif
     else   ! search only in a part of VB:
        if (.not. allocated(Mat_DOS%E)) then  ! VB is treated as atomic energy level:
            dE_cur = 0.0d0  ! [eV]
        else    ! it is a band, there is DOS, use it:
            N = size(Mat_DOS%int_DOS)   ! the last element tells the number of electrons in VB
            if (E .LT. Mat_DOS%E(N)) then
                call Find_in_array_monoton(Mat_DOS%E, E, M_temp) ! find the closes value in the precalculated array of energy losses
            else
                M_temp = N+1  ! the last element is possible, all of them are
            endif
            if (M_temp .GT. 1) then
                Sum_DOS = Mat_DOS%int_DOS(M_temp-1)   ! only these electrons of VB are available
                call random_number(RN)
                Tot_N = RN*Sum_DOS   ! chose this electron
                call Find_in_array_monoton(Mat_DOS%int_DOS, Tot_N, N_temmp) ! make increasing array
                dE_cur = Mat_DOS%E(N_temmp)    ! [eV]
            else    ! if it's too close to the Ip, we can't resolve it in out grid for VB
                dE_cur = 0.0d0 ! [eV]
            endif
        endif
    endif
end subroutine From_where_in_VB


subroutine count_for_Auger_shells(Target_atoms, NRG, coun)
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects
   real(8), intent(in) :: NRG   ! [eV] energy
   integer, intent(out) :: coun ! count how many shells can participate in the Auger
   integer i, j, N, M
   N = size(Target_atoms) ! number of elements
   coun = 0
   do i = 1, N
      M = size(Target_atoms(i)%Ip)
      do j = 1, M
        !if (Target_atoms(i)%Ip(j) .LT. NRG) then
         if ((Target_atoms(i)%Ip(j) .LT. NRG) .AND. &
         (Target_atoms(i)%Ip(j) .GE. Target_atoms(1)%Ip(size(Target_atoms(1)%Ip))) ) then ! it is possible
            coun = coun + Target_atoms(i)%Nel(j)
         endif
      enddo
   enddo
end subroutine count_for_Auger_shells

subroutine Choose_for_Auger_shell(Target_atoms, NRG, Shel, Sh1, KOA1)   ! chose the second participant (hole) of the Auger-decay
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects
   real(8), intent(in) :: NRG   ! 
   integer, intent(in) :: Shel  ! the decaying shell
   integer, intent(out) :: Sh1, KOA1    ! which atom's which shell 
   integer i, j, N, M, coun_sh
   N = size(Target_atoms) ! number of elements
   coun_sh = 0
   do i = 1, N
      M = size(Target_atoms(i)%Ip) ! number of shells
      do j = 1, M
         !if (Target_atoms(i)%Ip(j) .LT. NRG) then ! it is possible
         if ((Target_atoms(i)%Ip(j) .LT. NRG) .AND. &
          (Target_atoms(i)%Ip(j) .GE. Target_atoms(1)%Ip(size(Target_atoms(1)%Ip))) ) then ! it is possible
            coun_sh = coun_sh + Target_atoms(i)%Nel(j)
            if (coun_sh .GE. Shel) then 
               Sh1 = j ! this is the shell
               KOA1 = i
            endif ! (coun_sh .EQ. Shel) then
         endif ! (Ip_elements(i,j) .LT. Ip_elements(KOA, SHL)) then ! it is possible
         if (coun_sh .GE. Shel) exit
      enddo ! j
      if (coun_sh .GE. Shel) exit
   enddo ! i
end subroutine Choose_for_Auger_shell


! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
! Dealing with energies:

subroutine Electron_recieves_E(dE, Nat_cur, Nshl_cur, Target_atoms, Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, dE_cur, Error_message)
! => dE [eV] transferred energy to electron
    real(8), intent(in) :: dE   ! given energy [eV]
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
    type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
    integer, intent(in) :: Nat_cur, Nshl_cur, Lowest_Ip_At, Lowest_Ip_Shl   ! # of atom, shell, and # of atom and shell for the VB
    real(8), intent(out) :: dE_cur  ! [eV] kinetic energy of an electron
    type(Error_handling), intent(inout) :: Error_message	! error messages are dealed with as objects
    real(8) RN, Tot_N, Sum_DOS, E
    integer N, N_temmp, M_temp
    character(100) Error_descript
    
    E = dE - Target_atoms(Nat_cur)%Ip(Nshl_cur)    ! [eV] energy that electron might recieve
    N_temmp = 1
    if ((Nat_cur .EQ. Lowest_Ip_At) .AND. (Nshl_cur .EQ. Lowest_Ip_Shl)) then   ! it's valence band:
        if ((.not. allocated(Mat_DOS%E)) .OR. (dE .LE. Target_atoms(Nat_cur)%Ip(Nshl_cur))) then  ! VB is treated as atomic energy level:
            dE_cur = E  ! [eV]
        else    ! it is a band, there is DOS, use it:
            N = size(Mat_DOS%int_DOS)   ! the last element tells the number of electrons in VB
            if (E .LT. Mat_DOS%E(N)) then
                call Find_in_array_monoton(Mat_DOS%E, E, M_temp) ! find the closes value in the precalculated array of energy losses
            else
                M_temp = N+1  ! the last element is possible, all of them are
            endif
            if (M_temp .GT. 1) then
                Sum_DOS = Mat_DOS%int_DOS(M_temp-1)   ! only these electrons of VB are available
                call random_number(RN)
                Tot_N = RN*Sum_DOS   ! chose this electron
                call Find_in_array_monoton(Mat_DOS%int_DOS, Tot_N, N_temmp) ! make increasing array
                dE_cur = E - Mat_DOS%E(N_temmp)    ! [eV]
            else    ! if it's too close to the Ip, we can't resolve it in out grid for VB
                dE_cur = E ! [eV]
            endif
        endif
    else    ! it's just a shell:
        dE_cur = E  ! [eV] energy that electron recieves
    endif
    if (dE_cur .LT. 0.0d0) then
        Error_descript = 'Transferred energy to electron is negative!'    ! description of an error
        call Save_error_details(Error_message, 10, Error_descript) ! write it into the error-log file
        print*, trim(adjustl(Error_descript)) ! print it also on the sreen
        write(*,'(a,e,e,i3,i3)') 'dE, dE_cur:', dE, dE_cur, Nat_cur, Nshl_cur
        write(*,'(a,e)') 'Ionization potential ', Target_atoms(Nat_cur)%Ip(Nshl_cur)
        write(*,'(a,e,e,i3,i3)') 'Band:', N_temmp, Mat_DOS%E(N_temmp)
        !write(*, '(a,e,e)') 'SumDOS: ', Sum_DOS, Tot_N
        !pause 'STOPPED WORKING...'
        if (N_temmp .GT. 1) then
            print*, Mat_DOS%E(N_temmp), Mat_DOS%E(N_temmp-1)
        endif
    endif
end subroutine Electron_recieves_E


subroutine SHI_energy_transfer(SHI_loc, MFP_Object, Target_atoms, Matter, Mat_DOS, Nat_cur, Nshl_cur, NumPar, dE)
    type(Ion), intent(inout) :: SHI_loc ! SHI parameters if needed
    type(All_MFP), dimension(:), intent(in) :: MFP_Object ! integrated SHI differential MFP
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
    type(Solid), intent(in) :: Matter ! material parameters
    type(Density_of_states), intent(in) :: Mat_DOS ! DOS
    integer, intent(in) :: Nat_cur, Nshl_cur    ! number of atom and number of shell
    type(Flag), intent(in) :: NumPar
    real(8), intent(out) :: dE  ! [eV] transfer energy
    real(8) RN, Tot_N, E_cur, dL, dEdx
    integer i, N_temmp, M_temp, N, M
    
    select case (Target_atoms(Nat_cur)%KOCS_SHI(Nshl_cur)) ! which inelastic cross section to use (BEB vs CDF):
    case (1) ! CDF cross section
        N = size(MFP_Object(Nat_cur)%ELMFP(Nshl_cur)%L)
        call random_number(RN)
        
        E_cur = Target_atoms(Nat_cur)%Ip(Nshl_cur)  ! exactly ionization potential [eV]
        ! find the closest value in the precalculated array of energy losses:
        call Find_in_array_monoton(MFP_Object(Nat_cur)%ELMFP(Nshl_cur)%E, Target_atoms(Nat_cur)%Ip(Nshl_cur), M_temp)
        if (M_temp .GT. 1) then
            call Interpolate(5, MFP_Object(Nat_cur)%ELMFP(Nshl_cur)%E(M_temp-1), &
                MFP_Object(Nat_cur)%ELMFP(Nshl_cur)%E(M_temp), MFP_Object(Nat_cur)%ELMFP(Nshl_cur)%L(M_temp-1), &
                MFP_Object(Nat_cur)%ELMFP(Nshl_cur)%L(M_temp), E_cur, dL)  ! interpolate to find exact value
        else
            dL = MFP_Object(Nat_cur)%ELMFP(Nshl_cur)%L(1)
        endif
        
        Tot_N = 1.0d0/dL + RN*(1.0d0/MFP_Object(Nat_cur)%ELMFP(Nshl_cur)%L(N) - 1.0d0/dL)

        if (Tot_N .LT. 1d20) then
            ! find the closest value in the precalculated array of energy losses:
            call Find_in_array(1.0d0/MFP_Object(Nat_cur)%ELMFP(Nshl_cur)%L, Tot_N, N_temmp)
        else
            N_temmp = M_temp
        endif
        if (N_temmp .GT. M_temp) then
            call Interpolate(5, 1.0d0/MFP_Object(Nat_cur)%ELMFP(Nshl_cur)%L(N_temmp-1), &
                1.0d0/MFP_Object(Nat_cur)%ELMFP(Nshl_cur)%L(N_temmp), MFP_Object(Nat_cur)%ELMFP(Nshl_cur)%E(N_temmp-1), &
                MFP_Object(Nat_cur)%ELMFP(Nshl_cur)%E(N_temmp), Tot_N, dE)  ! interpolate to find exact value
        else
            dE = Target_atoms(Nat_cur)%Ip(Nshl_cur)  ! [eV]
        endif
        
    case default ! BEB cross section:
        call SHI_TotIMFP(SHI_loc, Target_atoms, Nat_cur, Nshl_cur, dL, dEdx, Matter, Mat_DOS, NumPar) ! get total MFP; from module "Cross_sections"
        call random_number(RN)
        dL = dL/RN   ! [A] sampled MFP
        call SHI_NRG_transfer_BEB(SHI_loc, Target_atoms, Nat_cur, Nshl_cur, dL, Matter, dE) ! from module "Cross_sections"
    end select

end subroutine SHI_energy_transfer


! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
! Dealing with shells:

subroutine Which_shell(MFP_Object, MFP_full, E, Nat, Nshl)
    type(All_MFP), dimension(:), intent(in) :: MFP_Object   ! SHI or electron mean free paths for all shells
    real(8), dimension(:,:), intent(in) :: MFP_full ! full MFP of the particle
    real(8), intent(in) :: E    ! [eV] energy
    integer, intent(out) :: Nat, Nshl   ! number of atom and shell chosen
    real(8) RN, MFP, MFP_tot, MFP_sum
    integer N_temmp, i, j
    
    real(8), dimension(size(MFP_Object), size(MFP_Object(1)%ELMFP)) :: Temp_MFPs
    Temp_MFPs = 0.0d0
    
    MFP_tot = 0.0d0
    n_at:do i = 1, size(MFP_Object)  ! number of atoms
       Nat = i
       do j = 1, size(MFP_Object(i)%ELMFP) ! number of shells
          Nshl = j
          call Find_in_array_monoton(MFP_Object(i)%ELMFP(j)%E, E, N_temmp) ! find the closest value in the precalculated array of energy losses
            if (N_temmp .EQ. 1) then
                MFP = 1.0d20
            else
                if ((MFP_Object(i)%ELMFP(j)%L(N_temmp-1) .EQ. MFP_Object(i)%ELMFP(j)%L(N_temmp)) .OR. &
                 (MFP_Object(i)%ELMFP(j)%L(N_temmp-1) .GT. 1d20)) then
                    MFP = MFP_Object(i)%ELMFP(j)%L(N_temmp-1)
                else
                    call Interpolate(5, MFP_Object(i)%ELMFP(j)%E(N_temmp-1), MFP_Object(i)%ELMFP(j)%E(N_temmp), &
                      MFP_Object(i)%ELMFP(j)%L(N_temmp-1), MFP_Object(i)%ELMFP(j)%L(N_temmp), E, MFP)  ! interpolate to find exact value
                endif
            endif
            Temp_MFPs(i,j) = 1.0d0/MFP
            MFP_tot = MFP_tot + Temp_MFPs(i,j)
       enddo
   enddo n_at

    call random_number(RN)
    MFP_tot = RN*MFP_tot
    MFP_sum = 0.0d0
    do i = 1, size(MFP_Object)  ! number of atoms
        Nat = i
        do j = 1, size(MFP_Object(i)%ELMFP) ! number of shells
            Nshl = j
            MFP_sum = MFP_sum + Temp_MFPs(i,j)
            if (MFP_sum .GE. MFP_tot) exit
        enddo   ! j
        if (MFP_sum .GE. MFP_tot) exit
    enddo   ! i
end subroutine Which_shell


subroutine Next_free_path_1d(E, MFP_E_array, MFP_L_array, MFP) ! temp = MFP [A] for this array, whatever it is
    real(8), intent(in) :: E    ! [eV] energy of particle
    real(8), dimension(:), intent(in) :: MFP_E_array    ! [eV] array with mean free paths
    real(8), dimension(:), intent(in) :: MFP_L_array    ! [A] array with mean free paths
    real(8), intent(out) :: MFP    ! [A] interpolated mean free path
    integer N_temmp, N_last ! number of element in the array closest to the one we are looking for
    call Find_in_array_monoton(MFP_E_array, E, N_temmp) ! find the closes value in the precalculated array of energy losses
    if (N_temmp .EQ. 1) then
        MFP = 1.5d21    ! [A] infinity
    else
        N_last = N_temmp-1  ! the last point
        if (MFP_L_array(N_last) .GE. 1.0d16) then
            MFP = MFP_L_array(N_last)
        else
            ! interpolate to find exact value:
            call Interpolate(5, MFP_E_array(N_last), MFP_E_array(N_temmp), MFP_L_array(N_last), MFP_L_array(N_temmp), E, MFP)
        endif
    endif
end subroutine Next_free_path_1d

subroutine Next_free_path_2d(E, MFP_array, MFP) ! temp = MFP [A] for this array, whatever it is
    real(8), intent(in) :: E    ! [eV] energy of particle
    real(8), dimension(:,:), intent(in) :: MFP_array    ! [eV, A] array with mean free paths
    real(8), intent(out) :: MFP    ! [A] interpolated mean free path
    integer N_temmp, N_last ! number of element in the array closest to the one we are looking for
    call Find_in_array_monoton(MFP_array, E, 1, N_temmp) ! find the closes value in the precalculated array of energy losses
    if (N_temmp .EQ. 1) then
        MFP = 1.5d21    ! [A] infinity
    else
        N_last = N_temmp-1  ! the last point
        if (MFP_array(2,N_last) .GE. 1.0d16) then
            MFP = MFP_array(2,N_last)
        else
            ! interpolate to find exact value:
            call Interpolate(5, MFP_array(1,N_last), MFP_array(1,N_temmp), MFP_array(2,N_last), MFP_array(2,N_temmp), E, MFP)
        endif
    endif
end subroutine Next_free_path_2d


subroutine How_many_electrons(Nat, El_IMFP, El_EMFP, Hole_IMFP, Hole_EMFP, Phot_IMFP, SHI_path, SHI_loss, target_atoms, &
                Matter, Lowest_Ip_At, Lowest_Ip_Shl, SHI, SHI_MFP, Total_el_MFPs, Elastic_MFP, Total_Hole_MFPs, &
                Elastic_Hole_MFP, Total_Photon_MFPs, All_electrons, Em_electrons, All_holes, All_photons, Nel, NumPar)
    integer, intent(in) :: Nat, Lowest_Ip_At, Lowest_Ip_Shl  ! number of atoms, and VB numbers
    real(8), dimension(:,:), allocatable, intent(inout) :: El_IMFP     ! to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), allocatable, intent(inout) :: El_EMFP     ! to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), allocatable, intent(inout) :: Hole_IMFP   ! to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), allocatable, intent(inout) :: Hole_EMFP   ! to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), allocatable, intent(inout) :: Phot_IMFP   ! to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), allocatable, intent(inout) :: SHI_path    ! to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), allocatable, intent(inout) :: SHI_loss    ! to use in a subroutine, we need this shape of an array
    type(Ion), intent(in) :: SHI   ! declare SHI as an object with atributes "Ion"
    type(All_MFP), dimension(:), allocatable, intent(in) :: SHI_MFP         ! SHI mean free paths for all shells
    type(Solid), intent(in) :: Matter   ! all material parameters
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
    type(All_MFP), dimension(:), intent(in) :: Total_el_MFPs        ! electron mean free paths for all shells
    type(MFP), intent(in) :: Elastic_MFP                    ! elastic mean free path
    type(All_MFP), dimension(:), intent(in) :: Total_Hole_MFPs      ! hole mean free paths for all shells
    type(MFP), intent(in) :: Elastic_Hole_MFP               ! hole elastic mean free path
    type(All_MFP), dimension(:), intent(in) :: Total_Photon_MFPs    ! photon mean free paths for all shells
    type(Electron), dimension(:), allocatable, intent(inout) :: All_electrons ! define array of electrons
    real(8), dimension(:), allocatable, intent(inout) :: Em_electrons ! emitted electrons energies
    type(Hole), dimension(:), allocatable, intent(inout) :: All_holes ! define array of holes
    type(Photon), dimension(:), allocatable, intent(inout) :: All_photons ! define array of photons
    integer, intent(out) :: Nel ! that's the size of Electron and Hole arrays
    type(Flag), intent(in) :: NumPar ! numerical parameters
    
    real(8) dEdx, SHI_dEdx
    integer N_temmp,i, j, Nshl
    
    !sssssssssssssssssssssssssssssssssss
    ! SHIs arrays:
    if(.not. allocated(SHI_path)) allocate(SHI_path(2,size(SHI_MFP(1)%ELMFP(1)%E)))   ! SHI Mean Free Path [eV, A]
    if(.not. allocated(SHI_loss)) allocate(SHI_loss(2,size(SHI_MFP(1)%ELMFP(1)%E)))   ! SHI losses [eV, eV/A]
    SHI_loss = 0.0d0
    SHI_path = 0.0d0
    SHI_loss(1,:) = SHI_MFP(1)%ELMFP(1)%E(:)    ! to use later we need an array of this shape
    SHI_path(1,:) = SHI_MFP(1)%ELMFP(1)%E(:)    ! to use later we need an array of this shape
    do i = 1, Nat
       Nshl = size(target_atoms(i)%Ip)  ! number of shells
       do j = 1, Nshl
           SHI_loss(2,:) = SHI_loss(2,:) + SHI_MFP(i)%ELMFP(j)%dEdx(:) ! to use later we need an array of this shape
           SHI_path(2,:) = SHI_path(2,:) + 1.0d0/SHI_MFP(i)%ELMFP(j)%L(:) ! to use later we need an array of this shape
       enddo
    enddo
    !SHI_path(2,:) = 1.0d0/SHI_path(2,:) ! total MFP [A]
    where(SHI_path(2,:) < 1.0d-10) ! too slow electron/hole
       SHI_path(2,:) = 1.0d30
    elsewhere   ! regular energy particle
       SHI_path(2,:) = 1.0d0/SHI_path(2,:) ! total MFP [A]
    endwhere

    call Find_in_array_monoton(SHI_loss, SHI%E, 1, N_temmp) ! find the closes value in the precalculated array of energy losses
    ! interpolate to find exact value:
    call Interpolate(5, SHI_loss(1,N_temmp-1), SHI_loss(1,N_temmp), SHI_loss(2,N_temmp-1), SHI_loss(2,N_temmp), SHI%E, SHI_dEdx)
    dEdx = SHI_dEdx*Matter%Layer ![eV] total energy loss by SHI within the given layer
    
    
    Nel = CEILING(dEdx/(Target_atoms(Lowest_Ip_At)%Ip(Lowest_Ip_Shl))) ! mean number of electrons estimated by the bandgap
    if (Nel .GT. 5000*Matter%Layer) Nel = 5000*Matter%Layer
    if (Nel .LT. 1000) Nel = 1000
    
    if (.not. allocated(All_electrons)) allocate(All_electrons(INT(Nel))) ! define arrays of electrons
    if (.not. allocated(Em_electrons)) allocate(Em_electrons(INT(Nel)))   ! define array of emitted electrons
    if (.not. allocated(All_holes)) allocate(All_holes(INT(Nel))) ! define arrays of holes
    if (.not. allocated(Hole_EMFP)) allocate(Hole_EMFP(2,size(Elastic_Hole_MFP%E)))
    if (.not. allocated(Hole_IMFP)) allocate(Hole_IMFP(2,size(Total_Hole_MFPs(1)%ELMFP(1)%E)))

    !eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
    ! Electrons (and some holes) arrays:
    do i = 1, size(All_electrons)
       ! set initial electron data:
       call Particle_event(All_electrons(i), E=0.0d0, t0=SHI%t0, tn=1d20, X=0.0d0, Y=0.0d0, Z=0.0d0, theta=0.0d0, phi=0.0d0)
       call Particle_event(All_holes(i), E=0.0d0, Ehkin=0.0d0, t0=SHI%t0, tn=1d21, X=0.0d0, Y=0.0d0, Z=0.0d0, L=1.0d30, &
                KOA=0, Shl=0, Mass=1.0d30, theta=0.0d0, phi=0.0d0) ! set initial hole data
    enddo
    
    Em_electrons = 0.0d0
    
    if (.not. allocated(El_EMFP)) allocate(El_EMFP(2,size(Elastic_MFP%E)))
    if (.not. allocated(El_IMFP)) allocate(El_IMFP(2,size(Total_el_MFPs(1)%ELMFP(1)%E)))
    El_EMFP = 0.0d0
    El_IMFP = 0.0d0
    El_EMFP(1,:) = Elastic_MFP%E(:)    ! to use later we need an array of this shape
    El_EMFP(2,:) = Elastic_MFP%L(:)    ! to use later we need an array of this shape
    
    El_IMFP(1,:) = Total_el_MFPs(1)%ELMFP(1)%E(:)    ! to use later we need an array of this shape
    do i = 1, Nat
       Nshl = size(target_atoms(i)%Ip)  ! number of shells
       do j = 1, Nshl
           El_IMFP(2,:) = El_IMFP(2,:) + 1.0d0/Total_el_MFPs(i)%ELMFP(j)%L(:) ! to use later we need an array of this shape
       enddo
    enddo
    !El_IMFP(2,:) = 1.0d0/El_IMFP(2,:)   ![A]
    where(El_IMFP(2,:) < 1.0d-10) ! too short
       El_IMFP(2,:) = 1.0d30
    elsewhere   ! regular energy particle
       El_IMFP(2,:) = 1.0d0/El_IMFP(2,:)   ![A]
    endwhere
    
    !hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
    ! Holes arrays:
    Hole_EMFP = 0.0d0
    Hole_IMFP = 0.0d0
    Hole_EMFP(1,:) = Elastic_Hole_MFP%E(:)    ! to use later we need an array of this shape
    Hole_EMFP(2,:) = Elastic_Hole_MFP%L(:)    ! to use later we need an array of this shape
    Hole_IMFP(1,:) = Total_Hole_MFPs(1)%ELMFP(1)%E(:)    ! to use later we need an array of this shape
    do i = 1, Nat
       Nshl = size(target_atoms(i)%Ip)  ! number of shells
       do j = 1, Nshl
           Hole_IMFP(2,:) = Hole_IMFP(2,:) + 1.0d0/Total_Hole_MFPs(i)%ELMFP(j)%L(:) ! to use later we need an array of this shape
       enddo
    enddo
    !Hole_IMFP(2,:) = 1.0d0/Hole_IMFP(2,:)   ![A]
    where(Hole_IMFP(2,:) < 1.0d-10) ! too short
       Hole_IMFP(2,:) = 1.0d30
    elsewhere   ! regular energy particle
       Hole_IMFP(2,:) = 1.0d0/Hole_IMFP(2,:)   ![A]
    endwhere
    
    !hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
    ! Photons arrays:
    if (NumPar%include_photons) then ! only if we include photons:
       if (.not. allocated(All_photons)) allocate(All_photons(INT(Nel))) ! define array of photons
       do i = 1, size(All_photons)
          ! set initial photon data:
          call Particle_event(All_photons(i), E=0.0d0, t0=SHI%t0, tn=1d20, X=0.0d0, Y=0.0d0, Z=0.0d0, theta=0.0d0, phi=0.0d0)
       enddo
       if (.not. allocated(Phot_IMFP)) allocate(Phot_IMFP(2,size(Total_Photon_MFPs(1)%ELMFP(1)%E)))
       Phot_IMFP = 0.0d0
       Phot_IMFP(1,:) = Total_Photon_MFPs(1)%ELMFP(1)%E(:)    ! to use later we need an array of this shape
       do i = 1, Nat
          Nshl = size(target_atoms(i)%Ip)  ! number of shells
          do j = 1, Nshl
             Phot_IMFP(2,:) = Phot_IMFP(2,:) + 1.0d0/Total_Photon_MFPs(i)%ELMFP(j)%L(:) ! to use later we need an array of this shape
          enddo
       enddo
       !Phot_IMFP(2,:) = 1.0d0/Phot_IMFP(2,:) ! [A]
       where(Phot_IMFP(2,:) < 1.0d-10) ! too short
          Phot_IMFP(2,:) = 1.0d30
       elsewhere   ! regular energy particle
          Phot_IMFP(2,:) = 1.0d0/Phot_IMFP(2,:) ! [A]
       endwhere
    endif
end subroutine How_many_electrons


subroutine Find_min_time_particle(SHI, Electrons, Holes, Photons, KOP, NOP, t_cur, NumPar)
   class(Ion), intent(in) :: SHI
   class(Electron), DIMENSION(:), intent(in) :: Electrons
   class(Hole), DIMENSION(:), intent(in) :: Holes
   type(Photon), dimension(:), intent(in) :: Photons
   integer, intent(out) :: KOP	! kind of particle
   integer, intent(out) :: NOP	! number of particle of this kind
   real(8), intent(out) :: t_cur	! time of collision of this particle
   type(Flag), intent(in) :: NumPar ! numerical parameters
   integer i
   
   if (NumPar%include_photons) then ! only if we include photons:
    ! kind of particle determines the next select case: order of arrays must coinside with select-case!!!:
    KOP = transfer( minloc((/SHI%tn,minval(Electrons%tn),minval(Holes%tn),minval(Photons%tn)/)) , i)
   else ! no photons:
    ! kind of particle determines the next select case: order of arrays must coinside with select-case!!!:
    KOP = transfer( minloc((/SHI%tn,minval(Electrons%tn),minval(Holes%tn)/)) , i)
   endif
   select case (KOP)
    case (1) ! SHI
        NOP =  1        ! we have only 1 SHI
        t_cur = SHI%tn  ! [fs] next collision of SHI
    case (2) ! electron
        NOP =  transfer(minloc(Electrons%tn), i)	 ! which electron from array
        t_cur = minval(Electrons%tn)   ! [fs] time of electron next collision
    case (3) ! hole
        NOP =  transfer(minloc(Holes%tn), i) ! which hole from array
        t_cur = minval(Holes%tn)    ! [fs] time hole decay
    case (4) ! photon
        NOP =  transfer(minloc(Photons%tn), i)	 ! which photon from array
        t_cur = minval(Photons%tn)   ! [fs] time of photon absorbtion
   endselect
end subroutine Find_min_time_particle

subroutine barrier_parameters (Matter, Em_E1, Em_gamma)
    real(8) EM_B, Em_L, Em_W1, Em_bb, Em_delta, Em_gam1, Em_gam2, Em_ksi
    type(Solid), intent(in) :: Matter
    real(8) :: work_function, bar_height
    real(8), intent(out) :: Em_gamma, Em_E1
    
    work_function = Matter%work_function
    bar_height = Matter%bar_height
    
    Em_L = Matter%bar_length*1.0d-10     !Barrier length     [m]
    Em_B = 2.0d0*bar_height - work_function + 2.0d0*sqrt(bar_height*bar_height - bar_height*work_function)
    
    Em_ksi = 0.5d0*sqrt(8.0d0*g_me*Em_L*Em_L*Em_B*g_e/((2.0d0*g_Pi*g_h)*(2.0d0*g_Pi*g_h))-1.0d0)
    Em_bb = cosh(2.0d0*g_pi*Em_ksi)
    Em_delta = 2.0d0*g_pi*Em_L*sqrt(2.0d0*g_me*g_e)/(2.0d0*g_pi*g_h)
    
    Em_E1 = bar_height + 2.0d0*sqrt(bar_height*(bar_height - work_function))*(acosh(Em_bb)/(Em_delta*(sqrt(bar_height) &
                + sqrt(bar_height - work_function)))-1.0d0)
    Em_gam1 = Em_delta*(sqrt(Em_E1) + sqrt(Em_E1 - work_function))
    Em_gam2 = Em_delta*(sqrt(Em_E1) - sqrt(Em_E1 - work_function))
    Em_gamma = (Em_gam1*sinh(Em_gam1) + 2.0d0*Em_gam2*sinh(Em_gam2))/(sqrt(Em_E1*(Em_E1 - work_function))*(Em_bb + cosh(Em_gam1)))
end subroutine barrier_parameters


subroutine set_time_grid (Tim, dt, dt_flag, time_grid)
    real(8), intent(in) :: Tim ! [fs] total time
    real(8), intent(in) :: dt  ! [fs] timestep
    integer, intent(in) :: dt_flag  ! kind of time-grid: 0=linear, 1=logarithmic
    real(8), dimension(:), allocatable, intent(inout) :: time_grid ! grid_points in time
    real(8) tim_glob
    integer i

    select case (dt_flag)   ! what kind of time-grid to use:
        case (:0)   ! linear time-grid
            allocate(time_grid(CEILING(Tim/dt)+1))
            time_grid = 0.0d0
            time_grid(1) = dt
            do i = 1, size(time_grid)-2
                time_grid(i+1) = time_grid(i) + dt
            enddo
            time_grid(size(time_grid)) = Tim + dt
        case (1:)   ! logarithmic time-grid
            i = 0
            tim_glob = 0.01d0
            do while (tim_glob .LE. Tim)
                i = i + 1
                tim_glob = tim_glob*dt	! [fs]
            enddo
            allocate(time_grid(i+1))
            time_grid(1) = 0.01d0   ! [fs] first point
            do i = 1, size(time_grid)-2
                time_grid(i+1) = time_grid(i)*dt	! [fs]
            enddo
            time_grid(size(time_grid)) = Tim + dt
    endselect
end subroutine set_time_grid 


! Monte-Carlo of the SHI passage
subroutine SHI_Monte_Carlo(SHI_MFP, SHI_path, SHI_loc, diff_SHI_MFP, Target_atoms, All_electrons, All_holes, Tot_Nel, &
                Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, Error_message, El_IMFP, El_EMFP, Hole_IMFP, Hole_EMFP, Matter, &
                NumPar)
    type(All_MFP), dimension(:), intent(in) :: SHI_MFP         ! SHI mean free paths for all shells
    type(All_MFP), dimension(:), intent(in) :: diff_SHI_MFP    ! SHI differential mean free paths for all shells
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
    integer, intent(in) :: Lowest_Ip_At, Lowest_Ip_Shl ! number of atom and of shell which correspond to the lowest ionization potential
    type(Density_of_states), intent(in) :: Mat_DOS  ! material DOS
    type(Error_handling), intent(inout) :: Error_message	! error messages are dealed with as objects
    type(Ion), intent(inout) :: SHI_loc ! define SHI for local use
    type(Electron), dimension(:), allocatable, intent(inout) :: All_electrons ! define array of electrons
    type(Hole), dimension(:), allocatable, intent(inout) :: All_holes ! define array of holes
    real(8), dimension(:,:), intent (inout) :: SHI_path    ! to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), intent(in) :: El_IMFP    ! total IMFP for electron, to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), intent(in) :: El_EMFP    ! to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), intent(in) :: Hole_IMFP    ! total IMFP for electron, to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), intent(in) :: Hole_EMFP    ! to use in a subroutine, we need this shape of an array
    type(Solid), intent(in) :: Matter   ! all material parameters
    integer, intent(inout) :: Tot_Nel
    type(Flag), intent(in) :: NumPar
    
    real(8) dE, X, Y, Z, L, dE_cur, mh
    real(8) phi, theta, hphi, htheta
    integer Nat_cur, Nshl_cur, i
    real(8) IMFP, EMFP, SHI_IMFP, MFP_tot
    real(8) RN
    character(200) :: Error_descript
        
    ! Find which shell of which atom is being ionized:
    call Which_shell(SHI_MFP, SHI_path, SHI_loc%E, Nat_cur, Nshl_cur)   ! => Nat_cur, Nshl_cur
    ! Find transferred energy:
    call SHI_energy_transfer(SHI_loc, diff_SHI_MFP, Target_atoms, Matter, Mat_DOS, Nat_cur, Nshl_cur, NumPar, dE)   ! => dE [eV]
    
    ! Next SHI collision parameters:
    call Next_free_path(SHI_loc%E, SHI_path, SHI_IMFP) ! => SHI_MFP [A]
    call random_number(RN)
    SHI_IMFP = -SHI_IMFP*log(RN)    ! sample next SHI free path
    ! SHI loses this amount of energy, change the current time, and Z-coordinate:
    Z = SHI_loc%Z+SHI_loc%L ! [A] new Z coordinate
    call Particle_event(SHI_loc, E=SHI_loc%E-dE, t0=SHI_loc%tn, Z=Z, L=SHI_IMFP) ! these SHI's parameters are updated
    call Get_time_of_next_event(SHI_loc, MFP=SHI_IMFP)  ! => SHI_loc%tn is updated [fs]
    call Equilibrium_charge_SHI(SHI_loc, Target_atoms)  ! new Barcas' equilibrium charge for the new energy
    
    Tot_Nel = Tot_Nel + 1   ! we have ionized an electron!
    call Check_size(All_electrons, All_holes, N=Tot_Nel)    ! check if electrons are too many and the size of arrays must be increased
    
    !eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
    
    ! How much energy an electron recieves ( ! => dE_cur [eV] kinetic energy of the electron):
    call Electron_recieves_E(dE, Nat_cur, Nshl_cur, Target_atoms, Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, dE_cur, Error_message)
    
    ! Find the correct angles of electron emmition/scattering:
    call Update_electron_angles(SHI_loc, dE, theta, phi)    ! => theta, phi

    
    ! Get electron mean free path:
    call Next_free_path(dE_cur, El_IMFP, IMFP) ! => IMFP of electron [A]
    call Next_free_path(dE_cur, El_EMFP, EMFP) ! => EMFP of electron [A]
    call random_number(RN)
    MFP_tot = -log(RN)/(1.0d0/IMFP + 1.0d0/EMFP)    ! [A] sample total electron free path (inelastic + elastic)
    
    ! Give initial parameters to this electron:
    L = Impact_parameter(SHI_loc, dE)   ! [A] SHI impact parameter
    X = SHI_loc%X+L*sin(phi)    ! [A] new X coordinate
    Y = SHI_loc%Y+L*cos(phi)    ! [A] new Y coordinate
    ! new electron parameters:
    call Particle_event(All_electrons(Tot_Nel), E=dE_cur, t0=SHI_loc%t0, X=X, Y=Y, Z=Z, L=MFP_tot, theta=theta, phi=phi)
    call Get_time_of_next_event(All_electrons(Tot_Nel), MFP=MFP_tot)   ! => All_electrons(Tot_Nel)%tn is updated for next collision [fs]
    
    call cut_off(Matter%cut_off, All_electrons=All_electrons(Tot_Nel)) ! compare electron energy with cut-off and update it's time
    
    if ((All_electrons(Tot_Nel)%E .LT. -1.0d-9) .OR. isnan(All_electrons(Tot_Nel)%E)) then  ! Error, electron got negative energy!
        ! description of an error:
        write(Error_descript, '(a,i,a,e)') 'SHI electron #', Tot_Nel, ' got negative energy ', All_electrons(Tot_Nel)%E
        call Save_error_details(Error_message, 20, Error_descript) ! write it into the error-log file
        print*, trim(adjustl(Error_descript)) ! print it also on the sreen
    endif
    !hhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
    ! We have also produced a hole:

    call Update_holes_angles_SHI((dE-dE_cur), htheta, hphi)         !Random angles of hole
    call Particle_event(All_holes(Tot_Nel), t0=SHI_loc%t0, X=X, Y=Y, Z=Z, KOA=Nat_cur, Shl=Nshl_cur, theta=htheta, phi=hphi)
    call Hole_parameters(All_holes(Tot_Nel), Matter, Mat_DOS, target_atoms, Hole_IMFP, Hole_EMFP, (dE-dE_cur), &
            Lowest_Ip_At, Lowest_Ip_Shl)
    call cut_off(Matter%cut_off, All_holes=All_holes(Tot_Nel)) ! compare hole energy with cut-off and update it's time
    
    if ((All_holes(Tot_Nel)%Ehkin .LT. -1.0d-9) .OR. isnan(All_holes(Tot_Nel)%Ehkin)) then  ! Error, hole got negative energy!
        ! description of an error:
        write(Error_descript, '(a,i,a,e)') 'SHI hole #', Tot_Nel, ' got negative energy ', All_holes(Tot_Nel)%Ehkin
        call Save_error_details(Error_message, 20, Error_descript) ! write it into the error-log file
        print*, trim(adjustl(Error_descript)) ! print it also on the sreen
    endif

end subroutine SHI_Monte_Carlo


! Monte-carlo of an electron
subroutine Electron_Monte_Carlo(All_electrons, All_holes, El_IMFP, El_EMFP, Hole_IMFP, Hole_EMFP, CDF_Phonon, Matter, target_atoms, &
            Total_el_MFPs, Elastic_MFP, Tot_Nel, NOP, Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, Error_message, Em_electrons, &
            Em_Nel, Em_gamma, Em_E1, At_NRG, Out_R, Out_Elat, Out_V, i, DSF_DEMFP, NumPar)
    integer, intent(in) :: Lowest_Ip_At, Lowest_Ip_Shl ! number of atom and of shell which correspond to the lowest ionization potential
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
    type(CDF), intent(in) :: CDF_Phonon ! declare CDF for phonons
    type(Solid), intent(in) :: Matter   ! all material parameters
    type(Density_of_states), intent(in) :: Mat_DOS  ! material DOS
    type(Error_handling), intent(inout) :: Error_message	! error messages are dealed with as objects
    type(All_MFP), dimension(:), intent(in), target :: Total_el_MFPs    ! electron mean free paths for all shells
    type(MFP), intent(in), target :: Elastic_MFP                        ! elastic mean free path
    type(Electron), dimension(:), allocatable :: All_electrons ! define array of electrons
    type(Hole), dimension(:), intent(inout), allocatable :: All_holes ! define array of holes
    real(8), dimension(:), intent(inout) :: Em_electrons                     ! define array of energies of emitted electrons
    real(8), dimension(:,:), intent(in) :: El_IMFP    ! total IMFP for electron, to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), intent(in) :: El_EMFP    ! to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), intent(in) :: Hole_IMFP    ! total IMFP for electron, to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), intent(in) :: Hole_EMFP    ! to use in a subroutine, we need this shape of an array
    real(8), intent(in) :: Em_gamma, Em_E1
    integer, intent(inout) :: Tot_Nel, Em_Nel
    integer, intent(in) :: NOP, i
    real(8), intent(inout) :: At_NRG
    real(8), dimension(:), intent(inout) :: Out_R   ! [A] radius for distributions
    real(8), dimension(:), intent(inout) :: Out_V  ! inverse volume of cilinder layers [1/A^3]
    real(8), dimension(:,:), intent(inout) :: Out_Elat  ! [eV/A^3] lattuce energy density vs time vs R
    type(Differential_MFP), dimension(:), intent(in) :: DSF_DEMFP
    type(Flag), intent(inout) :: NumPar
        
    real(8) Eel, RN, dE, dE_cur, mh, R, Em_Penetr, MFP_tot
    real(8) theta0, phi0, theta2, phi2, phi1, theta1, theta, phi, htheta, hphi
    real(8) IMFP, EMFP, L, X, Y, Z, dE_loc
    integer Nat_cur, Nshl_cur, j, ii
    character(200) :: Error_descript
    character(8) kind_of_particle
    
    kind_of_particle  = 'Electron'
    Eel = All_electrons(NOP)%E ! to use as shorter name; [eV] electron energy
    
    ! Find is it elastic or inelastic collision:
    call Next_free_path(Eel, El_IMFP, IMFP) ! => IMFP of electron [A]
    call Next_free_path(Eel, El_EMFP, EMFP) ! => EMFP of electron [A]
    
    call random_number(RN)
    L = All_electrons(NOP)%L            ! [A] old mean free path
    theta0 = All_electrons(NOP)%theta   ! old angle
    phi0 = All_electrons(NOP)%phi       ! old angle
    X = All_electrons(NOP)%X + L*sin(theta0)*sin(phi0)    ! [A] new X coordinate
    Y = All_electrons(NOP)%Y + L*sin(theta0)*cos(phi0)    ! [A] new Y coordinate
    Z = All_electrons(NOP)%Z + L*cos(theta0)             ! [A] new Z coordinate
    
    el_vs_inel:if (RN*(1.0d0/IMFP + 1.0d0/EMFP) .LT. 1.0d0/IMFP) then  ! inelastic
        ! Find which shell of which atom is being ionized:
        call Which_shell(Total_el_MFPs, El_IMFP, Eel, Nat_cur, Nshl_cur)   ! => Nat_cur, Nshl_cur
        
        Tot_Nel = Tot_Nel + 1   ! we have ionized a new electron!
        call Check_size(All_electrons, All_holes, N=Tot_Nel)    ! check if electrons are too many and the size of arrays must be increased
        
        ! => find IMFP of electron [A] needed for calculation of transferred energy:
        call Next_free_path(Eel, Total_el_MFPs(Nat_cur)%ELMFP(Nshl_cur)%E, Total_el_MFPs(Nat_cur)%ELMFP(Nshl_cur)%L, IMFP)
        ! => dE [eV] transferred energy:
        call Electron_energy_transfer(Eel, Target_atoms, Nat_cur, Nshl_cur, IMFP, dE, Matter, Mat_DOS, NumPar, kind_of_particle)
        call Update_electron_angles(All_electrons(NOP)%E, dE, theta, phi)    ! => theta, phi
        
        ! 222222222222222222222222222222222
        ! New ionized electrons parameters:
        ! How much energy an electron recieves (! => dE_cur [eV] kinetic energy of the electron):
        call Electron_recieves_E(dE, Nat_cur, Nshl_cur, Target_atoms, Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, dE_cur, Error_message)
        
        ! Its parameters:
        call Next_free_path(dE_cur, El_IMFP, IMFP) ! => IMFP of electron [A]
        call Next_free_path(dE_cur, El_EMFP, EMFP) ! => EMFP of electron [A]
        call random_number(RN)
        MFP_tot = -log(RN)/(1.0d0/IMFP + 1.0d0/EMFP)    ! [A] sample total electron free path (inelastic + elastic)
        theta2 = g_Pi/2.0d0 - theta ! scattering of particles with equal masses
        phi2 = phi + g_Pi
        call New_Angles_both(phi0, theta0, theta2, phi2, phi1, theta1)    ! => phi1, theta1
        ! new electron parameters:
        call Particle_event(All_electrons(Tot_Nel), E=dE_cur, t0=All_electrons(NOP)%tn, X=X, Y=Y, Z=Z, L=MFP_tot, theta=theta1, phi=phi1)
        call Get_time_of_next_event(All_electrons(Tot_Nel), MFP=MFP_tot)   ! => All_electrons(NOP)%tn is updated for next collision [fs]
        
        call cut_off(Matter%cut_off, All_electrons=All_electrons(Tot_Nel)) ! compare electron energy with cut-off        
        
        if ((All_electrons(Tot_Nel)%E .LT. -1.0d-9) .OR. isnan(All_electrons(Tot_Nel)%E)) then  ! Error, electron got negative energy!
            ! description of an error:
            write(Error_descript, '(a,i,a,e)') 'Impact electron #', Tot_Nel, ' got negative energy ', All_electrons(Tot_Nel)%E
            print*, 'Incident electron parameters: ', NOP, All_electrons(NOP)%E
            print*, kind_of_particle
            print*, All_electrons(NOP)%L, All_electrons(NOP)%tn, All_electrons(NOP)%t0
            print*, All_electrons(NOP)%theta, All_electrons(NOP)%phi 
            
            call Save_error_details(Error_message, 21, Error_descript) ! write it into the error-log file
            print*, trim(adjustl(Error_descript)) ! print it also on the sreen
        endif
        ! hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
        ! Get parameters of the created hole:
        
        call Update_holes_angles_SHI((dE-dE_cur), htheta, hphi)     !Random angles of hole
        call Particle_event(All_holes(Tot_Nel), t0=All_electrons(NOP)%t0, X=X, Y=Y, Z=Z, KOA=Nat_cur, &
                Shl=Nshl_cur, theta=htheta, phi=hphi) ! new hole parameters
        call Hole_parameters(All_holes(Tot_Nel), Matter, Mat_DOS, Target_atoms, Hole_IMFP, Hole_EMFP, &
                (dE-dE_cur), Lowest_Ip_At, Lowest_Ip_Shl)
        call cut_off(Matter%cut_off, All_holes=All_holes(Tot_Nel)) ! compare hole energy with cut-off and update it's time

        if ((All_holes(Tot_Nel)%Ehkin .LT. -1.0d-9) .OR. isnan(All_holes(Tot_Nel)%Ehkin)) then  ! Error, electron got negative energy!
            ! description of an error:
            write(Error_descript, '(a,i,a,e)') 'El-el impact hole #', Tot_Nel, ' got negative energy ', All_holes(Tot_Nel)%Ehkin
            call Save_error_details(Error_message, 20, Error_descript) ! write it into the error-log file
            print*, trim(adjustl(Error_descript)) ! print it also on the sreen
        endif
        
    else el_vs_inel ! elastic

        ! => find IMFP of electron [A] needed for calculation of transferred energy:
        call Next_free_path(Eel, Elastic_MFP%E, Elastic_MFP%L, EMFP)
        
        if (NumPar%kind_of_EMFP .EQ. 2) then           ! DSF cross-sections
            call NRG_transfer_elastic_DSF(DSF_DEMFP, Eel, EMFP, dE) ! module "Cross_sections"
        else if (NumPar%kind_of_EMFP .EQ. 1) then      ! CDF phonon peaks
            ! => dE [eV] transferred energy:
            call Electron_energy_transfer(Eel, EMFP, Target_atoms, CDF_Phonon, Matter, dE, NumPar, Mat_DOS, kind_of_particle)
        else                                    ! Atomic cross-sections of Mott 
            dE = 0.0d0                          ! [eV] transferred energy
            do ii = 1, size(Target_atoms)        ! for all atomic spicies:
               call NRG_transfer_elastic_atomic(Target_atoms, ii, Eel, dE_loc)
               !dE = dE + dE_loc                 ! [eV] sum of transferred energies during collision witn all kind of atoms
               dE = dE + dE_loc*Target_atoms(ii)%Pers ! [eV] sum of transferred energies during collision witn all kind of atoms
            enddo
            !dE = dE/size(Target_atoms)          ! [eV] Average dE over kinds of atoms
            dE = dE/dble(SUM(target_atoms(:)%Pers))          ! [eV} Average dE over kinds of atoms
        endif
        
        call Update_particle_angles_lat(target_atoms, Eel, dE, theta, phi)
	    !call Update_electron_angles_el(Eel, dE, theta, phi)
	
        At_NRG = At_NRG + dE    ! [eV] lattice energy
!         if (dE >= 0.0d0) then   ! lattice heating
! 	        lat_inc =  lat_inc + dE
! 	    else   ! latice cooling
! 	        lat_decr = lat_decr + dE
! 	    endif
        
        if (isnan(theta) .OR. isnan(phi)) then
            print*, 'ERROR in Electron_Monte_Carlo:'
            print*, 'Elastic', Eel, dE, theta, phi
            pause
        endif

        ! save this energy in the array of radial distributions of atomic energies:
        R = SQRT(X*X + Y*Y) ! [A] radius of this electron at the moment of scattering
        call Find_in_array_monoton(Out_R, R, j) ! find where in the distribution array it is
        Out_Elat(i,j) = Out_Elat(i,j) + dE*Out_V(j)   ! [eV/A^3] here is the scattering event happend
    endif el_vs_inel ! "named if"-statement
    ! 111111111111111111111111111111111
    ! Incident electron new parameters:
    call Next_free_path(Eel-dE, El_IMFP, IMFP) ! => IMFP of electron [A]
    call Next_free_path(Eel-dE, El_EMFP, EMFP) ! => EMFP of electron [A]
    call random_number(RN)
    MFP_tot = -log(RN)/(1.0d0/IMFP + 1.0d0/EMFP)    ! [A] sample total electron free path (inelastic + elastic)
    call New_Angles_both(phi0, theta0, theta, phi, phi1, theta1)    ! => phi1, theta1
    ! new electron parameters:
    call Particle_event(All_electrons(NOP), E=Eel-dE, t0=All_electrons(NOP)%tn, X=X, Y=Y, Z=Z, L=MFP_tot, theta=theta1, phi=phi1)
    call Get_time_of_next_event(All_electrons(NOP), MFP=MFP_tot)   ! => All_electrons(NOP)%tn is updated for next collision [fs]
    
    call cut_off(Matter%cut_off, All_electrons=All_electrons(NOP)) ! compare electron energy with cut-off    

    if (Matter%work_function .GT. 0) then      !Electron emission
        call calculate_emission (All_electrons, Em_electrons, NOP, Matter, Em_Nel, Em_gamma, Em_E1)                          
    endif
    
    if ((All_electrons(NOP)%E .LT. -1.0d-9) .OR. isnan(All_electrons(NOP)%E)) then  ! Error, electron got negative energy!
        ! description of an error:
        write(Error_descript, '(a,i,a,e)') 'Incident electron #', NOP, ' got negative energy ', All_electrons(NOP)%E
        call Save_error_details(Error_message, 22, Error_descript) ! write it into the error-log file
        print*, trim(adjustl(Error_descript)) ! print it also on the sreen
    endif
 end subroutine Electron_Monte_Carlo


subroutine calculate_emission (All_electrons, Em_electrons, NOP, Matter, Em_Nel, Em_gamma, Em_E1)
    type(Electron), dimension(:), intent(inout) :: All_electrons ! define array of electrons
    real(8), dimension(:), intent(inout) :: Em_electrons
    integer, intent(inout) :: Em_Nel
    integer, intent(in) :: NOP
    real(8), intent(in) :: Em_gamma, Em_E1
    type(Solid), intent(in) :: Matter   ! all material parameters
    
    real(8) Em_Penetr, RN
    
    out_z:if (All_electrons(NOP)%Z .LT. 0.0d0) then
        high_E:if (All_electrons(NOP)%E .GE. 1.5d0*Matter%bar_height) then        !Energy of an electron is greater than the barrier
            Em_Nel = Em_Nel + 1
            All_electrons(NOP)%tn = 1d30         ! to exclude further interaction of emitted electrons
            All_electrons(NOP)%L = 1d30
            Em_electrons(Em_Nel) = All_electrons(NOP)%E - Matter%work_function
        else high_E                                                       !Energy of an electron is less than the barrier
            call random_number(RN)                            
            Em_Penetr = 1.0d0/(1.0d0+exp(Em_gamma*(Em_E1 - All_electrons(NOP)%E)))      ! Calculation of penetrability
            if (RN .LT. Em_Penetr) then
                Em_Nel = Em_Nel + 1
                All_electrons(NOP)%tn = 1.0d30         ! to exclude further interaction of emitted electrons
                All_electrons(NOP)%L = 1.0d30
                Em_electrons(Em_Nel) = All_electrons(NOP)%E - Matter%work_function
            else
                if (cos(All_electrons(NOP)%theta) .LT. 0) then
                    All_electrons(NOP)%theta = g_pi - All_electrons(NOP)%theta      !Electron is reflected
                endif  
            endif
        endif high_E
    endif out_z
endsubroutine calculate_emission


! Monte-Carlo of a hole
subroutine Hole_Monte_Carlo(All_electrons, All_holes, All_photons, El_IMFP, El_EMFP, &
            Hole_IMFP, Hole_EMFP, Phot_IMFP, CDF_Phonon, Matter, target_atoms, &
            Total_Hole_MFPs, Elastic_Hole_MFP, Tot_Nel, Tot_Nphot, NOP, Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, Error_message, &
            At_NRG, Out_R, Out_Elat, Out_V, i, t_cur, DSF_DEMFP_H, NumPar)
    type(Electron), dimension(:), intent(inout), allocatable :: All_electrons   ! define array of electrons
    type(Hole), dimension(:), intent(inout), allocatable :: All_holes           ! define array of holes
    type(Photon), dimension(:), intent(inout), allocatable :: All_photons       ! define array of photons
    integer, intent(in) :: Lowest_Ip_At, Lowest_Ip_Shl ! number of atom and of shell which correspond to the lowest ionization potential
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
    type(CDF), intent(in) :: CDF_Phonon ! declare CDF for phonons
    type(Solid), intent(in) :: Matter   ! all material parameters
    type(Density_of_states), intent(in) :: Mat_DOS  ! material DOS
    type(Error_handling), intent(inout) :: Error_message	! error messages are dealed with as objects
    type(All_MFP), dimension(:), intent(in), target :: Total_Hole_MFPs    ! electron mean free paths for all shells
    type(MFP), intent(in), target :: Elastic_Hole_MFP                        ! elastic mean free path
    real(8), intent(in) :: t_cur
    type(Flag), intent(inout) :: NumPar
    real(8), dimension(:,:), intent(in) :: El_IMFP    ! total IMFP for electron, to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), intent(in) :: El_EMFP    ! to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), intent(in) :: Hole_IMFP    ! total IMFP for electron, to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), intent(in) :: Hole_EMFP    ! to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), intent(in) :: Phot_IMFP    ! total IMFP for photons, to use in a subroutine, we need this shape of an array
    integer, intent(inout) :: Tot_Nel, Tot_Nphot    ! Total numbers of electrons and photons
    integer, intent(in) :: NOP, i
    real(8), intent(inout) :: At_NRG
    type(Differential_MFP), dimension(:), intent(in) :: DSF_DEMFP_H     !DSF differential elastic cross-sections for holes
    real(8), dimension(:), intent(inout) :: Out_R   ! [A] radius for distributions
    real(8), dimension(:), intent(inout) :: Out_V  ! inverse volume of cilinder layers [1/A^3]
    real(8), dimension(:,:), intent(inout) :: Out_Elat  ! [eV/A^3] lattuce energy density vs time vs R
    !-------------------------
    real(8) Eel, Egap, RN, dE, dE_cur, mh, R, MFP_tot, Ehole, t_Auger, t_Radiat
    real(8) theta0, phi0, theta2, phi2, theta, phi, phi1, theta1, htheta, hphi, hphi1, htheta1, hphi2, htheta2, Etest
    real(8) IMFP, EMFP, HEMFP, HIMFP, L, X, Y, Z
    real(8) E_new1, E_new2, dE_loc
    integer Nat_cur, Nshl_cur, j, Sh1, KOA1, Sh2, KOA2, ii
    character(200) :: Error_descript
    character(8) kind_of_particle
    
    kind_of_particle = 'Hole'
    Egap = Target_atoms(1)%Ip(size(Target_atoms(1)%Ip))
    VB_vs_deep:if ((All_holes(NOP)%KOA .EQ. Lowest_Ip_At) .AND. (All_holes(NOP)%Shl .EQ. Lowest_Ip_Shl)) then   ! it's VB:
         
        Eel = All_holes(NOP)%Ehkin ! [eV] hole kinetic energy (total - Egap)

        ! Find is it elastic or inelastic collision:
        call Next_free_path(Eel, Hole_IMFP, HIMFP) ! => IMFP of hole [A]
        call Next_free_path(Eel, Hole_EMFP, HEMFP) ! => EMFP of hole [A]
        call random_number(RN)
        L = All_holes(NOP)%L            ! [A] old mean free path
        theta0 = All_holes(NOP)%theta   ! old angle
        phi0 = All_holes(NOP)%phi       ! old angle
        X = All_holes(NOP)%X + L*sin(theta0)*sin(phi0)    ! [A] new X coordinate
        Y = All_holes(NOP)%Y + L*sin(theta0)*cos(phi0)    ! [A] new Y coordinate
        Z = All_holes(NOP)%Z + L*cos(theta0)              ! [A] new Z coordinate
        
        !inel_vs_el:if ((RN*(1.0d0/HIMFP + 1.0d0/HEMFP) .LT. 1.0d0/HIMFP) .AND. (HEMFP .LT. 1d15)) then  ! inelastic
        !inel_vs_el:if (RN*(1.0d0/HIMFP + 1.0d0/HEMFP) .LT. 1.0d0/HIMFP) then  ! inelastic
        inel_vs_el:if ((RN*(1.0d0/HIMFP + 1.0d0/HEMFP) .LT. 1.0d0/HIMFP) .AND. (HIMFP .LT. 1d15)) then  ! inelastic
            ! Find which shell of which atom is being ionized:
            call Which_shell(Total_Hole_MFPs, Hole_IMFP, Eel, Nat_cur, Nshl_cur)   ! => Nat_cur, Nshl_cur
            
            Tot_Nel = Tot_Nel + 1   ! we have ionized a new electron!
            ! check if electrons are too many and the size of arrays must be increased:
            call Check_size(All_electrons, All_holes, N=Tot_Nel)
            
            ! => find IMFP of incident hole [A] needed for calculation of transferred energy                        
            call Next_free_path(Eel, Total_Hole_MFPs(Nat_cur)%ELMFP(Nshl_cur)%E, Total_Hole_MFPs(Nat_cur)%ELMFP(Nshl_cur)%L, HIMFP)
            ! => dE [eV] transferred energy:
            call Electron_energy_transfer(Eel, Target_atoms, Nat_cur, Nshl_cur, HIMFP, dE, Matter, Mat_DOS, NumPar, kind_of_particle)
            ! => htheta, hphi of ionized electron, htheta1, hphi1 angles of incident hole:
            call Update_holes_angles_el(All_holes(NOP), Eel, dE, htheta, hphi, htheta1, hphi1)
                                    
            ! New ionized electrons parameters:
            ! How much energy an electron recieves:
            call Electron_recieves_E(dE, Nat_cur, Nshl_cur, Target_atoms, Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, dE_cur, Error_message)
            
            call Next_free_path(dE_cur, El_IMFP, IMFP) ! => IMFP of electron [A]
            call Next_free_path(dE_cur, El_EMFP, EMFP) ! => EMFP of electron [A]
            call random_number(RN)
            MFP_tot = -log(RN)/(1.0d0/IMFP + 1.0d0/EMFP)    ! [A] sample total electron free path (inelastic + elastic)
            theta2 = htheta
            call random_number(RN)
            phi2 = hphi
            call New_Angles_both(phi0, theta0, theta2, phi2, phi1, theta1)    ! => phi1, theta1. Angles measured from Z axis
            ! new electron parameters:
            call Particle_event(All_electrons(Tot_Nel), E=dE_cur, t0=All_holes(NOP)%t0, X=X, Y=Y, Z=Z, L=MFP_tot, theta=theta1, phi=phi1)
            ! => All_electrons(NOP)%tn is updated for next collision [fs]:
            call Get_time_of_next_event(All_electrons(Tot_Nel), MFP=MFP_tot)
            
            call cut_off(Matter%cut_off, All_electrons=All_electrons(Tot_Nel)) ! compare electron energy with cut-off    
            
            if (isnan(theta1))  then
                print*, 'ERROR #1 in Hole_Monte_Carlo'
                print*, 'El 1 = ', phi0, theta0, theta2, phi2, phi1, theta1
                pause
            endif    
            
            if ((All_electrons(Tot_Nel)%E .LT. -1.0d-9) .OR. isnan(All_electrons(Tot_Nel)%E)) then  ! Error, electron got negative energy!
                write(Error_descript, '(a,i,a,e)') 'Hole-impact electron #', Tot_Nel, ' got negative energy ', &
                    All_electrons(Tot_Nel)%E ! description of an error
                call Save_error_details(Error_message, 40, Error_descript) ! write it into the error-log file
                print*, trim(adjustl(Error_descript)) ! print it also on the sreen
            endif
            
            ! Get parameters of the created holes:
            call Update_holes_angles_SHI((dE-dE_cur), htheta, hphi)
            call Particle_event(All_holes(Tot_Nel), t0=All_holes(NOP)%t0, X=X, Y=Y, Z=Z, &
                KOA=Nat_cur, Shl=Nshl_cur, theta=htheta, phi=hphi) ! new hole parameters
            call Hole_parameters(All_holes(Tot_Nel), Matter, Mat_DOS, Target_atoms, Hole_IMFP, Hole_EMFP, &
                (dE-dE_cur), Lowest_Ip_At, Lowest_Ip_Shl)
            call cut_off(Matter%cut_off, All_holes=All_holes(Tot_Nel)) ! compare hole energy with cut-off and update it's time
                       
            if ((All_holes(Tot_Nel)%Ehkin .LT. -1.0d-9) .OR. isnan(All_holes(Tot_Nel)%Ehkin)) then  ! Error, hole got negative energy!
                ! description of an error:
                write(Error_descript, '(a,i,a,e)') 'Hole-impact hole #', Tot_Nel, ' got negative energy ', All_holes(Tot_Nel)%Ehkin
                call Save_error_details(Error_message, 41, Error_descript) ! write it into the error-log file
                print*, trim(adjustl(Error_descript)) ! print it also on the sreen
            endif
            
            if (((All_holes(Tot_Nel)%E + All_holes(Tot_Nel)%Ehkin) .LT. Egap-1.0d-12)) then  ! Error, hole in band gap!
                write(Error_descript, '(a,i,a,e,e,e)') 'Hole-impact hole #', Tot_Nel, ' is in the band gap ', &
                    All_holes(Tot_Nel)%Ehkin, All_holes(Tot_Nel)%E, Egap ! description of an error
                call Save_error_details(Error_message, 41, Error_descript) ! write it into the error-log file
                print*, trim(adjustl(Error_descript)) ! print it also on the sreen
            endif
            
            ! Check parameters of incident hole
            call check_hole_parameters(Mat_DOS, Eel, dE, Ehole, All_electrons(Tot_Nel))
                                                                                   
        else inel_vs_el  !Elastic
            ! => find EMFP of hole [A] needed for calculation of transferred energy:
            call Next_free_path(Eel, Elastic_Hole_MFP%E, Elastic_Hole_MFP%L, HEMFP)
            !call Electron_energy_transfer(Eel, HEMFP, Target_atoms, CDF_Phonon, Matter, dE, NumPar, Mat_DOS, kind_of_particle) ! => dE [eV] transferred energy
            if (NumPar%kind_of_EMFP .EQ. 2) then           ! DSF cross-sections
                call NRG_transfer_elastic_DSF(DSF_DEMFP_H, Eel, HEMFP, dE) ! module "Cross_sections"
            else if (NumPar%kind_of_EMFP .EQ. 1) then      ! CDF phonon peaks
                ! => dE [eV] transferred energy:
                call Electron_energy_transfer(Eel, HEMFP, Target_atoms, CDF_Phonon, Matter, dE, NumPar, Mat_DOS, kind_of_particle)
            else                                    ! Atomic cross-sections of Mott
                dE = 0.0d0                          ! [eV] transferred energy
                do ii = 1, size(Target_atoms)        ! for all atomic spicies:
                    call NRG_transfer_elastic_atomic(Target_atoms, ii, Eel, dE_loc, M_eff = All_holes(NOP)%Mass)
                    !dE = dE + dE_loc                 ! [eV] sum of transferred energies during collision witn all kind of atoms
                    dE = dE + dE_loc*Target_atoms(ii)%Pers ! [eV] sum of transferred energies during collision witn all kind of atoms
                enddo
                !dE = dE/size(Target_atoms)          ! [eV] Average dE over kinds of atoms
                dE = dE/dble(SUM(target_atoms(:)%Pers))          ! [eV} Average dE over kinds of atoms
            endif

            call Update_particle_angles_lat(target_atoms, Eel, dE, htheta1, hphi1, All_holes(NOP)%Mass)
            
            call check_hole_parameters(Mat_DOS, Eel, dE, Ehole)
                                                                                   
            ! save this energy in the array of radial distributions of atomic energies:
            At_NRG = At_NRG + dE    ! [eV] lattice energy
!             if (dE >= 0.0d0) then   ! lattice heating
! 	            lat_inc =  lat_inc + dE
! 	        else   ! lattice cooling
! 	            lat_decr = lat_decr + dE
! 	        endif
            R = SQRT(X*X + Y*Y) ! [A] radius of this hole at the moment of scattering
            if (isnan(R)) then 
               print*, 'ERROR #2 in Hole_Monte_Carlo'
               print*, 'hole_elastic', NOP, dE
               print*, 'NEW:', X, Y
               print*, 'OLD:', All_holes(NOP)%X, All_holes(NOP)%Y
               print*, 'Stuff:', L, theta0, phi0
               print*, 'NRGs:', All_holes(NOP)%Ehkin, HIMFP, HEMFP, RN
            endif

            call Find_in_array_monoton(Out_R, R, j) ! find where in the distribution array it is
            Out_Elat(i,j) = Out_Elat(i,j) + dE*Out_V(j)   ! [eV/A^3] here is the scattering event happend
        endif inel_vs_el
        
        ! Incident hole new parameters:
        call New_Angles_both(phi0, theta0, htheta1, hphi1, hphi2, htheta2)
        call Particle_event(All_holes(NOP), t0=All_holes(NOP)%tn, X=X, Y=Y, Z=Z, &
            KOA=All_holes(NOP)%KOA, Shl=All_holes(NOP)%Shl, theta=htheta2, phi=hphi2) ! new hole parameters
        call Hole_parameters(All_holes(NOP), Matter, Mat_DOS, Target_atoms, &
            Hole_IMFP, Hole_EMFP, (Ehole+Egap), Lowest_Ip_At, Lowest_Ip_Shl)
        call cut_off(Matter%cut_off, All_holes=All_holes(NOP)) ! compare hole energy with cut-off and update it's time

        if (isnan(All_holes(NOP)%X)) print*, 'inc_hole', NOP, All_holes(NOP)%KOA, All_holes(NOP)%Shl, Ehole+Egap, htheta2, phi2
        
        if ((All_holes(NOP)%Ehkin .LT. -1.0d-9) .OR. isnan(All_holes(NOP)%Ehkin)) then  ! Error, electron got negative energy!
            write(Error_descript, '(a,i,a,e)') 'Hole-el incident hole #', NOP, ' got negative energy ', &
                All_holes(NOP)%Ehkin ! description of an error
            call Save_error_details(Error_message, 20, Error_descript) ! write it into the error-log file
            print*, trim(adjustl(Error_descript)) ! print it also on the sreen
        endif
            
    else VB_vs_deep    ! it's deep shell: Only deep shells can decay via Auger
!        Consistency test
!        Etest = All_holes(NOP)%E ! save initial holes energy to test later
!        if ( abs(All_holes(NOP)%E -  Target_atoms(All_holes(NOP)%KOA)%Ip(All_holes(NOP)%Shl) ) .GT. 1d-12) then
!           print*, 'Holes energy is not equal to its IP:', All_holes(NOP)%E,  Target_atoms(All_holes(NOP)%KOA)%Ip(All_holes(NOP)%Shl)
!           print*, Lowest_Ip_Shl, All_holes(NOP)%Shl
!           print*, Lowest_Ip_At, All_holes(NOP)%KOA
!        endif
        call random_number(RN)
        t_Auger = Target_atoms(All_holes(NOP)%KOA)%Auger(All_holes(NOP)%Shl)   ! Auger characteristic time
        t_Radiat = Target_atoms(All_holes(NOP)%KOA)%Radiat(All_holes(NOP)%Shl) ! Radiative decay characteristic time
        ! Auger vs Radiative decay:
        auger_vs_rad:if (RN*(1.0d0/t_Auger + 1.0d0/t_Radiat) .LT. 1.0d0/t_Auger) then ! it's Auger decay:
            call Auger_decay(All_holes(NOP)%KOA, All_holes(NOP)%Shl, Target_atoms, Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, &
                             Sh1, KOA1, Sh2, KOA2, dE, E_new1, E_new2, Error_message)   ! => Sh1, KOA1, Sh2, KOA2, dE, E_new1, E_new2
            
            if (abs(All_holes(NOP)%E - (dE+E_new1+E_new2)) .GT. 1d-10) then
               print*, 'ERROR #3 in Hole_Monte_Carlo'
               print*, 'Auger does not conserve energy!', Etest - (dE+E_new1+E_new2)
               write(*,'(f,f,f,f)') Etest, dE,  E_new1, E_new2
            endif
            
            possible_dec:if (Sh2 .GT. 0) then ! decay really happened
                ! hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
                ! Change the parameters of the decayed hole:
                call Update_holes_angles_SHI(E_new1, htheta, hphi)
                call Particle_event(All_holes(NOP), t0=All_holes(NOP)%tn, KOA=KOA1, Shl=Sh1, theta=htheta, phi=hphi) ! new hole parameters
                call Hole_parameters(All_holes(NOP), Matter, Mat_DOS, Target_atoms, &
                        Hole_IMFP, Hole_EMFP, E_new1, Lowest_Ip_At, Lowest_Ip_Shl) ! Update times for the old hole
                call cut_off(Matter%cut_off, All_holes=All_holes(NOP)) ! compare hole energy with cut-off and update it's time
               
!               Consistency test
!               if ((KOA1 .NE. Lowest_Ip_At) .OR. (Sh1 .NE. Lowest_Ip_Shl)) then
!                if ( abs(All_holes(NOP)%E -  Target_atoms(KOA1)%Ip(Sh1) ) .GT. 1d-12) then
!                    print*, 'Holes energy is not equal to its IP:', All_holes(NOP)%E,  Target_atoms(KOA1)%Ip(Sh1)
!                    pause 'Hole 1'
!                endif
!               endif
                
                Tot_Nel = Tot_Nel + 1   ! we have ionized a new electron!
                ! check if electrons are too many and the size of arrays must be increased:
                call Check_size(All_electrons, All_holes, N=Tot_Nel)
                
                ! The new hole created:
                call Update_holes_angles_SHI(E_new2, htheta, hphi)
                call Particle_event(All_holes(Tot_Nel), t0=All_holes(NOP)%t0, X=All_holes(NOP)%X, Y=All_holes(NOP)%Y, &
                        Z=All_holes(NOP)%Z, KOA=KOA2, Shl=Sh2, theta=htheta, phi=hphi) ! new hole parameters
                call Hole_parameters(All_holes(Tot_Nel), Matter, Mat_DOS, Target_atoms, &
                        Hole_IMFP, Hole_EMFP, E_new2, Lowest_Ip_At, Lowest_Ip_Shl)
                call cut_off(Matter%cut_off, All_holes=All_holes(Tot_Nel)) ! compare hole energy with cut-off and update it's time
                     
!               Consistency test
!                if ((KOA2 .NE. Lowest_Ip_At) .OR. (Sh2 .NE. Lowest_Ip_Shl)) then
!                 if ( abs(All_holes(Tot_Nel)%E -  Target_atoms(KOA2)%Ip(Sh2) ) .GT. 1d-12) then
!                    print*, 'Holes energy is not equal to its IP:', All_holes(Tot_Nel)%E,  Target_atoms(KOA2)%Ip(Sh2)
!                    pause 'Hole 2'
!                 endif
!                endif     
                                   
                ! The new ionized electorns:
                call Next_free_path(dE, El_IMFP, IMFP) ! => IMFP of electron [A]
                call Next_free_path(dE, El_EMFP, EMFP) ! => EMFP of electron [A]
                call random_number(RN)
                MFP_tot = -log(RN)/(1.0d0/IMFP + 1.0d0/EMFP)    ! [A] sample total electron free path (inelastic + elastic)
                call random_number(RN)
                phi1 = 2.0d0*g_Pi*RN    ! random angle
                call random_number(RN)
                theta1 = g_Pi*RN        ! random angle
                call Particle_event(All_electrons(Tot_Nel), E=dE, t0=All_holes(NOP)%t0, X=All_holes(NOP)%X, Y=All_holes(NOP)%Y, &
                    Z=All_holes(NOP)%Z, L=MFP_tot, theta=theta1, phi=phi1) ! new electron parameters
                ! => All_electrons(NOP)%tn is updated for next collision [fs]:
                call Get_time_of_next_event(All_electrons(Tot_Nel), MFP=MFP_tot)
                
                call cut_off(Matter%cut_off, All_electrons=All_electrons(Tot_Nel)) ! compare electron energy with cut-off    
                
                ! Error, electron got negative energy!:
                if ((All_electrons(Tot_Nel)%E .LT. -1.0d-9) .OR. isnan(All_electrons(Tot_Nel)%E)) then
                    ! description of an error:
                    print*, 'ERROR #4 in Hole_Monte_Carlo'
                    write(Error_descript, '(a,i6,a,f9.3)') 'Auger electron #', Tot_Nel, ' got negative energy ', All_electrons(Tot_Nel)%E
                    call Save_error_details(Error_message, 23, Error_descript) ! write it into the error-log file
                    print*, trim(adjustl(Error_descript)) ! print it also on the sreen
                    print*, 'KOA=', All_holes(NOP)%KOA, 'Shl=', All_holes(NOP)%Shl
                    print*, 'dE=', dE, 'E_new1=', E_new1, 'E_new2=', E_new2
                    print*, 'Sh1=',Sh1, 'KOA1=', KOA1, 'Sh2=', Sh2, 'KOA2=', KOA2
                endif
            else possible_dec   ! decay is not possible for some reason, let this hole go...
                !print*, 'Impossible Auger tried to happen...'
                call Particle_event(All_holes(NOP), t0=t_cur, tn=1d21) ! set changes in the hole properties
            endif possible_dec ! (Sh2 .GT. 0)
        else auger_vs_rad ! it's Radiative decay
!            print*, 'Wow, radiative decay!'
!            pause 'Appreciate this moment...'
            call Radiative_decay(All_holes(NOP)%KOA, All_holes(NOP)%Shl, Target_atoms, Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, &
                                 Sh1, KOA1, dE, E_new1, Error_message)   ! => Sh1, KOA1, dE
                                 
            ! Change the parameters of the decayed hole:
            call Update_holes_angles_SHI(E_new1, htheta, hphi)
            call Particle_event(All_holes(NOP), t0=All_holes(NOP)%tn, X=All_holes(NOP)%X, Y=All_holes(NOP)%Y, &
                Z=All_holes(NOP)%Z, KOA=KOA1, Shl=Sh1, theta=htheta, phi=hphi) ! new hole parameters
            call Hole_parameters(All_holes(NOP), Matter, Mat_DOS, Target_atoms, Hole_IMFP, Hole_EMFP, E_new1, Lowest_Ip_At, Lowest_Ip_Shl)
            call cut_off(Matter%cut_off, All_holes=All_holes(NOP)) ! compare hole energy with cut-off and update it's time
            
            ! We created a photon:
            Tot_Nphot = Tot_Nphot + 1
            call Check_size(Photons=All_photons, N=Tot_Nphot)    ! check if photons are too many and the size of arrays must be increased
            call Next_free_path(dE, Phot_IMFP, IMFP) ! => IMFP of a photon [A]
            call random_number(RN)
            MFP_tot = -log(RN)*IMFP ! [A] sample total photon free path (inelastic only)
            call random_number(RN)
            phi1 = 2.0d0*g_Pi*RN    ! random angle
            call random_number(RN)
            theta1 = g_Pi*RN        ! random angle
            call Particle_event(All_photons(Tot_Nphot), E=dE, t0=All_holes(NOP)%t0, X=All_holes(NOP)%X, Y=All_holes(NOP)%Y, &
                Z=All_holes(NOP)%Z, L=MFP_tot, theta=theta1, phi=phi1) ! new photon's parameters
            call Get_time_of_next_event(All_photons(Tot_Nphot), MFP=MFP_tot)   ! => All_photons(NOP)%tn is updated for next collision [fs]    
            if ((All_photons(Tot_Nphot)%E .LT. -1.0d-9) .OR. isnan(All_photons(Tot_Nphot)%E)) then  ! Error, electron got negative energy!
                ! description of an error:
                print*, 'ERROR #5 in Hole_Monte_Carlo'
                write(Error_descript, '(a,i6,a,f9.3)') 'Auger photon #', Tot_Nphot, ' got negative energy ', All_photons(Tot_Nphot)%E
                call Save_error_details(Error_message, 30, Error_descript) ! write it into the error-log file
                print*, trim(adjustl(Error_descript)) ! print it also on the sreen
                print*, 'KOA=', All_holes(NOP)%KOA, 'Shl=', All_holes(NOP)%Shl
                print*, 'dE=', dE, 'E_new1=', E_new1
                print*, 'Sh1=', Sh1, 'KOA1=', KOA1
            endif
        endif auger_vs_rad
    endif VB_vs_deep
end subroutine Hole_Monte_Carlo



subroutine Photon_Monte_Carlo(All_electrons, All_holes, All_photons, El_IMFP, El_EMFP, Hole_IMFP, Hole_EMFP, &
        Phot_IMFP, CDF_Phonon, Matter, target_atoms, &
        Total_Photon_MFPs, Tot_Nel, Tot_Nphot, NOP, Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, Error_message, &
        t_cur, NumPar)
    type(Electron), dimension(:), intent(inout), allocatable :: All_electrons   ! define array of electrons
    type(Hole), dimension(:), intent(inout), allocatable :: All_holes           ! define array of holes
    type(Photon), dimension(:), intent(inout), allocatable :: All_photons       ! define array of photons
    integer, intent(in) :: Lowest_Ip_At, Lowest_Ip_Shl ! number of atom and of shell which correspond to the lowest ionization potential
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
    type(CDF), intent(in) :: CDF_Phonon ! declare CDF for phonons
    type(Solid), intent(in) :: Matter   ! all material parameters
    type(Density_of_states), intent(in) :: Mat_DOS  ! material DOS
    type(Error_handling), intent(inout) :: Error_message	! error messages are dealed with as objects
    type(All_MFP), dimension(:), intent(in), target :: Total_Photon_MFPs     ! photon mean free paths for all shells
    real(8), intent(in) :: t_cur
    type(Flag), intent(inout) :: NumPar
    real(8), dimension(:,:), intent(in) :: El_IMFP    ! total IMFP for electron, to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), intent(in) :: El_EMFP    ! to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), intent(in) :: Hole_IMFP    ! total IMFP for electron, to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), intent(in) :: Hole_EMFP    ! to use in a subroutine, we need this shape of an array
    real(8), dimension(:,:), intent(in) :: Phot_IMFP    ! total IMFP for photons, to use in a subroutine, we need this shape of an array
    integer, intent(inout) :: Tot_Nel, Tot_Nphot    ! Total numbers of electrons and photons
    integer, intent(in) :: NOP
    real(8) Eel, dE_cur, mh, htheta, hphi, Egap, IMFP, EMFP, MFP_tot, RN, theta1, phi1, Z, Y, X, theta0, phi0, L
    integer Nat_cur, Nshl_cur
    character(200) :: Error_descript
    
    ! All a photon can de so far is photoabsorbtion:
    Eel = All_photons(NOP)%E ! [eV] photon energy
    L = All_photons(NOP)%L                ! [A] old mean free path
    theta0 = All_photons(NOP)%theta     ! old angle
    phi0 = All_photons(NOP)%phi         ! old angle
    X = All_photons(NOP)%X + L*sin(theta0)*sin(phi0)    ! [A]  X coordinate of the photoabsorbtion
    Y = All_photons(NOP)%Y + L*sin(theta0)*cos(phi0)    ! [A]  Y coordinate of the photoabsorbtion
    Z = All_photons(NOP)%Z + L*cos(theta0)              ! [A]  Z coordinate of the photoabsorbtion
    call Which_shell(Total_Photon_MFPs, Phot_IMFP, Eel, Nat_cur, Nshl_cur)   ! => Nat_cur, Nshl_cur
    
    Egap = Target_atoms(1)%Ip(size(Target_atoms(1)%Ip)) ! Egap [eV]
    !eeeeeeeeeeeeeeeeeeeeeeeeeeeee
    ! New electron:
    Tot_Nel = Tot_Nel + 1   ! we have ionized a new electron!
    call Check_size(All_electrons, All_holes, N=Tot_Nel)    ! check if electrons are too many and the size of arrays must be increased
    ! New ionized electrons parameters:
    ! How much energy an electron recieves (!=> dE_cur [eV] kinetic energy of the electron):
    call Electron_recieves_E(Eel, Nat_cur, Nshl_cur, Target_atoms, Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, dE_cur, Error_message)

    call Next_free_path(dE_cur, El_IMFP, IMFP) ! => IMFP of electron [A]
    call Next_free_path(dE_cur, El_EMFP, EMFP) ! => EMFP of electron [A]
    call random_number(RN)
    MFP_tot = -log(RN)/(1.0d0/IMFP + 1.0d0/EMFP)    ! [A] sample total electron free path (inelastic + elastic)
    call New_Angles_both(All_photons(NOP)%phi, All_photons(NOP)%theta, g_Pi/2.0d0, 0.0d0, phi1, theta1)    ! => phi1, theta1
    ! new electron parameters:
    call Particle_event(All_electrons(Tot_Nel), E=dE_cur, t0=All_photons(NOP)%t0, X=X, Y=Y, Z=Z, L=MFP_tot, theta=theta1, phi=phi1)
    call Get_time_of_next_event(All_electrons(Tot_Nel), MFP=MFP_tot)   ! => All_electrons(NOP)%tn is updated for next collision [fs]
    
    call cut_off(Matter%cut_off, All_electrons=All_electrons(Tot_Nel)) ! compare electron energy with cut-off    

    if ((All_electrons(Tot_Nel)%E .LT. -1.0d-9) .OR. isnan(All_electrons(Tot_Nel)%E)) then  ! Error, electron got negative energy!
        ! description of an error:
        write(Error_descript, '(a,i,a,e)') 'Photo-electron #', Tot_Nel, ' got negative energy ', All_electrons(Tot_Nel)%E
        call Save_error_details(Error_message, 50, Error_descript) ! write it into the error-log file
        print*, trim(adjustl(Error_descript)) ! print it also on the sreen
    endif
    !hhhhhhhhhhhhhhhhhhhhhhhhhhhhh
    ! New hole:
    ! Get parameters of the created holes:
    call Update_holes_angles_SHI((Eel-dE_cur), htheta, hphi)
    ! new hole parameters:
    call Particle_event(All_holes(Tot_Nel), t0=All_photons(NOP)%t0, X=X, Y=Y, Z=Z, KOA=Nat_cur, Shl=Nshl_cur, theta=htheta, phi=hphi)
    call Hole_parameters(All_holes(Tot_Nel), Matter, Mat_DOS, Target_atoms, Hole_IMFP, Hole_EMFP, (Eel-dE_cur), Lowest_Ip_At, Lowest_Ip_Shl)
    call cut_off(Matter%cut_off, All_holes=All_holes(Tot_Nel)) ! compare hole energy with cut-off and update it's time

    if ((All_holes(Tot_Nel)%Ehkin .LT. -1.0d-9) .OR. isnan(All_holes(Tot_Nel)%Ehkin)) then  ! Error, hole got negative energy!
        ! description of an error:
        write(Error_descript, '(a,i,a,e)') 'Photo-hole #', Tot_Nel, ' got negative energy ', All_holes(Tot_Nel)%Ehkin
        call Save_error_details(Error_message, 51, Error_descript) ! write it into the error-log file
        print*, trim(adjustl(Error_descript)) ! print it also on the sreen
    endif
    if (((All_holes(Tot_Nel)%E + All_holes(Tot_Nel)%Ehkin) .LT. Egap)) then  ! Error, hole in band gap!
        ! description of an error:
        write(Error_descript, '(a,i,a,e)') 'Photo-hole #', Tot_Nel, ' is in the band gap ', All_holes(Tot_Nel)%Ehkin
        call Save_error_details(Error_message, 52, Error_descript) ! write it into the error-log file
        print*, trim(adjustl(Error_descript)) ! print it also on the sreen
    endif
    !ppppppppppppppppppppppppppppp
    ! Photon is absorbed:
    Tot_Nphot = Tot_Nphot - 1   ! we have lost one photon
    call Particle_event(All_photons(NOP), E=0.0d0, t0=1.0d27, tn=1.0d27, X=0.0d0, Y=0.0d0, Z=0.0d0, L=1.0d26, theta=0.0d0, phi=0.0d0) ! photon is gone
end subroutine Photon_Monte_Carlo



subroutine Radiative_decay(KOA, SHL, Target_atoms, Lowest_Ip_At, Lowest_Ip_Shl, Mat_DOS, &
            Sh1, KOA1, dE, E_new1, Error_message)   ! => Sh1, KOA1, dE
   integer, intent(in) :: KOA, SHL ! atomic species and shell numbers
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects
   integer, intent(in) :: Lowest_Ip_At, Lowest_Ip_Shl ! number of atom and of shell which correspond to the lowest ionization potential
   type(Density_of_states), intent(in) :: Mat_DOS   ! material DOS
   integer, intent(out) :: KOA1, Sh1	! shell and atom to which deeper hole jumps up
   real(8), intent(out) :: dE, E_new1 ! [eV] energies of 1) ejected photon, 2) decayed hole
   type(Error_handling), optional, intent(inout) :: Error_message ! deals with errors, if any
   character(200) Writing_var
   real(8) RN, Energy_diff, dE_cur
   integer coun, Shel

   ! count how many electrons on the shells that can participate in Radiative-decay:
   call count_for_Auger_shells(Target_atoms, Target_atoms(KOA)%Ip(SHL), coun)
   call random_number(RN)
   Shel = 1 + nint(RN*(coun-1)) ! this shell is where the hole jumps up
   call Choose_for_Auger_shell(Target_atoms, Target_atoms(KOA)%Ip(SHL), Shel, Sh1, KOA1)
   dE_cur = 0.0d0
   if (allocated(Mat_DOS%E)) then    ! it is VB instead of just a shell
      if ((KOA1 .EQ. Lowest_Ip_At) .AND. (Sh1 .EQ. Lowest_Ip_Shl)) then ! it is VB
         call From_where_in_VB(Mat_DOS, dE_cur)
      endif
   endif
   E_new1 = dE_cur + Target_atoms(KOA1)%Ip(Sh1)  ! new energy of the hole [eV]
   dE = Target_atoms(KOA)%Ip(SHL) - E_new1   ! [eV] energy difference between the 2 levels goes to photon                        
end subroutine Radiative_decay


subroutine cut_off(cut_off1, All_electrons, All_holes)
    real(8), intent(in) :: cut_off1
    type(Electron), optional, intent(inout) :: All_electrons
    type(Hole), optional, intent(inout) :: All_holes
    if (present(All_electrons)) then
        if (All_electrons%E .LT. cut_off1) All_electrons%tn = 1.0d20
        !print*, 'Electron cut-off', All_electrons%E, cut_off1, All_electrons%tn
    endif
    if (present(All_holes)) then
        !if ((All_holes%Ehkin .GT. 0) .AND. (All_holes%Ehkin .LT. cut_off1)) All_holes%tn = 1.0d20
        ! condition for mass to make sure it is a VB hole:
        if ((All_holes%Mass .LT. 1d15) .AND. (All_holes%Ehkin .LT. cut_off1)) All_holes%tn = 1.0d20
        !print*, 'Hole cut-off', All_holes%Mass, All_holes%Ehkin, cut_off1, All_holes%tn
    endif
end subroutine cut_off



!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! Fields subroutines - unfinished and should not be used!

subroutine find_field(X, Y, Out_R, Out_field, Efield)
    real(8), intent(in) :: X, Y
    real(8), dimension(:), intent(in) :: Out_R, Out_field
    real(8), intent(out) :: Efield
    real(8) R
    integer i
    R = sqrt(X*X + Y*Y)
    call Find_in_array_monoton(Out_R, R, i) ! find where in the distribution array it is
    if (i .EQ. 1) then
        Efield = 0.0d0
    else
        Efield = Out_field(i-1)
    endif
end subroutine find_field


subroutine update_fields(NumPar, Matter, target_atoms, Out_R, All_electrons, All_holes, tim, Out_field, Tot_Nel)
    type(Flag), intent(in) :: NumPar                            ! parameter for feilds
    type(Solid), intent(in) :: Matter                           ! all material parameters
    type(Atom), dimension(:), intent(in) :: target_atoms        ! define target atoms as objects, we don't know yet how many they are
    real(8), dimension(:), intent(in) :: Out_R                  ! [A] radius for distributions
    type(Electron), dimension(:), intent(inout) :: All_electrons   ! define array of electrons
    type(Hole), dimension(:), intent(inout) :: All_holes           ! define array of holes
    real(8), intent(in) :: tim                                  ! [fs] present time instance
    real(8), dimension(:), intent(inout) :: Out_field           ! fields vs radius [V/m (or N/C)]
    integer, intent(in) :: Tot_Nel
    
    real(8), dimension(size(Out_field)) :: Out_field_loc        ! fields vs radius [V/m (or N/C)]
    real(8) V, L0, theta0, phi0, sinthet, X, Y, Z, R
    integer i, j, k
    
        Out_field = 0.0d0                                       ! renew the array
        Out_field_loc = 0.0d0
        
        ! total number of electrons
        all_el:do k = 1, Tot_Nel   ! for all electrons
            !eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
            if (All_electrons(k)%E .GT. 0.0d0) then
                call Get_velosity(All_electrons(k), V)
                L0 = V*(tim - All_electrons(k)%t0)*1.0d-5       ! [A] mean free path travelled until this time instance
                if (L0 .LT. 0.0d0) L0 = 0.0d0                   ! just for case...
                theta0 = All_electrons(k)%theta                 ! old angle
                phi0 = All_electrons(k)%phi                     ! old angle
            else
                L0 = 0.0d0
                theta0 = 0.0d0   ! old angle
                phi0 = 0.0d0       ! old angle
            endif
            sinthet = sin(theta0)
            X = All_electrons(k)%X + L0*sinthet*sin(phi0)       ! [A] new X coordinate
            Y = All_electrons(k)%Y + L0*sinthet*cos(phi0)       ! [A] new Y coordinate
            Z = All_electrons(k)%Z + L0*cos(theta0)             ! [A] new Z coordinate
            R = SQRT(X*X + Y*Y) ! [A] radius of this el
            if (isnan(R)) then
                print*, 'ERROR in update_fields'
                print*, 'In field calculation of radius NaN occured:'
                print*, 'R', R, X, Y, L0, theta0, phi0
                print*, All_electrons(k)%X, All_electrons(k)%Y, All_electrons(k)%E, All_electrons(k)%t0, All_electrons(k)%tn
            endif
            call Find_in_array_monoton(Out_R, R, j)             ! find where in the distribution array it is
            
            call Particle_event(All_electrons(k), X=X, Y=Y, Z=Z, t0=tim) ! Update electron coordinates and time to current value
            
            Out_field_loc(j) = Out_field_loc(j) - 1.0d0         ! here is the electron
            
            !hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
            if (All_holes(k)%Mass .LT. 1.0d6) then
                call Get_velosity(All_holes(k), V, target_atoms)    ! holes
                L0 = V*(tim - All_holes(k)%t0)*1.0d-5 ! [A] mean free path travelled until this time instance
                if (L0 .LT. 0.0d0) L0 = 0.0d0  ! just for case...
                theta0 = All_holes(k)%theta   ! old angle
                phi0 = All_holes(k)%phi       ! old angle
            else
                L0 = 0.0d0
                theta0 = 0.0d0   ! old angle
                phi0 = 0.0d0       ! old angle
            endif
            sinthet = sin(theta0)
            X = All_holes(k)%X + L0*sinthet*sin(phi0)    ! [A] new X coordinate
            Y = All_holes(k)%Y + L0*sinthet*cos(phi0)    ! [A] new Y coordinate
            Z = All_holes(k)%Z + L0*cos(theta0)
            R = SQRT(X*X + Y*Y) ! [A] radius of this el
            if (isnan(R)) then
                print*, 'For some reason, there is NaN occured during field calculation for hole # ', k
                print*, 'R', R, X, Y, Z, L0, theta0, phi0
                print*, 'H', All_holes(k)%X, All_holes(k)%Y
                print*, All_holes(k)%E, All_holes(k)%Ehkin, All_holes(k)%t0, All_holes(k)%tn
                print*, All_holes(k)%L, All_holes(k)%KOA, All_holes(k)%Shl, All_holes(k)%Mass, V
            endif
            call Find_in_array_monoton(Out_R, R, j) ! find where in the distribution array it is
            
            call Particle_event(All_holes(k), X=X, Y=Y, Z=Z, t0=tim)  ! Update hole coordinates and time to current value
            Out_field_loc(j) = Out_field_loc(j) + 1.0d0   ! here is the hole
        enddo all_el ! do k = 1, Tot_Nel
               
        ! Make it right dimension:
        do i = 1, size(Out_field)
            if (i .EQ. 1) then
                Out_field(i) = Out_field_loc(i)
            else
                Out_field(i) = Out_field(i-1) + Out_field_loc(i)
            endif
        enddo
        
        open(1234, file = 'Test.txt')
        do i = 1, size(Out_field)
            write(1234, '(e,e,e)') Out_R(i), Out_field_loc(i), Out_field(i) 
        enddo
        close(1234)
        
        Out_field(:) = Out_field(:)/Out_R(:)*g_e/(2.0d0*g_Pi*g_e0*Matter%Layer)*1d20  ! [V/m]
    
end subroutine update_fields


subroutine update_particle_velocities(All_electrons, All_holes, Out_field, Out_R, Tot_Nel, NumPar, field_time, Matter, Tot_field, &
                                      Error_message, Mat_DOS, t_cur)
    real(8), dimension(:), intent(in) :: Out_R                  ! [A] radius for distributions
    type(Electron), dimension(:), intent(inout) :: All_electrons   ! define array of electrons
    type(Hole), dimension(:), intent(inout) :: All_holes           ! define array of holes
    real(8), intent(in) :: field_time                           
    real(8), intent(in) :: t_cur                  
    real(8), dimension(:), intent(inout) :: Out_field           ! fields vs radius [V/m (or N/C)]
    integer, intent(in) :: Tot_Nel                              ! Number of electron-hole pairs
    real(8), intent(out) :: Tot_field                           ! Total energy of Electric field
    type(Flag), intent(in) :: NumPar
    type(Solid), intent(in) :: Matter
    type(Error_handling), intent(inout) :: Error_message	! error messages are dealed with as objects
    type(Density_of_states), intent(in) :: Mat_DOS
    
    character(200) :: Error_descript
    real(8) Vsq, Vx, Vy, Vz, V0, Vx0, Vy0, Vz0, dt
    real(8) R, X, Y, Z, theta, phi, sintheta, sinphi, cosphi
    real(8) Efield, Ex, Ey, Ez, dE, Enrg, ch, temp_c
    integer k, j
    real(8) F_accel_el, F_decel_el, F_accel_h, F_decel_h
    F_accel_el = 0.0d0
    F_decel_el = 0.0d0
    F_accel_h = 0.0d0
    F_decel_h = 0.0d0
        
        if (field_time .LT. NumPar%field_dt) then       ! for first moment of time
            dt = field_time*1.0d-15                     ! [s]                             
        else 
            dt = NumPar%field_dt*1.0d-15                ! [s] Time grid
        endif
    !eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
        all_el:do k = 1, Tot_Nel   ! for all electrons
                    
            high_E:if (All_electrons(k)%E .GT. max(0.0d0, Matter%cut_off)) then
            
                X = All_electrons(k)%X
                Y = All_electrons(k)%Y
                Z = All_electrons(k)%Z
                R = SQRT(X*X + Y*Y)                         ! [A] radius of this el
                call Find_in_array_monoton(Out_R, R, j)     ! find where in the distribution array it is
                
                theta = All_electrons(k)%theta              ! angles of electron
                phi = All_electrons(k)%phi
                sintheta = sin(theta)
                sinphi = sin(phi)
                cosphi = cos(phi)
                
                Efield = Out_field(j)               ! [V/m] = [N/C]
                Ex = Efield*sinphi                  ! [N/C] Field x-component
                Ey = Efield*cosphi                  ! [N/C] Field y-component
                Ez = 0                              ! Field is zero along z-axis (cyllindrical symmetry)
            
                call Get_velosity(All_electrons(k), V0)
                Vx0 = V0*sintheta*sinphi            ! [m/s] Old x-component of particle velosity
                Vy0 = V0*sintheta*cosphi            ! [m/s] Old y-component of particle velosity
                Vz0 = V0*cos(theta)                 ! [m/s] Old z-component of particle velosity
                
                ch = -1
                Temp_c = ch*g_e*dt/g_me
                Vx = Vx0 + Ex*Temp_c                ! [m/s] New x-component of particle velosity
                Vy = Vy0 + Ey*Temp_c                ! [m/s] New y-component of particle velosity
                Vz = Vz0 + Ez*Temp_c                ! [m/s] New z-component of particle velosity
                
                Vsq = Vx*Vx + Vy*Vy + Vz*Vz         ! [m/s] Nev velocity modulus
                Enrg = g_me*Vsq/2.0d0/g_e
                
                ! Change of particle energy during time dt. Negative - electron get energy, Positive - Electron lose energy:
                dE = All_electrons(k)%E - Enrg
                Tot_field = Tot_field + dE          ! Virtual vault of field energy (to check conservation of energy)
                
                if (dE .LT. 0) then
                    F_accel_el = F_accel_el + dE
                else
                    F_decel_el = F_decel_el + dE
                endif
                
                All_electrons(k)%E = Enrg           ! New energy of electron
                
                if ((All_electrons(k)%E .LT. -1.0d-9) .OR. isnan(All_electrons(k)%E)) then  ! Error, electron got negative energy!
                    write(Error_descript, '(a,i,a,e, a)') 'Electron #', k, ' got negative energy ', &
                        All_electrons(k)%E, ' during interaction with electric field.' ! description of an error
                    call Save_error_details(Error_message, 22, Error_descript) ! write it into the error-log file
                    print*, trim(adjustl(Error_descript)) ! print it also on the sreen
                endif
                
                !!! New direction of movement - new angles
                
                All_electrons(k)%theta = acos(Vz/sqrt(Vsq))             ! New theta-angle
                
                if ((Vx .GE. 0.0d0) .AND. (Vy .GT. 0.0d0)) then     
                    All_electrons(k)%phi = atan(Vx/Vy)                  ! New phi-angle
                else if ((Vx .LT. 0.0d0) .AND. (Vy .GT. 0.0d0)) then
                    All_electrons(k)%phi = atan(Vx/Vy) + 2.0d0*g_pi
                else if ((Vy .LT. 0.0d0)) then
                    All_electrons(k)%phi = atan(Vx/Vy) + g_pi
                else if (Vy .EQ. 0.0d0) then
                    if (Vx .GT. 0.0d0) then
                        All_electrons(k)%phi = g_pi/2.0d0
                    else
                        All_electrons(k)%phi = 3.0d0*g_pi/2.0d0
                    endif    
                endif   
                    
            endif high_E
        enddo all_el
    !hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh   
end subroutine update_particle_velocities

END MODULE Monte_Carlo

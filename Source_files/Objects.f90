!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains subroutines to deal with objects/types,
! in the framework of the object oriented programming:

MODULE Objects
  implicit none 

!==============================================
! Particles as objects:
! Dummy particle type, all particles will have at least these atributes:
type :: Basic_particle  
   real(8) E    ! [eV] energy
   real(8) t0   ! [fs] current time
   real(8) tn   ! [fs] time of the next collision
   real(8) X, Y, Z  ! [A] coordinates
   real(8) vX, vY, vZ  ! [A] velocities
   real(8) aX, aY, aZ  ! [A] accelerations
   real(8) L        ! [A] sampled mean free path
   real(8) theta    ! current theta angle
   real(8) phi      ! current phi angle
end type Basic_particle
! Extend this type for particular particles:

type, EXTENDS (Basic_particle) :: Ion ! SHI as an object contains:
   integer Zat      ! atomic number of the SHI
   real(8) Zeff     ! effective charge in units of electron charge
   real(8) Mass     ! ion mass in the units of proton mass
   character(3) Name    ! abbreviation of the ion
   character(30) Full_Name    ! full name of the ion
   integer Kind_Zeff    ! 0=Barkas; 1=Bohr; 2=Nikolaev-Dmitriev; 3=Schiwietz-Grande;
   integer Kind_ion     ! 0=point-like charge; 1=Brandt-Kitagawa
   real(8) fixed_Zeff
end type Ion

type, EXTENDS (Basic_particle) :: Electron ! electron as an object contains:
   integer emission_flag
end type Electron

type, EXTENDS (Basic_particle) :: Hole ! hole as an object contains:
   integer KOA   ! kind of atom which this hole is sitting in
   integer Shl   ! number of shell which this hole is in
   real(8) Mass     ! [me] effective hole mass in units of electron mass
   real(8) Ehkin     ! [eV] Valence hole kinetic energy (Etot - Egap), Holes%E is used to save holes potential energy (Egap or IP)
end type Hole

type, EXTENDS (Basic_particle) :: Photon ! photon as an object
    ! No additional attributes required
end type Photon
!==============================================

!==============================================
!R_z distributions
type :: all_R_Z_distr
    real(8), dimension(:,:,:), allocatable :: Electron_Ee
    real(8), dimension(:,:,:), allocatable :: Electron_Ne
    real(8), dimension(:,:,:), allocatable :: Lattice_En
    real(8), dimension(:,:,:), allocatable :: ValHoles_RZ_ee, ValHoles_RZ_ee_pot, ValHoles_RZ_ne
    real(8), dimension(:,:), allocatable :: El_ee_Z
endtype all_R_Z_distr
!==============================================



!==============================================
! Material parameters are all in here:
type :: Solid
    character(100) Target_name
    character(100) Chem ! chemical formula
    real(8) Dens    ! [g/cm^3] material density
    real(8) Vsound  ! [m/s] speed of sound in the material
    real(8) At_Dens ! [1/cm^3] atomic material density
    real(8) Layer   ! [A] thickness of the material layer analyzed
    real(8) N_VB_el ! number of electrons in the valence band per molecule
    real(8) work_function   ![eV] work function of the material
    real(8) bar_length      ![A] length of the surface potential barrier
    real(8) bar_height      ![eV] height of the surface potential barrier
    real(8) hole_mass       !effective mass of a hole in units of electron mass
    real(8) El_eff_mass     !effective mass of an electron
    real(8) E_F             ![eV] Fermi energy
    real(8) v_f             ![m/s] fermi velocity
    real(8) cut_off         ![eV] cut_off energy
    real(8) temp            ![K] temperature
    real(8) :: Egap         ! [eV] bandgap (used for single-pole CDF automatic fit)
end type solid
!==============================================


!==============================================
! Density of states of the material:
type :: Density_of_states
   real(8), dimension(:), allocatable :: E  ! [eV] energy
   !DOS from top of valence band (for dielectrics and semiconductors)
   real(8), dimension(:), allocatable :: DOS        ! [1/eV] DOS itself
   real(8), dimension(:), allocatable :: int_DOS    ! Integral of DOS, used for calculations of probabilities
   real(8), dimension(:), allocatable :: k          ! Wave vector [1/m] 
   real(8), dimension(:), allocatable :: Eff_m      ! [me] Corresponding effective mass of particle in free particle dispersion relation
   !DOS from bottom of conduction band (inverted, for metals)
   real(8), dimension(:), allocatable :: DOS_inv        ! [1/eV] DOS itself
   real(8), dimension(:), allocatable :: int_DOS_inv    ! Integral of DOS, used for calculations of probabilities
   real(8), dimension(:), allocatable :: k_inv      ! wave vector [1/m]
   real(8), dimension(:), allocatable :: Eff_m_inv  ! [me] Corresponding effective mass of particle in free particle dispersion relation
end type Density_of_states
!==============================================


!==============================================
! Flags of options
type :: Flag
    integer :: kind_of_EMFP  ! kind of inelastic mean free path (0=atomic; 1=CDF, 2=DSF)
    integer :: kind_of_CDF   ! kid of CDF used for inelastic CS: 0=Ritchie; 1=Single=pole
    integer :: kind_of_CDF_ph   ! kid of CDF used for elastic CS: 0=Ritchie; 1=Single=pole
    integer :: kind_of_DR    ! target electron dispersion relation used in CDF calculations
    integer :: field_include ! include fields (1) or not (0)
    real(8) :: field_dt      ! time-grid for fields update
    character(1) :: path_sep ! path separator
    integer :: dt_flag       ! kind of time-grid (0=linear;1=log)
    logical :: include_photons  ! to include or not radiative decays of holes and further photons propagation
    logical :: plasmon_Emax     ! to use maximal plasmon energy as upper integration limit in cross-sections calculations
    integer :: CDF_elast_Zeff ! kind of effective charge of target atoms (1=1, 0=Barkas-like Zeff)
    real(8) :: MD_dt
    integer :: MD_grid
    integer :: Num_Z_points
    real(8) :: Zout_min
    real(8) :: Zout_max
    real(8) :: Zout_dz
    ! Flags for automatic recalcultion of MFPs:
    logical :: redo_IMFP, redo_EMFP ! do we have to?
    integer :: Last_mod_time_CDF, Last_mod_time_DSF ! Time when the CDF- and DSF-files were last modified
    ! Printout for testing:
    logical :: verbose, very_verbose
end type Flag
!==============================================


!==============================================
!DSF cross-sections
type :: Differential_MFP
    real(8) :: E           ! Incedent electron energy [eV]
    real(8), dimension(:), allocatable :: dE    ! Transferred energy in the interaction event [eV]
    real(8), dimension(:), allocatable :: dL    ! Differential over transferred energy mean free path [A]
    real(8), dimension(:), allocatable :: dL_absorb, dL_emit ! Differential MFP for absorption (dE<0) and emission of energy (dE>0)
end type Differential_MFP
!==============================================


!==============================================
! Electron mean free path for one shell, object with three entities:
type :: MFP
   real(8), dimension(:), allocatable :: E  ! [eV] energy
   real(8), dimension(:), allocatable :: L  ! [A] MFP itself
   real(8), dimension(:), allocatable :: dEdx   ! [eV/A] mean energy loss
end type MFP
!==============================================


!==============================================
! All electron mean free paths for all shells combined together:
type :: MFP_elastic
   type(MFP) :: Total   ! total MFP for scatterign on phonons
   type(MFP) :: Absorb  ! MFP for phonon absorption (inrease of electron energy)
   type(MFP) :: Emit    ! MFP for phonon emission (derease of electron energy)
end type MFP_elastic
!==============================================

!==============================================
! All electron mean free paths for all shells combined together:
type :: All_MFP
   type(MFP), dimension(:), allocatable :: ELMFP   ! contains all data about mean free path for all shells
end type All_MFP
!==============================================


!==============================================
! Complex dielectric function as an object:
type :: CDF ! parameters entering Ritchi-CDF:
   real(8), DIMENSION(:), ALLOCATABLE :: E0       ! parameter E0
   real(8), DIMENSION(:), ALLOCATABLE :: A        ! parameter A
   real(8), DIMENSION(:), ALLOCATABLE :: Gamma    ! parameter gamma
   ! For Delta-CDF:
   real(8), DIMENSION(:), ALLOCATABLE :: alpha   ! parameter alpha
end type CDF
!==============================================


!==============================================
! Target atom as an object:
type :: Atom ! atom as an object contains the following info:
   integer Zat ! atomic number
   real(8) Mass ! atomic mass in [proton mass] units
   real(8) Pers  ! persentage of this kind of atoms in the compound
   character(3) Name    ! abbreviation
   character(30) Full_Name  ! full name
   integer N_shl ! number of shells for this kind of atom
   character(11), dimension(:), allocatable :: Shell_name   ! names of the shells
   integer, dimension(:), allocatable :: Shl_num    ! shell designator within EADL-standard: http://www-nds.iaea.org/epdl97/
   real(8), dimension(:), allocatable :: Nel    ! number of electrons in this shell for all shells
   real(8), dimension(:), allocatable :: Ip ! [eV] ionization potentials of all shells
   real(8), dimension(:), allocatable :: Ek ! [eV] average kinetic energy of all shells (from EADL)
   real(8), dimension(:), allocatable :: Auger  ! [fs] auger-decay time for all shells
   real(8), dimension(:), allocatable :: Radiat ! [fs] radiative-decay time for all shells
   class(CDF), dimension(:), allocatable :: Ritchi  ! Ritchi CDF defined for each shell
   integer, dimension(:), allocatable :: PQN   ! principal quantum number of this shell
   integer, dimension(:), allocatable :: KOCS  ! Kind of cross section for electrons and holes (0=CDF, 1=BEB, ...)
   integer, dimension(:), allocatable :: KOCS_SHI  ! Kind of cross section for SHI (0=CDF, 1=BEB, ...)
end type Atom
!==============================================

!==============================================
! For reading atomic data from our periodic table:
type Atomic_data    ! our internal atomic database "INPUT_atomic_data.dat"
   integer :: Z     ! atomic number
   real(8) :: Mass  ! atomic mass [atomic units]
   character(15) :: Full_Name   ! Full atomic name
   character(3) :: Name         ! element name 
   real(8) :: Nvb   ! number of valence electrons
endtype Atomic_data
!==============================================


!==============================================
! All distributions in the cylindrical layers:
type Distr_shells
    real(8), dimension(:,:), allocatable :: nh    ! density of holes for different shells vs R [1/cm^3] vs time [fs]
    real(8), dimension(:,:), allocatable :: Eh    ! energy of holes for different shells vs R [eV/A^3] vs time [fs]
    real(8), dimension(:,:), allocatable :: Ehkin    ! energy of holes for different shells vs R [eV/A^3] vs time [fs]
end type Distr_shells
! Now combine shells into one type:
type Distr_on_shells
    type(Distr_shells), dimension(:), allocatable :: Shl    ! which shell of the atom
end type Distr_on_shells
! Now combine it all into one type:
type Cylinder_distr
    real(8), dimension(:), allocatable :: R    ! cylindrical radii [A]
    real(8), dimension(:,:), allocatable :: ne   ! electron densities [1/cm^3] vs time [fs]
    real(8), dimension(:,:), allocatable :: Ee   ! electron energies [eV/A^3] vs time [fs]
    class(Distr_on_shells), dimension(:), allocatable :: Atom  ! which sort of atom
end type Cylinder_distr
!==============================================


!==============================================
! File names for printing out stuff:
type All_names
   character(300), dimension(:), allocatable :: F
end type All_names
!==============================================


!==============================================
! Error handling as an object:
type Error_handling
   LOGICAL Err		! indicates that an error occured
   integer Err_Num	! assign a number to an error
   integer File_Num		! number of the file with error log
   character(100) Err_descript	! describes more details about the error
end type
!==============================================


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 contains

! How to write the log about an error:
subroutine Save_error_details(Err_name, Err_num, Err_data)
   class(Error_handling) :: Err_name    ! object containing all details
   integer, intent(in) :: Err_num       ! number of error asigned
   character(100), intent(in) :: Err_data   ! description of the error
   integer FN   ! number of file where to write error log
   FN = Err_name%File_Num   ! this number is provided in the Err_name object
   Err_name%Err = .true.    ! error occured, mark it as "true"
   Err_name%Err_Num = Err_num   ! number of error we asign to it
   Err_name%Err_descript = Err_data ! descriptino of an error
   write(FN, '(a,i2,1x,a)') 'Error #', Err_name%Err_Num, trim(adjustl(Err_name%Err_descript))   ! write it all into the file
end subroutine Save_error_details


! To change parameters of a certain type:
subroutine Particle_event(Particle, E, t0, tn, X, Y, Z, L, & ! for each particle
                          Zat, Zeff, Mass, Name, Full_name, Kind_Zeff, &    ! for ion
                          vx, vy, vz, ax, ay, az, &
                          theta, phi, &   ! for electron or photon
                          KOA, Shl, Ehkin) ! for hole
   class(Basic_particle), intent(inout) :: Particle ! can be Ion, or Electron, or Hole
   real(8), intent(in), optional :: E    ! [eV] energy
   real(8), intent(in), optional :: Ehkin    ! [eV] energy
   real(8), intent(in), optional :: t0   ! [fs] current time
   real(8), intent(in), optional :: tn   ! [fs] time of the next collision
   real(8), intent(in), optional :: X, Y, Z  ! [A] coordinates
   real(8), intent(in), optional :: vx, vy, vz  ! [m/s] velosities
   real(8), intent(in), optional :: ax, ay, az  ! [m/s^2] accelerations
   real(8), intent(in), optional :: L        ! [A] neaxt mean free path
   integer, intent(in), optional :: Zat      ! atomic number of the SHI
   real(8), intent(in), optional :: Zeff     ! effective charge in units of electron charge
   real(8), intent(in), optional :: Mass     ! ion mass in the units of proton mass
   character(3), intent(in), optional :: Name    ! abbreviation of the ion
   character(30), intent(in), optional :: Full_Name    ! full name of the ion
   integer, intent(in), optional :: Kind_Zeff    ! 0=point-like charge, 1=Bradt-Kitagawa charge
   real(8), intent(in), optional :: theta    ! current theta angle
   real(8), intent(in), optional :: phi      ! current phi angle
   integer, intent(in), optional :: KOA   ! kind of atom which this hole is sitting in
   integer, intent(in), optional :: Shl   ! number of shell which this hole is in

   if (present(E)) Particle%E = E       ! [eV] change energy
   if (present(t0)) Particle%t0 = t0    ! [fs] change current time instance
   if (present(tn)) Particle%tn = tn    ! [fs] change time of next event
   if (present(X)) Particle%X = X       ! [A] change X-coordinate
   if (present(Y)) Particle%Y = Y       ! [A] change Y-coordinate
   if (present(Z)) Particle%Z = Z       ! [A] change Z-coordinate
   if (present(L)) Particle%L = L       ! next mean free path
   if (present(vx)) Particle%vx = vx       ! [m/s] change X-velosity component
   if (present(vy)) Particle%vy = vy       ! [m/s] change y-velosity component
   if (present(vz)) Particle%vz = vz       ! [m/s] change z-velosity component
   if (present(ax)) Particle%ax = ax       ! [m/s^2] change X-acceleration component
   if (present(ay)) Particle%ay = ay       ! [m/s^2] change y-acceleration component
   if (present(az)) Particle%az = az       ! [m/s^2] change z-acceleration component


   select type(Particle)
    class is (Ion)
        if (present(Zat)) Particle%Zat = Zat        ! atomic number of the SHI
        if (present(Zeff)) Particle%Zeff = Zeff     ! effective charge in units of electron charge
        if (present(Mass)) Particle%Mass = Mass     ! ion mass in the units of proton mass
        if (present(Name)) Particle%Name = Name     ! abbreviation of the ion
        if (present(Full_Name)) Particle%Full_Name = Full_Name    ! full name of the ion
        if (present(Kind_Zeff)) Particle%Kind_Zeff = Kind_Zeff    ! 0=point-like charge, 1=Bradt-Kitagawa charge
    class is (Electron)
        if (present(theta)) Particle%theta = theta  ! current theta angle
        if (present(phi)) Particle%phi = phi        ! current phi angle
    class is (Hole)
        if (present(KOA)) Particle%KOA = KOA        ! kind of atom which this hole is sitting in
        if (present(Shl)) Particle%Shl = Shl        ! number of shell which this hole is in
        if (present(theta)) Particle%theta = theta  ! current theta angle
        if (present(phi)) Particle%phi = phi        ! current phi angle
        if (present(Mass)) Particle%Mass = Mass     ! hole mass in the units of proton mass
        if (present(Ehkin)) Particle%Ehkin = Ehkin  ! hole kinetic energy
    class is (Photon)
        if (present(theta)) Particle%theta = theta  ! current theta angle
        if (present(phi)) Particle%phi = phi        ! current phi angle
   end select
end subroutine Particle_event


subroutine Check_size(Electrons, Holes, Photons, N)
    type(Electron), dimension(:), allocatable, intent(inout), optional :: Electrons
    type(Hole), dimension(:), allocatable, intent(inout), optional :: Holes
    type(Photon), dimension(:), allocatable, intent(inout), optional :: Photons
    integer, intent(in) :: N    ! the number of particle from the array
    integer i, j, M
    if (present(Electrons) .AND. present(Holes)) then
        i = size(Electrons)
        j = size(Holes)
        if (N .GT. i) then    ! increase the array size
            print*, 'Size of arrays All_electrons and All_holes was increased by 1000'
            M = i + 1000
            !M = N
            call resize_array(Electrons, Holes, i=i, M=M)
            Print*, 'and equals to ', M
        endif
    endif
    if (present(Photons)) then
        i = size(Photons)
        if (N .GT. i) then    ! increase the array size
            print*, 'Size of arrays All_photons was increased by 1000'
            M = i + 1000
            call resize_array(Photons=Photons, i=i, M=M)
            Print*, 'and equals to ', M
        endif
    endif
end subroutine Check_size

subroutine resize_array(Electrons, Holes, Photons, i, M)
    type(Electron), dimension(:), allocatable, intent(inout), optional :: Electrons
    type(Hole), dimension(:), allocatable, intent(inout), optional :: Holes
    type(Photon), dimension(:), allocatable, intent(inout), optional :: Photons
    integer, intent(in) :: i, M ! last size, needed size
    type(Electron), dimension(:), allocatable :: e_trans
    type(Hole), dimension(:), allocatable :: h_trans
    type(Photon), dimension(:), allocatable :: p_trans
    integer j
    if (present(Electrons) .AND. present(Holes)) then
        allocate(e_trans(i))    ! define arrays of electrons
        allocate(h_trans(i))    ! define arrays of holes
        e_trans = Electrons     ! save arrays of electrons
        h_trans = Holes         ! save arrays of holes
        deallocate(Electrons)   ! resize it
        allocate(Electrons(M))  ! resize it
        deallocate(Holes)       ! resize it
        allocate(Holes(M))      ! resize it
        Electrons(1:i) = e_trans(1:i)
        Holes(1:i) = h_trans(1:i)
        do j = i+1, M
            call Particle_event(Electrons(j), E=0.0d0, t0=Electrons(1)%t0, tn=1d20, X=0.0d0, Y=0.0d0, Z=0.0d0, L=1d20, theta=0.0d0, phi=0.0d0) ! set initial electron data
            call Particle_event(Holes(j), E=0.0d0, Ehkin=0.0d0, t0=Electrons(1)%t0, tn=1d21, X=0.0d0, Y=0.0d0, Z=0.0d0, L=1.0d30, KOA=0, Shl=0, Mass=1.0d30, theta=0.0d0, phi=0.0d0) ! set initial hole data
        enddo
        deallocate(e_trans)
        deallocate(h_trans)
    endif
    if (present(Photons)) then
        allocate(p_trans(i))    ! define arrays of photons
        p_trans = Photons     ! save arrays of photons
        deallocate(Photons)   ! resize it
        allocate(Photons(M))  ! resize it
        Photons(1:i) = p_trans(1:i)
        do j = i+1, M
            call Particle_event(Photons(j), E=0.0d0, t0=Photons(1)%t0, tn=1d20, X=0.0d0, Y=0.0d0, Z=0.0d0, L=1d20, theta=0.0d0, phi=0.0d0) ! set initial electron data
        enddo
        deallocate(p_trans)
    endif
end subroutine resize_array


subroutine Convert_Obj_to_arrays(Out_R_Z, Out_Electron_RZ_ee, Out_Electron_RZ_ne, Out_lattice_RZ, Out_ValHoles_RZ_ee, &
                                  Out_ValHoles_RZ_ee_pot, Out_ValHoles_RZ_ne, Out_Electron_Z_ne)
    type(all_R_Z_distr), intent(in) :: Out_R_Z
    real(8), dimension(:,:,:), intent(inout) :: Out_Electron_RZ_ee, Out_Electron_RZ_ne, Out_Lattice_RZ
    real(8), dimension(:,:,:), intent(inout) :: Out_ValHoles_RZ_ee, Out_ValHoles_RZ_ee_pot, Out_ValHoles_RZ_ne
    real(8), dimension(:,:), intent(inout) :: Out_Electron_Z_ne

    Out_Electron_RZ_ee = Out_R_Z%Electron_Ee
    Out_Electron_RZ_ne = Out_R_Z%Electron_ne
    Out_Lattice_RZ = Out_R_Z%Lattice_En
    Out_Electron_Z_ne = Out_R_Z%El_ee_Z
    Out_ValHoles_RZ_ee = Out_R_Z%ValHoles_RZ_ee
    Out_ValHoles_RZ_ee_pot = Out_R_Z%ValHoles_RZ_ee_pot
    Out_ValHoles_RZ_ne = Out_R_Z%ValHoles_RZ_ne

end subroutine Convert_Obj_to_arrays
 
end MODULE Objects

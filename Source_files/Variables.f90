!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains all global variables

module Variables
use Objects                 ! Objects.f90

implicit none

!-----------------------------------------------
! Define variables used in the main program:
! All these types are defined in the module "Objects":
type(Error_handling) Error_message	! error messages are dealed with as objects
type(All_names) :: File_names   ! all file names for printing out stuff
type(Atom), dimension(:), allocatable :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
type(Ion) :: SHI   ! declare SHI as an object with atributes "Ion"
type(All_MFP), dimension(:), allocatable, target :: Total_el_MFPs   ! electron mean free paths for all shells
type(All_MFP), dimension(:), allocatable, target :: Total_Hole_MFPs ! hole mean free paths for all shells
type(All_MFP), dimension(:), allocatable, target :: Total_Photon_MFPs   ! photon mean free paths for all shells

type(All_MFP), dimension(:), allocatable, target :: SHI_MFP         ! SHI mean free paths for all shells
type(All_MFP), dimension(:), allocatable, target :: diff_SHI_MFP    ! SHI differential mean free paths for all shells
type(CDF) :: CDF_Phonon ! declare CDF for phonons
!type(MFP), target :: Elastic_MFP         ! elastic mean free path of an electron
!type(MFP), target :: Elastic_Hole_MFP    ! elastic mean free path of a hole
type(MFP_elastic), target :: Elastic_MFP         ! elastic mean free path of an electron
type(MFP_elastic), target :: Elastic_Hole_MFP    ! elastic mean free path of a hole

type(Solid) :: Matter   ! all material parameters
type(Density_of_states) :: Mat_DOS  ! material DOS
!type(Density_of_states) :: Mat_DOS_inv  ! inverted material DOS (for dispersion relation of metals)
type(Cylinder_distr) :: Out_Distr   ! OUTPUT radial distributions

type(Differential_MFP), dimension(:), allocatable :: DSF_DEMFP, DSF_DEMFP_H
type(Flag) :: NumPar ! numerical parameters and flags are here
    
!-----------------------------------------------
! Those below are normal standard type variables:
real(8) Tim ! [fs] total time
real(8) dt  ! [fs] timestep
real(8) as1, SHI_E, SHI_IMFP, IMFP, SHI_dEdx, dEdx, Sigma, E
integer Nshtot
integer i, j, k, ctim(8), c1(8), MC_stat, Nit, FN, Nel, N_temmp
integer NMC ! number of MC iterations
integer Lowest_Ip_At, Lowest_Ip_Shl ! number of atom and shell corresponding to valence band
character(100) Input_files, Material_name, Output_path, Output_path_SHI
character(200), dimension(10) :: Output_file    ! for number of different files
character(100) Text_var
logical file_exist  ! for checking existance of files
logical read_well ! there was an error in the reading file, true/false?
real(8) lat_inc, lat_decr
! For OPENMP:
integer Num_th, my_id
!For linux
CHARACTER(len = 100) :: path

character(8) kind_of_particle

!-----------------------------------------------
! Since derived types do not work with openmp, we have to temporary decompose the 'Cylinder_distr' into arrays:
real(8), dimension(:), allocatable :: Out_tot_Ne
real(8), dimension(:), allocatable :: Out_tot_Nphot
real(8), dimension(:), allocatable :: Out_tot_E
real(8), dimension(:), allocatable :: Out_R
real(8), dimension(:), allocatable :: Out_V
real(8), dimension(:), allocatable :: Out_E_e
real(8), dimension(:), allocatable :: Out_E_phot
real(8), dimension(:), allocatable :: Out_E_at, Out_E_at_heating, Out_E_at_cooling
real(8), dimension(:,:,:), allocatable :: Out_E_h
real(8), dimension(:,:), allocatable :: Out_ne
real(8), dimension(:,:), allocatable :: Out_Ee
real(8), dimension(:,:), allocatable :: Out_nphot
real(8), dimension(:,:), allocatable :: Out_Ephot
real(8), dimension(:,:), allocatable :: Out_Ee_vs_E ! electron spectrum
real(8), dimension(:,:), allocatable :: Out_Eh_vs_E ! valence hole spectrum
real(8), dimension(:,:), allocatable :: Out_Elat, Out_Elat_heating, Out_Elat_cooling
real(8), dimension(:,:), allocatable :: Out_Eat_dens  ! [eV/A^3] atom's energy
real(8), dimension(:,:,:,:), allocatable :: Out_nh
real(8), dimension(:,:,:,:), allocatable :: Out_Eh
real(8), dimension(:,:,:,:), allocatable :: Out_Ehkin
real(8), dimension(:,:), allocatable :: Out_theta, Out_theta_h  ! velosity angular distribution for electrons and VB holes
real(8), dimension(:), allocatable :: Out_theta1
real(8), dimension(:,:), allocatable :: Out_Ee_vs_E_Em
real(8), dimension(:), allocatable:: Out_Ne_Em
real(8), dimension(:), allocatable :: Out_E_Em 
real(8), dimension(:), allocatable :: Out_E_field
real(8), dimension(:,:), allocatable :: Out_field_all ! [V/m] electrical fields vs time vs R

real(8), dimension(:), allocatable :: Out_diff_coeff


!-----------------------------------------------
character(100), parameter :: dashline = '--------------------------------------------------------'
character(100), parameter :: starline = '********************************************************'
!-----------------------------------------------
contains

subroutine get_path_separator(path_sep, Error_message, read_well)
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    CHARACTER(len=1), intent(out) :: path_sep   ! path separator
    logical, intent(inout) :: read_well ! did the data read well?
    CHARACTER(len = 100) :: path, Err_data
    CALL get_environment_variable("PATH",path)
    if (path(1:1) .EQ. '/') then        !unix based OS
        path_sep = '/'
    else if (path(3:3) .EQ. '\') then   !Windows OS
        path_sep = '\'
    else
        Err_data = 'Path separator is not defined' ! unknown OS
        call Save_error_details(Error_message, 0, Err_data)
        read_well = .false. ! didn't read well this data
    endif
end subroutine

endmodule

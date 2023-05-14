!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains all subroutines needed to calculate cross-sections:

MODULE Cross_sections
  use Universal_Constants               ! let it use universal constants
  use Objects                           ! since it uses derived types, it must know about them from module 'Objects'
  use Reading_files_and_parameters      ! use some useful subroutines from this module
implicit none
! this interface finds by itself which of the two subroutine to use depending on the parameters passed:
interface Electron_energy_transfer ! finds energy transfer in electron collision
    module procedure Electron_energy_transfer_inelastic
    module procedure Electron_energy_transfer_elastic
end interface Electron_energy_transfer

!private  ! hides items not listed on public statement 
public :: Electron_energy_transfer


contains

subroutine Tot_Phot_IMFP(Ele, Target_atoms, Nat, Nshl, Sigma, dEdx, Matter, NumPar)
    real(8), intent(in) :: Ele  ! electron energy [eV]
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    integer, intent(in) :: Nat, Nshl    ! number of atom, number of shell
    real(8), intent(out) :: Sigma, dEdx   ! calculated inverse mean free path (cross-section) [1/A], and the energy losses [eV/A]
    type(Solid), intent(in) :: Matter ! properties of material
    type(Flag), intent(in) :: NumPar
    integer i, j, n
    real(8) temp1, ImE, Sigma1

    select case (Target_atoms(Nat)%KOCS_SHI(Nshl)) ! which inelastic cross section to use (EPDL97 vs CDF):
    case (1) ! CDF cross section
       if (Ele .LT. Target_atoms(Nat)%Ip(Nshl)) then ! no ionization possible, IMFP = infinity
          Sigma = 1d30      ! [A] IMFP
          dEdx = Ele/Sigma  ! energy losses [eV/A]
       else
          call Imewq(Ele, 0.0d0, Target_atoms, Nat, Nshl, ImE, Matter, NumPar=NumPar, photon=.true.) ! constructs full Im(-1/e(w,q)) as a sum of Drude-like functions
          Sigma = g_cvel*g_h/(ImE*Ele*g_e)*1d10 ! [A] IMFP
          dEdx = Ele/Sigma    ! energy losses [eV/A]
       endif
    case default ! extract the cross section from EPDL97:
       if (Target_atoms(Nat)%Shl_num(Nshl) .GE. 63) then ! it's VB
          Sigma1 = 0.0d0
          do i = Nat, size(Target_atoms) ! for all atoms, find outermost shells that form the VB:
             if (i .NE. Nat) then
                j = Target_atoms(i)%Shl_num(Target_atoms(i)%N_shl)
                call next_designator(n, j) ! find the designator for the VB (next shell after last one)
                call get_photon_cross_section_EPDL(Ele, Target_atoms, i, Shl_dsgntr=n, Matter=Matter, Sigma=Sigma)
             else
                call get_photon_cross_section_EPDL(Ele, Target_atoms, i, Nshl=Nshl, Matter=Matter, Sigma=Sigma)
             endif
             Sigma1 = Sigma1 + Sigma*(Target_atoms(i)%Pers)/SUM(Target_atoms(:)%Pers)
          enddo
          temp1 = Matter%At_Dens*1d-24
          Sigma = 1.0d0/(temp1*Sigma1) ! IMFP [A]
       else
          call get_photon_cross_section_EPDL(Ele, Target_atoms, Nat, Nshl=Nshl, Matter=Matter, Sigma=Sigma)
          temp1 = Matter%At_Dens*1d-24*(Target_atoms(Nat)%Pers)/SUM(Target_atoms(:)%Pers)
          Sigma = 1.0d0/(temp1*Sigma) ! IMFP [A]
       endif
    end select
end subroutine Tot_Phot_IMFP


subroutine TotIMFP(Ele, Target_atoms, Nat, Nshl, Sigma, dEdx, Matter, Mat_DOS, NumPar, kind_of_particle)
    real(8), intent(in) :: Ele  ! electron energy [eV]
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    integer, intent(in) :: Nat, Nshl    ! number of atom, number of shell
    real(8), intent(out) :: Sigma, dEdx   ! calculated inverse mean free path (cross-section) [1/A], and the energy losses [eV/A]
    type(Density_of_states), intent(in) :: Mat_DOS
    type(Solid), intent(in) :: Matter ! properties of material
    type(Flag), intent(in) :: NumPar
    character(8), intent(in) :: kind_of_particle
       
    integer i, j, n, Mnum
    real(8) Emin, Emax, E, dE, dL, Ltot1, Ltot0, ddEdx, a, b, temp1, temp2, Eplasmon, Egap, Mass
    
    Emin = Target_atoms(Nat)%Ip(Nshl)   ! [eV] ionization potential of the shell is minimum possible transferred energy
    Egap = Target_atoms(1)%Ip(size(Target_atoms(1)%Ip))   ! band gap [eV]
    if (Emin .LE. 1.0d-3) Emin = 1.0d-3 ! for metals there is no band gap
    
    if (trim(adjustl(kind_of_particle)) .EQ. 'Electron') then 
        Emax = (Ele + Emin)/2.0e0 ! [eV] maximum energy transfer, accounting for equality of electrons
        Mass = 1.0d0
    else if (trim(adjustl(kind_of_particle)) .EQ. 'Hole') then
        if(Matter%hole_mass .GE. 0) then
            Mass = Matter%hole_mass
        else                    ! Define mass from DOS
            call find_in_array_monoton(Mat_DOS%E, Ele, Mnum)
            Mass = Mat_DOS%Eff_m(Mnum)
            !print*, 'Mass h', Ele, Mat_DOS%DOS(Mnum), Mass
        endif
        Emax = 4.0d0*Ele*Mass/((Mass+1.0d0)*(Mass+1.0d0))
    else
        Emax = (Ele + Emin)/2.0e0 ! [eV] maximum energy transfer, accounting for equality of electrons
        Mass = 1.0d0
    endif

    ! Use maximal plasmon energy as upper limit of integration
    if (Numpar%plasmon_Emax) then       ! If included
        if (Emin .EQ. Egap) then        ! For VB only
            Eplasmon = sqrt(Matter%N_VB_el*Matter%At_Dens*1d6*g_h*g_h/(g_me*g_e0) + Egap*Egap)    ! Maximal energy of plasmons
            if (Eplasmon .GE. Emax) Emax = Eplasmon ! single atom vs plasmon
            if (Ele .LT. Emax) Emax = Ele ! no more than the total electron energy
        endif
    endif
    
    select case (Target_atoms(Nat)%KOCS(Nshl)) ! which inelastic cross section to use (BEB vs CDF):
    case (1) ! CDF cross section
      DLTA:if (NumPar%kind_of_DR == 4) then  ! Delta-CDF
       Sigma = Integral_CDF_delta_CS(g_me, g_me, Ele, Target_atoms(Nat)%Ritchi(Nshl), Emin, Matter%At_Dens, .true.) ! [A^2] Below
       ! Convert into MFP:
       Sigma = MFP_from_sigma(Sigma, Matter%At_Dens)   ! [A] below
       dEdx = energy_loss_delta(Ele, g_me, 1.0d0, Emin, Matter%At_Dens, g_me, Target_atoms(Nat)%Ritchi(Nshl), .true.) ! [eV/A] below
      else DLTA  ! Ritchie-CDF:
       n = 20*(MAX(INT(Emin),10))    ! number of integration steps
       dE = (Emax - Emin)/(real(n)) ! differential of transferred momentum [kg*m/s]
       i = 1       ! to start integration
       E = Emin    ! to start integration
       Ltot1 = 0.0d0
       Ltot0 = 0.0d0
       call Diff_cross_section(Ele, E, Target_atoms, Nat, Nshl, Ltot0, Mass, Matter, Mat_DOS, NumPar)
       ddEdx = 0.0d0
       do while (E .LE. Emax) ! integration
          dE = (1.0d0/(E+1.0d0) + E)/real(n)
          ! If it's Simpson integration:
          a =  E + dE/2.0d0
          call Diff_cross_section(Ele, a, Target_atoms, Nat, Nshl, dL, Mass, Matter, Mat_DOS, NumPar)
          temp1 = dL
          b = E + dE
          call Diff_cross_section(Ele, b, Target_atoms, Nat, Nshl, dL, Mass, Matter, Mat_DOS, NumPar)
          
!           if (Ele == 1.0d4) PRINT*, Ele, E, dL
          
          temp2 = dE/6.0d0*(Ltot0 + 4.0d0*temp1 + dL)
          Ltot1 = Ltot1 + temp2
          ddEdx = ddEdx + E*temp2
          Ltot0 = dL
          E = E + dE  ! [eV]
       enddo
!        if (Ele == 1.0d4) then
!           print*, 'Ele', Ele
!           pause 
!        endif
       
       if (Ltot1 < 1.0d-10) then
           Sigma = 1.0d20
           dEdx = 0.0d0 !*dE ! energy losses [eV/A]
       else
           Sigma = 1.0d0/(Mass*Ltot1) !*dE ! [A]
           dEdx = Mass*ddEdx !*dE ! energy losses [eV/A]
       endif
      endif DLTA  ! Ritchie-CDF
    case default ! BEB cross section
       Sigma = Sigma_BEB(Ele,Target_atoms(Nat)%Ip(Nshl),Target_atoms(Nat)%Ek(Nshl),Target_atoms(Nat)%Nel(Nshl)) ! [A^2] cross section
       temp1 = Matter%At_Dens*1d-24*(Target_atoms(Nat)%Pers)/SUM(Target_atoms(:)%Pers)
       Sigma = 1.0d0/(Mass*temp1*Sigma) ! IMFP [A]
       ddEdx = dSigma_w_int_BEB(Ele, (Ele-1.0d0)/2.0d0, Target_atoms(Nat)%Ip(Nshl), Target_atoms(Nat)%Ek(Nshl), Target_atoms(Nat)%Nel(Nshl)) ! cross section integrated with energy [A^2*eV]
       dEdx = Mass*temp1*ddEdx ! energy losses [eV/A]
    end select
end subroutine TotIMFP


function Integral_CDF_delta_CS(M, mt, E, CDF_coefs, Ip, nat, identical, hw_phonon, Emax_in) result(Sigma)
   real(8) Sigma    ! [A^2]
   real(8), intent(in) :: M, mt ! [kg]
   real(8), intent(in) :: E ! [eV]
   type(CDF), intent(in), target :: CDF_coefs   ! CDF coefficients
   real(8), intent(in) :: Ip    ! [eV] ionization potential
   real(8), intent(in) :: nat   ! [1/cm^3] atomic density
   logical, intent(in) :: identical   ! identical particles, true or not
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency [eV]
   real(8), intent(in), optional :: Emax_in ! [eV] integration limit, used for defining transfered energy
   !---------------------------------------
   real(8) :: Mc2, mtc2, CS, P, Mc22, Emin, Emax, Eeq, dEed, Estart
   real(8) :: a, b  ! coefficients of the linear extrapolation function
   real(8), dimension(:), pointer :: alpha, E0
   integer :: i, Nosc
   
   ! How many delta-functions:
   Nosc = size(CDF_coefs%E0(:))
   
   Mc2 = rest_energy(M)   ! module "Relativity"
   mtc2 = rest_energy(mt)   ! module "Relativity"
   ! Parameters of the delta-functions CDF:
   E0 => CDF_coefs%E0
   alpha => CDF_coefs%alpha
   
   CS = 0.0d0   ! to start with
   do i = 1, Nosc
      Emin = W_min(Ip, Mc2, mtc2, E, E0(i)) ! module "CS_integration_limits"

      if (abs(M - mt)/M < 1.0d-6) then  ! identical particles:
         Estart = Ip
      else  ! non-identical
         Estart = minimal_sufficient_E(Ip, Mc2, mtc2) ! module "CS_integration_limits"
      endif
      
      ! Where to switch from the approximate formula to linear extrapolation:
      call  find_Wmax_equal_Wmin(Mc2, mtc2, identical, E, Ip, E0(i), Eeq)   ! module "CS_integration_limits"
      
      if (E <= Estart) then ! cannot ionize
         CS = 0.0d0
         P = 0.0d0
      else
         dEed = Eeq/100.0d0    ! shift a little the point at which we switch to linear extrapolation
         
         if (E <= Eeq+dEed) then   ! delta-CDF model does not apply here, replace the CDF with the linear extrapolation cross section
            ! Find the coefficients of the linear approximation:
            call Find_linear_a_b(alpha(i), M, Mc2, mtc2, E0(i), Ip, identical, nat, Eeq+dEed, a, b) ! below
            CS = a*E + b    ! linearly approximated cross section
            P = 1.0d0   ! no prefactor needed in this case, it alreaxy was included into the linear extrapolation
         else   ! delta-model works
            if (present(hw_phonon)) then  ! scattering on phonons:
               Emax = W_max(Mc2, mtc2, identical, E, Ip, hw_phonon) ! module "CS_integration_limits"
            else  ! scattering on electrons or atoms:
               Emax = W_max(Mc2, mtc2, identical, E, Ip) ! module "CS_integration_limits"
            endif
            if (present(Emax_in)) then
               if (Emax_in < Emin) then
                  Emax = Emin
               elseif (Emax_in < Emax) then
                  Emax = Emax_in
               endif
            endif
            CS = CS - (integrated_delta_CDF_CS(alpha(i), Mc2, E0(i), mtc2, Emax, Ip, E) - integrated_delta_CDF_CS(alpha(i), Mc2, E0(i), mtc2, Emin, Ip, E))   ! below
            ! Prefactor:
            P = P_prefactor(M, E, nat)   ! module "CDF_delta"
         endif ! (E <= Eeq) 
      endif ! (E <= Ip)
   enddo ! i = 1, Nosc
   
   ! Total cross section:
   Sigma = abs(CS)*P   ! [A^2]

   nullify(alpha, E0)
end function Integral_CDF_delta_CS


pure function P_prefactor(M, E, nat) result(P)
   real(8) :: P
   real(8), intent(in) :: M, E  ! mass [kg] and energy [eV] of the projectile
   real(8), intent(in) :: nat   ! [1/cm^3] atomic density
   real(8) :: Mc2, beta, v
   ! beta-factor: v/c
   v = velosity_from_kinetic_energy(E, M) ! module "Relativity"
   beta = beta_factor(v)    ! module "Relativity"
   P = 1.0d24/(g_Pi*g_a0*nat*g_me_eV*(beta*beta)) ! converting into [A^2/eV]
end function P_prefactor


pure function velosity_from_kinetic_energy(Ekin, M0) result(v)
   real(8) v		! velosity  [m/s]
   real(8), intent(in) :: Ekin	! kinetic energy [eV]
   real(8), intent(in) :: M0	! rest mass [kg]
   real(8) :: fact, Erest
   if (M0 < 1.0d-10*g_me) then   ! assume massless particle
      v = g_cvel    ! [m/s]
   else
      Erest = rest_energy(M0)	! [eV] rest energy
      fact = Ekin/Erest + 1.0d0
      v = g_cvel*sqrt(1.0d0 - 1.0d0/(fact*fact))	! [m/s]
   endif
end function velosity_from_kinetic_energy

pure function beta_factor(v) result(beta)
   real(8) beta
   real(8), intent(in) :: v	! [m/s] velosity
   beta = v/g_cvel
end function beta_factor

pure function rest_energy(M0) result(E)
   real(8) E	! [eV] total energy
   real(8), intent(in) :: M0	! [kg] rest mass of the particle
   E = M0*g_cvel*g_cvel/g_e	! [eV]
end function rest_energy

! Limits of transferred energy:
pure function W_min(Ip, Mc2, mtc2, E, E0) result(Wmin)
   real(8) :: Wmin  ! minimal transfered energy [eV]
   real(8), intent(in) :: Ip    ! [eV] ionization potential, or band gap, or the limiting low energy (which may also be zero)
   real(8), intent(in) :: Mc2, mtc2 ! [eV] mass of the incident and target particle
   real(8), intent(in) :: E ! [eV] incident particle energy
   real(8), intent(in), optional :: E0 ! [eV] position of maximum of delta-model CDF
   !-----------------------
   real(8) :: E0min
   real(8) :: mtc2_2
   
   Wmin = Ip

   if (present(E0)) then
      if (abs(Mc2 - mtc2)/Mc2 < 1.0d-6) then   ! e.g. electron-electron
         E0min = E0*(1.0d0-0.25d0*E0/E)    ! delta-model CDF
      else ! e.g. ion-electron
           if (E < 1.0d9*E0) then ! exact expression:
              mtc2_2 = mtc2*mtc2
              E0min = (-2.0d0*Mc2*sqrt(E)*(Mc2*sqrt(E)-sqrt(Mc2*E0*mtc2-E0*mtc2_2+mtc2_2*E))/(Mc2-mtc2)+2.0d0*Mc2*E-E0*mtc2)/(Mc2-mtc2)
           else   ! approximation would do:
              E0min = E0*(1.0d0-0.25d0*E0/E*Mc2/mtc2)    ! delta-model CDF
           endif
      endif
      Wmin = max(Wmin,E0min)
   endif
end function W_min


pure function W_max(Mc2, mtc2, identical, E, Ip, hw_phonon) result(Wmax)
   real(8) :: Wmax  ! maximal transfered energy [eV]
   real(8), intent(in) :: Mc2, mtc2 ! [eV] mass of the incident particles, and the scattering center
   logical, intent(in) :: identical   ! identical particles, true or not
   real(8), intent(in) :: E ! [eV] incident particle energy
   real(8), intent(in) :: Ip    ! [eV] minimum energy (ionization potential)
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency [eV]
   real(8) :: eps, Emtc22, EMmt, Mm, Wmax1
!    eps = 1.0d-6 ! within what margin we consider particles masses equal
   if (identical) then   ! elecron-electron, etc.
      Wmax = (E+Ip)*0.50d0
   else ! not identical particles
      Emtc22 =  2.0d0*mtc2*E
      Mm = Mc2+mtc2
      Wmax = Emtc22*(E + Mc2*2.0d0) / (Emtc22 + Mm*Mm)    ! Eq.(A.25) [1]
      if (present(hw_phonon)) then  ! it's scattering on atoms/phonons, so may be it's a phonon:
         if (Wmax < hw_phonon) Wmax = hw_phonon
      endif
      ! Classical limit:
!       Wmax1 = 4.0d0*Mc2*mtc2/(Mm*Mm)*E
!       write(*,'(a,e,e,e,e,e)') 'Wmax', E, Wmax1, Wmax, Mc2, mtc2
   endif
end function W_max


 subroutine find_Wmax_equal_Wmin(Mc2, mtc2, identical, E, Ip, E0, Eeq, hw_phonon)
   real(8), intent(in) :: Mc2, mtc2 ! [eV] mass of the incident particles, and the scattering center
   logical, intent(in) :: identical   ! identical particles, true or not
   real(8), intent(in) :: E ! [eV] incident particle energy
   real(8), intent(in) :: Ip    ! [eV] minimum energy (ionization potential)
   real(8), intent(in) :: E0    ! [eV]
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency [eV]
   real(8), intent(out) :: Eeq  ! point where Wmax=Wmin
   !----------------------------
   real(8) :: mtc2_2, Mc2_2, Eeq1, Mm, Mm4, adj_coef
  
   if (.not.identical) then ! scattering between two different particles (ion-electron, electron-ion, positron-electron, etc.)
      if (present(hw_phonon)) then  ! scattering on phonon:
         Eeq =0.25d0*E0*(-2.0d0*hw_phonon+E0)/(-hw_phonon+E0)
      else ! scattering on a particle
         mtc2_2 = mtc2*mtc2
         Mc2_2 = Mc2*Mc2         
         Mm = Mc2+mtc2
         Mm4 = Mm*Mm*Mm*Mm
         adj_coef = 1.28d0  ! shift the crossing point by this much
         
         !Eeq = adj_coef*E0*(Mc2/mtc2)*(Mm4/abs(4.0d0*Mc2*Mc2*Mc2*Mc2 - Mm4))
         Eeq = E0*(Mc2/mtc2)*(Mm4/abs(4.0d0*Mc2*Mc2 - Mm4))
         
!            Eeq1 = 0.125d0*E0*Mm*Mm/(Mc2*mtc2)   ! just for comparison
!            write(*,'(a,e,e,e,f)') 'Eeq', E, Eeq, Eeq1, E0
      endif
   else ! scattering between identical particles (electron-electron etc.)
      Eeq = 5.0d0/4.0d0*E0 - Ip/2.0d0 + 0.25d0*sqrt(17.0d0*E0*E0 - 12.0d0*Ip*E0 + 4.0d0*Ip*Ip)
   endif
end subroutine find_Wmax_equal_Wmin
 

! For a given minimal transferred energy, what is ion energy that is sufficient to transfer such an amount:
pure function minimal_sufficient_E(Ip, Mc2, mtc2) result(Emin)
   real(8) :: Emin
   real(8), intent(in) :: Ip    ! [eV] transferred energy
   real(8), intent(in) :: Mc2 ! [eV] ion mass
   real(8), intent(in) :: mtc2  ! [eV] target particle mass (electron)
   real(8) :: Mc2_2, mtc2_2
   
   if (abs(Mc2-mtc2)/Mc2 < 1.0d-6) then ! identical
      Emin = Ip
   else ! not identical
      Mc2_2 = Mc2*Mc2
      mtc2_2 = mtc2*mtc2

      Emin = 0.50d0*(Ip*Mc2 - 2.0d0*mtc2*Mc2 + 2.0d0*Ip*mtc2 + sqrt(-Ip*Ip*Mc2_2 + 4.0d0*mtc2_2*Mc2_2 - &
               6.0d0*mtc2_2*Mc2*Ip + 2.0d0*Ip*Ip*mtc2_2 + 2.0d0*Ip*Mc2*Mc2_2))/(Mc2 - Ip)
   endif
end function minimal_sufficient_E


! Find the coefficients of the linear approximation:
 subroutine Find_linear_a_b(alpha, M, Mc2, mtc2, E0, Ip, identical, nat, Eeq, a, b)
   real(8), intent(in) :: alpha, M, Mc2, mtc2, E0, Ip
   logical, intent(in) :: identical
   real(8), intent(in) :: nat   ! [1/cm^3] atomic density
   real(8), intent(in) :: Eeq   ! [eV] point where W_min = W_max - there we switch from delta-model to the linear extrapolation
   real(8), intent(out) :: a, b  ! coefficients of the linear extrapolation function
   real(8) :: Wmin_lim, Wmax_lim, CS, P, IpMm
   
   Wmin_lim = W_min(Ip, Mc2, mtc2, Eeq, E0) ! module "CS_integration_limits"
   Wmax_lim = W_max(Mc2, mtc2, identical, Eeq, Ip) ! module "CS_integration_limits"
   P = P_prefactor(M, Eeq, nat)   ! below
   CS = -P * (integral_CS(alpha, Mc2, mtc2, E0, Wmax_lim) - integral_CS(alpha, Mc2, mtc2, E0, Wmin_lim))  ! below
   
   IpMm = minimal_sufficient_E(Ip, Mc2, mtc2)    ! module "CS_integration_limits"
   
   a = CS/(Eeq-IpMm)
   b = -CS*Ip/(Eeq-IpMm)
end subroutine Find_linear_a_b

function integrated_delta_CDF_CS(alpha, Mc2, E0, mtc2, W, Ip, E) result(CS)
   real(8) CS   ! a part of the cross section
   real(8), intent(in) :: alpha, Mc2, E0, mtc2, W, Ip, E
   !------------------------------------------
   real(8) :: Mc22, Wmin
   ! Low integration limit:
   Wmin = W_min(Ip, Mc2, mtc2, E, E0) ! module "CS_integration_limits"

   if ((W < Wmin) .or. (E <= Ip)) then ! cannot ionize or transfer energy
      CS = 0.0d0
   else
      !CS = alpha/(g_me_eV * (Mc22 - E0)) * ( (Mc22 - mtc2)*log(Mc22 + W - E0) + Mc22*(mtc2/E0 - 1.0d0)*log(W) - mtc2*(Mc22/E0 - 1.0d0)*log(abs(W - E0)) )
      CS =  integral_CS(alpha, Mc2, mtc2, E0, W)    ! below
   endif
end function integrated_delta_CDF_CS

pure function integral_CS(alpha, Mc2, mtc2, E0, W) result(CS)
   real(8) CS   ! part of the cross section of scattering (not normalized, prefactor missing)
   real(8), intent(in) :: alpha, Mc2, E0, mtc2, W
   real(8) :: Mc22
   Mc22 = 2.0d0*Mc2
   CS = alpha/(g_me_eV * (Mc22 - E0)) * ( (Mc22 - mtc2)*log(Mc22 + W - E0) + Mc22*mtc2/E0*( log(W) - log(abs(W - E0))) )
end function integral_CS


pure function int_energy_loss(alpha, Mc2, mtc2, E0, Wmax, Wmin) result (WCS)
   real(8) :: WCS
   real(8), intent(in) :: alpha, Mc2, mtc2, Wmin, Wmax, E0
   WCS =-alpha*(mtc2*log(abs((Wmax-E0)/(Wmin-E0))) + (2.0d0*Mc2-mtc2)*log(abs((-Wmax+E0-2.0d0*Mc2)/(-Wmin+E0-2.0d0*Mc2))))/g_me_eV
end function int_energy_loss


function energy_loss_delta(E, M, Zeff, Ip, nat, Mt, CDF_coefs, identical, hw_phonon) result(Se)
   real(8) :: Se    ! [eV/A] energy loss
   real(8), intent(in) :: E, M, Zeff, Ip, nat, Mt
   type(CDF), intent(in), target :: CDF_coefs	! CDF coefficients
   logical, intent(in) :: identical   ! identical particles, true or not
   real(8), intent(in), optional :: hw_phonon   ! [eV]
   !---------------------------------------
   real(8) :: Mc2, mtc2, S_cur, P, Mc22, Emin, Emax, Eeq, dEed, Estart
   real(8) :: a, b  ! coefficients of the linear extrapolation function
   real(8), dimension(:), pointer :: alpha, E0
   integer :: i, Nosc
   
   ! How many delta-functions:
   Nosc = size(CDF_coefs%E0(:))
   
   Mc2 = rest_energy(M)   ! module "Relativity"
   mtc2 = rest_energy(mt)   ! module "Relativity"
   ! Parameters of the delta-functions CDF:
   E0 => CDF_coefs%E0
   alpha => CDF_coefs%alpha

   S_cur = 0.0d0   ! to start with
   do i = 1, Nosc
      Emin = W_min(Ip, Mc2, mtc2, E, E0(i)) ! module "CS_integration_limits"

      if (abs(M - mt)/M < 1.0d-6) then  ! identical particles:
         Estart = Ip
      else  ! non-identical
         Estart = minimal_sufficient_E(Ip, Mc2, mtc2) ! module "CS_integration_limits"
      endif
      
      ! Where to switch from the approximate formula to linear extrapolation:
      call  find_Wmax_equal_Wmin(Mc2, mtc2, identical, E, Ip, E0(i), Eeq)   ! module "CS_integration_limits"
      
      if (E <= Estart) then ! cannot ionize
         S_cur = 0.0d0
         P = 0.0d0
      else
         dEed = Eeq/100.0d0    ! shift a little the point at which we switch to linear extrapolation
         
         if (present(hw_phonon)) then  ! scattering on phonons:
            Emax = W_max(Mc2, mtc2, identical, E, Ip, hw_phonon) ! module "CS_integration_limits"   
         else  ! scattering on electrons or atoms:
            Emax = W_max(Mc2, mtc2, identical, E, Ip) ! module "CS_integration_limits"
         endif

         if (E <= Eeq+dEed) then   ! delta-CDF model does not apply here, replace the CDF with the linear extrapolation cross section
            ! Find the coefficients of the linear approximation:
            call Find_linear_a_b(alpha(i), M, Mc2, mtc2, E0(i), Ip, identical, nat, Eeq+dEed, a, b) ! below
            S_cur = 0.5d0*a*(Emax - Estart)    ! linearly approximated cross section
            P = 1.0d0   ! no prefactor needed in this case, it alreaxy was included into the linear extrapolation
         else   ! delta-model works
            !S_cur = S_cur - (integral_CS_x_W(alpha(i), Mc2, E0(i), Estart, mtc2, Emax) - integral_CS_x_W(alpha(i), Mc2, E0(i), Estart, mtc2, Emin)) ! below
            S_cur = S_cur -  int_energy_loss(alpha(i), Mc2, mtc2, E0(i), Emax, Emin)  ! below
            ! Prefactor:
            P = P_prefactor(M, E, nat)   ! below
         endif ! (E <= Eeq) 
      endif ! (E <= Ip)
   enddo ! i = 1, Nosc
   
   ! Energy loss including effective charge:
   Se = abs(S_cur)*(Zeff*Zeff)*P*nat*1.0d-24   ! [eV/A]
   
   nullify(alpha, E0)   
end function energy_loss_delta


pure function MFP_from_sigma(sigma, nat) result(lambda)
   real(8) lambda   ! [1/A]
   real(8), intent(in) :: sigma ! [A^2]
   real(8), intent(in) :: nat   ! [cm^-3]
   real(8) :: na
   if (sigma > 1.0d-24) then  ! finite MFP
      na = nat * 1.0d-24    ! [A^-3] converted from [cm^-3]
      lambda = 1.0d0 / (sigma * na)    ! [A] mean free path
   else
      lambda = 1.0d30   ! infinity
   endif
end function MFP_from_sigma




subroutine Electron_energy_transfer_inelastic(Ele, Target_atoms, Nat, Nshl, L_tot, dE_out, Matter, Mat_DOS, NumPar, kind_of_particle)
    real(8), intent(in) :: Ele  ! electron energy [eV]
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    integer, intent(in) :: Nat, Nshl    ! number of atom, number of shell
    real(8), intent(in) :: L_tot    ! [A] total mean free path
    real(8), intent(out) :: dE_out   ! the transferred energy [eV]
    type(Solid), intent(in) :: Matter
    type(Density_of_states), intent(in) :: Mat_DOS
    type(Flag), intent(in) :: NumPar
    character(8), intent(in) :: kind_of_particle

    real(8) Emin, Emax, E, RN, L_need, L_cur, Eplasmon, Egap
    real(8) Mass
    integer Mnum
    
    call random_number(RN)    
    L_need = L_tot/RN   ! [A] we need to reach
    
    Emin = Target_atoms(Nat)%Ip(Nshl)   ! [eV] ionization potential of the shell is minimum possible transferred energy
    Egap = Target_atoms(1)%Ip(size(Target_atoms(1)%Ip))   ! band gap [eV]
    if (Emin .LE. 1.0d-3) Emin = 1.0d-3 ! for metals there is no band gap
      
    if (trim(adjustl(kind_of_particle)) .EQ. 'Electron') then 
        Emax = (Ele + Emin)/2.0e0 ! [eV] maximum energy transfer, accounting for equality of electrons
        Mass = 1.0d0
    else if (trim(adjustl(kind_of_particle)) .EQ. 'Hole') then 
        if(Matter%hole_mass .GE. 0) then
            Mass = Matter%hole_mass
        else                    ! Define mass from DOS
            call find_in_array_monoton(Mat_DOS%E, Ele, Mnum)
            Mass = Mat_DOS%Eff_m(Mnum)
        endif
        Emax = 4.0d0*Ele*Mass/((Mass+1.0d0)*(Mass+1.0d0))
    else ! just for case, default value:
        Emax = (Ele + Emin)/2.0e0 ! [eV] maximum energy transfer, accounting for equality of electrons
        Mass = 1.0d0
    endif
    
!   if (trim(adjustl(NumPar%kind_of_particle)) .EQ. 'Electron') then 
!        if (Mass .GT. 1.0d0) then
!            print*, "ooooooooooooooooooooooooooooooooooooooo"
!            print*, Mass, NumPar%kind_of_particle
!            print*, "ooooooooooooooooooooooooooooooooooooooo"
!            pause
!        endif
!   endif
    
    
    ! Use maximal plasmon energy as upper limit of integration
    if (Numpar%plasmon_Emax) then       ! If included
        if (Emin .EQ. Egap) then        ! For VB only
            Eplasmon = sqrt(Matter%N_VB_el*Matter%At_Dens*1d6*g_h*g_h/(g_me*g_e0) + Egap*Egap)    ! Maximal energy of plasmons
            if (Eplasmon .GE. Emax) Emax = Eplasmon ! single atom vs plasmon
            if (Ele .LT. Emax) Emax = Ele ! no more than the total electron energy
        endif
    endif
    
    select case (Target_atoms(Nat)%KOCS(Nshl)) ! which inelastic cross section to use (BEB vs CDF):
    case (1) ! CDF cross section energy transfer:
        call Electron_NRG_transfer_CDF(Ele, Target_atoms, Nat, Nshl, L_need, E, Matter, Mat_DOS, NumPar, Mass, Emin, Emax)
    case default ! do BEB energy transfer:
        call Electron_NRG_transfer_BEB(Ele, Target_atoms, Nat, Nshl, L_need, E, Matter, Mass, Emin, Emax)
    end select
    if (E .GE. Emax) E = Emax
    dE_out = E  ! energy transfer [eV]
    
    !if (trim(adjustl(NumPar%kind_of_particle)) .EQ. 'Electron') then 
        if (dE_out .LT. Egap) then
             print*, "###########################################"
             print*, 'ERROR in Electron_energy_transfer_inelastic:'
             print*, dE_out, Egap, Emin, Emax, Ele
             print*, Mass, kind_of_particle
             print*, "###########################################"
             pause
        endif
    !endif
    
end subroutine Electron_energy_transfer_inelastic


subroutine Electron_NRG_transfer_CDF(Ele, Target_atoms, Nat, Nshl, L_need, E, Matter, Mat_DOS, NumPar, Mass, Emin, Emax)
    real(8), intent(in) :: Ele  ! electron energy [eV]
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    integer, intent(in) :: Nat, Nshl    ! number of atom, number of shell
    real(8), intent(in) :: L_need    ! [A] total mean free path
    real(8), intent(out) :: E   ! the transferred energy [eV]
    type(Solid), intent(in) :: Matter
    type(Density_of_states), intent(in) :: Mat_DOS
    type(Flag), intent(in) :: NumPar
    real(8), intent(in) :: Mass, Emin, Emax ! electron/hole mass; min and max energies for integration
    real(8) L_cur, Ltot0, Ltot1, a, b, temp1, temp2, dE, dL
    integer n, i

    
 if (NumPar%kind_of_DR .EQ. 4) then    ! Delta-CDF
    E = get_inelastic_energy_transfer(Ele, Matter, Target_atoms, numpar, Nat, Nshl, Emin)  ! below
 else
    n = 10*(MAX(INT(Emin),10))    ! number of integration steps
    dE = (Emax - Emin)/(real(n)) ! differential of transferred momentum [kg*m/s]
    i = 1       ! to start integration
    E = Emin    ! to start integration
    Ltot1 = 0.0d0
    Ltot0 = 0.0d0
    L_cur = 1.0d10
    call Diff_cross_section(Ele, E, Target_atoms, Nat, Nshl, Ltot0, Mass, Matter, Mat_DOS, NumPar)
    !ddEdx = 0.0d0

    do while (L_cur .GT. L_need)
        dE = (1.0d0/(E+1.0d0) + E)/real(n)
        ! If it's Simpson integration:
        a =  E + dE/2.0d0
        call Diff_cross_section(Ele, a, Target_atoms, Nat, Nshl, dL, Mass, Matter, Mat_DOS, NumPar)
        temp1 = dL
        b = E + dE
        call Diff_cross_section(Ele, b, Target_atoms, Nat, Nshl, dL, Mass, Matter, Mat_DOS, NumPar)
        temp2 = dE/6.0d0*(Ltot0 + 4.0d0*temp1 + dL)
        Ltot1 = Ltot1 + temp2
        !ddEdx = ddEdx + E*temp2
        Ltot0 = dL
        L_cur = 1.0d0/(Mass*Ltot1)
        E = E + dE  ! [eV]
        if (E .GE. Emax) exit
    enddo
 endif
end subroutine Electron_NRG_transfer_CDF




! Interface to select the model of electron inelastic scattering cross section:
function get_inelastic_energy_transfer(Ee, Matter, Target_atoms, numpar, j, k, Ip, CDF_dispers, CDF_m_eff, hw_phonon) result(dE)
   real(8) :: dE   ! [eV] sampled transferred energy
   real(8), intent(in) :: Ee	! [eV] electron kinetic energy
   type(Solid), intent(in) :: Matter
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
   type(Flag), intent(in) :: NumPar
   integer, intent(in) :: j, k	! number or element and its shell
   real(8), intent(in) :: Ip    ! [eV] ionization potential or band gap
   integer, intent(in), optional :: CDF_dispers, CDF_m_eff	! dispersion relation and effective mass
   real(8), intent(in), optional :: hw_phonon   ! maximal phonon frequency, for particle-phonon scattering
   !---------------------------------------
   real(8) :: eps, Zeff, max_E0, Eeq, RN, CS_sampled, CS_cur, E_left, E_right, E_cur
   real(8) :: CS_tot	! [A^2] precalculated total cross section
   integer :: dispers, m_eff, El_inelast
   ! Set accepteble margin of precision for the angle:
   eps = 1.0d-3
   
   ! Set the model parameters:
   Zeff = 1.0d0	! electron charge
!    El_inelast = numpar%El_inelast   ! to chose the model of CS calculations below
!    if (present(CDF_dispers)) then	! user provided alternative (e.g. for all but VB)
!       dispers = CDF_dispers
!    else	! for the VB
!       dispers = numpar%CDF_dispers
!    endif
!    if (present(CDF_m_eff)) then	! user provided alternative (e.g. for all but VB)
!       m_eff = CDF_m_eff
!    else
!       m_eff = numpar%hole_mass
!    endif

!    VAL0:if ( (Material%Elements(j)%valent(k)) .and. (allocated(Material%CDF_valence%A)) ) then    ! Valence band
!       max_E0 = maxval(Material%CDF_valence%E0(:))
!    else VAL0
!       max_E0 = maxval(Material%Elements(j)%CDF(k)%E0(:))
!    endif VAL0
   max_E0 = maxval(Target_atoms(j)%Ritchi(k)%E0(:))
   call find_Wmax_equal_Wmin(0.0d0, 0.0d0, .true., Ee, Ip, max_E0, Eeq)   ! module "CS_integration_limits"
   
   E_left = Ip ! [eV] minimal transferred energy
   E_right = (Ip + Ee)*0.5d0    ! [eV] maximal transferred energy
   ! Use of the maximal plasmon energy as the upper limit of integration is not included yet !!!         
   
   ! Sample the cross section:
   call random_number(RN)
!    CS_tot = get_integral_inelastic_CS(Ee, Material, DOS, numpar, j, k, Ip, El_inelast, Zeff, dispers, m_eff, E_right)   ! below
!    CS_tot = CDF_total_CS_delta(Ee, g_me, Zeff, Ip, Material%At_Dens, Material%CDF_valence, g_me, .true., Emax_in = E_right)    ! module "CDF_delta"
   CS_tot = Integral_CDF_delta_CS(g_me, g_me, Ee, Target_atoms(j)%Ritchi(k), Ip, Matter%At_Dens, .true., Emax_in = E_right) ! [A^2] Below
   CS_sampled = RN*CS_tot
   
   ! Start finding CS:
   E_cur = (E_left + E_right)*0.5d0
!    CS_cur = get_integral_inelastic_CS(Ee, Material, DOS, numpar, j, k, Ip, El_inelast, Zeff, dispers, m_eff, E_cur)   ! below
!    CS_cur = CDF_total_CS_delta(Ee, g_me, Zeff, Ip, Material%At_Dens, Material%CDF_valence, g_me, .true., Emax_in = E_cur)    ! module "CDF_delta"
   CS_cur = Integral_CDF_delta_CS(g_me, g_me, Ee, Target_atoms(j)%Ritchi(k), Ip, Matter%At_Dens, .true., Emax_in = E_cur) ! [A^2] Below
   
   ! Search by bisection method:
   do while (ABS(CS_cur - CS_sampled)/CS_sampled > eps)
      if (CS_cur > CS_sampled) then
         E_right = E_cur
      else
         E_left = E_cur
      endif
      E_cur = (E_left + E_right)/2.0d0
      if (abs(E_left - E_right) < eps) exit  ! precise enough
!       CS_cur = get_integral_inelastic_CS(Ee, Material, DOS, numpar, j, k, Ip, El_inelast, Zeff, dispers, m_eff, E_cur)   ! below
!       CS_cur = CDF_total_CS_delta(Ee, g_me, Zeff, Ip, Material%At_Dens, Material%CDF_valence, g_me, .true., Emax_in = E_cur)    ! module "CDF_delta"
      CS_cur = Integral_CDF_delta_CS(g_me, g_me, Ee, Target_atoms(j)%Ritchi(k), Ip, Matter%At_Dens, .true., Emax_in = E_cur) ! [A^2] Below
   enddo
   
   ! Output: sampled transferred energy:
   dE = E_cur
end function get_inelastic_energy_transfer




subroutine Electron_NRG_transfer_BEB(Ele, Target_atoms, Nat, Nshl, L_need, E, Matter, Mass, Emin, Emax)
    real(8), intent(in) :: Ele  ! electron energy [eV]
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    integer, intent(in) :: Nat, Nshl    ! number of atom, number of shell
    real(8), intent(in) :: L_need    ! [A] total mean free path
    real(8), intent(out) :: E   ! the transferred energy [eV]
    type(Solid), intent(in) :: Matter
    real(8), intent(in) :: Mass, Emin, Emax ! electron/hole mass; min and max energies for integration
    real(8) L_cur, Emin1, Emax1, Sigma_cur, temp1
    integer coun
    
   Emin1 = Emin
   Emax1 = (Ele - Target_atoms(Nat)%Ip(Nshl))/2.0d0
   coun = 0
   E = 0.0d0    ! to start iterations
   L_cur = 1.0d10
   Sigma_cur = dSigma_int_BEB(Ele, E, Target_atoms(Nat)%Ip(Nshl), Target_atoms(Nat)%Ek(Nshl), Target_atoms(Nat)%Nel(Nshl))
   temp1 = Matter%At_Dens*1d-24*(Target_atoms(Nat)%Pers)/SUM(Target_atoms(:)%Pers)
   L_cur = 1.0d0/(Mass*temp1*Sigma_cur) ! IMFP [A]
   do while (ABS(L_cur-L_need)/L_need .GT. 0.001d0) ! search the transferred energy by bisection:
     coun = coun + 1	! just count the loops
     Sigma_cur = dSigma_int_BEB(Ele, E, Target_atoms(Nat)%Ip(Nshl), Target_atoms(Nat)%Ek(Nshl), Target_atoms(Nat)%Nel(Nshl))
     temp1 = Matter%At_Dens*1d-24*(Target_atoms(Nat)%Pers)/SUM(Target_atoms(:)%Pers)
     L_cur = 1.0d0/(Mass*temp1*Sigma_cur) ! IMFP [A]
     if (L_cur .GT. L_need) then
        Emin1 = E
     else
        Emax1 = E
     endif
     E = (Emax1 + Emin1)/2.0d0
     !write(*,'(a,f,f,f,f,f)') 'T:', Emax1, Emin1, L_need, L_cur, L_tot
     if (coun .GE. 1d3) then 
        write(*,'(a,e,e,e,e,e,e)') 'TOO MANY LOOPS IN BEB-BISECTION:', Ele, Emax1, Emin1, L_need, L_cur
        exit
     endif
   enddo
   E = E + Target_atoms(Nat)%Ip(Nshl) ! For BEB it's counted like that
end subroutine Electron_NRG_transfer_BEB


subroutine SHI_NRG_transfer_BEB(SHI, Target_atoms, Nat, Nshl, L_need, Matter, E)
    type(Ion), intent(inout) :: SHI ! all SHI parameters
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    integer, intent(in) :: Nat, Nshl    ! number of atom, number of shell
    real(8), intent(in) :: L_need    ! [A] total mean free path
    real(8), intent(out) :: E   ! the transferred energy [eV]
    type(Solid), intent(in) :: Matter
    real(8) L_cur, Emin1, Emax1, Sigma_cur, temp1, Ele, MSHI, Zeff, Emin, Emax
    integer coun
    
    Ele = SHI%E  ! energy of the particle [eV]
    MSHI = g_Mp*SHI%Mass    ! SHI mass [kg]
    Zeff = SHI%Zeff   ! SHI effective charge (must be precalculated)
    Emin = Target_atoms(Nat)%Ip(Nshl) ! min transferred energy [eV]
    if (Emin .LE. 1.0d-3) Emin = 1.0d-3 ! for metals there is no band gap, but we need a finite value
    Emax = 4.0d0*Ele*g_me*MSHI/((MSHI+g_me)*(MSHI+g_me)) ! [eV] maximum energy transfer of SHI

   Emin1 = Emin
   !Emax1 = (Ele - Target_atoms(Nat)%Ip(Nshl))/2.0d0
   Emax1 = Emax
   coun = 0
   E = 0.0d0    ! to start iterations
   L_cur = 1.0d10
   !Sigma_cur = dSigma_int_BEB(Ele, E, Target_atoms(Nat)%Ip(Nshl), Target_atoms(Nat)%Ek(Nshl), Target_atoms(Nat)%Nel(Nshl))
   Sigma_cur = dSigma_int_BEB_SHI(Ele, E, Target_atoms(Nat)%Ip(Nshl), Target_atoms(Nat)%Ek(Nshl), Target_atoms(Nat)%Nel(Nshl), MSHI, Zeff)
   temp1 = Matter%At_Dens*1d-24*(Target_atoms(Nat)%Pers)/SUM(Target_atoms(:)%Pers)
   L_cur = 1.0d0/(temp1*Sigma_cur) ! IMFP [A]
   do while (ABS(L_cur-L_need)/L_need .GT. 0.001d0) ! search the transferred energy by bisection:
     coun = coun + 1	! just count the loops
     !Sigma_cur = dSigma_int_BEB(Ele, E, Target_atoms(Nat)%Ip(Nshl), Target_atoms(Nat)%Ek(Nshl), Target_atoms(Nat)%Nel(Nshl))
     Sigma_cur = dSigma_int_BEB_SHI(Ele, E, Target_atoms(Nat)%Ip(Nshl), Target_atoms(Nat)%Ek(Nshl), Target_atoms(Nat)%Nel(Nshl), MSHI, Zeff)
     temp1 = Matter%At_Dens*1d-24*(Target_atoms(Nat)%Pers)/SUM(Target_atoms(:)%Pers)
     L_cur = 1.0d0/(temp1*Sigma_cur) ! IMFP [A]
     if (L_cur .GT. L_need) then
        Emin1 = E
     else
        Emax1 = E
     endif
     E = (Emax1 + Emin1)/2.0d0
     !write(*,'(a,f,f,f,f,f)') 'T:', Emax1, Emin1, L_need, L_cur, L_tot
     if (coun .GE. 1d3) then 
        write(*,'(a,e,e,e,e,e,e)') 'TOO MANY LOOPS IN BEB-BISECTION:', Ele, Emax1, Emin1, L_need, L_cur
        exit
     endif
   enddo
   E = E + Target_atoms(Nat)%Ip(Nshl) ! For BEB it's counted like that
end subroutine SHI_NRG_transfer_BEB


subroutine Diff_cross_section(Ele, hw, Target_atoms, Nat, Nshl, Diff_IMFP, Mass, Matter, Mat_DOS, NumPar)
    real(8), intent(in), target :: Ele  ! SHI energy [eV]
    REAL(8), INTENT(in), target :: hw   ! transferred energy [eV]
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    integer, intent(in) :: Nat, Nshl    ! number of atom, and number of shell    
    real(8), intent(out) :: Diff_IMFP   ! differential inverse mean free path 1/lambda(Ele,dE)
    real(8), intent(in) :: Mass
    type(Solid), intent(in) :: Matter
    type(Density_of_states), intent(in) :: Mat_DOS
    type(Flag), intent(in) :: NumPar
        
    integer i, n
    real(8), pointer :: Ee, dE
    real(8) dLs, qmin, qmax, hq, ddq, dq, Ime, dLs0, dL, hq0, dq_save
    real(8) a, b, x, temp1, temp2, T_fact
    Ee => Ele        ! energy [eV]
    dE => hw         ! transferred energy [eV]
    
    if (dE > Ee) then
       qmin = sqrt(2.0d0*Mass*g_me)/g_h*(sqrt(Ee))
       qmax = sqrt(2.0d0*Mass*g_me)/g_h*(sqrt(Ee))
    else
       qmin = sqrt(2.0d0*Mass*g_me)/g_h*(sqrt(Ee) - sqrt((Ee - dE))) ! min transferred momentum [sqrt(eV/J) /m] (not [kg*m/s])
       qmax = sqrt(2.0d0*Mass*g_me)/g_h*(sqrt(Ee) + sqrt((Ee - dE))) ! max transferred momentum [sqrt(eV/J) /m] (not [kg*m/s])
    endif
!     print*, 'Mass', Mass, Ee, dE
!     print*, g_h*g_h*qmin*qmin/(2.0d0*Mass*g_me),  g_h*g_h*qmax*qmax/(2.0d0*Mass*g_me)          
!     pause 
    
    dLs = 0.0d0 ! starting integration, mean free path per energy [A/eV]^(-1)
    hq = qmin    ! transient transferred momentum for integration [kg*m/s]
    n = 100
    dq = (qmax - qmin)/real(n) ! differential of transferred momentum [kg*m/s]
    dLs0 = 0.0d0
    do while (hq .LT. qmax) ! no matter how many points, go till the end
        dq = hq/real(n)
        
        ! If it's Simpson integration:
        a = hq + dq/2.0d0
        call Imewq(hw, a, Target_atoms, Nat, Nshl, ImE, Matter, Mat_DOS, NumPar)
        temp1 = ImE
        b = hq + dq
        call Imewq(hw, b, Target_atoms, Nat, Nshl, ImE, Matter, Mat_DOS, NumPar)
        dL = ImE
        dLs = dLs + dq/6.0d0*(dLs0 + 4.0d0*temp1 + dL)/hq
        dLs0 = dL
!         if ((abs(Ee - 10.0) < 1.0d-2) .and. (abs(dE - 4.0d0) < 1.0d-2)) then
!            write(*,'(f,f,f,f,f)') g_h*g_h*hq*hq/(2.0d0*Mass*g_me), dL/( g_h*g_h*hq*hq/(2.0d0*Mass*g_me)), dLs
!         endif
        hq = hq + dq
    enddo
!     if (abs(Ee - 10.0) < 1.0d-2) then
!        write(*,'(f,f,f,f,f)') Ee, dE , dLs
!     endif
   
    if (Matter%temp > 0.0d0) then   ! non-zero temperature
       T_fact = 1.0d0/(1.0d0-exp(-hw/Matter%temp*g_kb))
    else    ! zero temeprature
       T_fact = 1.0d0
    endif
    Diff_IMFP = 1.0d0/(g_Pi*g_a0*Ele)*dLs*T_fact
    nullify(Ee)
    nullify(dE)
end subroutine Diff_cross_section


subroutine Electron_energy_transfer_elastic(Ele, L_tot, Target_atoms, CDF_Phonon, Matter, dE_out, NumPar, Mat_DOS, kind_of_particle)
    real(8), intent(in), target :: Ele  ! electron energy [eV]
    real(8), intent(in) :: L_tot    ! total mean free path [A]
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    type(CDF), intent(in) :: CDF_Phonon ! phonon CDF parameters
    type(Solid), intent(in) :: Matter   ! all material parameters
    real(8), intent(out) :: dE_out   ! calculated inverse mean free path (cross-section) [A], and the energy losses [eV/A]
    type(Flag), intent(in) :: NumPar
    type(Density_of_states), intent(in) :: Mat_DOS
    character(8), intent(in) :: kind_of_particle
    
    integer i, j, n, Mnum
    real(8) Emin, Emax, E, dE, dL, Ltot1, Ltot0, ddEdx, a, b, RN, temp1, temp2, qdebay, Edebay, Mtarget, L_cur, L_need, Mass
    real(8) :: Zt, Zeff
    real(8), pointer :: Ee
    
    call random_number(RN)
    L_need = L_tot/RN   ! [A] we need to reach
    
    qdebay = (6.0d0*g_Pi*g_Pi*(Matter%At_Dens*1e6))**(0.33333333d0)   ! maximum energy of phonons analyzed [1/m], debay energy
    Edebay = 2.0d0*g_h*Matter%Vsound*qdebay/g_e  ! maximum energy of phonons analyzed [eV], debay energy // 2 is artificial...
    Ee => Ele   ! just to use this name later
    Emin = 0.1d-8    !Target_atoms(Nat)%Ip(Nshl)   ! [eV] ionization potential of the shell is minimum possible transferred energy

    Mtarget = g_Mp*SUM(target_atoms(:)%Mass*dble(target_atoms(:)%Pers))/dble(SUM(target_atoms(:)%Pers)) ! average mass of target atom [kg]
    
    if (trim(adjustl(kind_of_particle)) .EQ. 'Electron') then
        Mass = 1.0d0 
    else if (trim(adjustl(kind_of_particle)) .EQ. 'Hole') then
        if(Matter%hole_mass .GE. 0) then
            Mass = Matter%hole_mass
        else                    ! Define mass from DOS
            call find_in_array_monoton(Mat_DOS%E, Ele, Mnum)
            Mass = Mat_DOS%Eff_m(Mnum)
        endif
    else    ! default value just in case
        Mass = 1.0d0 
    endif

    ! Target mean atomic charge:
    if (NumPar%CDF_elast_Zeff /= 1) then   ! Barkas-like charge
        Zt = SUM(target_atoms(:)%Zat*dble(target_atoms(:)%Pers))/dble(SUM(target_atoms(:)%Pers)) ! mean atomic number of target atoms
        Zeff = 1.0d0 + Equilibrium_charge_Target(Ee, g_me, Zt, (Zt-1.0e0), 0, 1.0e0) ! Equilibrium charge, see below
    else  ! one, as used in old CDF expression
        Zeff = 1.0d0    ! electron charge
    endif
        
    Emax =  4.0e0*Ee*Mass*g_me*Mtarget/((Mtarget+Mass*g_me)*(Mtarget+Mass*g_me))    ! [eV] maximum energy transfer
    
    if (Edebay .GE. Emax) Emax = Edebay ! single atom vs phonon
    if (Ee .LT. Emax) Emax = Ee ! no more than the total electron energy

    n = 10*(MAX(INT(Emin),10))    ! number of integration steps
    dE = (Emax - Emin)/(real(n)) ! differential of transferred energy
    i = 1       ! to start integration
    E = Emin    ! to start integration
    Ltot1 = 0.0d0
    Ltot0 = 0.0d0
    L_cur = 1.0d10
    call Diff_cross_section_phonon(Ele, E, CDF_Phonon, Ltot0, Mtarget, Mass, Matter%temp, 1.0d0)
    ddEdx = 0.0e0
    do while (L_cur .GT. L_need)
        dE = (1.0d0/(E+1.0d0) + E)/real(n)
        ! If it's Simpson integration:
        a =  E + dE/2.0d0
        call Diff_cross_section_phonon(Ele, a, CDF_Phonon, dL, Mtarget, Mass, Matter%temp, 1.0d0)
        temp1 = dL
        b = E + dE
        call Diff_cross_section_phonon(Ele, b, CDF_Phonon, dL, Mtarget, Mass, Matter%temp, 1.0d0)
        temp2 = dE/6.0d0*(Ltot0 + 4.0d0*temp1 + dL)
        Ltot1 = Ltot1 + temp2
        ddEdx = ddEdx + E*temp2
        Ltot0 = dL
        !L_cur = 1.0d0/Mass/Ltot1
        L_cur = 1.0d0/(Zeff*Zeff*Mass*Ltot1) ! include effective charge of target atoms
        E = E + dE  ! [eV]
        if (E .GE. Emax) exit
    enddo
    if (E .GE. Emax) E = Emax
    dE_out = E ! energy transferred [eV]
    nullify(Ee)
end subroutine Electron_energy_transfer_elastic


subroutine SHI_Total_IMFP(SHI, Target_atoms, Nat, Nshl, Sigma, dEdx, Matter, Mat_DOS, NumPar)
    class(Ion), intent (inout) :: SHI  ! all the data for the SHI
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    integer, intent(in) :: Nat, Nshl    ! number of atom, number of shell
    real(8), intent(out) :: Sigma, dEdx   ! calculated inverse mean free path (cross-section) [1/A], and the energy losses [eV/A]
    type(Solid), intent(in) :: Matter
    type(Density_of_states), intent(in) :: Mat_DOS
    type(Flag), intent(in) :: NumPar
    
    select case (SHI%Kind_ion)
    case (1)     ! use Brandt-Kitagawa charge for SHI
       call SHI_TotIMFP_BK(SHI, Target_atoms, Nat, Nshl, Sigma, dEdx, Matter, Mat_DOS, NumPar)
    case default ! use point-like charge for SHI
       call SHI_TotIMFP(SHI, Target_atoms, Nat, Nshl, Sigma, dEdx, Matter, Mat_DOS, NumPar)
    endselect
end subroutine SHI_Total_IMFP


subroutine SHI_TotIMFP(SHI, Target_atoms, Nat, Nshl, Sigma, dEdx, Matter, Mat_DOS, NumPar, dSedE)
    class(Ion), intent (inout) :: SHI  ! all the data for the SHI
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    integer, intent(in) :: Nat, Nshl    ! number of atom, number of shell
    real(8), intent(out) :: Sigma, dEdx   ! calculated inverse mean free path (cross-section) [1/A], and the energy losses [eV/A]
    type(All_MFP), dimension(:), allocatable, optional :: dSedE ! SHI integral of d(MFP)/d(dE) [eV, A/eV]
    real(8) :: Ele      ! energy [eV]
    real(8) :: M_SHI    ! atomic mass of SHI
    real(8) :: Z_SHI    ! nuclear charge of SHI
    type(Solid), intent(in) :: Matter
    type(Density_of_states), intent(in) :: Mat_DOS
    type(Flag), intent(in) :: NumPar
    
    integer i, j, i0, j0, k0, n, n0, Nmax
    real(8) Emin, Emax, E, E0, dE, dL, Ltot1, Ltot0, ddEdx, a, b, temp1, temp2, MSHI, Zeff
    real(8) Egap, Eplasmon
    
    call Equilibrium_charge_SHI(SHI, Target_atoms)  ! get Barcas' equilibrium charge from module Cross_sections
    Ele = SHI%E  ! energy of the particle [eV]
    MSHI = g_Mp*SHI%Mass    ! SHI mass [kg]
    Zeff = SHI%Zeff   ! SHI effective charge (must be precalculated)
    
    Egap = Target_atoms(1)%Ip(size(Target_atoms(1)%Ip))   ! band gap [eV]
    Emin = Target_atoms(Nat)%Ip(Nshl) ! min transferred energy [eV]
    if (Emin .LE. 1.0d-3) Emin = 1.0d-3 ! for metals there is no band gap, but we need a finite value
    Emax = 4.0d0*Ele*g_me*MSHI/((MSHI+g_me)*(MSHI+g_me)) ! [eV] maximum energy transfer of SHI

    if (Numpar%plasmon_Emax) then       ! Use maximal plasmon energy as upper limit of integration
        if (Emin .EQ. Egap) then
            Eplasmon = sqrt(Matter%N_VB_el*Matter%At_Dens*1d6*g_h*g_h/(g_me*g_e0) + Egap*Egap)    ! [g/cm^3] material density
            if (Eplasmon .GE. Emax) Emax = Eplasmon ! single atom vs plasmon
        endif
    endif 
    
    Nmax = 50
    n = 10*(MAX(INT(Emin),Nmax))    ! number of integration steps
    i = 0       ! to start integration
    E = Emin    ! to start integration
    
    loss_f:if (present(dSedE)) then    ! first, count how many lines:
        if (.not. allocated(dSedE)) then
            allocate(dSedE(size(Target_atoms))) ! nmumber of atoms
            do i0 = 1, size(Target_atoms)   ! number of atoms
                allocate(dSedE(i0)%ELMFP(size(Target_atoms(i0)%Ip))) ! number of atoms
                do j0 = 1, size(Target_atoms(i0)%Ip) ! number of shells
                    n0 = 10*(MAX(INT(Target_atoms(i0)%Ip(j0)),Nmax))
                    E0 = Target_atoms(i0)%Ip(j0)    ! [eV] to start counting
                    dE = (Emax - E0)/real(n0) ! differential of transferred momentum [kg*m/s]
                    k0 = 0
                    do while (E0 .LE. Emax) ! integration
                        k0 = k0 + 1
                        dE = (1.0d0/(E0+1.0d0) + E0)/real(n0)
                        E0 = E0 + dE  ! [eV]
                    enddo
                    allocate(dSedE(i0)%ELMFP(j0)%E(k0))
                    dSedE(i0)%ELMFP(j0)%E = 0.0d0
                    allocate(dSedE(i0)%ELMFP(j0)%L(k0))
                    dSedE(i0)%ELMFP(j0)%L = 0.0d0
                    allocate(dSedE(i0)%ELMFP(j0)%dEdx(k0))
                    dSedE(i0)%ELMFP(j0)%dEdx = 0.0d0
                enddo
            enddo
            dE = (Emax - Emin)/(real(n)) ! differential of transferred momentum [kg*m/s]
        endif
    endif loss_f

    select case (Target_atoms(Nat)%KOCS_SHI(Nshl)) ! which inelastic cross section to use (BEB vs CDF):
    case (1) ! CDF cross section
       Ltot1 = 0.0d0
       Ltot0 = 0.0d0
       call SHI_Diff_cross_section(Ele, MSHI, Emax, E, Target_atoms, Nat, Nshl, Ltot0, Matter, Mat_DOS, NumPar)
       ddEdx = 0.0d0
       do while (E <= Emax) ! integration
          i = i + 1
          if (present(dSedE)) then
             if (i > size(dSedE(Nat)%ELMFP(Nshl)%E)) exit
          endif
          dE = (1.0d0/(E+1.0d0) + E)/real(n)
          ! If it's Simpson integration:
          a =  E + dE/2.0d0
          call SHI_Diff_cross_section(Ele, MSHI, Emax, a, Target_atoms, Nat, Nshl, dL, Matter, Mat_DOS, NumPar)
          temp1 = dL
          b = E + dE
          call SHI_Diff_cross_section(Ele, MSHI, Emax, b, Target_atoms, Nat, Nshl, dL, Matter, Mat_DOS, NumPar)
          temp2 = dE/6.0d0*(Ltot0 + 4.0d0*temp1 + dL)
          Ltot1 = Ltot1 + temp2
          ddEdx = ddEdx + dE/6.0d0*(E*Ltot0 + a*4.0d0*temp1 + b*dL)
          Ltot0 = dL
          if (present(dSedE)) then
            dSedE(Nat)%ELMFP(Nshl)%E(i) = E         ! [eV]
!             print*, 'dE=', E, Emin, Emax
            dSedE(Nat)%ELMFP(Nshl)%L(i) = (g_Pi*g_a0*Ele*g_me)/(MSHI*Zeff*Zeff*Ltot1)     ! [A]
            dSedE(Nat)%ELMFP(Nshl)%dEdx(i) = (g_Pi*g_a0*Ele*g_me)/(MSHI*Zeff*Zeff*ddEdx)  ! [eV/A]
          endif
          E = E + dE  ! [eV]
       enddo
       Sigma = 1.0d0/(g_Pi*g_a0*Ele)*MSHI/g_me*Zeff*Zeff*Ltot1     ! [1/A]
       if (Sigma > 1d30) Sigma = 1d30 ! get rid of infinities
       dEdx =  1.0d0/(g_Pi*g_a0*Ele)*MSHI/g_me*Zeff*Zeff*ddEdx     ! energy losses [eV/A]
    case default ! BEB cross section:
       Sigma = Sigma_BEB_SHI(Ele, Emax, Target_atoms(Nat)%Ip(Nshl),Target_atoms(Nat)%Ek(Nshl),Target_atoms(Nat)%Nel(Nshl), MSHI, Zeff) ! [A^2] cross section
       temp1 = Matter%At_Dens*1d-24*(Target_atoms(Nat)%Pers)/SUM(Target_atoms(:)%Pers)
       Sigma = 1.0d0/(temp1)*Sigma ! IMFP [A]
       if (Sigma > 1d30) Sigma = 1d30 ! get rid of infinities
       ddEdx = dSigma_w_int_BEB_SHI(Ele, Emax, Target_atoms(Nat)%Ip(Nshl), Target_atoms(Nat)%Ek(Nshl), Target_atoms(Nat)%Nel(Nshl), MSHI, Zeff) ! cross section integrated with energy [A^2*eV]
       dEdx = temp1*ddEdx !*MSHI/g_me ! energy losses [eV/A]
    end select
end subroutine SHI_TotIMFP



pure function Equilibrium_charge_Target(Ekin, Mass, ZSHI, Zmean, Kind_Zeff, fixed_Zeff) result (Zeff)  ! Equilibrium charge
   real(8) Zeff	! effective SHI state
   real(8), intent(in) :: Ekin	! [eV] kinetic energy of SHI
   real(8), intent(in) :: Mass  ! [kg] SHI mass
   real(8), intent(in) :: ZSHI	! SHI atomic number
   real(8), intent(in) :: Zmean	! mean atomic number of elements in the target
   integer, intent(in) :: Kind_Zeff	! model for effective charge
   real(8), intent(in) :: fixed_Zeff	! for the case of user-provided fixed charge of SHI
   !--------------------------
   real(8) x, x2, x4, Zt, Zp, vp, vpvo, c1, c2
   vp = sqrt(2.0d0*Ekin*g_e/Mass)  ! incident particle velocity
   !vp = velosity_from_kinetic_energy(Ekin, Mass, afs=.false.)     ! module "Relativity"


   !Zp = dble(ZSHI) ! SHI atomic number
   Zp = (ZSHI) ! SHI atomic number
   select case (Kind_Zeff)   ! 0=Barkas; 1=Bohr; 2=Nikolaev-Dmitriev; 3=Schiwietz-Grande;
      case (1)   ! Original Bohr:
         Zeff = Zp*(1.0d0-exp(-(vp/g_v0/(Zp**(0.66666666d0)))))
      case (2)   ! Nikolaev, Dmitriev, Phys. Lett. 28A, 277 (1968):
         c1 = 0.6d0      ! k
         c2 = 0.45d0     ! alpha
         Zeff = Zp*(1.0d0 + (vp/(Zp**c2*g_v0*4.0d0/3.0d0))**(-1.0d0/c1))**(-c1)
      case (3) ! Schiwietz et al, NIMB 225, 4 (2004):
         Zt = Zmean	! mean atomic number of target atoms
         c1 = 1.0d0 - 0.26d0*exp(-Zt/11.0d0 - (Zt-Zp)*(Zt-Zp)/9.0d0)
         vpvo = Zp**(-0.543d0)*vp/g_v0
         c2 = 1.0d0 + 0.03d0*vpvo*dlog(Zt)
         x = c1*(vpvo/c2/1.54d0)**(1.0d0 + 1.83d0/Zp)
         x2 = x*x
         x4 = x2*x2
         Zeff = Zp*(8.29d0*x + x4)/(0.06d0/x + 4.0d0 + 7.4d0*x + x4)
      case (4)
         Zeff = fixed_Zeff
      case default ! Barkas:
         Zeff = Zp*(1.0d0-exp(-(vp*125.0d0/g_cvel/(Zp**(0.66666666d0)))))
   endselect
end function Equilibrium_charge_Target


subroutine Equilibrium_charge_SHI(SHI, Target_atoms)  ! Equilibrium charge
   class(Ion), intent (inout) :: SHI  ! all the data for the SHI
   type(Atom), dimension(:), intent(in), optional :: Target_atoms  ! all data for target atoms
   real(8) x, x2, x4, Zt, Zp, vp, y, Zp052, vpvo, c1, c2
   vp = dsqrt(2.0d0*SHI%E*g_e/(SHI%Mass*g_Mp))  ! SHI velocity
   if (present(Target_atoms)) then
      !Zt = SUM(Target_atoms(:)%Zat)/size(Target_atoms)  ! mean atomic number of target atoms
      Zt = SUM(target_atoms(:)%Zat*dble(target_atoms(:)%Pers))/dble(SUM(target_atoms(:)%Pers)) ! mean atomic number of target atoms
      Zp = dble(SHI%Zat) ! SHI atomic number
      select case (SHI%Kind_Zeff)   ! 0=Barkas; 1=Bohr; 2=Nikolaev-Dmitriev; 3=Schiwietz-Grande;
         case (1)   ! Original Bohr:
            SHI%Zeff = Zp*(1.0d0-dexp(-(vp/g_v0/(Zp**(0.66666666d0)))))
         case (2)   ! Nikolaev, Dmitriev, Phys. Lett. 28A, 277 (1968):
            c1 = 0.6d0      ! k
            c2 = 0.45d0     ! alpha
            SHI%Zeff = Zp*(1.0d0 + (vp/(Zp**c2*g_v0*4.0d0/3.0d0))**(-1.0d0/c1))**(-c1)
         case (3) ! Schiwietz et al, NIMB 225, 4 (2004):
            c1 = 1.0d0 - 0.26d0*dexp(-Zt/11.0d0 - (Zt-Zp)*(Zt-Zp)/9.0d0)
            vpvo = Zp**(-0.543d0)*vp/g_v0
            c2 = 1.0d0 + 0.03d0*vpvo*dlog(Zt)
            x = c1*(vpvo/c2/1.54d0)**(1.0d0 + 1.83d0/Zp)
            x2 = x*x
            x4 = x2*x2
            SHI%Zeff = Zp*(8.29d0*x + x4)/(0.06d0/x + 4.0d0 + 7.4d0*x + x4)
         case (4)
            SHI%Zeff = SHI%fixed_Zeff
         case default ! Barkas:
            SHI%Zeff = Zp*(1.0d0-dexp(-(vp*125.0d0/g_cvel/(Zp**(0.66666666d0)))))
      endselect
   else ! only Barkas is possibe:
      SHI%Zeff = Zp*(1.0d0-dexp(-(vp*125.0d0/g_cvel/(Zp**(0.66666666d0)))))
   endif
end subroutine Equilibrium_charge_SHI


subroutine SHI_Diff_cross_section(Ele, MSHI, Emax, hw, Target_atoms, Nat, Nshl, Diff_IMFP, Matter, Mat_DOS, NumPar)
    real(8), intent(in), target :: Ele  ! SHI energy [eV]
    real(8), intent(in) :: MSHI    ! atomic mass of SHI [kg]
    real(8), intent(in) :: Emax    ! maximum transferred energy [eV]
    REAL(8), INTENT(in) :: hw   ! transferred energy [eV]
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    integer, intent(in) :: Nat, Nshl    ! number of atom, and number of shell    
    real(8), intent(out) :: Diff_IMFP   ! differential inverse mean free path 1/lambda(Ele,dE)
    type(Solid), intent(in) :: Matter
    type(Density_of_states), intent(in) :: Mat_DOS
    type(Flag), intent(in) :: NumPar
    
    integer i, n
    real(8) dLs, qmin, qmax, hq, ddq, dq, dE, Ime, dLs0, dL, hq0, dq_save
    real(8), pointer :: Ee
    real(8) a, b, x, temp1, temp2, T_fact
    Ee => Ele        ! energy [eV] -- don't copy, just point onto it
    
    !qmin = sqrt(2.0e0*g_me)/g_h*(sqrt(Ee) - sqrt((Ee - dE))) ! min transferred momentum [kg*m/s] electron
    !qmax = sqrt(2.0e0*g_me)/g_h*(sqrt(Ee) + sqrt((Ee - dE))) ! max transferred momentum [kg*m/s]
    qmin = hw/g_h/sqrt(2.0e0*Ee/MSHI)       ! min transferred momentum [kg*m/s] SHI
    qmax = sqrt(2.0e0*g_me)/g_h*sqrt(Emax)  ! min transferred momentum [kg*m/s] SHI

    dLs = 0.0e0 ! starting integration, mean free path per energy [A/eV]^(-1)
    hq = qmin    ! transient transferred momentum for integration [sqrt(eV/J) /m] (not [kg*m/s])
    n = 100
    !n = 200
    dq = (qmax - qmin)/real(n) ! differential of transferred momentum [kg*m/s]
    dLs0 = 0.0e0
    do while (hq .LT. qmax) ! no matter how many points, go till the end
        dq = hq/real(n)
        ! If it's Simpson integration:
        a = hq + dq/2.0e0
        call Imewq(hw, a, Target_atoms, Nat, Nshl, ImE, Matter, Mat_DOS, NumPar)
        temp1 = ImE
        b = hq + dq
        call Imewq(hw, b, Target_atoms, Nat, Nshl, ImE, Matter, Mat_DOS, NumPar)
        dL = ImE
        dLs = dLs + dq/6.0e0*(dLs0 + 4.0e0*temp1 + dL)/hq
        dLs0 = dL
        hq = hq + dq
    enddo    
    if (Matter%temp > 0.0d0) then   ! non-zero temperature
       T_fact = 1.0d0/(1.0d0-exp(-hw/Matter%temp*g_kb))
    else    ! zero temeprature
       T_fact = 1.0d0
    endif
    Diff_IMFP = dLs*T_fact
    
    nullify(Ee)
end subroutine SHI_Diff_cross_section


subroutine Imewq(hw, hq, Target_atoms, Nat, Nshl, ImE, Matter, Mat_DOS, NumPar, photon) ! constructs full Im(-1/e(w,q)) as a sum of Drude-like functions
    REAL(8), INTENT(in) ::  hw    ! transferred energy [eV]
    REAL(8), INTENT(in) ::  hq    ! transferred momentum [sqrt(eV/J) /m] (not [kg*m/s] !)
    REAL(8), INTENT(out) :: ImE   ! loss function Im(-1/e(w,q))
    type(Atom), dimension(:), intent(in), target :: Target_atoms  ! all data for target atoms
    integer, intent(in) :: Nat, Nshl    ! number of atom, and number of shell
    type(Density_of_states), intent(in), optional :: Mat_DOS
    type(Solid), intent(in) :: Matter
    type(Flag), intent(in) :: NumPar
    logical, optional, intent(in) :: photon
    
    real(8), pointer :: A, Gamma, E ! no need to copy variable, just point onto it!
    real(8) sumf
    integer i, N
    N = size(Target_atoms(Nat)%Ritchi(Nshl)%A)  ! that's how many fit functions we have
    ImE = 0.0e0
    do i = 1, N !Nff_esh(Nosh)
        A => Target_atoms(Nat)%Ritchi(Nshl)%A(i)            !A0(Nosh, i)
        E => Target_atoms(Nat)%Ritchi(Nshl)%E0(i)           !E0(Nosh, i)
        Gamma => Target_atoms(Nat)%Ritchi(Nshl)%Gamma(i)    !Gamma0(Nosh, i)
        if (present(photon)) then ! it's a photon
            sumf = Diel_func(A,E,Gamma, hw, hq, photon=.true.) ! Diel_func function, see below
        else ! it's a particale
            if (Matter%El_eff_Mass .EQ. 0) then
                 if (Target_atoms(1)%Ip(size(Target_atoms(1)%Ip)) .LT. 0.2d0) then  ! metallic material: DOS measured from bottom of CB
                    !Use inverted DOS
                    sumf = Diel_func(A,E,Gamma, hw, hq, NumPar=NumPar, Matter=Matter, k=Mat_DOS%k_inv, Eff_m=Mat_DOS%Eff_m_inv) ! Diel_func function, see below
                else         ! Insulator or semiconductor: DOS measured from top of VB
                    sumf = Diel_func(A,E,Gamma, hw, hq, NumPar=NumPar, Matter=Matter, k=Mat_DOS%k, Eff_m=Mat_DOS%Eff_m) ! Diel_func function, see below
                endif    
            else
                sumf = Diel_func(A,E,Gamma, hw, hq, NumPar=NumPar, Matter=Matter) ! Diel_func function, see below
            endif
        endif
        ImE = ImE + sumf
        nullify(A)
        nullify(E)
        nullify(Gamma)
    enddo
end subroutine Imewq


function Diel_func(A,E,Gamma,dE,dq, NumPar, Matter, Mtarget, photon, k, Eff_m)     ! fit functions in Ritchi algorithm
    real(8) A, E, Gamma, dE, dq     ! parameters and variable
    real(8) Diel_func               ! function itself
    real(8) E0                      ! temporary parameter
    type(Flag), optional :: NumPar
    type(Solid), optional :: Matter
    real(8), dimension(:), optional :: k, Eff_m
    real(8), optional :: Mtarget    ! [kg] target atoms mass
    logical, optional, intent(in) :: photon  ! for photon always q=0
        
    real(8) qlim, EEE, sqq, Mass, Ef_m, Gamma1, hq2, dE2, E02
    integer i, j

    phot:if (present(photon)) then ! it's a photon:
        dE2 = dE*dE
        E02 = E*E
        Diel_func = A*Gamma*dE/((dE2 - E02)*(dE2 - E02) + Gamma*Gamma*dE2)
    else phot ! an electron or a hole:
        hq2 = g_h*g_h*dq*dq

        mtar:if (present(Mtarget)) then   ! if scattering center is atom
            E0 = E + hq2/(2.0d0*Mtarget)  ! for phonons
            Gamma1 = Gamma
        else mtar ! scattering center is electron
            if (present(k) .AND. present(Eff_m)) then ! effective mass from VB or CB
                qlim = abs(dq)*sqrt(g_e)
                if (qlim .LE. k(size(k))) then
                    call find_in_array_monoton(k, qlim, j)
                    Mass = Eff_m(j)
                else
                    Mass = 1.0d0
                endif
            else if (Matter%El_eff_Mass .GT. 0) then  ! effective mass as a constant
                Mass = Matter%El_eff_Mass
            else  ! free-electron mass
                Mass = 1.0d0
            endif
            
            sqq = hq2/(2.0d0*Mass*g_me)            
            
            select case(NumPar%kind_of_DR)
                case(1) ! free electron dispersion relation
                    E0 = E + sqq
                    Gamma1 = Gamma
                case(2) ! Plasmon-pole approximation                         
                    E0 = sqrt(E*E + Matter%v_f*Matter%v_f*hq2*0.3333333333333d0 + sqq*sqq)
                    Gamma1 = Gamma
                case(3) ! Extended dielectric model of Ritchie                     
                    !E0 = (E**(2.0d0/3.0d0) + (sqq)**(2.0d0/3.0d0))**(3.0d0/2.0d0)
                    E0 = (E**(0.666666666666d0) + (sqq)**(0.666666666666d0))**1.5d0
                    Gamma1 = sqrt(Gamma*Gamma + sqq*sqq)
                case default ! free electron by default
                    E0 = E + sqq
                    Gamma1 = Gamma
            end select    
        endif mtar
        dE2 = dE*dE
        E02 = E0*E0
!         print*, 'E0', dE, sqq, E0
        Diel_func = A*Gamma1*dE/((dE2 - E02)*(dE2 - E02) + Gamma1*Gamma1*dE2)
!         write(*,'(a,f,f,f,f,f)') 'E0', dE, sqq, Diel_func
    endif phot
end function Diel_func



subroutine SHI_TotIMFP_BK(SHI, Target_atoms, Nat, Nshl, Sigma, dEdx, Matter, Mat_DOS, NumPar)
    class(Ion), intent (inout) :: SHI  ! all the data for the SHI
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    integer, intent(in) :: Nat, Nshl    ! number of atom, number of shell
    real(8), intent(out) :: Sigma, dEdx   ! calculated inverse mean free path (cross-section) [1/A], and the energy losses [eV/A]
    type(Solid), intent(in) :: Matter
    type(Density_of_states), intent(in) :: Mat_DOS
    type(Flag), intent(in) :: NumPar
    
    real(8) :: Ele      ! energy [eV]    
    real(8) :: M_SHI    ! atomic mass of SHI
    real(8) :: Z_SHI    ! nuclear charge of SHI
    
    integer i, j, n
    real(8) Emin, Emax, E, dE, dL, Ltot1, Ltot0, ddEdx, a, b, temp1, temp2, MSHI, Zeff
    real(8) Egap, Eplasmon
    
    call Equilibrium_charge_SHI(SHI, Target_atoms)  ! get Barcas' equilibrium charge from module Cross_sections
    Ele = SHI%E  ! energy of the particle [eV]
    MSHI = g_Mp*SHI%Mass    ! SHI mass [kg]
    Zeff = SHI%Zeff   ! SHI effective charge (must be precalculated)
    Egap = Target_atoms(1)%Ip(size(Target_atoms(1)%Ip))   ! band gap [eV]
    Emin = Target_atoms(Nat)%Ip(Nshl) ! min transferred energy [eV]
    
    if (Emin .LE. 1.0d-3) Emin = 1.0d-3 ! for metals there is no band gap, but we need a finite value
    Emax = 4.0d0*Ele*g_me*MSHI/((MSHI+g_me)*(MSHI+g_me)) ! [eV] maximum energy transfer of SHI

    if (Numpar%plasmon_Emax) then       ! Use maximal plasmon energy as upper limit of integration
        if (Emin .EQ. Egap) then
            Eplasmon = sqrt(Matter%N_VB_el*Matter%At_Dens*1d6*g_h*g_h/(g_me*g_e0) + Egap*Egap)    ! [g/cm^3] material density
            if (Eplasmon .GE. Emax) Emax = Eplasmon ! single atom vs plasmon
        endif
    endif 
    
    n = 10*(MAX(INT(Emin),10))    ! number of integration steps
    dE = (Emax - Emin)/(real(n)) ! differential of transferred momentum [kg*m/s]
    i = 1       ! to start integration
    E = Emin    ! to start integration
    Ltot1 = 0.0d0
    Ltot0 = 0.0d0
    call SHI_Diff_cross_section_BK(SHI, Emax, E, Target_atoms, Nat, Nshl, Ltot0, Matter, Mat_DOS, NumPar)
    ddEdx = 0.0e0
    do while (E .LE. Emax) ! integration
        dE = (1.0d0/(E+1.0d0) + E)/real(n)
        ! If it's Simpson integration:
        a =  E + dE/2.0e0
        call SHI_Diff_cross_section_BK(SHI, Emax, a, Target_atoms, Nat, Nshl, dL, Matter, Mat_DOS, NumPar)
        temp1 = dL
        b = E + dE
        call SHI_Diff_cross_section_BK(SHI, Emax, b, Target_atoms, Nat, Nshl, dL, Matter, Mat_DOS, NumPar)
        temp2 = dE/6.0d0*(Ltot0 + 4.0d0*temp1 + dL)
        Ltot1 = Ltot1 + temp2
        ddEdx = ddEdx + E*temp2
        Ltot0 = dL
        E = E + dE  ! [eV]
    enddo
    Sigma = 1.0d0/(g_Pi*g_a0*Ele)*MSHI/g_me*Ltot1     ! [1/A]
    if (Sigma > 1d30) Sigma = 1d30 ! get rid of infinities
    dEdx =  1.0d0/(g_Pi*g_a0*Ele)*MSHI/g_me*ddEdx     ! energy losses [eV/A]
end subroutine SHI_TotIMFP_BK


subroutine SHI_Diff_cross_section_BK(SHI, Emax, hw, Target_atoms, Nat, Nshl, Diff_IMFP, Matter, Mat_DOS, NumPar)
    class(Ion), intent (in) :: SHI  ! all the data for the SHI
    real(8), intent(in) :: Emax    ! maximum transferred energy [eV]
    REAL(8), INTENT(in) :: hw   ! transferred energy [eV]
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    integer, intent(in) :: Nat, Nshl    ! number of atom, and number of shell    
    real(8), intent(out) :: Diff_IMFP   ! differential inverse mean free path 1/lambda(Ele,dE)
    type(Solid), intent(in) :: Matter
    type(Density_of_states), intent(in) :: Mat_DOS
    type(Flag), intent(in) :: NumPar
    
    integer i, n
    real(8) dLs, qmin, qmax, hq, ddq, dq, Ee, dE, Ime, dLs0, dL, hq0, dq_save
    real(8) a, b, x, temp1, temp2, MSHI, Zeff, Z_SHI, rho
    
    Ee = SHI%E  ! energy of the particle [eV]
    MSHI = g_Mp*SHI%Mass    ! SHI mass [kg]
    Zeff = SHI%Zeff   ! SHI effective charge (must be precalculated)
    Z_SHI = SHI%Zat ! atomic number of SHI
    dE = hw         ! transferred energy [eV]
    
    !qmin = sqrt(2.0e0*g_me)/g_h*(sqrt(Ee) - sqrt((Ee - dE))) ! min transferred momentum [kg*m/s] electron
    !qmax = sqrt(2.0e0*g_me)/g_h*(sqrt(Ee) + sqrt((Ee - dE))) ! max transferred momentum [kg*m/s]
    qmin = hw/g_h/sqrt(2.0d0*Ee/MSHI)       ! min transferred momentum [kg*m/s] SHI
    qmax = sqrt(2.0e0*g_me)/g_h*sqrt(Emax)  ! min transferred momentum [kg*m/s] SHI

    dLs = 0.0d0 ! starting integration, mean free path per energy [A/eV]^(-1)
    hq = qmin    ! transient transferred momentum for integration [kg*m/s]
    n = 100
    dq = (qmax - qmin)/real(n) ! differential of transferred momentum [kg*m/s]
    dLs0 = 0.0d0
    do while (hq .LT. qmax) ! no matter how many points, go till the end
        dq = hq/real(n)
        ! If it's Simpson integration:
        a = hq + dq/2.0d0
        call Imewq(hw, a, Target_atoms, Nat, Nshl, ImE, Matter, Mat_DOS, NumPar)
        call Brand_Kitagawa(a, Z_SHI, Zeff, rho) 
        temp1 = ImE*rho*rho
        b = hq + dq
        call Imewq(hw, b, Target_atoms, Nat, Nshl, ImE, Matter, Mat_DOS, NumPar)
        call Brand_Kitagawa(b, Z_SHI, Zeff, rho) 
        dL = ImE*rho*rho
        dLs = dLs + dq/6.0d0*(dLs0 + 4.0d0*temp1 + dL)/hq
        dLs0 = dL
        hq = hq + dq
    enddo    
    Diff_IMFP = dLs/(1-exp(-hw/Matter%temp*g_kb))
end subroutine SHI_Diff_cross_section_BK


subroutine Brand_Kitagawa(hq, Z_SHI, Zeff, rho)
    real(8), intent(in) :: hq
    real(8), intent(in) :: Z_SHI
    real(8), intent(in) :: Zeff
    real(8), intent(out) :: rho
    real(8) kl2, a, Z, kl
    a = 0.2400519147d0
    Z = (Z_SHI - Zeff)/Z_SHI
    kl = hq*(g_a0*1d-10*sqrt(g_e))*2.0d0*a*Z**(2.0d0/3.0d0)/(Z_SHI**(2.0d0/3.0d0)*(1.0d0-Z/7.0d0))
    kl2 = kl*kl
    rho = Z_SHI*(1.0d0 - Z + kl2)/(1.0d0 + kl2)
end subroutine Brand_Kitagawa


subroutine Elastic_cross_section(Ee, CDF_Phonon, Target_atoms, Matter, EMFP, dEdx, NumPar, Mat_DOS, kind_of_particle, prefact)
   real(8), intent(in) :: Ee    ! [eV] electron energy
   type(CDF), intent(in) :: CDF_Phonon
   type(Atom), dimension(:), intent(in) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   type(Solid), intent(inout) :: Matter   ! all material parameters
   real(8), intent(out) :: EMFP ! [A] elastic mean free path
   real(8), intent(out) :: dEdx   ! [eV/A] energy loss
   type(Flag), intent(in) :: NumPar
   type(Density_of_states) :: Mat_DOS
   character(8), intent(in) :: kind_of_particle
   real(8), intent(in), optional :: prefact  ! prefactor if required (e.g. for negative frequencies)
   !--------------------------
         
   real(8) Sigma_Tot    ! [cm^2] total elastic cross-section
   real(8) Sigma_el, Zeff, Zt, Mass
   integer i, Mnum

   dEdx = 0.0d0
   ! Get the total cross-section:
   if (NumPar%kind_of_EMFP .EQ. 1) then   ! CDF cross section

      ! Target mean atomic charge:
      if (numpar%CDF_elast_Zeff /= 1) then   ! Barkas-like charge
         Zt = SUM(target_atoms(:)%Zat*dble(target_atoms(:)%Pers))/dble(SUM(target_atoms(:)%Pers)) ! mean atomic number of target atoms
         Zeff = 1.0d0 + Equilibrium_charge_Target(Ee, g_me, Zt, (Zt-1.0e0), 0, 1.0d0) ! Equilibrium charge, see below
      else  ! one, as used in old CDF expression
         Zeff = 1.0d0    ! electron charge
      endif

      if (present(prefact)) then
         call Tot_EMFP(Ee, Target_atoms, CDF_Phonon, Matter, EMFP, dEdx, NumPar, Mat_DOS, kind_of_particle, Zeff, prefact) ! below
      else  ! no prefactor needed, default option
         call Tot_EMFP(Ee, Target_atoms, CDF_Phonon, Matter, EMFP, dEdx, NumPar, Mat_DOS, kind_of_particle, Zeff) ! below
      endif
   else  ! Mott cross section
      if (kind_of_particle .EQ. 'Electron') then
         Mass = 1.0d0
      else if (kind_of_particle .EQ. 'Hole') then
         if(Matter%hole_mass .GE. 0) then
               Mass = Matter%hole_mass
         else  ! Define mass from DOS
               call find_in_array_monoton(Mat_DOS%E, Ee, Mnum)
               Mass = Mat_DOS%Eff_m(Mnum)
         endif
      endif
      Sigma_Tot = 0.0d0 ! [cm^2] cross-section
      do i = 1, size(Target_atoms)  ! for all atomic spicies:
         call Atomic_elastic_sigma(Target_atoms, i, Ee, Sigma_el) ! total cross-section of elastic scattering on an atom:
         ! [cm^2] total cross-section as sum of atomic for all spicies:
         !Sigma_Tot = Sigma_Tot + Sigma_el
         Sigma_Tot = Sigma_Tot + Sigma_el*dble(target_atoms(i)%Pers) ! account for stoichiometry
      enddo
      !Sigma_Tot = Sigma_Tot/size(Target_atoms)
      Sigma_Tot = Sigma_Tot/dble(SUM(target_atoms(:)%Pers))*Mass ! account for stoichiometry

      EMFP = 1.0d8/(Sigma_Tot*Matter%At_Dens) ! [A] elastic mean free path
   endif
end subroutine Elastic_cross_section


subroutine Tot_EMFP(Ele, Target_atoms, CDF_Phonon, Matter, Sigma, dEdx, NumPar, Mat_DOS, kind_of_particle, Zeff, prefact)
    real(8), intent(in), target :: Ele  ! electron energy [eV]
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    type(CDF), intent(in) :: CDF_Phonon ! phonon CDF parameters
    type(Solid), intent(inout) :: Matter   ! all material parameters
    real(8), intent(out) :: Sigma, dEdx   ! calculated inverse mean free path (cross-section) [A], and the energy losses [eV/A]
    type(Flag), intent(in) :: NumPar
    type(Density_of_states), intent(in) :: Mat_DOS     
    character(8), intent(in) :: kind_of_particle
    real(8), intent(in) :: Zeff  ! effective charge of target atoms [electron charge]
    real(8), intent(in), optional :: prefact ! prefactor for Emax, if required (e.g., use for negative frequencies)
    !--------------------
    integer :: i, j, n, Mnum
    real(8) :: Emin, Emax, E, dE, dL, Ltot1, Ltot0, ddEdx
    real(8) :: a, b, temp1, temp2, qdebay, Edebay, Mtarget, Mass, pref
    real(8), pointer :: Ee

    if (present(prefact)) then
      pref = prefact ! defined by the user
    else
      pref = 1.0d0   ! no prefactor required
    endif
    qdebay = (6.0d0*g_Pi*g_Pi*(Matter%At_Dens*1e6))**(0.33333333d0)   ! maximum energy of phonons analyzed [1/m], debay energy
    Edebay = 2.0d0*g_h*Matter%Vsound*qdebay/g_e  ! maximum energy of phonons analyzed [eV], debay energy // 2 is artificial...
    Ee => Ele   ! just to use this name later
    Emin = pref * 0.1d-8    !Target_atoms(Nat)%Ip(Nshl)   ! [eV] ionization potential of the shell is minimum possible transferred energy
    Mtarget = g_Mp*SUM(target_atoms(:)%Mass*dble(target_atoms(:)%Pers))/dble(SUM(target_atoms(:)%Pers)) ! average mass of a target atom [kg]
    
    if (trim(adjustl(kind_of_particle)) .EQ. 'Electron') then
        Mass = 1.0d0
    else if (trim(adjustl(kind_of_particle)) .EQ. 'Hole') then
        if(Matter%hole_mass .GE. 0) then
            Mass = Matter%hole_mass
        else                    ! Define mass from DOS
            call find_in_array_monoton(Mat_DOS%E, Ele, Mnum)
            Mass = Mat_DOS%Eff_m(Mnum)
        endif 
    else
        Mass = 1.0d0
    endif    
        
    Emax = 4.0e0*Ee*Mass*g_me*Mtarget/((Mtarget+Mass*g_me)*(Mtarget+Mass*g_me))    ! [eV] maximum energy transfer
      
    if (Edebay .GE. Emax) Emax = Edebay ! single atom vs phonon
    if (Ee .LT. Emax) Emax = Ee ! no more than the total electron energy
    Emax = pref * Emax  ! if user defined some prefactor (e.g. negatve for phonon absorption)

    n = 20*(MAX(INT(Emin),10))    ! number of integration steps
    dE = (Emax - Emin)/(dble(n)) ! differential of transferred energy
    i = 1       ! to start integration
    E = Emin    ! to start integration
    Ltot1 = 0.0d0
    Ltot0 = 0.0d0
    call Diff_cross_section_phonon(Ele, E, CDF_Phonon, Ltot0, Mtarget, Mass, Matter%temp, 1.0d0)
    ddEdx = 0.0e0
    do while (abs(E) .LE. abs(Emax)) ! integration
        dE = (1.0d0/(E+1.0d0) + E)/dble(n)
        dE = pref * dE  ! if user requested a prefactor
        ! If it's Simpson integration:
        a =  E + dE/2.0d0
        call Diff_cross_section_phonon(Ele, a, CDF_Phonon, dL, Mtarget, Mass, Matter%temp, 1.0d0)
        temp1 = dL
        b = E + dE
        call Diff_cross_section_phonon(Ele, b, CDF_Phonon, dL, Mtarget, Mass, Matter%temp, 1.0d0)
        temp2 = abs(dE)/6.0d0*(Ltot0 + 4.0d0*temp1 + dL)
        Ltot1 = Ltot1 + temp2
        ddEdx = ddEdx + E*temp2
        Ltot0 = dL
        E = E + dE  ! [eV]
        
        ! test:
        !if ((Ele-1.0d0) < 1.0d-2) then
        !   print*, 'Tot_EMFP', Ele, E, 1.0d0/Mass/Ltot1
        !endif
    enddo
    !Sigma = 1.0d0/Mass/Ltot1 !*dE ! [A]
    Sigma = 1.0d0/(Zeff*Zeff*Mass*Ltot1) !*dE ! [A]
    if (Sigma > 1d30) Sigma = 1d30 ! get rid of infinities

    !dEdx = Mass*ddEdx !*dE ! energy losses [eV/A]
    dEdx = (Zeff*Zeff) * Mass*ddEdx !*dE ! energy losses [eV/A]

    nullify(Ee)
end subroutine Tot_EMFP


subroutine Diff_cross_section_phonon(Ele, hw, CDF_Phonon, Diff_IMFP, Mtarget, Mass, Ttarget, pref)
    real(8), intent(in), target :: Ele  ! SHI energy [eV]
    REAL(8), INTENT(in), target :: hw   ! transferred energy [eV]
    type(CDF), intent(in) :: CDF_Phonon ! phonon CDF parameters
    real(8), intent(out) :: Diff_IMFP   ! differential inverse mean free path 1/lambda(Ele,dE)
    real(8), intent(in) :: Mtarget, Mass  ! average mass of atoms of the target [kg]
    real(8), intent(in) :: Ttarget      ! temperature of the sample [K]
    real(8), intent(in) :: pref  ! prefactor defined by the user (e.g., for negative energies)
    !---------------------------
    integer i, n
    real(8), pointer :: Ee, dE
    real(8) dLs, qmin, qmax, hq, ddq, dq, Ime, dLs0, dL, hq0, dq_save
    real(8) a, b, x, temp1, temp2
    Ee => Ele        ! energy [eV]
    dE => hw         ! transferred energy [eV]
    
    qmin = sqrt(2.0d0*Mass*g_me)/g_h*(sqrt(Ee) - sqrt((Ee - dE)))        ! min transferred momentum [kg*m/s]
    qmax = pref * sqrt(2.0d0*Mass*g_me)/g_h*(sqrt(Ee) + sqrt((Ee - dE))) ! max transferred momentum [kg*m/s]
    
    dLs = 0.0d0  ! starting integration, mean free path per energy [A/eV]^(-1)
    hq = qmin    ! transient transferred momentum for integration [kg*m/s]
    n = 100
    dq = (qmax - qmin)/dble(n) ! differential of transferred momentum [kg*m/s]
    dLs0 = 0.0d0
    do while (abs(hq) .LT. abs(qmax)) ! no matter how many points, go till the end
        dq = hq/dble(n)
        ! If it's Simpson integration:
        a = hq + dq/2.0d0
        call Imewq_phonon(hw, a, CDF_Phonon, ImE, Mtarget)
        temp1 = ImE
        b = hq + dq
        call Imewq_phonon(hw, b, CDF_Phonon, ImE, Mtarget)
        dL = ImE
        dLs = dLs + dq/6.0d0*(dLs0 + 4.0d0*temp1 + dL)/hq
        dLs0 = dL
        hq = hq + dq
        
!         if ( (abs(Ele-1.0d0) < 1.0d-2) .and. ( abs(hw - 0.1d0) < 0.01d0) ) &
!            print*, 'Diff_cross_section_phonon', Ele, hw, g_h*g_h*hq*hq/(2.0d0*Mtarget), dL
    enddo    
    Diff_IMFP = 1.0d0/(g_Pi*g_a0*Ele)*dLs/(1-exp(-hw/Ttarget*g_kb))
    nullify(Ee)
    nullify(dE)
end subroutine Diff_cross_section_phonon


subroutine Imewq_phonon(hw, hq, CDF_Phonon, ImE, Mtarget) ! constructs full Im(-1/e(w,q)) as a sum of Drude-like functions
    REAL(8), INTENT(in) ::  hw    ! transferred energy [eV]
    REAL(8), INTENT(in) ::  hq    ! transferred momentum [kg*m/s]
    REAL(8), INTENT(out) :: ImE
    real(8), intent(in) :: Mtarget
    type(CDF), intent(in), target :: CDF_Phonon ! phonon CDF parameters
    
    real(8), pointer :: A, Gamma, E ! no need to copy variable, just point onto it!
    real(8) sumf
    integer i, N
    N = size(CDF_Phonon%A)  ! that's how many fit functions we have
    ImE = 0.0d0
    do i = 1, N !Nff_esh(Nosh)
        A => CDF_Phonon%A(i)            !A0(Nosh, i)
        E => CDF_Phonon%E0(i)           !E0(Nosh, i)
        Gamma => CDF_Phonon%Gamma(i)    !Gamma0(Nosh, i)
        sumf = Diel_func(A,E,Gamma, hw, hq, Mtarget=Mtarget) ! Diel_func function, see below
        ImE = ImE + sumf
        nullify(A)
        nullify(E)
        nullify(Gamma)
    enddo
end subroutine Imewq_phonon


subroutine Atomic_elastic_sigma(Target_atoms, KOA, Ee, sigma_el) ! total cross-section of elastic scattering on an atom:
   type(Atom), dimension(:), intent(in), target :: Target_atoms  ! all the target atoms parameters
   integer, intent(in) :: KOA   ! kind of atom
   real(8), intent(in) :: Ee    ! [eV] incident electron energy
   REAL(8), INTENT(out) :: sigma_el   ! cross-section of elastic scteering [cm^2]
   
   real(8) :: Zat   ! atomic number of an atom of the media
   real(8) nc, pc, mec2e, Zat137, RyEe, beta2
   
   Zat = real(Target_atoms(KOA)%Zat)   ! atom atomic number, point onto it to use below
   mec2e = g_me*g_cvel*g_cvel/g_e   ! a parameter enterring eq. below:
   Zat137 = Zat/137.0d0   ! a parameter enterring eq. below:
   RyEe = g_Ry/Ee   ! a parameter enterring eq. below:
   beta2 = 2.0d0*Ee/(mec2e) ! v/c, a parameter enterring eq. below:

   pc = 1.7d-5*(Zat**(2.0d0/3.0d0))*(1.0d0-beta2)/(beta2)
   nc = pc*(1.13d0 + 3.76d0*(Zat137*Zat137)/(beta2)*sqrt(Ee/(Ee+mec2e)))
   sigma_el = g_Pi*g_a0*g_a0*Zat*(Zat+1.0d0)/(nc*(nc+1.0d0))*RyEe*RyEe*1d-16  ! Mott cross-section [cm^2]

end subroutine Atomic_elastic_sigma


subroutine NRG_transfer_elastic_atomic(Target_atoms, KOA, Ee, dE, M_eff) ! total cross-section of elastic scattering on an atom:
   type(Atom), dimension(:), intent(in), target :: Target_atoms  ! all the target atoms parameters
   integer, intent(in) :: KOA   ! kind of atom
   real(8), intent(in) :: Ee    ! [eV] incident electron energy
   REAL(8), INTENT(out) :: dE   ! [eV] transferred energy
   real(8), intent(in), optional :: M_eff ! factor to rescale incident article mass (use for holes scattering)s
 
   REAL(8) :: Mat   ! mass of an atom of the media [kg]
   real(8) :: M_effective
   real(8) :: Zat   ! atomic number of an atom of the media
   real(8) nc, pc, mec2e, Zat137, RyEe, Masses, x
   
   call random_number(x)

   if (present(M_eff)) then
      M_effective = M_eff
   else
      M_effective = 1.0d0
   endif
   
   Mat = Target_atoms(KOA)%Mass*g_Mp  ! [kg] atom mass
   Zat = dble(Target_atoms(KOA)%Zat)   ! atom atomic number
   mec2e = g_me*g_cvel*g_cvel/g_e   ! a parameter enterring eq. below:
   Zat137 = Zat/137.0d0   ! a parameter enterring eq. below:
   RyEe = g_Ry/Ee   ! a parameter enterring eq. below:
   Masses = Mat + g_me*M_effective    ! a parameter enterring eq. below:
   pc = 1.7d-5*(Zat**(2.0d0/3.0d0))*(1.0d0-2.0d0*Ee/(mec2e))/(2.0d0*Ee/(mec2e))
   nc = pc*(1.13d0 + 3.76d0*(Zat137*Zat137)/(2.0d0*Ee/mec2e)*sqrt(Ee/(Ee+mec2e)))
   dE = 4.0d0*Ee*M_effective*g_me*Mat/(Masses*Masses)*(nc*(1.0d0-x)/(nc+x))
end subroutine NRG_transfer_elastic_atomic 




subroutine NRG_transfer_elastic_DSF(Elastic_MFP, DSF_DEMFP, Eel, dE)
   type(MFP_elastic), intent(in) :: Elastic_MFP    ! Total elastic mean free paths
   type(Differential_MFP), dimension(:), intent(in) :: DSF_DEMFP  ! Differential EMFPs
   real(8), intent(in) :: Eel               ! Incident electron (or hole) energy [eV]
   real(8), intent(out) :: dE               ! Transferred energy in this collision
   !------------------------
   integer :: NumE, NumL, i, kk, j, i_MFP
   real(8) :: L_need, RN, Value1, EMFP_tot, EMFP_absorb, EMFP_emit
   real(8), dimension(:,:), allocatable :: dLdE
   logical :: it_is_emission


   ! 0) Get the index of the array according to the total energy of the particle:
   call Find_in_array_monoton(DSF_DEMFP%E, Eel, NumE) ! module "Reading_files_and_parameters"
   if (NumE > 1) then
      NumE = NumE-1
   endif
   if ((Eel < DSF_DEMFP(NumE)%E) .or. (Eel > DSF_DEMFP(NumE+1)%E)) then
      if (NumE > 1) print*, 'Problem in NRG_transfer_elastic_DSF: energy transfer above DSF limits', Eel, DSF_DEMFP(NumE)%E, DSF_DEMFP(NumE+1)%E
   endif

   ! To interpolate the data between two energy points:
   allocate(dLdE(2,size(DSF_DEMFP(NumE)%dL)), source = 0.0d0)


   ! 1) Compare absorption vs emission of energy into atomic system:
   ! Find energy in the array:
   call Find_in_array_monoton(Elastic_MFP%Total%E, Eel, i_MFP) ! module "Reading_files_and_parameters"
   if (i_MFP > 1) then
      i_MFP = i_MFP - 1
   endif
   !print*, Elastic_MFP%Total%E(i_MFP), Elastic_MFP%Absorb%L(i_MFP), Eel
   !print*, Elastic_MFP%Total%E(i_MFP+1), Elastic_MFP%Absorb%L(i_MFP+1)

   ! Interpolate MFPs:
   EMFP_tot = linear_interpolation(Elastic_MFP%Total%L(i_MFP), Elastic_MFP%Total%L(i_MFP+1), &
                                 Elastic_MFP%Total%E(i_MFP), Elastic_MFP%Total%E(i_MFP+1), Eel) ! below
   EMFP_emit = linear_interpolation(Elastic_MFP%Emit%L(i_MFP), Elastic_MFP%Emit%L(i_MFP+1), &
                                 Elastic_MFP%Total%E(i_MFP), Elastic_MFP%Total%E(i_MFP+1), Eel) ! below
   EMFP_absorb = linear_interpolation(Elastic_MFP%Absorb%L(i_MFP), Elastic_MFP%Absorb%L(i_MFP+1), &
                                 Elastic_MFP%Total%E(i_MFP), Elastic_MFP%Total%E(i_MFP+1), Eel) ! below

   !write(*,'(f,f,es,es,es)') Eel, EMFP_tot/EMFP_emit, EMFP_tot, EMFP_emit, EMFP_absorb

   ! Sample which process: emission vs absorption:
   call random_number(RN)
   if (RN < EMFP_tot/EMFP_emit) then
      it_is_emission = .true.    ! emission
   else
      it_is_emission = .false.   ! absorption
   endif
!    print*, 'EMFP', RN, it_is_emission, EMFP_tot/EMFP_emit
   !pause

   ! 2) Get the transferred energy:

   ! Set the interpolated values:
   if (it_is_emission) then ! emission (dE>0)
      if (NumE == size(DSF_DEMFP%E)) then
         dLdE(1,:) = DSF_DEMFP(NumE)%dL_emit(:)
         dLdE(2,:) = DSF_DEMFP(NumE)%dE(:)
      else
         Value1 = (Eel - DSF_DEMFP(NumE)%E)/(DSF_DEMFP(NumE+1)%E - DSF_DEMFP(NumE)%E)
         dLdE(1,:) = DSF_DEMFP(NumE)%dL_emit(:) + (DSF_DEMFP(NumE+1)%dL_emit(:) - DSF_DEMFP(NumE)%dL_emit(:)) * Value1
         dLdE(2,:) = DSF_DEMFP(NumE)%dE(:) + (DSF_DEMFP(NumE+1)%dE(:) - DSF_DEMFP(NumE)%dE(:)) * Value1
      endif
      ! Make sure interpolation didn't go wrong:
      do i = 1, size(dLdE,2)
         if (dLdE(1,i) < 0.0d0) then
            dLdE(1,i) = 0.0d0
            !print*, 'Trouble #1 in interpolation in NRG_transfer_elastic_DSF:', dLdE(1,i), DSF_DEMFP(NumE+1)%dL_emit(i), DSF_DEMFP(NumE)%dL_emit(i)
         endif
      enddo
   else ! absorption (dE<0)
      if (NumE == size(DSF_DEMFP%E)) then
         dLdE(1,:) = DSF_DEMFP(NumE)%dL_absorb(:)
         dLdE(2,:) = DSF_DEMFP(NumE)%dE(:)
      else
         Value1 = (Eel - DSF_DEMFP(NumE)%E)/(DSF_DEMFP(NumE+1)%E - DSF_DEMFP(NumE)%E)
         dLdE(1,:) = DSF_DEMFP(NumE)%dL_absorb(:) + (DSF_DEMFP(NumE+1)%dL_absorb(:) - DSF_DEMFP(NumE)%dL_absorb(:)) * Value1
         dLdE(2,:) = DSF_DEMFP(NumE)%dE(:) + (DSF_DEMFP(NumE+1)%dE(:) - DSF_DEMFP(NumE)%dE(:)) * Value1
      endif
      ! Make sure interpolation didn't go wrong:
      do i = 1, size(dLdE,2)
         if (dLdE(1,i) < 0.0d0) then
            dLdE(1,i) = 0.0d0
            !print*, 'Trouble #2 in interpolation in NRG_transfer_elastic_DSF:', i, dLdE(1,i), DSF_DEMFP(NumE+1)%dL_absorb(i), DSF_DEMFP(NumE)%dL_absorb(i)
         endif
      enddo
   endif

   ! Define the sampled MFP:
   call random_number(RN)

!   do j = 1, 100   ! Testing
!       RN = dble(j)/100.0d0    ! Testing
   if (it_is_emission) then ! emission (dE>0)
      L_need = EMFP_emit/RN   ! [A] we need to reach
   else
      L_need = EMFP_absorb/RN   ! [A] we need to reach
   endif

   call Linear_approx_2x1d_DSF(dLdE(1,:), dLdE(2,:), L_need, dE) ! module "Reading_files_and_parameters"

   if (dE .GT. 1.0d0) then ! Potentially unphysically large energy transfer
      if (it_is_emission) then ! emission (dE>0)
         print*, "Problem in NRG_transfer_elastic_DSF: unphyically large energy transfer (emission):"
         print*, dE, Eel, EMFP_emit, L_need
      else
         print*, "Problem in NRG_transfer_elastic_DSF: unphyically large energy transfer (absorption):"
         print*, dE, Eel, EMFP_absorb, L_need
      endif
   endif

   !--------------
!    ! Testing:
!    if (it_is_emission) then ! emission (dE>0)
!       kk = size(DSF_DEMFP(NumE)%dL_emit)
!       i = 1
!       do while (DSF_DEMFP(NumE)%dL_emit(i) >= L_need)
!         if (i .EQ. kk) exit
!         i = i+1
!       enddo
!       NumL = i
!       write(*,'(a,f,f,f,f,f,f,f,f)') 'dE', DSF_DEMFP(NumE)%E, Eel, RN, L_need, EMFP_emit, DSF_DEMFP(NumE)%dL_emit(NumL), DSF_DEMFP(NumE)%dE(NumL)
!    else
!       kk = size(DSF_DEMFP(NumE)%dL_absorb)
!       i = 1
!       do while (DSF_DEMFP(NumE)%dL_absorb(i) >= L_need)
!         if (i .EQ. kk) exit
!         i = i+1
!       enddo
!       NumL = i
!       write(*,'(a,f,f,f,f,f,f,f,f)') 'dE', DSF_DEMFP(NumE)%E, Eel, RN, L_need, EMFP_absorb, DSF_DEMFP(NumE)%dL_absorb(NumL), DSF_DEMFP(NumE)%dE(NumL)
!    endif
!  enddo  !testing
  !dE = DSF_DEMFP(NumE)%dE(NumL)
  deallocate(dLdE)
!   pause 'NRG_transfer_elastic_DSF'
endsubroutine NRG_transfer_elastic_DSF


pure function linear_interpolation(y1, y2, x1, x2, x_exact) result (y_out)
   real(8) y_out
   real(8), intent(in) :: y1, y2, x1, x2, x_exact
   y_out = y1 + (y2 - y1)/(x2 - x1)*(x_exact - x1)
end function linear_interpolation


subroutine NRG_transfer_elastic_DSF_OLD(DSF_DEMFP, Eel, EMFP, dE)
   type(Differential_MFP), dimension(:), intent(in) :: DSF_DEMFP
   real(8), intent(in) :: Eel               ! Incident electron energy [eV]
   real(8), intent(in) :: EMFP              ! Total EMFP of electron with energy Eel [A]
   real(8), intent(out) :: dE               ! Transferred energy in this collision
   !------------------------
   integer :: NumE, NumL, i, kk, j
   real(8) :: L_need, RN, Value1
   real(8), dimension(:,:), allocatable :: dLdE

!    i=1
!    do while (DSF_DEMFP(i)%E .LT. Eel)
!       i = i + 1
!    enddo
!    NumE = i-1

   ! Get the index of the array according to the total energy of the particle:
   call Find_in_array_monoton(DSF_DEMFP%E, Eel, NumE) ! module "Reading_files_and_parameters"
   NumE = NumE-1
   if ((Eel < DSF_DEMFP(NumE)%E) .or. (Eel > DSF_DEMFP(NumE+1)%E)) then
      print*, 'Problem in NRG_transfer_elastic_DSF: energy transfer above DSF limits', Eel, DSF_DEMFP(NumE)%E, DSF_DEMFP(NumE+1)%E
   endif

   ! Interpolate the data between two energy points:
   allocate(dLdE(2,size(DSF_DEMFP(NumE)%dL)), source = 0.0d0)
   ! Set the interpolated values:
   if (NumE == size(DSF_DEMFP%E)) then
      dLdE(1,:) = DSF_DEMFP(NumE)%dL(:)
      dLdE(2,:) = DSF_DEMFP(NumE)%dE(:)
   else
      Value1 = (Eel - DSF_DEMFP(NumE)%E)/(DSF_DEMFP(NumE+1)%E - DSF_DEMFP(NumE)%E)
      dLdE(1,:) = DSF_DEMFP(NumE)%dL(:) + (DSF_DEMFP(NumE+1)%dL(:) - DSF_DEMFP(NumE)%dL(:)) * Value1

      dLdE(2,:) = DSF_DEMFP(NumE)%dE(:) + (DSF_DEMFP(NumE+1)%dE(:) - DSF_DEMFP(NumE)%dE(:)) * Value1
   endif
   ! Make sure interpolation didn't go wrong:
   do i = 1, size(dLdE,2)
      if (dLdE(1,i) < 0.0d0) print*, 'Trouble #1 in interpolation in NRG_transfer_elastic_DSF:', dLdE(1,i), DSF_DEMFP(NumE+1)%dL(i), DSF_DEMFP(NumE)%dL(i)
   enddo


   ! Define the sampled MFP:
   call random_number(RN)

  do j = 1, 100   ! Testing
      RN = dble(j)/100.0d0    ! Testing
   L_need = EMFP/RN   ! [A] we need to reach

   !call Linear_approx_2x1d_DSF(DSF_DEMFP(NumE)%dL, DSF_DEMFP(NumE)%dE, L_need, dE) ! module "Reading_files_and_parameters"
   call Linear_approx_2x1d_DSF(dLdE(1,:), dLdE(2,:), L_need, dE) ! module "Reading_files_and_parameters"

   if (dE .GT. 1.0d0) then ! maybe unphysically large energy transfer
      print*, "Problem in NRG_transfer_elastic_DSF: unphyically large energy transfer", dE, Eel, EMFP, L_need
   endif

    kk = size(DSF_DEMFP(NumE)%dL)
    i = 1
    do while (DSF_DEMFP(NumE)%dL(i) >= L_need)
 !         print*, 'dL', i,  DSF_DEMFP(NumE)%dE(i), DSF_DEMFP(NumE)%dL(i), L_need
        if (i .EQ. kk) exit
        i = i+1
    enddo
    NumL = i
    if (NumL > 1) then
        Value1 = DSF_DEMFP(NumE)%dE(NumL-1)+(DSF_DEMFP(NumE)%dE(NumL)-DSF_DEMFP(NumE)%dE(NumL-1))/(DSF_DEMFP(NumE)%dL(NumL)-DSF_DEMFP(NumE)%dL(NumL-1))*(L_need - DSF_DEMFP(NumE)%dL(NumL-1))
    else
        Value1 = DSF_DEMFP(NumE)%dE(NumL)
    endif

!    write(*,'(a,f,i,f)') ' 1 :', Eel, NumE, DSF_DEMFP(NumE)%E
!    write(*,'(a,f,f)') ' 2 :', EMFP, RN
!    write(*,'(a,f,f)') ' 3 :', L_need, dE
!    write(*,'(a,f,f)') ' 4 :', DSF_DEMFP(NumE)%dL(NumL), DSF_DEMFP(NumE)%dE(NumL)
   write(*,'(a,f,f,f,f,f,f,f,f)') 'dE', DSF_DEMFP(NumE)%E, Eel, RN, L_need, EMFP, DSF_DEMFP(NumE)%dL(NumL), Value1, DSF_DEMFP(NumE)%dE(NumL)
 enddo  !testing
  !dE = DSF_DEMFP(NumE)%dE(NumL)
  deallocate(dLdE)
  pause 'NRG_transfer_elastic_DSF'
endsubroutine NRG_transfer_elastic_DSF_OLD


!BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB
! Ionization of an atom: BEB cross-section
! Eq.(57) from [Y.K.Kim, M.E.Rudd, Phys.Rev.A 50 (1994) 3954]
function Sigma_BEB(T,B,U,N)
   real(8) Sigma_BEB ! cross-section [A^2]
   real(8) T	! kinetic energy of incident electron [eV]
   real(8) B	! bindning energy of atomic electron [eV]
   real(8) U	! mean kinetic energy of electron in sub-shell [eV]
   real(8) N	! ocupation number of electrons in this shell
   real(8) t0, u0, S
   if (T .LE. B) then
      Sigma_BEB = 0.0d0 ! [A^2] cross section for energies lower than the Ip
   else
      S = 4.0d0*g_Pi*g_a0*g_a0*N*(g_Ry/B)*(g_Ry/B) ! Eq.(4)
      t0 = T/B	! energy normalized to the Rydberg constant ! Eq.(4)
      u0 = U/B
      Sigma_BEB = S/(t0+u0+1.0d0)*(log(t0)*0.5d0*(1.0d0-1.0d0/(t0*t0)) + (1.0d0-1.0d0/t0) -log(t0)/(t0+1.0d0)) ! [A^2]
   endif
end function Sigma_BEB


function Sigma_BEB_SHI(T,Emax,B,U,N,MSHI,ZSHI)
   real(8) Sigma_BEB_SHI ! cross-section [A^2]
   real(8) T	! kinetic energy of incident electron [eV]
   real(8) B	! bindning energy of atomic electron [eV]
   real(8) U	! mean kinetic energy of electron in sub-shell [eV]
   real(8) N	! ocupation number of electrons in this shell
   real(8) MSHI,ZSHI ! mass and charge of the SHI
   real(8) Emax, dSigma0, dSigma, w0, t0, u0, S
   !Emax = 4.0d0*T*g_me*MSHI/((MSHI+g_me)*(MSHI+g_me)) ! [eV] maximum energy transfer of SHI
   if (Emax .LE. B) then
      Sigma_BEB_SHI = 0.0d0 ! [A^2] cross section for energies lower than the Ip
   else
      !S = 4.0d0*g_Pi*g_a0*g_a0*N*(g_Ry/B)*(g_Ry/B) ! Eq.(4)
      !t0 = T/B/MSHI*g_Mp	! energy normalized for SHI (~Eq.(2)) to fit the data
      !u0 = U/B
      !Sigma_BEB_SHI = S/(t0+u0+1.0d0)*(log(t0)*0.5d0*(1.0d0-1.0d0/(t0*t0)) + (1.0d0-1.0d0/t0) -log(t0)/(t0+1.0d0)) ! [A^2]
      dSigma0 = dSigma_int_BEB_SHI(T, 0.0d0, B, U, N, MSHI, ZSHI)
      dSigma = dSigma_int_BEB_SHI(T, Emax, B, U, N, MSHI, ZSHI)
      Sigma_BEB_SHI = (dSigma - dSigma0) ! renormalize for SHI
   endif
end function Sigma_BEB_SHI

function dSigma_int_BEB_SHI(T, w, B, U, N, MSHI, ZSHI)
   real(8), intent(in) :: w 	! [eV] transferred energy
   real(8), intent(in) :: T 	! [eV] kinetic energy of the electron
   real(8), intent(in) ::  B	! bindning energy of atomic electron [eV]
   real(8), intent(in) ::  U	! mean kinetic energy of electron in the sub-shell [eV]
   real(8), intent(in) ::  N	! ocupation number of electrons in this shell
   real(8) MSHI,ZSHI ! mass and charge of the SHI
   real(8) :: dSigma_int_BEB_SHI    ! differential cross-section

   real(8) t0, u0, w0, S, dSigma, dSigma0, Emax

   !Emax = 4.0d0*T*g_me*MSHI/((MSHI+g_me)*(MSHI+g_me)) ! [eV] maximum energy transfer of SHI
!   if (w .LE. B) then ! for the case if incident electron kinetic energy is lower than the Ip
!      dSigma_int_BEB_SHI = 0.0d0
!   else
      S = 4.0e0*g_Pi*g_a0*g_a0*N*(g_Ry/B)*(g_Ry/B) ! Eq.(4)
      t0 = T/B*g_Mp/MSHI	! energy normalized to fit the data
      u0 = U/B  	! kinetic energy normalized
      w0 = w/B	! transferred energy normalized

      dSigma0 = dSigma_dw_int(S, t0, u0, 0.0d0) ! function see below
      dSigma = dSigma_dw_int(S, t0, u0, w0) - dSigma0
      dSigma_int_BEB_SHI = ZSHI*ZSHI*dSigma
!   endif
end function dSigma_int_BEB_SHI


function dSigma_int_BEB(T, w, B, U, N)
   real(8), intent(in) :: w 	! [eV] transferred energy
   real(8), intent(in) :: T 	! [eV] kinetic energy of the electron
   real(8), intent(in) ::  B	! bindning energy of atomic electron [eV]
   real(8), intent(in) ::  U	! mean kinetic energy of electron in the sub-shell [eV]
   real(8), intent(in) ::  N	! ocupation number of electrons in this shell
   real(8) :: dSigma_int_BEB    ! differential cross-section

   real(8) t0, u0, w0, S, Sigma, dSigma0

!   if (w .LE. B) then ! for the case if transferred energy is lower than the Ip
!      dSigma_int_BEB = 0.0d0
!   else
      S = 4.0e0*g_Pi*g_a0*g_a0*N*(g_Ry/B)*(g_Ry/B) ! Eq.(4)
      t0 = T/B	! energy normalized to the Rydberg constant ! Eq.(4)
      u0 = U/B  	! kinetic energy normalized
      w0 = w/B	! transferred energy normalized

      dSigma0 = dSigma_dw_int(S, t0, u0, 0.0d0) ! function see below
      dSigma_int_BEB = dSigma_dw_int(S, t0, u0, w0) - dSigma0
!   endif
end function dSigma_int_BEB

function dSigma_dw_int(S, t0, u0, w0)
   real(8) t0, w0, u0, S, dSigma_dw_int
   dSigma_dw_int = S/(t0+u0+1.0e0)*(-(log(w0+1.0e0)-log(abs(t0-w0)))/(t0+1.0e0) + (1.0e0/(t0-w0)-1.0e0/(w0+1.0e0)) + log(t0)*0.5e0*(1.0e0/((t0-w0)*(t0-w0)) - 1.0e0/((w0+1.0e0)*(w0+1.0e0))))
end function dSigma_dw_int

function dSigma_w_int_BEB_SHI(T, w, B, U, N, MSHI, ZSHI)
   real(8), intent(in) :: w 	! [eV] transferred energy
   real(8), intent(in) :: T 	! [eV] kinetic energy of the electron
   real(8), intent(in) ::  B	! bindning energy of atomic electron [eV]
   real(8), intent(in) ::  U	! mean kinetic energy of electron in the sub-shell [eV]
   real(8), intent(in) ::  N	! ocupation number of electrons in this shell
   real(8) :: MSHI, ZSHI    ! mass of charge of the SHI
   real(8) :: dSigma_w_int_BEB_SHI ! differential cross-section
   real(8) t0, u0, w0, S, Sigma, dSigma0

!   if (w .LE. B) then ! for the case if incident electron kinetic energy is lower than the Ip
!      dSigma_w_int_BEB_SHI = 0.0d0
!   else
      S = 4.0e0*g_Pi*g_a0*g_a0*N*(g_Ry/B)*(g_Ry/B) ! Eq.(4)
      t0 = T/B*g_me/MSHI	! energy normalized to the Rydberg constant ! Eq.(4)
      u0 = U/B  	! kinetic energy normalized
      w0 = w/B	    ! transferred energy normalized

      dSigma0 = dSigma_dw_w_int(S, t0, u0, 0.0d0) ! function see below
      dSigma_w_int_BEB_SHI = ZSHI*ZSHI*B*(dSigma_dw_w_int(S, t0, u0, w0) - dSigma0)
!   endif
end function dSigma_w_int_BEB_SHI


function dSigma_w_int_BEB(T, w, B, U, N)
   real(8), intent(in) :: w 	! [eV] transferred energy
   real(8), intent(in) :: T 	! [eV] kinetic energy of the electron
   real(8), intent(in) ::  B	! bindning energy of atomic electron [eV]
   real(8), intent(in) ::  U	! mean kinetic energy of electron in the sub-shell [eV]
   real(8), intent(in) ::  N	! ocupation number of electrons in this shell
   real(8) :: dSigma_w_int_BEB ! differential cross-section
   real(8) t0, u0, w0, S, Sigma, dSigma0

!   if (w .LE. B) then ! for the case if incident electron kinetic energy is lower than the Ip
!      dSigma_w_int_BEB = 0.0d0
!   else
      S = 4.0e0*g_Pi*g_a0*g_a0*N*(g_Ry/B)*(g_Ry/B) ! Eq.(4)
      t0 = T/B	! energy normalized to the Rydberg constant ! Eq.(4)
      u0 = U/B  ! kinetic energy normalized
      w0 = w/B	! transferred energy normalized

      dSigma0 = dSigma_dw_w_int(S, t0, u0, 0.0d0) ! function see below
      dSigma_w_int_BEB = B*(dSigma_dw_w_int(S, t0, u0, w0) - dSigma0)
!   endif
end function dSigma_w_int_BEB

function dSigma_dw_w_int(S, t0, u0, w0)
   real(8) t0, w0, u0, S, dSigma_dw_w_int
   real(8) A, B, C, tw, logwt, overwt, overw1, overt1
   tw = 1.0d0/(t0+u0+1.0d0)
   logwt = log((w0+1.0d0)*abs(w0-t0))
   overwt = 1.0d0/(w0-t0)
   overw1 = 1.0d0/(w0+1.0d0)
   overt1 = 1.0d0/(t0+1.0d0)
   A = tw*overt1*overt1*(t0*t0*log(abs(w0-t0)) + t0*logwt + log(w0+1.0d0))
   B = tw*( -t0*overwt + overw1 + logwt)
   C = tw*log(t0)*(0.5d0*(overw1*overw1 + t0*overwt*overwt) + overwt - overw1)
   dSigma_dw_w_int = S*(A+B+C)
end function dSigma_dw_w_int
!BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB_BEB

 
END MODULE Cross_sections

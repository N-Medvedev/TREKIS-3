!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains subroutines for dealing with EADL and EPDL97 databases

MODULE Dealing_with_EADL
use Universal_Constants   ! let it use universal constants
use Objects   ! since it uses derived types, it must know about them from module 'Objects'

implicit none


 contains

! Find which elements constitute given chemical formula:
subroutine Decompose_compound(Name, path_sep, Target_atoms, read_well)
   type(Atom), dimension(:), allocatable, intent(inout), optional :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   character(*), intent(in) :: Name ! compound name (SiO2, Al2O3, etc.)
   character(*), intent(in) :: path_sep ! path separator
   logical, intent(inout) :: read_well  ! did we read the file well?
   type(atomic_data), dimension(:), allocatable :: Periodic_table ! this is an internal module variable
   integer FN, Reason, leng, C, i, cur, num, N
   character(3), dimension(100) :: ElNames  ! all elements in the compound, assuming they are not more then 100
   real(8), dimension(100) :: ElPersent     ! persentage of this element in the compound
   integer, dimension(100) :: ElNumbers     ! numbers of elements in the compound
   character(100) :: Folder_name, File_name
   character(15) :: Full_name
   real(8) :: M
   character(3) :: El
   character(*), parameter :: numbers = '0123456789'
   CHARACTER(*), PARAMETER :: LowCase = 'abcdefghijklmnopqrstuvwxyz'
   CHARACTER(*), PARAMETER :: UpCase  = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   logical :: file_exists, num_vs_char, found_atom, file_opened
   read_well = .true.
   Folder_name = 'INPUT_EADL'  ! here we keep databases
   File_name = trim(adjustl(Folder_name))//trim(adjustl(path_sep))//'INPUT_atomic_data.dat' ! fixed name of the database
   inquire(file=trim(adjustl(File_name)),exist=file_exists) ! check if input file is there
   exists:if (file_exists) then
      open (newunit=FN, file=trim(adjustl(File_name)))
      call Count_lines_in_file(FN, N)
      if (.not.allocated(Periodic_table)) allocate(Periodic_table(N-1))
      read(FN,*,IOSTAT=Reason) ! skip first line with names
      do i = 1, N-1
         read(FN,*,IOSTAT=Reason) Periodic_table(i)%Z, Periodic_table(i)%Full_name, Periodic_table(i)%Name, Periodic_table(i)%Mass, Periodic_table(i)%Nvb
         !print*, Periodic_table(i)%Z, Periodic_table(i)%Full_name, Periodic_table(i)%Name, Periodic_table(i)%Mass, Periodic_table(i)%Nvb
         if (Reason /= 0) exit
      enddo
      if (Reason /= 0) goto 911
      
      leng = LEN(trim(adjustl(Name))) ! how many characters are in the name
      
      num = 0 ! how many different elements are in this compound
      El = '   ' ! start a new name
      do i = 1,leng ! compare all name character by character
         if (verify(trim(adjustl(Name(i:i))),trim(adjustl(numbers))) == 0) then ! it's an integer number
         ! it tells you how many of these atoms are in the compound
            read(Name(i:i),*) cur ! read it into an integer
            C = C*10 + cur
         else if (verify(trim(adjustl(Name(i:i))),trim(adjustl(UpCase))) == 0) then ! it's an upper-case latine letter
         ! the element name starts
            if (num > 0) then ! at least one element was found:
                if (C <= 0) then
                   ElPersent(num) = 1  ! there is at least one
                else
                   ElPersent(num) = C  ! that's how many of this element in the compound
                endif
                ElNames(num) = trim(adjustl(El))
            endif
            num = num + 1   ! new element found in the compound
            ! Start a new element:
            C = 0           ! to start
            El = '   '      ! start a new name
            El = trim(adjustl(El))//Name(i:i)   ! write it's name 
         else if (verify(trim(adjustl(Name(i:i))),trim(adjustl(LowCase))) == 0) then ! it's a lower-case latine letter
         ! the element name still goes on
            El = trim(adjustl(El))//Name(i:i)   ! continue writing it's name
         else ! it's another symbole
         ! no idea what that might be...
            print*, 'Symbol ', trim(adjustl(Name(i:i))), ' in the compound formula could not be identified'
         endif
      enddo
      if (num > 0) then ! the last element:
        if (C <= 0) then
           ElPersent(num) = 1  ! there is at least one
        else
           ElPersent(num) = C  ! that's how many of this element in the compound
        endif
        ElNames(num) = trim(adjustl(El))
      endif
       
      if (present(Target_atoms)) then
         if (.not. allocated(Target_atoms)) allocate(Target_atoms(num)) ! that's how many atom kinds we have
         do i = 1, num
            call find_atomic_number(ElNames(i), ElNumbers(i), Full_name, M, Periodic_table, read_well)
            Target_atoms(i)%Zat = ElNumbers(i)
            Target_atoms(i)%Pers = ElPersent(i)
            Target_atoms(i)%Name = ElNames(i) ! name of the element
            Target_atoms(i)%Full_Name = Full_Name  ! full name of the element
            Target_atoms(i)%Mass = M  ! mass of the element in the proton-mass units
         enddo
      else
         do i = 1, num
            call find_atomic_number(ElNames(i), ElNumbers(i), Periodic_table=Periodic_table, read_well=read_well)
            print*, ElNames(i), ElPersent(i), ElNumbers(i)
         enddo
      endif
   else exists
      print*, 'Could not find the file: ', trim(adjustl(File_name))
      read_well = .false.
   endif exists
   
911 if (allocated(Periodic_table)) deallocate(Periodic_table)
   if (Reason /= 0) then
      print*, 'Could not read input file: ', trim(adjustl(File_name))
      read_well = .false.
   endif
   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it
end subroutine Decompose_compound


subroutine find_atomic_number(Name, Z, Full_name, M, Periodic_table, read_well)
   character(*), intent(in) :: Name ! element abbreviation
   integer, intent(out) :: Z        ! atomic number
   type(atomic_data), dimension(:) :: Periodic_table ! this is an internal module variable
   character(15), intent(out), optional :: Full_name ! element name
   real(8), intent(out), optional  :: M ! mass
   logical, intent(inout) :: read_well  ! did we read the file well?
   integer N, i
   read_well = .true.
   N = size(Periodic_table) ! that's how many elements we have in our table
   do i = 1, N ! compare names
      if (trim(adjustl(Name)) == trim(adjustl(Periodic_table(i)%Name))) then
         Z = Periodic_table(i)%Z
         if (present(Full_name)) Full_name = Periodic_table(i)%Full_name
         if (present(M)) M = Periodic_table(i)%Mass
         exit
      endif
      if (i == N) then
         print*, trim(adjustl(Name)), ' - such an element was not found in our database...'
         read_well = .false.
      endif
   enddo
end subroutine find_atomic_number



subroutine Count_lines_in_file(File_num, N)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of lines in this file
    integer i
    i = 0
    do
        read(File_num, *, end=603)
        i = i + 1
    enddo
    603 continue
    rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
    N = i
end subroutine Count_lines_in_file


subroutine check_atomic_parameters(NumPar, Target_atoms, N_at, cur_shl, shl, Error_message, read_well) ! from module 'Dealing_with_EADL'
   type(Flag), intent(inout) :: NumPar ! numerical parameters
   type(Atom), dimension(:), allocatable, intent(inout) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   integer, intent(in), optional :: N_at
   integer, intent(in), optional :: cur_shl, shl  ! current atoms, shell, and total number of shells
   type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
   logical, intent(inout) :: read_well  ! did we read the file well?
   
   integer INFO, i, j, k, FN, FN1, Z, Shl_dsgtr, Nat
   character(100) :: File_name_EADL, Folder_name, Error_descript
   character(100) :: File_name_EPDL97, File_name, File_name2
   logical :: File_opened
   
   call Check_EADP_EPDL97(NumPar%path_sep, File_name_EADL, File_name_EPDL97, INFO) ! check if these database-files exist
   if (INFO .EQ. 0) then ! files are there, reading them:
      if (present(cur_shl)) then ! for this shell only:
         Z = Target_atoms(N_at)%Zat ! atomic number

         if (Target_atoms(N_at)%Nel(cur_shl) .LE. 0) then ! nondefined in the input file, find the default value from EADL:
            call READ_EADL_TYPE_FILE_int (FN1, File_name_EADL, Z, 912, Target_atoms, N_at, cur_shl, shl) ! read PQNs and shell-names
         endif
         
         if (Target_atoms(N_at)%Ip(cur_shl) .LE. -1.0d-14) then ! nondefined in the input file, find the default value from EADL:
            call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 913, Target_atoms(N_at)%Ip, cur_shl, shl, Target_atoms(N_at)%Shl_num(cur_shl)) ! read binding energies
         endif
         
         if (Target_atoms(N_at)%Ek(cur_shl) .LE. 0.0d0) then ! nondefined in the input file, find the default value from EADL:
            if (Target_atoms(N_at)%Shl_num(cur_shl) .GE. 62) then ! find approximate Ekin for the VB:
               call next_designator(Shl_dsgtr, Target_atoms(N_at)%Shl_num(cur_shl-1)) ! find the designator for the VB (next shell after last one)
            else
               Shl_dsgtr = Target_atoms(N_at)%Shl_num(cur_shl)
            endif
            !print*, Target_atoms(N_at)%Ek(cur_shl), Target_atoms(N_at)%Ip(cur_shl), Shl_dsgtr
            call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 914, Target_atoms(N_at)%Ek, cur_shl, shl, Shl_dsgtr) ! read kinetic energies
            !print*, Target_atoms(N_at)%Ek(cur_shl), 'Ek'
         endif
         
         if ((Target_atoms(N_at)%Shl_num(cur_shl) .GE. 63) .OR. (.not.NumPar%include_photons)) then
            Target_atoms(N_at)%Radiat(cur_shl)=1d23 ! [fs] not possible process (or not included) => infinite time
         else
                 !READ_EADL_TYPE_FILE_real(FN, File_name, Z_needed, I_needed, Array1, cur_shl, shl_tot, Shl_dsgtr) 
            call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 921, Target_atoms(N_at)%Radiat, cur_shl, shl, Target_atoms(N_at)%Shl_num(cur_shl)) ! read radiative decay-times
            Target_atoms(N_at)%Radiat(cur_shl)=1d15*g_h/(g_e*Target_atoms(N_at)%Radiat(cur_shl)) ! [fs]
            print*, N_at, cur_shl, Target_atoms(N_at)%Radiat(cur_shl)
         endif
         
         if (Target_atoms(N_at)%Shl_num(cur_shl) .GE. 63) then
            Target_atoms(N_at)%Auger(cur_shl)=1.0d23 ! [fs] not possible process => infinite time
         else if ((Target_atoms(N_at)%Auger(cur_shl) .LE. 0.0d0) .OR. (Target_atoms(N_at)%Auger(cur_shl) .GT. 1.0d30)) then
            call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 922, Target_atoms(N_at)%Auger, cur_shl, shl, Target_atoms(N_at)%Shl_num(cur_shl)) ! read auger-times
            Target_atoms(N_at)%Auger(cur_shl)=1d15*g_h/(g_e*Target_atoms(N_at)%Auger(cur_shl)) ! [fs]
         endif
      else ! for all shells of this atom:
         Nat = size(Target_atoms)
         do i = 1, Nat
            Z = Target_atoms(i)%Zat ! atomic number
            call READ_EADL_TYPE_FILE_int (FN1, File_name_EADL, Z, 912, Target_atoms, i) ! read PQNs and shell-names
            call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 913, Target_atoms(i)%Ip) ! read binding energies
            call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 914, Target_atoms(i)%Ek) ! read kinetic energies
            call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 921, Target_atoms(i)%Radiat) ! read radiative decay-times
            call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 922, Target_atoms(i)%Auger) ! read auger-times
            !print*, 'Atom:',i, Nat, size(Target_atoms(i)%Ip), size(Target_atoms(i)%Ek), size(Target_atoms(i)%Radiat), size(Target_atoms(i)%Auger)
            do j = 1,size(Target_atoms(i)%Ip)
               !print*, 'Shell',j, size(Target_atoms(i)%Ip)
               Target_atoms(i)%Auger(j)=1d15*g_h/(g_e*Target_atoms(i)%Auger(j)) ! [fs]
               if ((Target_atoms(i)%Shl_num(j) .GE. 63) .OR. (.not.NumPar%include_photons)) then
                  Target_atoms(i)%Radiat(j)=1d23 ! [fs] not possible process (or not included) => infinite time
               else
                  Target_atoms(i)%Radiat(j)=1d15*g_h/(g_e*Target_atoms(i)%Radiat(j)) ! [fs]
               endif
            enddo ! j
         enddo ! i
      endif
   else ! couldn't find databases:
      Folder_name = 'INPUT_EADL'  ! here we keep databases
      if (INFO .EQ. 1) then
         File_name = trim(adjustl(Folder_name))//trim(adjustl(NumPar%path_sep))//'eadl.all'
         Error_descript = 'File '//trim(adjustl(File_name))//' is not found!'    ! description of an error
         call Save_error_details(Error_message, 21, Error_descript) ! write it into the error-log file
         print*, trim(adjustl(Error_descript)) ! print it also on the sreen
      endif
      if (INFO .EQ. 2) then
         File_name2 = trim(adjustl(Folder_name))//trim(adjustl(NumPar%path_sep))//'epdl97.all'
         Error_descript = 'File '//trim(adjustl(File_name2))//' is not found!'    ! description of an error
         call Save_error_details(Error_message, 22, Error_descript) ! write it into the error-log file
         print*, trim(adjustl(Error_descript)) ! print it also on the sreen
      endif
      read_well = .false.   ! it didn't read well the input file...
   endif
   inquire(file=trim(adjustl(File_name_EADL)),opened=File_opened)
   if (File_opened) close(FN1)   ! close the EADL database
end subroutine check_atomic_parameters


! Reading EADL data for case of integer arrays:
subroutine READ_EADL_TYPE_FILE_int(FN, File_name, Z_needed, I_needed, Target_atoms, N_at, Shl, Shl_tot)
   integer, intent (inout) :: FN
   character(100), intent(in) :: File_name
   integer, intent(in) :: Z_needed, I_needed 
   type(Atom), dimension(:), allocatable, intent(inout) :: Target_atoms  ! define target atoms as objects, we don't know yet how many they are
   integer, intent(in) :: N_at ! number of atom
   integer, intent(in), optional :: Shl, Shl_tot  ! shell number
   integer Z, A, Yi, Yo, Date, C, I, S, EndCheck, counter
   real(8) AW, X1
   real(8) READ1, READ2, READ3, READ4
   integer run_i, Reason, imin, imax
   logical File_opened

   inquire(file=trim(adjustl(File_name)),opened=File_opened)
   if (.not. File_opened) then
      !FN = 336
      open(newunit=FN, FILE = trim(adjustl(File_name)), status = 'old')
   endif

   Z = 0
   do while (Z .NE. Z_needed)
      call read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
      call How_many_lines(FN, counter, EndCheck) ! to know how many lines are there
      if (Z .GE. 100) print*, 'The element is not found in the EADL database...'
      if (Z .GE. 100) exit
   enddo

   if (Z .EQ. Z_needed) then ! this is the element we were looking for
      do while (I .NE. I_needed) ! find the value that we need
         call read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
         call How_many_lines(FN, counter, EndCheck) ! to know how many lines are there
         if (Z .GT. Z_needed) print*, 'INT: The I-value is not found in the EADL database...'
         if (Z .GT. Z_needed) exit
      enddo
      do run_i = 1, counter+1
         backspace(FN)
      enddo

      if (present(Shl)) then ! only for this shell:
         if (.not. allocated(Target_atoms(N_at)%Shell_name)) then
            allocate(Target_atoms(N_at)%Shell_name(Shl_tot)) ! now we know the dimensions...
         endif
         if (.not. allocated(Target_atoms(N_at)%PQN)) then
            allocate(Target_atoms(N_at)%PQN(Shl_tot)) ! now we know the dimensions...
         endif
         READ1 = 0
         if (Target_atoms(N_at)%Shl_num(Shl) .LT. 63) then
            do while (READ1 .NE. Target_atoms(N_at)%Shl_num(Shl))
               read(FN,9991,IOSTAT=Reason) READ1, READ2 ! Shell designator, No. of electrons in the shell
               if (READ1 .EQ. Target_atoms(N_at)%Shl_num(Shl)) exit
               IF (Reason .LT. 0) THEN ! end of file reached, nothing found
                  rewind(FN)
                  exit
               endif
            enddo
            if (Reason .LT. 0) then ! try with the next value:
                Z = 0
                do while (Z .NE. Z_needed)
                   call read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
                   call How_many_lines(FN, counter, EndCheck) ! to know how many lines are there
                   if (Z .GE. 100) print*, 'The element is not found in the EADL database...'
                   if (Z .GE. 100) exit
                enddo
                do while (I .NE. I_needed) ! find the value that we need
                   call read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
                   call How_many_lines(FN, counter, EndCheck) ! to know how many lines are there
                   if (Z .GT. Z_needed) print*, 'INT: The I-value is not found in the EADL database...'
                   if (Z .GT. Z_needed) exit
                enddo
                do run_i = 1, counter+1
                   backspace(FN)
                enddo
               
               call select_imin_imax(imin, imax, Target_atoms(N_at)%Shl_num(Shl))
               
               READ3 = 0.0d0
               READ4 = 0.0d0
               do run_i = imin, imax
                  read(FN,9991,IOSTAT=Reason) READ1, READ2 ! Shell designator, No. of electrons in the shell
                  if (READ1 .GE. imin) then
                     READ3 = READ3 + READ1 ! number of electrons in this shell
                     READ4 = READ4 + READ2 ! number of electrons in this shell
                     !print*, READ1, READ4, READ2
                  endif
                  if (READ1 .EQ. imax) exit
               enddo
               Target_atoms(N_at)%Nel(Shl) = READ4 ! number of electrons in this shell
               call define_PQN(Target_atoms(N_at)%Shl_num(Shl), Target_atoms(N_at)%Shell_name(Shl), Target_atoms(N_at)%PQN(Shl))
            else
               Target_atoms(N_at)%Nel(Shl) = READ2 ! number of electrons in this shell
               call define_PQN(INT(READ1), Target_atoms(N_at)%Shell_name(Shl), Target_atoms(N_at)%PQN(Shl))
            endif
         else
            call define_PQN(Target_atoms(N_at)%Shl_num(Shl), Target_atoms(N_at)%Shell_name(Shl), Target_atoms(N_at)%PQN(Shl))
         endif
      else ! for all shells:
         Target_atoms(N_at)%N_shl = counter ! that's how many shells we have in this atoms
         if (.not. allocated(Target_atoms(N_at)%Nel)) allocate(Target_atoms(N_at)%Nel(counter)) ! allocate number of electrons
         if (.not. allocated(Target_atoms(N_at)%Shell_name)) allocate(Target_atoms(N_at)%Shell_name(counter)) ! allocate shell names
         if (.not. allocated(Target_atoms(N_at)%PQN)) allocate(Target_atoms(N_at)%PQN(counter)) ! allocate principal quantum number
         if (.not. allocated(Target_atoms(N_at)%Shell_name)) allocate(Target_atoms(N_at)%Shell_name(counter)) ! allocate shell-names for each shell
         if (.not. allocated(Target_atoms(N_at)%Shl_num)) allocate(Target_atoms(N_at)%Shl_num(counter)) ! allocate shell disignator for each shell
         if (.not. allocated(Target_atoms(N_at)%Ip)) allocate(Target_atoms(N_at)%Ip(counter)) ! allocate ionization potentials
         if (.not. allocated(Target_atoms(N_at)%Ek)) allocate(Target_atoms(N_at)%Ek(counter)) ! allocate mean kinetic energies of the shells
         if (.not. allocated(Target_atoms(N_at)%KOCS)) allocate(Target_atoms(N_at)%KOCS(counter)) ! allocate mean kinetic energies of the shells
         if (.not. allocated(Target_atoms(N_at)%KOCS_SHI)) allocate(Target_atoms(N_at)%KOCS_SHI(counter)) ! allocate mean kinetic energies of the shells
         if (.not. allocated(Target_atoms(N_at)%Radiat)) then
            allocate(Target_atoms(N_at)%Radiat(counter)) ! allocate auger-times
            Target_atoms(N_at)%Radiat = 1d-24 ! to be inversed later
         endif
         if (.not. allocated(Target_atoms(N_at)%Auger)) then
            allocate(Target_atoms(N_at)%Auger(counter)) ! allocate radiative times
            Target_atoms(N_at)%Auger = 1d-24 ! to be inversed later
         endif
         do run_i = 1, counter
            read(FN,9991) READ1, READ2 ! Shell designator; No. of electrons in the shell
            Target_atoms(N_at)%Nel(run_i) = READ2 ! number of electrons in this shell
            Target_atoms(N_at)%Shl_num(run_i) = READ1 ! shell designator 
            call define_PQN(INT(READ1), Target_atoms(N_at)%Shell_name(run_i), Target_atoms(N_at)%PQN(run_i))
         enddo
      endif
   endif
   flush(FN)
   rewind(FN)
9999   format(I3,I3,I2,I2,E11.4,I6)	! 1st line
9998   format(I2,I3,I3,E11.4)		! 2st line
9997   format(71X,I1)			! last line
9991   format(6E11.4,6E11.4,6E11.4,6E11.4) ! 3d line if C=91,92
end subroutine READ_EADL_TYPE_FILE_int


! Reading EADL data for case of real arrays:
subroutine READ_EADL_TYPE_FILE_real(FN, File_name, Z_needed, I_needed, Array1, cur_shl, shl_tot, Shl_dsgtr) 
!    call READ_EADL_TYPE_FILE_real(FN1, File_name_EADL, Z, 922, Target_atoms(N_at)%Auger, N_at, cur_shl, shl) ! read auger-times
   integer, intent(inout) :: FN
   character(100), intent(in) :: File_name
   integer, intent(in) :: Z_needed, I_needed
   real(8), dimension(:), allocatable, intent(inout) :: Array1 ! array with variables that we need
   integer, intent(in), optional :: cur_shl, shl_tot, Shl_dsgtr  ! shell number, number of shells, shell designator
   
   integer Z, A, Yi, Yo, Date, C, I, S, EndCheck, counter
   real(8) AW, X1
   real(8) READ1, READ2, READ3, READ4
   integer run_i, imax, imin, Reason, icont
   logical File_opened, Exist_val

   Exist_val = .true.

   inquire(file=trim(adjustl(File_name)),opened=File_opened)
   if (.not. File_opened) then
      open(newunit=FN, FILE = trim(adjustl(File_name)), status = 'old') ! EADL
   endif

   Z = 0
   do while (Z .NE. Z_needed)
      call read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
      call How_many_lines(FN, counter, EndCheck) ! to know how many lines are there
      if (Z .GE. 100) print*, 'The element is not found in the EADL database...'
      if (Z .GE. 100) exit
   enddo

   if (Z .EQ. Z_needed) then ! this is the element we were looking for
      do while (I .NE. I_needed) ! find the value that we need
         call read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
         call How_many_lines(FN, counter, EndCheck) ! to know how many lines are there
         if (Z .GT. Z_needed) then
            print*, 'The I-value is not found in the EADL database:'
            print*, Z, Z_needed, I_needed
            Exist_val = .false.
         endif
         if (Z .GT. Z_needed) exit
      enddo
      do run_i = 1, counter+1
         backspace(FN)
      enddo

      if (present(cur_shl)) then ! this shell only:
         if (.not. allocated(Array1)) then
            allocate(Array1(shl_tot)) ! now we know the dimensions...
            Array1 = 0.0d0
         endif
         if (.not.Exist_val) then ! if there is no value in the database, nothing else to do
            Array1(cur_shl) = 1.0d-30
            goto 5555 ! just exit then...
         endif
         if (Shl_dsgtr .LT. 63) then
            READ1 = 0
            do while (READ1 .NE. Shl_dsgtr)
               read(FN,9991, IOSTAT=Reason) READ1, READ2
               if (READ1 .EQ. Shl_dsgtr) exit
               IF (Reason .LT. 0) THEN ! end of file reached, nothing found
                  rewind(FN)
                  exit
               endif
            enddo
             if (Reason .LT. 0) then ! If there is no such shell in the database, try to sum up sub-shells:
                Z = 0
                do while (Z .NE. Z_needed)
                   call read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
                   call How_many_lines(FN, counter, EndCheck) ! to know how many lines are there
                   if (Z .GE. 100) print*, 'The element is not found in the EADL database...'
                   if (Z .GE. 100) exit
                enddo
                do while (I .NE. I_needed) ! find the value that we need
                   call read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
                   call How_many_lines(FN, counter, EndCheck) ! to know how many lines are there
                   if (Z .GT. Z_needed) print*, 'INT: The I-value is not found in the EADL database...'
                   if (Z .GT. Z_needed) exit
                enddo
                do run_i = 1, counter+1
                   backspace(FN)
                enddo
               READ3 = 0.0d0
               call select_imin_imax(imin, imax, Shl_dsgtr) ! sum subshells:
               
               READ3 = 0.0d0
               icont = 0
               do run_i = imin, imax
                  read(FN,9991,IOSTAT=Reason) READ1, READ2 ! Shell designator, No. of electrons in the shell
                  if (READ1 .GE. imin) then
                     if (isnan(READ2)) then
                     else
                        READ3 = READ3 + READ2 ! number of electrons in this shell
                     endif
                     icont = icont + 1
                  endif
                  if (READ1 .EQ. imax) exit
               enddo
               Array1(cur_shl) = (READ3/real(icont))*1d6
            else
               Array1(cur_shl) = READ2*1d6
            endif
         else
            Array1(cur_shl) = 1d22
         endif
         if (isnan(Array1(cur_shl))) Array1(cur_shl) = 0.0d0 ! just in case...
      else ! for all shells:
         if (.not. allocated(Array1)) then
            allocate(Array1(counter)) ! now we know the dimensions...
            Array1 = 0.0d0
         endif
         if (Exist_val) then	! if the value is found in the database, use it, if not leave default
            do run_i = 1, counter
               read(FN,9991) READ1, READ2
               Array1(run_i) = READ2*1d6
            enddo
         else
            Array1 = 1.0d-30 ! if the value does not exists
         endif
      endif
   endif
   flush(FN)
5555   rewind(FN)
9999   format(I3,I3,I2,I2,E11.4,I6)	! 1st line
9998   format(I2,I3,I3,E11.4)		! 2st line
9997   format(71X,I1)			! last line
9991   format(6E11.4,6E11.4,6E11.4,6E11.4) ! 3d line if C=91,92
end subroutine READ_EADL_TYPE_FILE_real


subroutine select_imin_imax(imin, imax, shl_dsgntr)
   integer, intent(in) :: shl_dsgntr
   integer, intent(out) :: imin, imax ! according to EADL database, those are corresponding shells:
   select case (shl_dsgntr) ! sum subshells:
   case (2)
      imin = 3
      imax = 6
   case (4)
      imin = 5
      imax = 6
   case (7)
      imin = 8
      imax = 14
   case (9)
      imin = 10
      imax = 11
   case (12)
      imin = 13
      imax = 14
   case (15)
      imin = 16
      imax = 25
   case (26)
      imin = 27
      imax = 39
   case (40)
      imin = 41
      imax = 56
   case (57)
      imin = 58
      imax = 61
   case default
      imin = 0
      imax = 0
   endselect
end subroutine select_imin_imax


subroutine next_designator(Shl_dsgtr, Last_shl) ! find the designator for the VB (next shell after last one)
   integer, intent(in) :: Last_shl
   integer, intent(out) :: Shl_dsgtr
   select case (Last_shl) ! sum subshells:
   case (:1)
      Shl_dsgtr = 2
   case (2:6)
      Shl_dsgtr = 7
   case (7:14)
      Shl_dsgtr = 15
   case (15:25)
      Shl_dsgtr = 26
   case (26:39)
      Shl_dsgtr = 40
   case (40:)
      Shl_dsgtr = 57
   endselect
end subroutine next_designator



! EPDL97_EPDL97_EPDL97_EPDL97_EPDL97_EPDL97_EPDL97_EPDL97_EPDL97_EPDL97
subroutine get_photon_cross_section_EPDL(Ele, Target_atoms, Nat, Shl_dsgntr, Nshl, Matter, Sigma)
    real(8), intent(in) :: Ele  ! electron energy [eV]
    type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
    integer, intent(in) :: Nat ! number of atom
    integer, intent(in), optional :: Nshl, Shl_dsgntr    ! shell number, shell designator
    type(Solid), intent(in) :: Matter ! properties of material
    real(8), intent(out) :: Sigma ! calculated inverse mean free path (out of cross section) [1/A]
    integer Z
    integer FN2    
    character(100) File_name, Folder_name
    Folder_name = 'INPUT_EADL'  ! here we keep databases
    File_name = trim(adjustl(Folder_name))//'/epdl97.all'
    Z = Target_atoms(Nat)%Zat ! atomic number
    FN2 = 5558
    shl_sign:if (present(Shl_dsgntr)) then
       if (Shl_dsgntr .GE. 62) then ! VB
          if (present(Nshl)) then
             call READ_EPDL_TYPE_FILE_real(FN2,File_name, Z, 73, 0, 91, Sigma, Ele, Shl_dsgntr, Last_dsgtr=Target_atoms(Nat)%Shl_num(Nshl-1)) ! by subshells
          else
             call READ_EPDL_TYPE_FILE_real(FN2,File_name, Z, 73, 0, 91, Sigma, Ele, Shl_dsgntr) ! by subshells
          endif
       else ! correct shell designator for the shell forming VB has already been set in the parent subroutine:
          call READ_EPDL_TYPE_FILE_real(FN2,File_name, Z, 73, 0, 91, Sigma, Ele, Shl_dsgntr) ! by subshells
       endif
    else shl_sign ! Nshl is present
       if (present(Nshl)) then
          call READ_EPDL_TYPE_FILE_real(FN2,File_name, Z, 73, 0, 91, Sigma, Ele, Target_atoms(Nat)%Shl_num(Nshl)) ! by subshells
       else
          print*, 'Niether Nshl nor Shl_dsgntr is present in get_photon_cross_section_EPDL'
          pause 'What to do?'
       endif
    endif shl_sign
end subroutine get_photon_cross_section_EPDL


! Reading data from EPDL file: photoabsrotbtion cross-section for given photon energy
subroutine READ_EPDL_TYPE_FILE_real(FN2, File_name, Z_needed, C_needed, I_needed, S_needed, Photoabs_sigma, Eph, Shl_design, Last_dsgtr, Photoabs_sigma_tot)
   integer, intent (inout) :: FN2
   character(100), intent(in) :: File_name
   integer, intent(in) :: Z_needed, C_needed, I_needed, S_needed
   real(8), intent(in) :: Eph ! photon energy [eV]
   !type(Atom), dimension(:), intent(in) :: Target_atoms  ! all data for target atoms
   !integer, intent(in) :: Nat, Shl ! number of atom and a shell
   integer, intent(in) :: Shl_design ! shell designator - for this shell we are looking the cross section
   integer, intent(in), optional :: Last_dsgtr ! last shell which wasn't VB
   real(8), intent(inout) :: Photoabs_sigma ! Cross section for this shell [A^2]
   real(8), intent(inout), optional ::  Photoabs_sigma_tot  ! Total cross section [A^2]
   integer Z, A, Yi, Yo, Date, Iflag, C, I, S, EndCheck, counter, run_i
   real(8) AW, X1
   integer imin, imax, Shl_dsgtr
   real(8) READ1, READ2, READ3, READ4, E_photon
   logical File_opened

   inquire(file=trim(adjustl(File_name)),opened=File_opened)
   if (.not. File_opened) then
      !FN2 = 338
      open(FN2, FILE = trim(adjustl(File_name)), status = 'old') ! EPDL97
   endif

   E_photon = Eph/1e6 ! [MeV] units in the EPDL-database
   Photoabs_sigma = 0.0d0
   Shl_dsgtr = Shl_design
   if (present(Last_dsgtr)) then
       if (Shl_design .GE. 63) then
          call next_designator(Shl_dsgtr, Last_dsgtr) ! find the designator for the VB (next shell after last one)
       else
          Shl_dsgtr = Shl_design
       endif
   endif
   call select_imin_imax(imin, imax, Shl_dsgtr) ! find subshell designators
   
   Z = 0
   do while (Z .NE. Z_needed)
      !print*, 'Z:', Z, Z_needed
      call read_throu_EPDL(FN2, Z, A, Yi, Yo, AW, Date, Iflag, C, I, S, X1, EndCheck)
      call How_many_lines(FN2, counter, EndCheck) ! to know how many lines are there
      if (Z .GE. 100) print*, 'The element is not found in the EPDL97 database...'
      if (Z .GE. 100) exit
   enddo

   if (Z .EQ. Z_needed) then ! this is the element we were looking for
     run_i = 0.0e0
     do while (Z .LT. Z_needed + 1) ! do throughout the whole data for this element: all shells
        call read_throu_EPDL(FN2, Z, A, Yi, Yo, AW, Date, Iflag, C, I, S, X1, EndCheck)
        if ((C .EQ. C_needed) .AND. (I .EQ. I_needed) .AND. (S .EQ. S_needed)) then ! read it:
           SELECT CASE (S_needed) ! which data are these:
              CASE (0) ! total cross-section
                 if (present(Photoabs_sigma_tot)) then
                    call Find_value_while_reading(FN2, Iflag, counter, EndCheck, E_photon, READ2)
                    Photoabs_sigma_tot = READ2*1e-8 ! [A^2]
                 endif
              CASE (91) ! cross-section by shell
                 call Find_value_while_reading(FN2, Iflag, counter, EndCheck, E_photon, READ2)
                 run_i = run_i + 1
                 !write(*,'(a,f7.2,i,i,i)') 'run_i', X1, Shl_dsgtr, imin, imax
                 if ((X1 .EQ. Shl_dsgtr) .OR. ( (X1 .GE. imin) .AND. (X1 .LE. imax) ) )then
                     !print*, 'summed', Photoabs_sigma, READ2
                     Photoabs_sigma = Photoabs_sigma + READ2*1e-8 ! [A^2]
                 endif
              CASE DEFAULT
                 call How_many_lines(FN2, counter, EndCheck) ! to know how many lines are there
                 print*, 'God knows what is going on...'
           END SELECT
        else ! just skip all these lines:
           call How_many_lines(FN2, counter, EndCheck) ! to know how many lines are there
        endif ! ((C .EQ. C_needed) .AND. (I .EQ. I_needed) .AND. (S .EQ. S_needed))
     enddo !while (Z .LT. Z_needed + 1) ! do throughout the whole data for this element: all shells
   endif ! (Z .EQ. Z_needed)
   flush(FN2)
   rewind(FN2)
9999   format(I3,I3,I2,I2,E11.4,I6,I1)	! 1st line
9998   format(I2,I3,I3,E11.4)		! 2st line
9997   format(71X,I1)			! last line
9991   format(6E11.4,6E11.4,6E11.4,6E11.4) ! 3d line if C=91,92
end subroutine READ_EPDL_TYPE_FILE_real


subroutine Check_EADP_EPDL97(path_sep, File_name, File_name2, INFO)
    character(1), intent(in) :: path_sep
    character(100), intent(inout) :: File_name
    character(100), intent(inout) :: File_name2
    integer INFO
    character(100) :: Folder_name
    logical file_exist
    Folder_name = 'INPUT_EADL'  ! here we keep databases
    File_name = trim(adjustl(Folder_name))//trim(adjustl(path_sep))//'eadl.all'
    File_name2 = trim(adjustl(Folder_name))//trim(adjustl(path_sep))//'epdl97.all'
    inquire(file=trim(adjustl(File_name)),exist=file_exist) ! check if input file is there
    if (.not. file_exist) then
       write(*,'(a,a,a)') 'File ',trim(adjustl(File_name)),' is not found.'
       INFO = 1
    else
       INFO = 0
    endif
    inquire(file=trim(adjustl(File_name2)),exist=file_exist) ! check if input file is there
    if (.not. file_exist) then
       write(*,'(a,a,a)') 'File ',trim(adjustl(File_name2)),' is not found.'
       INFO = 2
    endif
end subroutine Check_EADP_EPDL97


! Reading first 2 lines with the descriptions in EPDL database: 
subroutine read_throu_EPDL(FN, Z, A, Yi, Yo, AW, Date1, Iflag, C, I, S, X1, EndCheck)
   integer, intent (inout) :: FN, Z, A, Yi, Yo, Date1, Iflag, C, I, S, EndCheck
   real(8), intent (inout) :: AW, X1
   real(8) skipread
   EndCheck = 0
   read(FN,'(i3, i3,i2, i2, E11.4, 2x, i8, i1)') Z, A, Yi, Yo, AW, Date1, Iflag  	! 1st line
   !read(FN,'(i2,i3,i3,E11.4,E11.4)') C, I, S, skipread, X1  			! 2st line
   read(FN,'(i2,i3,i3,13X,E11.4)') C, I, S, X1  			! 2st line
   read(FN,9997) EndCheck  			! 3st line, cheking...
9009   format(i3,i3,i2,i2,e11.4,i6)	! 1st line
9999   format(I3,I3,I2,I2,E11.4,I6,I1)	! 1st line
9998   format(I2,I3,I3,E11.4)		! 2st line
9997   format(71X,I1)			! last line
end subroutine read_throu_EPDL


! Reading first 2 lines with the descriptions in EADL database:
subroutine read_throu(FN, Z, A, Yi, Yo, AW, Date, C, I, S, X1, EndCheck)
   integer, intent (inout) :: FN, Z, A, Yi, Yo, Date, C, I, S, EndCheck
   real(8), intent (inout) :: AW, X1
   EndCheck = 0
   read(FN,9999) Z, A, Yi, Yo, AW, Date  	! 1st line
   read(FN,9998) C, I, S, X1  			! 2st line
   read(FN,9997) EndCheck  			! 3st line, cheking...
9999   format(I3,I3,I2,I2,E11.4,I6)	! 1st line
9998   format(I2,I3,I3,E11.4)		! 2st line
9997   format(71X,I1)			! last line
end subroutine read_throu 


! Counts how many lines contain the data in a block of EADL, EPDL-databases
subroutine How_many_lines(FN,counter, EndCheck)
   integer, intent(in) :: FN
   integer, intent(inout) :: counter, EndCheck
   counter = 0
   do while (EndCheck .NE. 1) ! until end of data for this sub-shell/input
      counter = counter + 1 ! how many lines are there
      read(FN,9997) EndCheck  			! 3st line, cheking...
      if (counter .GE. 5000) EndCheck = -1	! indicates some ERROR
      if (counter .GE. 5000) exit	! if something went wrong...
   enddo
9997   format(71X,I1)			! last line
end subroutine How_many_lines


! EADL database, that's how to find the cross-section for given photon energy for each subshell:
subroutine Find_value_while_reading(FN, Iflag, counter, EndCheck, E_photon, OUT_value)
   integer, intent(in) :: FN, Iflag
   integer, intent(inout) :: counter, EndCheck
   real(8), intent(in) :: E_photon
   real(8), intent(out) :: OUT_value
   real(8) READ1, READ2, E1, E2, Sigma1, Sigma2
   LOGICAL :: Firsttime
   Firsttime = .true.
   counter = 0
   do while (EndCheck .NE. 1) ! until end of data for this sub-shell/input
      counter = counter + 1 ! how many lines are there
      read(FN,9997) EndCheck  			! 3st line, cheking...
      if (EndCheck .NE. 1) then
         backspace(FN)
         read(FN,9991) READ1, READ2
      endif
      if ((READ1 .GT. E_photon) .AND. (Firsttime)) then
         Firsttime = .false.
         backspace(FN)
         backspace(FN)
         read(FN,9991) READ1, READ2
         if (counter .EQ. 1) then ! too small energy!
            OUT_value = 0.0e0
         else ! photon energy is above the ionization potential
            E1 = READ1
            Sigma1 = READ2
            read(FN,9991) READ1, READ2
            E2 = READ1
            Sigma2 = READ2
            call Interpolate_EPDL(Iflag, E1, E2, Sigma1, Sigma2, E_photon, OUT_value)
         endif
      endif
      if (counter .GE. 5000) EndCheck = -1	! indicates some ERROR
      if (counter .GE. 5000) exit	! if something went wrong...
   enddo
9991   format(6E11.4,6E11.4,6E11.4,6E11.4) ! 3d line if C=91,92
9997   format(71X,I1)			! last line
end subroutine Find_value_while_reading


! Interpolation of photoabsorbtion cross-section according to EADL database:
subroutine Interpolate_EPDL(Iflag, E1, E2, Sigma1, Sigma2, E_needed, OUT_value)
   integer, intent(in) :: Iflag
   real(8), intent(in) :: E1, E2, Sigma1, Sigma2, E_needed
   real(8), intent(out) :: OUT_value
   real(8) E2log, E1log, E_needed_log, Sigma1log, Sigma2log
   select case(Iflag) ! what interpolation to use:
      case(0,2) ! linear x and y
         OUT_value = Sigma1 + (Sigma2 - Sigma1)/(E2 - E1)*(E_needed - E1)
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
end subroutine Interpolate_EPDL


! According to EADL-data format, the atomic shells are numerated as:
subroutine define_PQN(READ1, Shell_name, PQN)
   integer, intent(in) :: READ1 ! shell designator
   integer, INTENT(inout), optional :: PQN !Principal quantum number
   character(11), intent(inout), optional :: Shell_name ! names
   SELECT CASE(INT(READ1))
   CASE ( : 1)
      if (present(PQN)) PQN = 1 ! K-shell
      if (present(Shell_name)) Shell_name = 'K-shell'
   CASE (2)
      if (present(PQN)) PQN = 2 ! L-shell
      if (present(Shell_name)) Shell_name = 'L-shell'
   CASE (3)
      if (present(PQN)) PQN = 2 ! L-shell
      if (present(Shell_name)) Shell_name = 'L1-shell'
   CASE (4)
      if (present(PQN)) PQN = 2 ! L-shell
      if (present(Shell_name)) Shell_name = 'L23-shell'
   CASE (5)
      if (present(PQN)) PQN = 2 ! L-shell
      if (present(Shell_name)) Shell_name = 'L2-shell'
   CASE (6)
      if (present(PQN)) PQN = 2 ! L-shell
      if (present(Shell_name)) Shell_name = 'L3-shell'
   CASE (7)
      if (present(PQN)) PQN = 3 ! M-shell
      if (present(Shell_name)) Shell_name = 'M-shell'
   CASE (8)
      if (present(PQN)) PQN = 3 ! M-shell
      if (present(Shell_name)) Shell_name = 'M1-shell'
   CASE (9)
      if (present(PQN)) PQN = 3 ! M-shell
      if (present(Shell_name)) Shell_name = 'M23-shell'
   CASE (10)
      if (present(PQN)) PQN = 3 ! M-shell
      if (present(Shell_name)) Shell_name = 'M2-shell'
   CASE (11)
      if (present(PQN)) PQN = 3 ! M-shell
      if (present(Shell_name)) Shell_name = 'M3-shell'
   CASE (12)
      if (present(PQN)) PQN = 3 ! M-shell
      if (present(Shell_name)) Shell_name = 'M45-shell'
   CASE (13)
      if (present(PQN)) PQN = 3 ! M-shell
      if (present(Shell_name)) Shell_name = 'M4-shell'
   CASE (14)
      if (present(PQN)) PQN = 3 ! M-shell
      if (present(Shell_name)) Shell_name = 'M5-shell'
   CASE (15)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N-shell'
   CASE (16)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N1-shell'
   CASE (17)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N23-shell'
   CASE (18)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N2-shell'
   CASE (19)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N3-shell'
   CASE (20)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N45-shell'
   CASE (21)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N4-shell'
   CASE (22)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N5-shell'
   CASE (23)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N67-shell'
   CASE (24)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N6-shell'
   CASE (25)
      if (present(PQN)) PQN = 4 ! N-shell
      if (present(Shell_name)) Shell_name = 'N7-shell'
   CASE (26)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O-shell'
   CASE (27)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O1-shell'
   CASE (28)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O23-shell'
   CASE (29)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O2-shell'
   CASE (30)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O3-shell'
   CASE (31)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O45-shell'
   CASE (32)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O4-shell'
   CASE (33)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O5-shell'
   CASE (34)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O67-shell'
   CASE (35)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O6-shell'
   CASE (36)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O7-shell'
   CASE (37)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O89-shell'
   CASE (38)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O8-shell'
   CASE (39)
      if (present(PQN)) PQN = 5 ! O-shell
      if (present(Shell_name)) Shell_name = 'O9-shell'
   CASE (40)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P-shell'
   CASE (41)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P1-shell'
   CASE (42)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P23-shell'
   CASE (43)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P2-shell'
   CASE (44)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P3-shell'
   CASE (45)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P45-shell'
   CASE (46)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P4-shell'
   CASE (47)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P5-shell'
   CASE (48)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P67-shell'
   CASE (49)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P6-shell'
   CASE (50)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P7-shell'
   CASE (51)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P89-shell'
   CASE (52)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P8-shell'
   CASE (53)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P9-shell'
   CASE (54)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P1011-shell'
   CASE (55)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P10-shell'
   CASE (56)
      if (present(PQN)) PQN = 6 ! P-shell
      if (present(Shell_name)) Shell_name = 'P11-shell'
   CASE (57)
      if (present(PQN)) PQN = 7 ! Q-shell
      if (present(Shell_name)) Shell_name = 'Q-shell'
   CASE (58)
      if (present(PQN)) PQN = 7 ! Q-shell
      if (present(Shell_name)) Shell_name = 'Q1-shell'
   CASE (59)
      if (present(PQN)) PQN = 7 ! Q-shell
      if (present(Shell_name)) Shell_name = 'Q23-shell'
   CASE (60)
      if (present(PQN)) PQN = 7 ! Q-shell
      if (present(Shell_name)) Shell_name = 'Q2-shell'
   CASE (61)
      if (present(PQN)) PQN = 7 ! Q-shell
      if (present(Shell_name)) Shell_name = 'Q3-shell'
   CASE (62)
      if (present(PQN)) PQN = 7 ! Q-shell
      if (present(Shell_name)) Shell_name = 'Q...-shell'
   CASE (63:)
      if (present(PQN)) PQN = 8 ! R-shell
      if (present(Shell_name)) Shell_name = 'Valence'
   CASE DEFAULT
      if (present(PQN)) PQN = 1 ! K-shell, if something else happened...
      if (present(Shell_name)) Shell_name = 'K-shell'
   END SELECT
end subroutine define_PQN


subroutine Find_element_name(Z, Name, Full_Name, M)
   integer, intent (in) :: Z	! atomic number
   real(8), intent(out) :: M	! mass [proton mass]
   Character(3), intent(out) :: Name ! short name of element
   Character(30), intent(out) :: Full_Name ! full name of element

   select case (Z)
   case(1)
 	Name = 'H'
 	Full_Name = 'Hydrogen' 	
 	M = 1.0082e0
   case(2) 	
	Name = 'He' 	
	Full_Name = 'Helium' 	
	M = 4.002602e0
   case(3)	
	Name = 'Li' 	
	Full_Name = 'Lithium' 	
	M = 6.942e0
   case(4)	
	Name = 'Be' 	
	Full_Name = 'Beryllium' 
	M = 9.012182e0
   case(5)
 	Name = 'B' 	
	Full_Name = 'Boron' 
	M = 10.812e0
   case(6)	
 	Name = 'C' 	
	Full_Name = 'Carbon' 
	M = 12.0112e0
   case(7)
	Name = 'N' 	
	Full_Name = 'Nitrogen'  	
 	M = 14.0072e0
   case(8)
	Name = 'O' 	
	Full_Name = 'Oxygen'  	
 	M = 15.9992e0
   case(9)
	Name = 'F' 	
	Full_Name = 'Fluorine'  	
 	M = 18.9984032e0
   case(10)
	Name = 'Ne' 	
	Full_Name = 'Neon'  	
 	M = 20.1797e0
   case(11)
	Name = 'Na' 	
	Full_Name = 'Sodium'  	
 	M = 22.98976928e0
   case(12)
	Name = 'Mg' 	
	Full_Name = 'Magnesium'  	
 	M = 24.3059e0
   case(13)
	Name = 'Al' 	
	Full_Name = 'Aluminium' 	
 	M = 26.9815386e0
   case(14)
	Name = 'Si' 	
	Full_Name = 'Silicon'  	
 	M = 28.0854e0
   case(15)
	Name = 'P' 	
	Full_Name = 'Phosphorus'  	
 	M = 30.973762e0
   case(16)
	Name = 'S' 	
	Full_Name = 'Sulfur'  	
 	M = 32.062e0
   case(17)
	Name = 'Cl' 	
	Full_Name = 'Chlorine' 	
 	M = 35.452e0
   case(18)
	Name = 'Ar' 	
	Full_Name = 'Argon'
	M = 39.948e0
   case(19)
	Name = 'K' 	
	Full_Name = 'Potassium'  	
 	M = 39.0983e0
   case(20)
	Name = 'Ca' 	
	Full_Name = 'Calcium'  	
 	M = 40.078e0
   case(21)
	Name = 'Sc' 	
	Full_Name = 'Scandium'  	
 	M = 44.955912e0
   case(22)
	Name = 'Ti' 	
	Full_Name = 'Titanium'  	
 	M = 47.867e0
   case(23)
	Name = 'V'	
	Full_Name = 'Vanadium'  	
 	M = 50.9415e0
   case(24)
	Name = 'Cr' 	
	Full_Name = 'Chromium'  	
 	M = 51.9961e0
   case(25)
	Name = 'Mn' 	
	Full_Name = 'Manganese' 	
 	M = 54.938045e0
   case(26)
	Name = 'Fe' 	
	Full_Name = 'Iron'  	
 	M = 55.845e0
   case(27)
	Name = 'Co' 	
	Full_Name = 'Cobalt' 	
 	M = 58.933195e0
   case(28)
	Name = 'Ni' 	
	Full_Name = 'Nickel'  	
 	M = 58.6934e0
   case(29)
	Name = 'Cu' 	
	Full_Name = 'Copper'  	
 	M = 63.546e0 
   case(30)
	Name = 'Zn' 	
	Full_Name = 'Zinc'  	
 	M = 65.38e0
   case(31)
	Name = 'Ga' 	
	Full_Name = 'Gallium'  	
 	M = 69.723e0
   case(32)
	Name = 'Ge' 	
	Full_Name = 'Germanium'  	
 	M = 72.630e0
   case(33)
	Name = 'As' 	
	Full_Name = 'Arsenic' 	
 	M = 74.92160e0
   case(34)
	Name = 'Se' 	
	Full_Name = 'Selenium' 	
 	M = 78.96e0
   case(35)
	Name = 'Br' 	
	Full_Name = 'Bromine'  	
 	M = 79.9049e0
   case(36)
	Name = 'Kr' 	
	Full_Name = 'Krypton'  	
 	M = 83.798e0
   case(37)
	Name = 'Rb' 	
	Full_Name = 'Rubidium'  	
 	M = 85.4678e0
   case(38)
	Name = 'Sr' 	
	Full_Name = 'Strontium'  	
 	M = 87.62e0
   case(39)
	Name = 'Y' 	
	Full_Name = 'Yttrium'  	
 	M = 88.90585e0
   case(40)
	Name = 'Zr' 	
	Full_Name = 'Zirconium'  	
 	M = 91.224e0
   case(41)
	Name = 'Nb' 	
	Full_Name = 'Niobium'  	
 	M = 92.90638e0
   case(42)
	Name = 'Mo' 	
	Full_Name = 'Molybdenum'  	
 	M = 95.96e0
   case(43)
	Name = 'Tc' 	
	Full_Name = 'Technetium' 
 	M = 98e0
   case(44)
	Name = 'Ru' 	
	Full_Name = 'Ruthenium'  	
 	M = 101.07e0
   case(45)
	Name = 'Rh' 	
	Full_Name = 'Rhodium'  	
 	M = 102.90550e0
   case(46)
	Name = 'Pd' 	
	Full_Name = 'Palladium'  	
 	M = 106.42e0
   case(47)
	Name = 'Ag' 	
	Full_Name = 'Silver'  	
 	M = 107.8682e0
   case(48)
	Name = 'Cd' 	
	Full_Name = 'Cadmium'  	
 	M = 112.411e0
   case(49)
	Name = 'In' 	
	Full_Name = 'Indium'  	
 	M = 114.818e0
   case(50)
	Name = 'Sn' 	
	Full_Name = 'Tin'
 	M = 118.710e0
   case(51)
	Name = 'Sb' 	
	Full_Name = 'Antimony'  	
 	M = 121.760e0
   case(52)
	Name = 'Te' 	
	Full_Name = 'Tellurium'  	
 	M = 127.60e0
   case(53)
	Name = 'I'	
	Full_Name = 'Iodine'  	
 	M = 126.90447e0
   case(54)
	Name = 'Xe' 	
	Full_Name = 'Xenon'  	
 	M = 131.293e0
   case(55)
	Name = 'Cs' 	
	Full_Name = 'Caesium'  	
 	M = 132.9054519e0
   case(56)
	Name = 'Ba' 	
	Full_Name = 'Barium'  	
 	M = 137.327e0
   case(57)
	Name = 'La' 	
	Full_Name = 'Lanthanum'  	
 	M = 138.90547e0
   case(58)
	Name = 'Ce' 	
	Full_Name = 'Cerium'  	
 	M = 140.116e0
   case(59)
	Name = 'Pr' 	
	Full_Name = 'Praseodymium'  	
 	M = 140.90765e0
   case(60)
	Name = 'Nd' 	
	Full_Name = 'Neodymium'  	
 	M = 144.242e0
   case(61)
	Name = 'Pm' 	
	Full_Name = 'Promethium' 
 	M = 145e0
   case(62)
	Name = 'Sm' 	
	Full_Name = 'Samarium'  	
 	M = 150.36e0
   case(63)
	Name = 'Eu' 	
	Full_Name = 'Europium' 	
 	M = 151.964e0
   case(64)
	Name = 'Gd' 	
	Full_Name = 'Gadolinium'  	
 	M = 157.25e0
   case(65)
	Name = 'Tb' 	
	Full_Name = 'Terbium'  	
 	M = 158.92535e0
   case(66)
	Name = 'Dy' 	
	Full_Name = 'Dysprosium'  	
 	M = 162.500e0
   case(67)
	Name = 'Ho' 	
	Full_Name = 'Holmium'  	
 	M = 164.93032e0
   case(68)
	Name = 'Er' 	
	Full_Name = 'Erbium'  	
 	M = 167.259e0
   case(69)
	Name = 'Tm' 	
	Full_Name = 'Thulium'  	
 	M = 168.93421e0
   case(70)
	Name = 'Yb' 	
	Full_Name = 'Ytterbium'  	
 	M = 173.054e0
   case(71)
	Name = 'Lu' 	
	Full_Name = 'Lutetium'  	
 	M = 174.9668e0
   case(72)
	Name = 'Hf' 	
	Full_Name = 'Hafnium'  	
 	M = 178.49e0
   case(73)
	Name = 'Ta' 	
	Full_Name = 'Tantalum'  	
 	M = 180.94788e0
   case(74)
	Name = 'W' 	
	Full_Name = 'Tungsten'  	
 	M = 183.84e0 
   case(75)
	Name = 'Re' 	
	Full_Name = 'Rhenium'  	
 	M = 186.207e0 
   case(76)
	Name = 'Os' 	
	Full_Name = 'Osmium' 	
 	M = 190.23e0
   case(77)
	Name = 'Ir' 	
	Full_Name = 'Iridium'  	
 	M = 192.217e0
   case(78)
	Name = 'Pt' 	
	Full_Name = 'Platinum'  	
 	M = 195.084e0
   case(79)
	Name = 'Au' 	
	Full_Name = 'Gold'  	
 	M = 196.966569e0
   case(80)
	Name = 'Hg' 	
	Full_Name = 'Mercury'  	
 	M = 200.592e0 
   case(81)
	Name = 'Tl' 	
	Full_Name = 'Thallium' 	
 	M = 204.389e0
   case(82)
	Name = 'Pb' 	
	Full_Name = 'Lead' 
	M = 207.2e0
   case(83)
	Name = 'Bi' 	
	Full_Name = 'Bismuth'  	
 	M = 208.98040e0
   case(84)
	Name = 'Po' 	
	Full_Name = 'Polonium' 
 	M = 209.0e0
   case(85)
	Name = 'At' 	
	Full_Name = 'Astatine' 
 	M = 210.0e0
   case(86)
	Name = 'Rn' 	
	Full_Name = 'Radon' 	
 	M = 222.0e0
   case(87)
	Name = 'Fr' 	
	Full_Name = 'Francium' 	
 	M = 223.0e0
   case(88)
	Name = 'Ra' 	
	Full_Name = 'Radium' 
 	M = 226.0e0
   case(89)
	Name = 'Ac' 	
	Full_Name = 'Actinium' 
 	M = 227.0e0
   case(90)
	Name = 'Th' 	
	Full_Name = 'Thorium'  	
 	M = 232.03806e0
   case(91)
	Name = 'Pa' 	
	Full_Name = 'Protactinium'
 	M = 231.035880e0
   case(92)
	Name = 'U' 	
	Full_Name = 'Uranium'  	
 	M = 238.02891e0
   case(93)
	Name = 'Np' 	
	Full_Name = 'Neptunium'
 	M = 237.0e0
   case(94)
	Name = 'Pu' 	
	Full_Name = 'Plutonium'
 	M = 244.0e0
   case(95)
	Name = 'Am' 	
	Full_Name = 'Americium' 
	M = 243.0e0
   case(96)
	Name = 'Cm' 	
	Full_Name = 'Curium' 
 	M = 247.0e0
   case(97)
	Name = 'Bk' 	
	Full_Name = 'Berkelium' 
 	M = 247.0e0
   case(98)
	Name = 'Cf' 	
	Full_Name = 'Californium' 
 	M = 251.0e0
   case(99)
	Name = 'Es' 	
	Full_Name = 'Einsteinium' 
 	M = 252.0e0
   case(100)
	Name = 'Fm' 	
	Full_Name = 'Fermium'
 	M = 257.0e0
   case(101)
	Name = 'Md' 	
	Full_Name = 'Mendelevium' 
 	M = 258.0e0
   case(102)
	Name = 'No' 	
	Full_Name = 'Nobelium'
 	M = 259.0e0
   case(103)
	Name = 'Lr' 	
	Full_Name = 'Lawrencium'
 	M = 262.0e0
   case(104)
	Name = 'Rf' 	
	Full_Name = 'Rutherfordium' 
 	M = 267.0e0
   case(105)
	Name = 'Db' 	
	Full_Name = 'Dubnium' 
 	M = 268.0e0
   case(106)
	Name = 'Sg' 	
	Full_Name = 'Seaborgium' 
 	M = 269.0e0
   case(107)
	Name = 'Bh' 	
	Full_Name = 'Bohrium'
 	M = 270.0e0
   case(108)
	Name = 'Hs' 	
	Full_Name = 'Hassium' 
 	M = 269.0e0
   case(109)
	Name = 'Mt' 	
	Full_Name = 'Meitnerium' 
 	M = 278.0e0
   case(110)
	Name = 'Ds' 	
	Full_Name = 'Darmstadtium' 	
 	M = 281.0e0
   case(111)
	Name = 'Rg' 	
	Full_Name = 'Roentgenium'
 	M = 281.0e0
   case(112)
	Name = 'Cn' 	
	Full_Name = 'Copernicium'
 	M = 285.0e0
   case(113)
	Name = 'Uut' 	
	Full_Name = 'Ununtrium'
 	M = 286.0e0
   case(114)
	Name = 'Fl' 	
	Full_Name = 'Flerovium' 
 	M = 289.0e0
   case(115)
	Name = 'Uup' 	
	Full_Name = 'Ununpentium' 
 	M = 288.0e0
   case(116)
	Name = 'Lv' 	
	Full_Name = 'Livermorium'
 	M = 293.0e0 
   case(117)
	Name = 'Uus' 	
	Full_Name = 'Ununseptium'	
 	M = 294.0e0
   case(118)
	Name = 'Uuo' 	
	Full_Name = 'Ununoctium'
 	M = 294.0e0
   case(119:)
	Name = 'UFO'
	Full_Name = 'Unknown'
	print*, 'Mass of this element is not in our database.'
	print*, 'Please, provide the value in the units of the proton mass:'
	read*, M
end select
end subroutine Find_element_name


END MODULE Dealing_with_EADL

!***************************************************************
! This file is part of TREKIS-3
!***************************************************************
! This module contains subroutines to deal with integration of a given function

MODULE Integration
use Objects

implicit none 

! this interface finds by itself which of the two subroutine to use depending on the dimensions of the array passed:
interface Integrate_function ! integrate function
    module procedure Integrate_function_one
    module procedure Integrate_function_save
end interface Integrate_function

interface Trapeziod
    module procedure Trapeziod_one
    module procedure Trapeziod_save
end interface Trapeziod


!private  ! hides items not listed on public statement 
public :: Integrate_function, Trapeziod



contains  

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
        
        select case(Int_type)   ! for different types of integration:
        case default
            call Trapeziod_save(x,f,temp_x0,temp_xn,i_0,i_n,res) ! here is the method of integration
        endselect
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

END MODULE Integration

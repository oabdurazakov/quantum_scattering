! rk4.f90:  4th order Runge-Kutta solution for harmonic oscillator
!
! From: "A SURVEY OF COMPUTATIONAL PHYSICS" 
! by RH Landau, MJ Paez, and CC BORDEIANU 
! Copyright Princeton University Press, Princeton, 2008.
! Electronic Materials copyright: R Landau, Oregon State Univ, 2008;
! MJ Paez, Univ Antioquia, 2008; and CC BORDEIANU, Univ Bucharest, 2008.
! Supported by the US National Science Foundation
!
!
! code cleaned up/modernized/modularized by RTC 09/2009
!			
module rk4_mod
  use deriv_mod
  implicit none

  private 
  public :: rk4

contains

! rk4: 4th order Runge-Kutta method
!
! x : independent variable. The subroutine does not modify this.
! xstep : stepsize. 
! y : dependent variables
! n : number of equations
!
  subroutine rk4(x, xstep, y,E,n)
    integer,intent(in) :: n
    real(kind=8),intent(in) :: x,xstep,E
    complex(kind=8),intent(inout),dimension(n) :: y

    complex(kind=8),dimension(n) :: k1, k2, k3, k4, t1, t2, t3
    real(kind=8) :: h
    integer :: i
    
    h=xstep/2.0_8
    do i = 1,n
       k1(i) = xstep * deriv(x, y, E,i)
       t1(i) = y(i) + 0.5_8*k1(i)
    enddo
    do i = 1,n
       k2(i) = xstep * deriv(x+h, t1,E,i)
       t2(i) = y(i) + 0.5_8*k2(i)
    enddo
    do i = 1,n
       k3(i) = xstep * deriv(x+h, t2,E,i)
       t3(i) = y(i) + k3(i)
    enddo
    do i = 1,n
       k4(i) = xstep * deriv(x+xstep, t3,E,i)
       y(i) = y(i) + (k1(i) + (2.0_8*(k2(i) + k3(i))) + k4(i))/6.0_8
       enddo
    
    return
  end subroutine rk4

end module rk4_mod

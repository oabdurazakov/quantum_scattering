! 
! derivatives for simple harmonic motion
!hbar=m=V0=0
module deriv_mod
  implicit none

  private
  public:: deriv

contains

  function deriv(x, y, E,i)
    real(kind=8)::E,V,V0
    complex(kind=8)::deriv
    real(kind=8),intent(in) :: x
    complex(kind=8),intent(in) :: y(:)
    integer,intent(in) :: i
    V0=0.3d0
    select case(i)
    case (1)
       deriv=y(2)
    case (2)
      if(x>=0.0d0.and.x<=5.0d0) then
         V=V0
      else if(x>=10.0d0.and.x<=15.0d0)then
         V=V0
         else
            V=0.0d0
         end if
       deriv=-26.25d0*0.067*(E-V)*y(1) ! 2*m/h^2=26.25 eV^-1*nm^-2
    end select

    return
  end function deriv

end module deriv_mod

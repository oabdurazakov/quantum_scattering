!
!1D quantum scattering through the double square barrier
!hbar=m=V0=1.0
program test
  use deriv_mod
  use rk4_mod
  implicit none

  integer,parameter :: n=2
  complex,parameter :: icmplx=dcmplx(0.0d0,1.0d0)
  complex(kind=8)::a,b
  complex(kind=8),dimension(n) :: y
  real(kind=8) :: dx,x,E,k,dE,r,t,exact,V0,a
  integer :: i,j
  
  !k=sqrt(2.0d0*E) !m=1,hbar=1


  V0=2.0d0
  a=2.0
  E=0.02d0 !E=0ev-0.5ev
  dE=0.01
  open(10,file='exact.dat',status='replace')

  do i=1,998    
     y(1)=dcmplx(1.0d0,0.0d0) !initial cond
     y(2)=dcmplx(0.0d0,k) !initial condition
 

  x=a !width of the barrier
  dx=-0.001d0 !step 

   
   do j=1,2000!integrate from right to left
      k=sqrt(2.0d0*E) !m=1,hbar=1
      call rk4(x,dx,y,E,n)
    
     ! if(i==150)then
      !   write(10,*)x,real(y(1)),E
   !end if
   x=x+dx
   end do
   x=0.0d0
   a=0.5*(y(1)+y(2)/(icmplx*k))
   b=0.5*(y(1)-y(2)/(icmplx*k))
   r=(abs(b)**2)/(abs(a)**2)!reflection
   t=1.0d0/(abs(a)**2)!transmission
   if(E<V0)then
      exact=1.0d0/(1.0d0+(V0**2)/(4*E*(V0-E))*sinh(2*a*sqrt(2*(V0-E)))**2)
  ! else if(E==V0)then
   !   exact=1.0d0/(1+2*E*a**2)
   else if(E>V0)then
      exact=1.0d0/(1+V0**2/(4*E*(E-V0))*sin(2*a*sqrt(2*(E-V0)))**2)
   end if
  write(10,*)E,t,exact
   E=E+dE
end do
  close(10)

end program test

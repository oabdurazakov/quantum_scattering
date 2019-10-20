!
!1D quantum scattering through the double square barrier

program test
  use deriv_mod
  use rk4_mod
  implicit none

  integer,parameter :: n=2
  complex,parameter :: icmplx=dcmplx(0.0d0,1.0d0)
  complex(kind=8)::a1,a2
  complex(kind=8),dimension(n) :: y
  real(kind=8) :: dx,x,E,k,dE,r,t
  integer :: i,j
  
  !k=sqrt(2.0d0*E) !m=1,hbar=1

 
  E=0.02d0 !E=0ev-0.5ev
  dE=0.01
  open(10,file='small.dat',status='replace')

  do i=1,998    
     y(1)=dcmplx(1.0d0,0.0d0) !initial cond
     y(2)=dcmplx(0.0d0,k) !initial condition
 

  x=15.0d0 !width of the barrier
  dx=-0.001d0 !step
 

 
   do j=1,15000 !integrate from right to left
      k=sqrt(26.250d0*0.067*E) !m=0.067m_e:effective mass of an electron, 
      call rk4(x,dx,y,E,n)
    
      if(i==15)then
         write(10,*)x,real(y(1))
   end if
   x=x+dx
   end do
   x=0.0d0
   a1=0.5*(y(1)+y(2)/(icmplx*k))
   a2=0.5*(y(1)-y(2)/(icmplx*k))
   r=(abs(a2)**2)/(abs(a1)**2)!reflection
   t=1.0d0/(abs(a1)**2)!transmission
   !write(10,*)E,log(t)
   E=E+dE
end do
  close(10)

end program test

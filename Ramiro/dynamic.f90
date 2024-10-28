!> The program solves classical Newton's Eq. for e -Ne+He (PES from Ramiro) \\
!> Non-adiabatic transitions take place according to Landau-Zener Surface Hopping algorithm (see J. Chem. Phys. 142, 104307 (2015)) \\
!> The PES must be given as follows R, x1, energies
program dynamic_icec
implicit none

integer, parameter :: ndof = 2
double precision, dimension(ndof) :: mass, q, p, v

double precision :: x1, r

double precision :: ti, tf
double precision :: dt, dr, dx1

!> allows to rescale the velocities after hopping (0: no rescale, 1: rescale)
integer :: rescal

!> fsta defines the current state 
integer :: fsta

!> fpot is the file providing the PES, finit provides masses, initial pos. and vel.
character(64) :: fpot, finit, arg

integer :: i, j, l

!> all inputs read from namelist input
namelist /input/ dt, dr, dx1, fpot, finit, tf, fsta, rescal

if(iargc()==0) then
   write(*,*) "Please provide the input file in command line"
   stop
else
 do i = 1, iargc()
   call getarg(i, arg)
!   write(*,*) arg
 enddo
endif

 write(*,'(A)')"Here is the input"
 open(unit=10,file=arg)
   read(10,nml=input)
   write(*,nml=input)
 close(10)
 write(*,*)
 tf = tf * 41d0

! reads the initial conditions
write(*,'(A)')"Reads the initial conditions"
write(*,'(A)')"!!!  Enter R first !!!"
 open(unit=10,file=finit)
 do i = 1, ndof
   read(10,*)mass(i),q(i),p(i)
   v(i) = p(i) / mass(i)
   write(100,'(6(f20.15,1X))')mass(i),q(i),p(i)
   mass(i)=mass(i)*1836.15d0
 enddo
 close(10)
 write(100,*)
 write(*,*)

! starts with the dynamics

 ti = 0d0
 call dyn(ndof,mass,q,v,ti,tf,fpot,fsta,dt,dr,dx1,rescal)

  write(100,*)fsta
  write(100,'(6(f20.10,1X))')(q(i),v(i),i=1, ndof)
  write(100,*)

end program dynamic_icec

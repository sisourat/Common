subroutine dyn(ndof,mass,q,v,ti,tf,fpot,fsta,dt,dr,dx1,rescal)
use bspline_module
use interpolation
use RDistributions
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Nico 29.10.2024 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! performs the dynamics until tf                        !!
! uses Verlet algorithm                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! for the system
integer, intent(in) :: ndof
double precision, dimension(ndof), intent(in) :: mass
double precision, dimension(ndof), intent(inout) :: q, v

! for the dynamics
double precision, intent(in) :: dt, dr, dx1 ! in atomic units
double precision, intent(inout) :: ti, tf
double precision, dimension(ndof) :: qm, qt, qnew, qdr, ddr
double precision :: time, grade
double precision, dimension(:), allocatable :: valdrp, valdrm, val

!for the electronic state
character(64) :: fpot

!for Landau-Zener 
double precision, dimension(:), allocatable :: epair, epair_t1, d_epair, d_epair_t1, d2_epair
double precision :: plz

double precision :: egap, scal, vscal, vunscal, ekin1, ekin2
integer, intent(in) :: rescal ! rescal if equals to 1

!for the interpolation
integer :: nsta
integer :: nr, nx1
double precision, dimension(:), allocatable :: r, x1
double precision, dimension(:,:,:), allocatable :: energy

integer, intent(inout) :: fsta
double precision, dimension(:,:), allocatable :: tr, tx1

integer,parameter :: kr = 4     !! order in r
integer,parameter :: kx1 = 4     !! order in x1
integer,parameter :: iknot = 0  !! automatically select the knots

double precision :: tol
logical :: fail, file_e

integer :: inbvx,inbvy,inbvz
integer :: iloy,iloz,iflag
integer :: idr, idx1

double precision :: newr, newx1

integer :: i, j, k, ista

integer :: nhop = 0

! Reads the data

inquire( file=trim(fpot), exist=file_e )
if ( file_e .eqv. .false. ) then
 write(*,*) trim(fpot), " does not exist"
 stop
endif

write(*,'(A,A)')'Read PESs in ',trim(fpot)
open(unit=10,file=trim(fpot))
 read(10,*)nr, nx1, nsta
 write(*,'(A,I4,A,I4,A,I4,A)')'There are',nr,' * ',nx1,' points'
 write(*,'(A,I4,A)')'There are', nsta, 'states'

allocate(r(nr),x1(nx1))
allocate(energy(nsta,nr,nx1))
!
do i = 1, nr
  do j = 1, nx1
    read(10,*)r(i),x1(j),(energy(ista,i,j),ista=1,nsta)
!    write(*,*)r(i),x1(j),(energy(ista,i,j),ista=1,nsta)
  enddo
enddo 
close(10)
write(*,'(A)')'Reading PESs done '
write(*,*)
write(*,'(A)')'Lets start the dynamics'
!!energy(:,:,:)=0d0

allocate(epair(nsta),epair_t1(nsta),d_epair(nsta),d_epair_t1(nsta),d2_epair(nsta))
epair(:) = 0d0
epair_t1(:) = 0d0
d_epair(:) = 0d0
d_epair_t1(:) = 0d0
d2_epair(:) = 0d0

!
! Set-up the interpolation
inbvx = 1
inbvy = 1
inbvz = 1
iloy  = 1
iloz  = 1

idr = 0
idx1 = 0

fail = .false.
tol = 1.0e-14

ddr(1) = dr
ddr(2) = dx1

allocate(tr(nsta,nr+kr),tx1(nsta,nx1+kx1))
allocate(val(nsta),valdrp(nsta),valdrm(nsta))

do ista = 1, nsta
 iflag=0
 call db2ink(r,nr,x1,nx1,energy(ista,:,:),kr,kx1,iknot,tr(ista,:),tx1(ista,:),energy(ista,:,:),iflag)

 if(iflag/=0) then
  write(*,*)"error in db2ink", iflag
  return
 endif
enddo

! landau-zenner surf. hopp. stuff
newr=q(1)
newx1=q(2)
do ista=1,nsta
 call db2val(newr,newx1,idr,idx1,tr(ista,:),tx1(ista,:),nr,nx1,kr,kx1,energy(ista,:,:),val(ista),iflag,inbvx,inbvy,iloy)
 if(iflag/=0) then
  write(*,*)"error in db2val at time (1)",time,newr,newx1
  stop
 endif
enddo

do ista=1,nsta
  epair(ista) = abs(val(ista)-val(fsta))
enddo

! first step of the Verlet algorithm
do i = 1, ndof
  qm(i) = q(i)
  qt(i) = q(i) + v(i)*dt
enddo

time=ti
call init_random_seed()
do while(time<tf)

! landau-zenner surf. hopp. stuff
 epair_t1(:) = epair(:)
 d_epair_t1(:) = d_epair(:)
 epair(:)=0d0
 
! computes the energy at position at time t
 newr=qt(1)
 newx1=qt(2)
 do ista=1,nsta
 call db2val(newr,newx1,idr,idx1,tr(ista,:),tx1(ista,:),nr,nx1,kr,kx1,energy(ista,:,:),val(ista),iflag,inbvx,inbvy,iloy)
 if(iflag/=0) then
  write(*,*)"error in db2val at time (2)",time, newr, newx1, ista
  stop
 endif
 enddo
 write(200,'(4(f20.10,1X),i3,10(f20.10,1X))')time*0.024,newr,newx1,val(fsta),fsta,v(2)
 !write(*,*)time*0.024,qt(2),val(fsta),(qt(2)-qm(2))/dt,fsta
 if(qt(2)>100.0 .or. qt(2)<-110.0) then
     exit
 endif

! apply LZ surface hopping here
do ista=1,nsta
  epair(ista) = abs(val(ista)-val(fsta))
enddo

 d_epair(:) = (epair(:)-epair_t1(:))/dt
 d2_epair(:) = (d_epair(:)-d_epair_t1(:))/dt

 do ista=1,nsta 
  if(d_epair(ista)*d_epair_t1(ista) < 0d0 .and. d2_epair(ista)>0d0) then
      plz = exp(-0.5d0*pi*sqrt(epair(ista)**3/d2_epair(ista)))

! is there enough kinetic energy to fill the energy gap, if not => frustrated hop
      ekin1 = 0d0
        do i = 1, ndof
             vunscal = (qt(i)-qm(i))/dt
             ekin1 = ekin1 + 0.5d0*mass(i)*vunscal**2
        enddo
    
      if(plz>rand_uniform(0d0,1d0) .and. epair(ista)<ekin1) then
        write(*,*)"HOP", ista
        if(rescal==1) then
          egap = val(fsta)-val(ista)
          write(*,*)"Egap,Ekin",egap,ekin1
          ekin1 = 0d0
          ekin2 = 0d0
          scal = egap ! only the electron velocity is rescaled
          do i = 2, ndof
             vunscal = (qt(i)-qm(i))/dt
             write(*,*)vunscal
             ekin1 = ekin1 + 0.5d0*mass(i)*vunscal**2
             if(1d0+2d0*scal/(mass(i)*vunscal**2)>0d0) then
                vscal = vunscal*sqrt(1d0+2d0*scal/(mass(i)*vunscal**2))
             write(*,*)vscal,'a'
                qt(i) = qt(i) + (vscal-vunscal)*dt
             else !! in this case the evenly distributed correction is larger than the veloc. in this coord=> we put the velocity to zero (energy is not strictly conserved)
                qt(i) = qm(i)
             endif
             vunscal = (qt(i)-qm(i))/dt
             ekin2 = ekin2 + 0.5d0*mass(i)*vunscal**2
          enddo
          write(*,*)"Ekin_unscal,Ekinscal,Energy conservation",ekin1, ekin2, ekin2-egap-ekin1
        endif
       fsta=ista
       epair_t1(:) = 0d0
       epair(:) = 0d0
       d_epair(:) = 0d0
       d_epair_t1(:) = 0d0
      exit ! only one hop allowed
      endif
  endif
 enddo

!! end surface hopping

! computes the new position
 
 do i = 1, ndof
   qdr(:) = qt(:) 
   qdr(i) = qt(i) + ddr(i)
   newr=qdr(1)
   newx1=qdr(2)
   do ista=1,nsta
    call db2val(newr,newx1,idr,idx1,tr(ista,:),tx1(ista,:),nr,nx1,kr,kx1,energy(ista,:,:),valdrp(ista),iflag,inbvx,inbvy,iloy)
   enddo

   qdr(:) = qt(:) 
   qdr(i) = qt(i) - ddr(i)
   newr=qdr(1)
   newx1=qdr(2)
   do ista=1,nsta
    call db2val(newr,newx1,idr,idx1,tr(ista,:),tx1(ista,:),nr,nx1,kr,kx1,energy(ista,:,:),valdrm(ista),iflag,inbvx,inbvy,iloy)
   enddo

    grade = 0.5*( (valdrp(fsta)-val(fsta))/ddr(i) - (valdrm(fsta)-val(fsta))/ddr(i) )
    !!if(i==1) grade = 0d0 ! fixed nuclei approx. 
    qnew(i) = 2*qt(i) - qm(i) - grade*dt**2/mass(i) 

 enddo
!! write(*,*)

 qm(:) = qt(:)
 qt(:) = qnew(:)
 time = time + dt

enddo
write(*,'(A)')'Dynamics done'
write(*,*)
!stop

q(:) = qnew(:)
v(:) = (qt(:)-qm(:))/dt

deallocate(r,x1,energy)
deallocate(tr,tx1)
deallocate(val,valdrp,valdrm)
deallocate(epair,epair_t1,d_epair,d2_epair)

end subroutine dyn

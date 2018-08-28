SUBROUTINE Initialize()
  use CONSTANTS
  use WAVE_FUNC
  use TD_CALC
  implicit none

  integer :: seedsize, is, ik, ix,ib
  integer, allocatable :: seed(:)
  real(8) :: tt, s
  real(8),dimension(:) :: fst(0:Nx-1)

  !Settings of randomnumber
  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  do is = 1, seedsize
    call system_clock(count=seed(is))
  end do
  call random_seed(put=seed(:))

  !Initial condition of wave function
  do ik = -LNk, RNk, 1
    do ib = 1, Nb, 1
      call random_number(fst)
      !do ix = Nx/2+1, Nx-1, 1
      !  fst(ix)=fst(-ix+Nx)
      !end do
      s=1.d0/sqrt(sum(abs(fst(:))**2)*dx)
      u(:,ik,ib)=fst(:)*s
    end do
  end do

  !Settings of k-points

  do ik = -LNk, RNk, 1
    k(ik)=dk*real(ik)
  end do

  !dk_1=dk
  !dk_2=dk*LNk_1/LNk_2

  !do ik = -LNk, RNk, 1
  !  if(ik <= LNk_1 .AND. ik >= -LNk_1) then
  !    k(ik)=dk_1*real(ik)
  !  else if ( (ik > LNk_1) .AND. (ik < LNk_1+LNk_2*(LNk_part-2)) ) then
  !    k(ik)=dk_2*real(ik-LNk_1)+dk_1*real(LNk_1)
  !  else if ( (ik < -LNk_1) .AND. (ik > -LNk_1-LNk_2*(LNk_part-2))) then
  !    k(ik)=dk_2*real(ik+LNk_1)-dk_1*real(LNk_1)
  !  else if ( (ik >= LNk_1+LNk_2*(LNk_part-2) )) then
  !    k(ik)=dk_1*real(ik-LNk_1-LNk_2*(LNk_part-2))+dk_1*real(LNk_1)+dk_2*real(LNk_2*(LNk_part-2))
  !  else if ( (ik <= -LNk_1-LNk_2*(LNk_part-2) )) then
  !    k(ik)=dk_1*real(ik+LNk_1+LNk_2*(LNk_part-2))-dk_1*real(LNk_1)-dk_2*real(LNk_2*(LNk_part-2))
  !  end if
  !end do


  ! Create the mask function for electric-field and fourie transformation
  do it = 0, Nt, 1
    tt=real(it)*dt
    mask(it)=sin(tt/tau)**4



    !Electric field
    !Et(it) = E0*cos(omega*tt)*mask(it)
    !At(it) = E0*sin(omega*tt)*mask(it)/omega
    !Et = E0*mask(it)
    Et(it)= E0*sin(omega*(tt-tau*pi/2.d0))*sin(tt/tau)**4

  end do
  do it = 0, Nt, 1
    if(it==0 .OR. it==Nt)then
      At(it)=0
      Et(it)=0
    else
      At(it)=-(sum(Et(0:it-1))+Et(it)*0.5d0)*dt
      !Et(it)=-(At(it+1)-At(it-1))/(2.d0*dt)
    end if
  end do


  !Create the potential
  do ix=0,Nx-1
    !pot(ix) = V0 * (exp(-((ix*dx-width/2.d0)/a0)**2) + exp(-((ix*dx-width*3.d0/2.d0)/a0)**2) + exp(-((ix*dx+width/2.d0)/a0)**2)  )
    pot(ix) = V0*(1+cos(2*pi*ix*dx/width))
  end do
!  pot(:)=0.d0

  allocate(modx(0:Nx*2+Nd-1))

  do ix=0,Nx*2+Nd-1
    modx(ix) = mod(ix,Nx)
  end do

  call fd_coef()

END SUBROUTINE Initialize
SUBROUTINE td_init()
  use CONSTANTS
  use WAVE_FUNC
  implicit none
  integer :: ib, ik
  do ib = 1, Nbt, 1
    do ik = -LNk, RNk, 1
      zu(:,ik,ib)=u(:,ik,ib)
    end do
  end do
END SUBROUTINE
Subroutine fd_coef
  use CONSTANTS
  implicit none
  integer :: i,j,k,n
  real(8) :: t,s,s0

  lap=0.d0
  nab=0.d0

  n=Nd
  s=0.d0
  do i=-n,n
    if ( i==0 ) cycle
    do j=-n,n
      if ( j==i .or. j==0 ) cycle
      s=s+1.d0/dble(i*j)
    end do
  end do
  lap(0)=s

  do j=1,n
    t=1.d0
    do i=-n,n
      if ( i==j ) cycle
      t=t*(j-i)
    end do
    s=0.d0
    do k=-n,n
      if ( k==j .or. k==0 ) cycle
      s0=1.d0
      do i=-n,n
        if ( i==k .or. i==j .or. i==0 ) cycle
        s0=s0*(-i)
      end do
      s=s+s0
    end do
    lap( j)=2.d0*s/t
    lap(-j)=lap(j)
  end do

  do j=1,n
    t=1.d0
    do i=-n,n
      if ( i==j ) cycle
      t=t*(j-i)
    end do
    s=1.d0
    do i=-n,n
      if ( i==j .or. i==0 ) cycle
      s=s*(-i)
    end do
    nab( j)=s/t
    nab(-j)=-nab(j)
  end do

  lap(-Nd:Nd)=lap(-Nd:Nd)/dx**2
  nab(-Nd:Nd)=nab(-Nd:Nd)/dx

  return
End Subroutine fd_coef
Subroutine fd_coef_k
  use CONSTANTS
  use FD_K
  implicit none
  integer :: i,j,k,n
  real(8) :: t,s,s0

  nab_k=0.d0

  n=Nd_k

  do j=1,n
    t=1.d0
    do i=-n,n
      if ( i==j ) cycle
      t=t*(j-i)
    end do
    s=1.d0
    do i=-n,n
      if ( i==j .or. i==0 ) cycle
      s=s*(-i)
    end do
    nab_k( j)=s/t
    nab_k(-j)=-nab_k(j)
  end do

  nab_k(-Nd_k:Nd_k)=nab_k(-Nd_k:Nd_k)/dk

  return
End Subroutine fd_coef_k

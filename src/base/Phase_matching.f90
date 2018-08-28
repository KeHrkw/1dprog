SUBROUTINE Phase_matching()
  use CONSTANTS
  use WAVE_FUNC
  integer :: ib, ik, ix, iphase, ipm
  real(8) :: uu
  complex(kind(0d0)) :: uix
  write(*,*) "phase matching"

  do ib = 1, Nb, 1
    uix=sum(u(0,:,ib))
    iphase=0
    do ix = 1, Nx-1, 1
      if ( abs(uix) < abs(sum(u(ix,:,ib))) ) then
        uix=sum(u(ix,:,ib))
        iphase=ix
      end if
    end do
    do ik = -LNk, RNk, 1
      uix=u(iphase,ik,ib)
      !write(8,'(<3>i,<2>e)') ib,iphase,ik,uix
      u(:,ik,ib)=u(:,ik,ib)/uix
      uu=sum(abs(u(:,ik,ib))**2)*dx
      u(:,ik,ib)=u(:,ik,ib)/sqrt(uu)
    end do
    do ix = 0, Nx-1, 1
      do ipm = 0, 1, 1
        if(ipm==0)then
          phase(ix,ib,ipm)=1.d0*exp(zI*2.d0*pi*real(ix-iphase)/real(Nx))
        else if(ipm==1)then
          phase(ix,ib,ipm)=1.d0*exp(-zI*2.d0*pi*real(ix-iphase)/real(Nx))
        end if
      end do
    end do
  end do


END SUBROUTINE

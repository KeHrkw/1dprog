SUBROUTINE base_plot()
  use CONSTANTS
  use WAVE_FUNC
  use COEFF
  implicit none
  integer :: ik, ix
  integer :: jb, ib
  real(8) :: uu
  complex(kind(0d0)) :: uix

  open(11,file='./output/base.dat')
  write(*,*) "Each Energy"
  write(*,*) "Energy < 0"
  do ib = 1, Nb, 1
    do ik = -LNk, RNk, 1
      k_in=k(ik)
      u_in(0:Nx-1)=u(0:Nx-1,ik,ib)
      call h0_operation()
      eps(ik,ib)=real(sum(conjg(u(:,ik,ib))*hu(:))*dx)
      if( eps(ik,ib)<0.d0 ) write(*,*) 'VB',ib,ik,eps(ik,ib)
      write(11,'(<2>i,<2>e)') ib,ik,k(ik),eps(ik,ib)
    end do
    write(11,*)
    write(11,*)
    hav_base(ib)=sum(eps(:,ib))/real(Nk)
    !hav_base(ib)=(sum(eps(-LNk_1:LNk_1,ib))+sum(eps(LNk_1+1:LNk_1+LNk_2*(LNk_part-2),ib))*real(LNk_1/LNk_2)+sum(eps(LNk_1+LNk_2*(LNk_part-2)+1:RNk,ib)) &
    !          & +sum(eps(-LNk_1-LNk_2*(LNk_part-2):-LNk_1-1,ib))*real(LNk_1/LNk_2)+sum(eps(-LNk:-LNk_1-LNk_2*(LNk_part-2)-1,ib)) )/real(LNk_1*LNk_part*2+1)
  end do
  close(11)
  do ib = 1, Nb, 1
    write(*,'(<1>i,<2>e)') ib, hav_base(ib),sum(hav_base(1:ib))
  end do
  open(11,file='./output/bphi.dat')
  do ix = 0, Nx-1, 1
    write(11,'(<2>i,<1>e)') 0, ix, pot(ix)
  end do
  write(11,*) ""
  write(11,*) ""
  do ib = 1, Nb, 1
    do ix = 0, Nx-1, 1
      write(11,'(<2>i,<1>e)') ib, ix, sum(abs(u(ix,:,ib))**2)/real(Nk)
    end do
    write(11,*) ""
    write(11,*) ""
  end do
  close(11)
  open(11,file='./output/bphi_k.dat')
  do ib = 1, Nb, 1
    do ik = -LNk, RNk, 1
      do ix = 0, Nx-1, 1
        write(11,'(<3>i,<2>e)') ib,ik,ix,u(ix,ik,ib)
      end do
      write(11,*) ""
      write(11,*) ""
    end do
  end do
  close(11)

  open(11,file='./output/sp_bphi.dat')
  do ib = 1, Nb, 1
    do ix = 0, Nx-1, 1
      do ik = -LNk-Nk, -LNk-1, 1
        uix=phase(ix,ib,0) * u(ix,ik+Nk,ib)
        write(11,'(<2>i,<10>e)') ix,ik,aimag(u(ix,ik+Nk,ib)),real(u(ix,ik+Nk,ib)),aimag(uix),real(uix),atan2(aimag(uix),real(uix)),aimag(log(uix))
      end do
      do ik = -LNk, RNk, 1
        write(11,'(<2>i,<10>e)') ix,ik,aimag(u(ix,ik,ib)),real(u(ix,ik,ib)),aimag(u(ix,ik,ib)),real(u(ix,ik,ib)),atan2(aimag(u(ix,ik,ib)),real(u(ix,ik,ib))),aimag(log(uix))
      end do
      do ik = RNk+1, RNk+Nk, 1
        uix=phase(ix,ib,1) * u(ix,ik-Nk,ib)
        write(11,'(<2>i,<10>e)') ix,ik,aimag(u(ix,ik-Nk,ib)),real(u(ix,ik-Nk,ib)),aimag(uix),real(uix),atan2(aimag(uix),real(uix)),aimag(log(uix))
      end do
      write(11,*)""
    end do
    write(11,*)
    write(11,*)
  end do
  close(11)
END SUBROUTINE

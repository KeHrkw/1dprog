SUBROUTINE base_plot()
  use CONSTANTS
  use WAVE_FUNC
  use COEFF
  implicit none
  integer :: ik, ix, ik_min, ik_max
  integer :: jb, ib
  real(8) :: uu, en_min, en_max
  complex(kind(0d0)) :: uix

  do ib=1,Nb
    en_min=eps(1,ib)
    ik_min=1
    en_max=eps(1,ib)
    ik_max=1
    do ik=2,Nk
      if(eps(ik,ib)>en_max)then
         en_max=eps(ik,ib)
         ik_max=ik
      else if(eps(ik,ib)<en_min)then
         en_min=eps(ik,ib)
         ik_min=ik
      end if
    end do
    if(en_max<0.d0) write(*,*) "Valence Band       ib =",ib
    if(en_min>0.d0) write(*,*) "Conduction Band    ib =",ib
      write(*,*) "  ik_max =",ik_max,"    en_max =",en_max
      write(*,*) "  ik_min =",ik_min,"    en_min =",en_min
  end do

  open(11,file='./output/base.dat')
  do ib = 1, Nb, 1
    do ik = 1, Nk, 1
      write(11,'(<2>i,<2>e)') ib,ik,k(ik),eps(ik,ib)
    end do
    write(11,*)
    write(11,*)
  end do
  close(11)
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
    do ik = 1, Nk, 1
      do ix = 0, Nx-1, 1
        write(11,'(<3>i,<3>e)') ib,ik,ix,k(ik),u(ix,ik,ib)
      end do
      write(11,*) ""
      write(11,*) ""
    end do
  end do
  close(11)

  open(11,file='./output/sp_bphi.dat')
  do ib = 1, Nb, 1
    do ix = 0, Nx-1, 1
      do ik = -Nk+1, 0, 1
        uix=phase(ix,ib,0) * u(ix,ik+Nk,ib)
        write(11,'(<2>i,<10>e)') ix,ik,aimag(u(ix,ik+Nk,ib)),real(u(ix,ik+Nk,ib)),aimag(uix),real(uix),atan2(aimag(uix),real(uix)),aimag(log(uix))
      end do
      do ik = 1, Nk, 1
        write(11,'(<2>i,<10>e)') ix,ik,aimag(u(ix,ik,ib)),real(u(ix,ik,ib)),aimag(u(ix,ik,ib)),real(u(ix,ik,ib)),atan2(aimag(u(ix,ik,ib)),real(u(ix,ik,ib))),aimag(log(uix))
      end do
      do ik = Nk+1, Nk+Nk, 1
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

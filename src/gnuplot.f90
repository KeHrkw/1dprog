SUBROUTINE gnuplot()
  use FD_K, only:isplot
  implicit none
open(10, file = 'ft_blo.plt', status = 'replace')
write(10,*) 'set terminal postscript eps color enhanced "Osaka" 25'
write(10,*) 'set output "./bphi.eps'
write(10,*) 'set ytics nomirror # 第一y軸の目盛は左側のみにする'
write(10,*) 'set y2tics         # 第二y軸を描画することを指定'
write(10,*) 'set grid           # 目盛線を描画することを指定'
write(10,*) 'set xlabel "ix"  # x軸の見出しを指定'
write(10,*) 'set ylabel "bphi"           # 第一y軸の見出しを指定'
write(10,*) 'set y2label "Potential" # 第二y軸の見出しを指定'
write(10,*) 'plot "./bphi.data" index 1 u 2:3 w l lw 4  title"n=0", \'
write(10,*) '     "./bphi.data" index 2 u 2:3 w l lw 4  title"n=1", \'
write(10,*) '     "./bphi.data" index 3 u 2:3 w l lw 4  title"n=2", \'
write(10,*) '      "./bphi.data" index 0 u 2:3 w l axes x1y2 lw 4  title"Pot"'
write(10,*) 'reset'
write(10,*) 'set output "./base.eps"'
write(10,*) 'set ylabel "Energy"'
write(10,*) 'set xlabel "k"'
write(10,*) 'plot "./base.data" u 3:4 w l lw 4  notitle'
write(10,*) 'reset'
write(10,*) 'set output "./td_norm.eps"'
write(10,*) 'set ytics nomirror'
write(10,*) 'set xlabel "time"'
write(10,*) 'set ylabel "norm"'
write(10,*) 'set grid'
write(10,*) 'plot "./td.data" u 2:3 w l lw 4  title "norm"'
write(10,*) 'reset'
write(10,*) 'set output "./td_energy.eps"'
write(10,*) 'set ytics nomirror'
write(10,*) 'set y2tics         # 第二y軸を描画することを指定'
write(10,*) 'set xlabel "time"'
write(10,*) 'set ylabel "Energy at n=0"'
write(10,*) 'set y2label "Energy at n=Nbt"'
write(10,*) 'set grid'
write(10,*) 'plot "./td.data" u 2:7 w l lw 4  title "n=0",\'
write(10,*) '     "./td.data" u 2:9 w l lw 4  title "n=Nbt" axes x1y2'
write(10,*) 'reset'
write(10,*) 'set output "./td_velocity.eps"'
write(10,*) 'set ytics nomirror'
write(10,*) 'set y2tics         # 第二y軸を描画することを指定'
write(10,*) 'set xlabel "time"'
write(10,*) 'set ylabel "Velocity at n=0"'
write(10,*) 'set y2label "Velocity at n=Nbt"'
write(10,*) 'set grid'
write(10,*) 'plot "./td.data" u 2:4 w l lw 4  title "n=0",\'
write(10,*) '     "./td.data" u 2:6 w l lw 4  title "n=Nbt" axes x1y2'
write(10,*) 'reset'
write(10,*) 'set output "./td_Elewave.eps"'
write(10,*) 'set ytics nomirror'
write(10,*) 'set y2tics         # 第二y軸を描画することを指定'
write(10,*) 'set xlabel "time"'
write(10,*) 'set ylabel "Et"'
write(10,*) 'set y2label "At"'
write(10,*) 'set grid'
write(10,*) 'plot "./td.data" u 2:10 w l lw 4  title "Et",\'
write(10,*) '     "./td.data" u 2:11 w l lw 4  title "At" axes x1y2'
write(10,*) 'reset'
write(10,*) 'set terminal postscript eps enhanced color'
write(10,*) 'set logscale y'
write(10,*) 'set format y "10^{%L}"'
write(10,*) 'set xrange[0:50]'
write(10,*) 'set mxtics 5'
write(10,*) 'set grid'
write(10,*) 'set grid mxtics'
write(10,*) 'set xlabel "harmonic energy"'
write(10,*) 'set ylabel "Intensity"'
write(10,*) 'set output "./ft0.eps"'
write(10,*) 'plot "./fts.data" u 3:4 lw 4  w l'
write(10,*) 'reset'
if( isplot==1 )then
  write(10,*) 'set output "./splt_bphi_ib1_Im.eps"'
  write(10,*) 'set pm3d map'
  write(10,*) 'set xlabel "x"'
  write(10,*) 'set ylabel "k-point"'
  write(10,*) 'splot "./sp_bphi.data" index 0 u 1:2:5 notitle'
  write(10,*) 'set output "./splt_bphi_ib1_Re.eps"'
  write(10,*) 'splot "./sp_bphi.data" index 0 u 1:2:6 notitle'
  write(10,*) 'set output "./splt_bphi_ib2_Im.eps"'
  write(10,*) 'splot "./sp_bphi.data" index 1 u 1:2:5 notitle'
  write(10,*) 'set output "./splt_bphi_ib2_Re.eps"'
  write(10,*) 'splot "./sp_bphi.data" index 1 u 1:2:6 notitle'
  write(10,*) 'set output "./splt_bphi_ib3_Im.eps"'
  write(10,*) 'splot "./sp_bphi.data" index 2 u 1:2:5 notitle'
  write(10,*) 'set output "./splt_bphi_ib3_Re.eps"'
  write(10,*) 'splot "./sp_bphi.data" index 2 u 1:2:6 notitle'
  write(10,*) 'set output "./splt_bphi_ib4_Im.eps"'
  write(10,*) 'splot "./sp_bphi.data" index 3 u 1:2:5 notitle'
  write(10,*) 'set output "./splt_bphi_ib4_Re.eps"'
  write(10,*) 'splot "./sp_bphi.data" index 3 u 1:2:6 notitle'
  write(10,*) 'set output "./splt_bphi_ib5_Im.eps"'
  write(10,*) 'splot "./sp_bphi.data" index 4 u 1:2:5 notitle'
  write(10,*) 'set output "./splt_bphi_ib5_Re.eps"'
  write(10,*) 'splot "./sp_bphi.data" index 4 u 1:2:6 notitle'
  write(10,*) 'set output "./splt_bphi_ib6_Im.eps"'
  write(10,*) 'splot "./sp_bphi.data" index 5 u 1:2:5 notitle'
  write(10,*) 'set output "./splt_bphi_ib6_Re.eps"'
  write(10,*) 'splot "./sp_bphi.data" index 5 u 1:2:6 notitle'
end if
close(10)

call system("gnuplot ./ft_blo.plt")
END SUBROUTINE
SUBROUTINE plot_wf(zu,cur_kf,hav_kf,Et,At)
  use CONSTANTS, only : Nx,Nk,Nbt,dx
  use WAVE_FUNC, only : k
  implicit none
  real(8),intent(in)  ::  Et,At
  complex(8),intent(in) :: zu(0:Nx-1,Nk,Nbt)
  complex(8),intent(in) :: cur_kf(1:Nk), hav_kf(1:Nk)
  complex(8)  ::  zu_ib(Nbt)
  real(8) ::  zu_abs(Nbt),ss(Nk)
  integer :: ix, ib, ik
  open(9,file='td_wf.data',status='old',position='append')
  do ix = 0, Nx-1, 1
    do ib = 1, Nbt, 1
      zu_ib(ib)=sum(zu(ix,1:Nk,ib))/real(Nk)
      ss(:)=abs(zu(ix,:,ib))**2
      zu_abs(ib)=sum(ss(:))/real(Nk)
    end do
    write(9,'(i,e,e,e,e)') ix,sum(zu_ib(:))/real(Nbt), sum(zu_abs(:))/real(Nbt),Et*real(ix)*dx
  end do
  write(9,*) ""
  write(9,*) ""
  close(9)

  open(9,file='td_kf.data',status='old',position='append')
  do ik = 1, Nk, 1
    do ib = 1, Nbt, 1
      zu_abs(ib)=sum(abs(zu(:,ik,ib))**2)*dx
    end do
    write(9,'(i,e,e,e,e,e)') ik,k(ik),sum(abs(zu(0,ik,:)))/real(Nbt),real(cur_kf(ik)),real(hav_kf(ik)),At*k(ik)
  end do
  write(9,*) ""
  write(9,*) ""
  close(9)
END SUBROUTINE

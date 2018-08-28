SUBROUTINE gnuplot()
  implicit none
  integer,parameter :: isplot=0
open(10, file = 'ft_blo.plt', status = 'replace')
write(10,*) 'set terminal postscript eps enhanced color'
write(10,*) 'set output "./output/bphi.eps'
write(10,*) 'set ytics nomirror # 第一y軸の目盛は左側のみにする'
write(10,*) 'set y2tics         # 第二y軸を描画することを指定'
write(10,*) 'set grid           # 目盛線を描画することを指定'
write(10,*) 'set xlabel "ix"  # x軸の見出しを指定'
write(10,*) 'set ylabel "bphi"           # 第一y軸の見出しを指定'
write(10,*) 'set y2label "Potential" # 第二y軸の見出しを指定'
write(10,*) 'plot "./output/bphi.dat" index 1 u 2:3 w l title"n=0", \'
write(10,*) '     "./output/bphi.dat" index 2 u 2:3 w l title"n=1", \'
write(10,*) '     "./output/bphi.dat" index 3 u 2:3 w l title"n=2", \'
write(10,*) '      "./output/bphi.dat" index 0 u 2:3 w l axes x1y2 title"Pot"'
write(10,*) 'reset'
write(10,*) 'set output "./output/base.eps"'
write(10,*) 'set ylabel "Energy"'
write(10,*) 'set xlabel "k"'
write(10,*) 'plot "./output/base.dat" u 3:4 w l notitle'
write(10,*) 'reset'
write(10,*) 'set output "./output/td_norm.eps"'
write(10,*) 'set ytics nomirror'
write(10,*) 'set y2tics         # 第二y軸を描画することを指定'
write(10,*) 'set xlabel "time"'
write(10,*) 'set ylabel "norm at n=0"'
write(10,*) 'set y2label "norm at n=1"'
write(10,*) 'set grid'
write(10,*) 'plot "./output/td.dat" u 1:2 w l title "n=0",\'
write(10,*) ' "./output/td.dat" u 1:3 w l title "n=Nbt" axes x1y2'
write(10,*) 'reset'
write(10,*) 'set output "./output/td_energy.eps"'
write(10,*) 'set ytics nomirror'
write(10,*) 'set y2tics         # 第二y軸を描画することを指定'
write(10,*) 'set xlabel "time"'
write(10,*) 'set ylabel "Energy at n=0"'
write(10,*) 'set y2label "Energy at n=1"'
write(10,*) 'set grid'
write(10,*) 'plot "./output/td_V.dat" index 0 u 2:6 w l title "n=0",\'
write(10,*) '     "./output/td_V.dat" index 1 u 2:6 w l title "n=1" axes x1y2'
write(10,*) 'reset'
write(10,*) 'set output "./output/td_velocity.eps"'
write(10,*) 'set ytics nomirror'
write(10,*) 'set y2tics         # 第二y軸を描画することを指定'
write(10,*) 'set xlabel "time"'
write(10,*) 'set ylabel "Velocity at n=0"'
write(10,*) 'set y2label "Velocity at n=1"'
write(10,*) 'set grid'
write(10,*) 'plot "./output/td_V.dat" index 0 u 2:4 w l title "n=0",\'
write(10,*) '     "./output/td_V.dat" index 1 u 2:4 w l title "n=1" axes x1y2'
write(10,*) 'reset'
write(10,*) 'set terminal postscript eps enhanced color'
write(10,*) 'set logscale y'
write(10,*) 'set format y "10^{%L}"'
write(10,*) 'set xrange[0:100]'
write(10,*) 'set mxtics 5'
write(10,*) 'set grid'
write(10,*) 'set grid mxtics'
write(10,*) 'set xlabel "harmonic energy"'
write(10,*) 'set ylabel "Intensity"'
write(10,*) 'set output "./output/ft0.eps"'
write(10,*) 'plot "./output/fts.dat" index 0 u 4:5 w l title "n=0",\'
write(10,*) '     "./output/fts.dat" index 1 u 4:5 w l title "n=1"'
write(10,*) 'set output "./output/integ.eps"'
write(10,*) 'set xrange[0:300]'
write(10,*) 'plot   "./output/integ.dat" index 0 u 2:3 w lp pt 2,\'
write(10,*) '	"./output/integ.dat" index 1 u 2:3 w lp pt 2'
write(10,*) 'reset'
if( isplot==1 )then
  write(10,*) 'set output "./output/splt_bphi_ib1_Im.eps"'
  write(10,*) 'set pm3d map'
  write(10,*) 'set xlabel "x"'
  write(10,*) 'set ylabel "k-point"'
  write(10,*) 'splot "./output/sp_bphi.dat" index 0 u 1:2:5 notitle'
  write(10,*) 'set output "./output/splt_bphi_ib1_Re.eps"'
  write(10,*) 'splot "./output/sp_bphi.dat" index 0 u 1:2:6 notitle'
  write(10,*) 'set output "./output/splt_bphi_ib2_Im.eps"'
  write(10,*) 'splot "./output/sp_bphi.dat" index 1 u 1:2:5 notitle'
  write(10,*) 'set output "./output/splt_bphi_ib2_Re.eps"'
  write(10,*) 'splot "./output/sp_bphi.dat" index 1 u 1:2:6 notitle'
  write(10,*) 'set output "./output/splt_bphi_ib3_Im.eps"'
  write(10,*) 'splot "./output/sp_bphi.dat" index 2 u 1:2:5 notitle'
  write(10,*) 'set output "./output/splt_bphi_ib3_Re.eps"'
  write(10,*) 'splot "./output/sp_bphi.dat" index 2 u 1:2:6 notitle'
  write(10,*) 'set output "./output/splt_bphi_ib4_Im.eps"'
  write(10,*) 'splot "./output/sp_bphi.dat" index 3 u 1:2:5 notitle'
  write(10,*) 'set output "./output/splt_bphi_ib4_Re.eps"'
  write(10,*) 'splot "./output/sp_bphi.dat" index 3 u 1:2:6 notitle'
  write(10,*) 'set output "./output/splt_bphi_ib5_Im.eps"'
  write(10,*) 'splot "./output/sp_bphi.dat" index 4 u 1:2:5 notitle'
  write(10,*) 'set output "./output/splt_bphi_ib5_Re.eps"'
  write(10,*) 'splot "./output/sp_bphi.dat" index 4 u 1:2:6 notitle'
  write(10,*) 'set output "./output/splt_bphi_ib6_Im.eps"'
  write(10,*) 'splot "./output/sp_bphi.dat" index 5 u 1:2:5 notitle'
  write(10,*) 'set output "./output/splt_bphi_ib6_Re.eps"'
  write(10,*) 'splot "./output/sp_bphi.dat" index 5 u 1:2:6 notitle'
end if
close(10)

call system("gnuplot ./ft_blo.plt")
END SUBROUTINE

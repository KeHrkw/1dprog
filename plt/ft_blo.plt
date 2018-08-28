set terminal postscript eps enhanced color
set output "./output/bphi.eps"
set ytics nomirror # 第一y軸の目盛は左側のみにする
set y2tics         # 第二y軸を描画することを指定
set grid           # 目盛線を描画することを指定
set xlabel "ix"  # x軸の見出しを指定
set ylabel "bphi"           # 第一y軸の見出しを指定
set y2label "Potential" # 第二y軸の見出しを指定
plot "./output/bphi.dat" index 1 u 2:3 w l title"n=0", \
     "./output/bphi.dat" index 2 u 2:3 w l title"n=1", \
     "./output/bphi.dat" index 3 u 2:3 w l title"n=2", \
      "./output/bphi.dat" index 0 u 2:3 w l axes x1y2 title"Pot"
reset
set output "./output/base.eps"
set ylabel "Energy"
set xlabel "k"
plot "./output/base.dat" u 3:4 w l notitle
reset
set output "./output/td_norm.eps"
set ytics nomirror
set y2tics         # 第二y軸を描画することを指定
set xlabel "time"
set ylabel "norm at n=0"
set y2label "norm at n=1"
set grid
plot "./output/td.dat" u 1:2 w l title "n=0", "./output/td.dat" u 1:3 w l title "n=Nbt" axes x1y2
reset
set output "./output/td_energy.eps"
set ytics nomirror
set y2tics         # 第二y軸を描画することを指定
set xlabel "time"
set ylabel "Energy at n=0"
set y2label "Energy at n=1"
set grid
plot "./output/td_V.dat" index 0 u 2:6 w l title "n=0",\
     "./output/td_V.dat" index 1 u 2:6 w l title "n=1" axes x1y2
reset
set output "./output/td_velocity.eps"
set ytics nomirror
set y2tics         # 第二y軸を描画することを指定
set xlabel "time"
set ylabel "Velocity at n=0"
set y2label "Velocity at n=1"
set grid
plot "./output/td_V.dat" index 0 u 2:4 w l title "n=0",\
     "./output/td_V.dat" index 1 u 2:4 w l title "n=1" axes x1y2
reset
set terminal postscript eps enhanced color
set logscale y
set format y "10^{%L}"
set xrange[0:300]
set mxtics 5
set grid
set grid mxtics
set xlabel "harmonic order"
set ylabel "Intensity"
set output "./output/ft0.eps"
plot './output/fts.dat' index 0 u 4:5 w l title "n=0",\
     './output/fts.dat' index 1 u 4:5 w l title "n=1"
set output "./output/integ.eps"
plot   "./output/integ.dat" index 0 u 2:3 w lp pt 2,\
	"./output/integ.dat" index 1 u 2:3 w lp pt 2
reset
set output "./output/splt_bphi_ib1_Im.eps"
set pm3d map
set xlabel "x"
set ylabel "k-point"
splot "./output/sp_bphi.dat" index 0 u 1:2:5 notitle
set output "./output/splt_bphi_ib1_Re.eps"
splot "./output/sp_bphi.dat" index 0 u 1:2:6 notitle
set output "./output/splt_bphi_ib2_Im.eps"
splot "./output/sp_bphi.dat" index 1 u 1:2:5 notitle
set output "./output/splt_bphi_ib2_Re.eps"
splot "./output/sp_bphi.dat" index 1 u 1:2:6 notitle
set output "./output/splt_bphi_ib3_Im.eps"
splot "./output/sp_bphi.dat" index 2 u 1:2:5 notitle
set output "./output/splt_bphi_ib3_Re.eps"
splot "./output/sp_bphi.dat" index 2 u 1:2:6 notitle
set output "./output/splt_bphi_ib4_Im.eps"
splot "./output/sp_bphi.dat" index 3 u 1:2:5 notitle
set output "./output/splt_bphi_ib4_Re.eps"
splot "./output/sp_bphi.dat" index 3 u 1:2:6 notitle
set output "./output/splt_bphi_ib5_Im.eps"
splot "./output/sp_bphi.dat" index 4 u 1:2:5 notitle
set output "./output/splt_bphi_ib5_Re.eps"
splot "./output/sp_bphi.dat" index 4 u 1:2:6 notitle
set output "./output/splt_bphi_ib6_Im.eps"
splot "./output/sp_bphi.dat" index 5 u 1:2:5 notitle
set output "./output/splt_bphi_ib6_Re.eps"
splot "./output/sp_bphi.dat" index 5 u 1:2:6 notitle

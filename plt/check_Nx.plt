set terminal postscript color enhanced
set xlabel "time"
set ylabel "Velocity"
set grid
set output "td_velocity_n0.eps"
plot "./Nx20/Nk500/output/td_V.dat" index 0 u 2:4 w l linewidth 2 title "Nx=20 n=0",\
     "./Nx30/Nk500/output/td_V.dat" index 0 u 2:4 w l linewidth 2 title "Nx=30 n=0",\
     "./Nx50/Nk500/output/td_V.dat" index 0 u 2:4 w l linewidth 2 title "Nx=50 n=0",\
     "./Nx70/Nk500/output/td_V.dat" index 0 u 2:4 w l linewidth 2 title "Nx=70 n=0",\
     "./Nx80/Nk500/output/td_V.dat" index 0 u 2:4 w l linewidth 2 title "Nx=80 n=0"
set output "td_velocity_n1.eps"
plot "./Nx20/Nk500/output/td_V.dat" index 1 u 2:4 w l linewidth 2 title "Nx=20 n=1",\
     "./Nx30/Nk500/output/td_V.dat" index 1 u 2:4 w l linewidth 2 title "Nx=30 n=1",\
     "./Nx50/Nk500/output/td_V.dat" index 1 u 2:4 w l linewidth 2 title "Nx=50 n=1",\
     "./Nx70/Nk500/output/td_V.dat" index 1 u 2:4 w l linewidth 2 title "Nx=70 n=1",\
     "./Nx80/Nk500/output/td_V.dat" index 1 u 2:4 w l linewidth 2 title "Nx=80 n=1"
reset
set terminal postscript color enhanced
set xlabel "time"
set ylabel "Energy"
set grid
set output "td_energy_n0.eps"
plot "./Nx20/Nk500/output/td_V.dat" index 0 u 2:6 w l linewidth 2 title "Nx=20 n=0",\
     "./Nx30/Nk500/output/td_V.dat" index 0 u 2:6 w l linewidth 2 title "Nx=30 n=0",\
     "./Nx50/Nk500/output/td_V.dat" index 0 u 2:6 w l linewidth 2 title "Nx=50 n=0",\
     "./Nx70/Nk500/output/td_V.dat" index 0 u 2:6 w l linewidth 2 title "Nx=70 n=0",\
     "./Nx80/Nk500/output/td_V.dat" index 0 u 2:6 w l linewidth 2 title "Nx=80 n=0"
set output "td_energy_n1.eps"
plot "./Nx20/Nk500/output/td_V.dat" index 1 u 2:6 w l linewidth 2 title "Nx=20 n=1",\
     "./Nx30/Nk500/output/td_V.dat" index 1 u 2:6 w l linewidth 2 title "Nx=30 n=1",\
     "./Nx50/Nk500/output/td_V.dat" index 1 u 2:6 w l linewidth 2 title "Nx=50 n=1",\
     "./Nx70/Nk500/output/td_V.dat" index 1 u 2:6 w l linewidth 2 title "Nx=70 n=1",\
     "./Nx80/Nk500/output/td_V.dat" index 1 u 2:6 w l linewidth 2 title "Nx=80 n=1"
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
set output "ft_n0.eps"
plot "./Nx20/Nk500/output/fts.dat" index 0 u 4:5 w l linewidth 2 title "Nx=20 n=0",\
     "./Nx30/Nk500/output/fts.dat" index 0 u 4:5 w l linewidth 2 title "Nx=30 n=0",\
     "./Nx50/Nk500/output/fts.dat" index 0 u 4:5 w l linewidth 2 title "Nx=50 n=0",\
     "./Nx70/Nk500/output/fts.dat" index 0 u 4:5 w l linewidth 2 title "Nx=70 n=0",\
     "./Nx80/Nk500/output/fts.dat" index 0 u 4:5 w l linewidth 2 title "Nx=80 n=0"
    set output "ft_n1.eps"
    plot "./Nx20/Nk500/output/fts.dat" index 1 u 4:5 w l linewidth 2 title "Nx=20 n=1",\
         "./Nx30/Nk500/output/fts.dat" index 1 u 4:5 w l linewidth 2 title "Nx=30 n=1",\
         "./Nx50/Nk500/output/fts.dat" index 1 u 4:5 w l linewidth 2 title "Nx=50 n=1",\
         "./Nx70/Nk500/output/fts.dat" index 1 u 4:5 w l linewidth 2 title "Nx=70 n=1",\
         "./Nx80/Nk500/output/fts.dat" index 1 u 4:5 w l linewidth 2 title "Nx=80 n=1"

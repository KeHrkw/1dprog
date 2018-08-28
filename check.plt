set terminal postscript color enhanced
set xlabel "time"
set ylabel "Velocity"
set grid
set output "td_velocity_n0.eps"
plot "~/Velocity/output/td_V.dat" index 0 u 2:4 w l linewidth 2 title "Velocity n=0",\
      "~/Length/output/td_V.dat" index 0 u 2:4 w l linewidth 2 title "Length n=0"
set output "td_velocity_n1.eps"
plot "~/Velocity/output/td_V.dat"  index 1 u 2:4 w l linewidth 2 title "Velocity n=1",\
      "~/Length/output/td_V.dat"  index 1 u 2:4 w l linewidth 2 title "Length n=1"
reset
set terminal postscript color enhanced
set xlabel "time"
set ylabel "Energy"
set grid
set output "td_energy_n0.eps"
plot "~/Velocity/output/td_V.dat" index 0 u 2:6 w l linewidth 2 title "Velocity n=0",\
      "~/Length/output/td_V.dat" index 0 u 2:6 w l linewidth 2 title "Length n=0"
set output "td_energy_n1.eps"
plot "~/Velocity/output/td_V.dat"  index 1 u 2:6 w l linewidth 2 title "Velocity n=1",\
      "~/Length/output/td_V.dat"  index 1 u 2:6 w l linewidth 2 title "Length n=1"
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
plot '~/Velocity/output/fts.dat' index 0 u 4:5 w l linewidth 2 title "Velocity n=0",\
    '~/Length/output/fts.dat' index 0 u 4:5 w l linewidth 2 title "Length n=0"
    set output "ft_n1.eps"
    plot '~/Velocity/output/fts.dat' index 1 u 4:5 w l linewidth 2 title "Velocity n=1",\
         '~/Length/output/fts.dat' index 1 u 4:5 w l linewidth 2 title "Length n=1"

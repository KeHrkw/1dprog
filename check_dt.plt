set terminal postscript eps enhanced color
set grid           # 目盛線を描画することを指定
set xlabel "time"  # x軸の見出しを指定
set ylabel "Velofity" # 第二y軸の見出しを指定
set key top left   # 凡例は左上に描画
set output "./output/check_velo0.eps"
plot "~/Bloch/test_dt/td.dat" u 1:4 w l title "Bloch",\
     "~/testsch/test_dt/td.dat" u 1:4 w l title "Periodic"
     set output "./output/check_velo1.eps"
     plot "~/Bloch/test_dt/td.dat" u 1:5 w l title "Bloch",\
          "~/testsch/test_dt/td.dat" u 1:5 w l title "Periodic"
set output "./output/check_Ene.eps"
set ylabel "Energy"           # 第一y軸の見出しを指定
plot "~/Bloch/test_dt/td.dat" u 1:3 w l title "Bloch",\
     "~/testsch/test_dt/td.dat" u 1:3 w l title "Periodic"
reset
set logscale y
set format y "10^{%L}"
set xrange[0:150]
set mxtics 5
set grid
set grid mxtics
set terminal postscript eps enhanced color
set output "./output/check_ft0.eps"
plot "~/Bloch/test_dt/fts.dat" index 0 u 4:5 w l ,\
     "~/testsch/test_dt/fts.dat" index 0 u 3:4 w l
set output "./output/check_ft1.eps"
plot "~/Bloch/test_dt/fts.dat" index 1 u 4:5 w l ,\
     "~/testsch/test_dt/fts.dat" index 1 u 3:4 w l

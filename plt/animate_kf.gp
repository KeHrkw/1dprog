#-------------------------------------------------------------------------------
# gnuplotの設定
#-------------------------------------------------------------------------------
reset
set nokey                 # 凡例の非表示
set xlabel "k"
set ylabel "abs(u(ik))"
set xrange [-0.4:0.4]
set yrange [2.2:2.7]       # y軸方向の範囲の設定
set size square           # 図を正方形にする

set term gif animate      # 出力をgifアニメに設定
set output "sample_kf.gif"  # 出力ファイル名の設定

#-------------------------------------------------------------------------------
# 変数の設定
#-------------------------------------------------------------------------------
n0 = 0    # ループ変数の初期値
n1 = 131   # ループ変数の最大値
dn = 1    # ループ変数の増加間隔

delay = 70

#-------------------------------------------------------------------------------
# ループの開始
#-------------------------------------------------------------------------------
load "animate_kf.plt"

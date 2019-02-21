#-------------------------------------------------------------------------------
# ループ処理
#-------------------------------------------------------------------------------
if(exist("n")==0 || n<0) n = n0  # ループ変数の初期化

#-------------------------------------------------------------------------------
# プロット
#-------------------------------------------------------------------------------
plot "td_kf.data"  index n using 2:3 with lines lw 4
#-------------------------------------------------------------------------------
# ループ処理
#-------------------------------------------------------------------------------
n = n + dn            # ループ変数の増加
if ( n < n1 ) reread  # ループの評価
undefine n            # ループ変数の削除